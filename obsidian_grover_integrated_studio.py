import datetime
import json
import os
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

try:
    from pypdf import PdfReader
    HAS_PYPDF = True
except ImportError:
    HAS_PYPDF = False

try:
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.metrics.pairwise import cosine_similarity
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

STOP = {
    "der", "die", "das", "und", "oder", "ein", "eine", "einer", "einem", "einen",
    "ist", "sind", "war", "wie", "was", "wer", "welche", "welcher", "welches",
    "nicht", "mit", "von", "zu", "im", "in", "auf", "für", "fuer", "den", "dem",
    "des", "dass", "daß", "auch", "als", "an", "es", "er", "sie", "ich", "wir",
    "man", "wird", "werden", "kann", "können", "koennen", "hat", "haben", "nach",
    "bei", "aus", "über", "ueber", "unter", "wenn", "dann", "nur", "noch", "schon",
    "sehr", "mehr", "diese", "dieser", "dieses", "zum", "zur", "vom", "am", "sich",
    "the", "a", "an", "and", "or", "of", "to", "in", "is", "are", "was", "what",
    "how", "why", "which", "with", "for", "from", "that", "this",
}

SUGGESTIONS = {
    "Was steht zur Perihelpräzession im Vault?": "Was steht zur Perihelpräzession im Vault?",
    "Welche Zettel behandeln E8?": "Welche Zettel behandeln E8?",
    "Was ist zum Quantenspeicher festgehalten?": "Was ist zum Quantenspeicher festgehalten?",
}

NOTE_EXTENSIONS = {".md", ".txt", ".pdf"}


def query_terms(query):
    return [t for t in re.findall(r"\w+", query.lower()) if t not in STOP and len(t) > 1]


def body_text(content):
    if content.startswith("---"):
        parts = content.split("---", 2)
        if len(parts) >= 3:
            return parts[2].strip()
    return content


def vault_fingerprint(path):
    root = Path(path)
    if not root.exists():
        return "missing"
    count = 0
    latest = 0.0
    for file_path in root.rglob("*"):
        if file_path.is_file() and file_path.suffix.lower() in NOTE_EXTENSIONS:
            count += 1
            latest = max(latest, file_path.stat().st_mtime)
    return f"{count}:{latest:.0f}"


def scan_obsidian_vault(path):
    notes = []
    if not os.path.exists(path):
        return notes
    for root, _, files in os.walk(path):
        for file in files:
            ext = os.path.splitext(file)[1].lower()
            full_path = os.path.join(root, file)
            content = ""
            if ext in [".md", ".txt"]:
                try:
                    with open(full_path, "r", encoding="utf-8", errors="ignore") as handle:
                        content = handle.read()
                except OSError:
                    continue
            elif ext == ".pdf" and HAS_PYPDF:
                try:
                    reader = PdfReader(full_path)
                    content = "\n".join([page.extract_text() or "" for page in reader.pages])
                except Exception:
                    continue

            if content.strip():
                stem = os.path.splitext(file)[0]
                wiki_links = set(re.findall(r"\[\[([^\]\|]+)(?:\|[^\]]+)?\]\]", content))
                tags = list(set(re.findall(r"#[\w/-]+", content)))
                notes.append({
                    "title": stem,
                    "filename": file,
                    "path": full_path,
                    "content": content,
                    "links": wiki_links,
                    "tags": tags,
                    "word_count": len(content.split()),
                })
    return notes


@st.cache_data(ttl="5m", show_spinner="Vault wird gelesen")
def load_vault(path, fingerprint):
    """fingerprint gehört zum Cache-Schlüssel und wechselt, wenn Dateien sich ändern."""
    if not fingerprint or fingerprint == "missing":
        return []
    return scan_obsidian_vault(path)


def grover_search_simulation(scores, threshold_val, num_iterations=None):
    n_items = len(scores)
    if n_items == 0:
        return np.array([]), 0, np.array([])

    targets = np.where(scores >= threshold_val)[0]
    marked = len(targets)
    if marked == 0:
        targets = np.array([int(np.argmax(scores))])
        marked = 1

    state = np.ones(n_items) / np.sqrt(n_items)
    opt_k = max(1, int(round((np.pi / 4.0) * np.sqrt(n_items / marked))))
    iterations = num_iterations if num_iterations is not None else opt_k

    oracle = np.ones(n_items)
    oracle[targets] = -1.0

    for _ in range(iterations):
        state = state * oracle
        mean_amp = np.mean(state)
        state = 2 * mean_amp - state

    return state ** 2, iterations, targets


def analyze_vault_gaps(notes, gap_threshold):
    if not notes or not HAS_SKLEARN:
        return {"missing_links": [], "orphans": []}

    corpus = [note["content"] for note in notes]
    try:
        vectorizer = TfidfVectorizer(stop_words=sorted(STOP), min_df=1)
        tfidf_matrix = vectorizer.fit_transform(corpus)
        sim_matrix = cosine_similarity(tfidf_matrix)
    except Exception:
        return {"missing_links": [], "orphans": []}

    missing_links = []
    orphans = []

    for i, note in enumerate(notes):
        title_a = note["title"]
        links_a = note["links"]
        incoming = sum(1 for other in notes if title_a in other["links"])
        if incoming == 0 and len(links_a) == 0:
            orphans.append(title_a)

        for j in range(i + 1, len(notes)):
            title_b = notes[j]["title"]
            similarity = float(sim_matrix[i, j])
            is_linked = (title_b in links_a) or (title_a in notes[j]["links"])
            if similarity >= gap_threshold and not is_linked:
                missing_links.append({
                    "note_a": title_a,
                    "note_b": title_b,
                    "similarity": similarity,
                })

    missing_links.sort(key=lambda item: item["similarity"], reverse=True)
    return {"missing_links": missing_links, "orphans": orphans}


@st.cache_data(ttl="5m", show_spinner=False)
def cached_gaps(path, fingerprint, gap_threshold):
    notes = load_vault(path, fingerprint)
    return analyze_vault_gaps(notes, gap_threshold)


def keyword_scores(query, notes):
    terms = query_terms(query)
    scores = np.zeros(len(notes))
    if not terms:
        return scores
    for index, note in enumerate(notes):
        haystack = f"{note['title']}\n{note['content']}".lower()
        scores[index] = sum(haystack.count(term) for term in terms)
    peak = float(scores.max())
    if peak > 0:
        scores = scores / peak
    return scores


def score_notes(query, notes, search_mode):
    if not notes:
        return np.array([])
    if search_mode.startswith("Semantisch") and HAS_SKLEARN:
        corpus = [f"{note['title']}\n{note['content']}" for note in notes]
        try:
            vectorizer = TfidfVectorizer(stop_words=sorted(STOP), min_df=1)
            matrix = vectorizer.fit_transform(corpus)
            sims = cosine_similarity(vectorizer.transform([query]), matrix).flatten()
            return np.asarray(sims, dtype=float)
        except Exception:
            return keyword_scores(query, notes)
    return keyword_scores(query, notes)


def best_excerpt(note, query, limit=700):
    text = body_text(note["content"]).strip()
    terms = query_terms(query)
    chunks = [chunk.strip() for chunk in re.split(r"\n\s*\n", text) if chunk.strip()]
    if not chunks:
        chunks = [text] if text else []

    def chunk_score(chunk):
        lowered = chunk.lower()
        return sum(lowered.count(term) for term in terms)

    ranked = sorted(chunks, key=chunk_score, reverse=True)
    parts = []
    size = 0
    for chunk in ranked:
        if terms and chunk_score(chunk) == 0 and parts:
            break
        piece = chunk[:500]
        parts.append(piece)
        size += len(piece)
        if len(parts) >= 2 or size >= limit:
            break
    excerpt = "\n\n".join(parts).strip()
    excerpt = re.sub(r"(?m)^#{1,6}\s+", "", excerpt or text)
    return excerpt.strip()[:limit]


def retrieve_grover(query, notes, search_mode, threshold, gap_data):
    scores = score_notes(query, notes, search_mode)
    if len(scores) == 0 or float(np.max(scores)) <= 0:
        return {
            "empty": True,
            "hits": [],
            "chart": [],
            "iterations": 0,
            "target_count": 0,
            "note_count": len(notes),
            "fallback": False,
        }

    raw_targets = np.where(scores >= threshold)[0]
    probabilities, iterations, targets = grover_search_simulation(scores, threshold)
    target_set = {int(index) for index in targets}
    order = [int(index) for index in np.argsort(scores)[::-1] if int(index) in target_set]
    if not order:
        order = [int(index) for index in np.argsort(probabilities)[::-1][:3]]

    orphans = set(gap_data.get("orphans", []))
    hits = []
    for index in order[:5]:
        note = notes[index]
        relatives = []
        for link in gap_data.get("missing_links", []):
            if link["note_a"] == note["title"]:
                relatives.append(link["note_b"])
            elif link["note_b"] == note["title"]:
                relatives.append(link["note_a"])
        hits.append({
            "title": note["title"],
            "path": note["path"],
            "probability": float(probabilities[index]),
            "score": float(scores[index]),
            "excerpt": best_excerpt(note, query),
            "tags": note["tags"],
            "orphan": note["title"] in orphans,
            "missing": relatives[:3],
        })

    chart_order = [int(index) for index in np.argsort(probabilities)[::-1][:12]]
    seen = set()
    chart = []
    for index in chart_order:
        title = notes[index]["title"]
        label = title if title not in seen else f"{title} ({index})"
        seen.add(label)
        chart.append({
            "Zettel": label[:48],
            "Wahrscheinlichkeit": float(probabilities[index]),
            "Farbe": "#ff7f0e" if index in target_set else "#4c78a8",
        })

    return {
        "empty": False,
        "hits": hits,
        "chart": chart,
        "iterations": int(iterations),
        "target_count": int(len(targets)),
        "note_count": len(notes),
        "fallback": len(raw_targets) == 0,
    }


def format_grover_answer(query, result):
    if result["empty"]:
        return (
            f"Zur Frage „{query}“ verstärkt Grover keinen Zettel. "
            "Im Vault gibt es dazu keinen Treffer."
        )
    lines = [f"Zur Frage „{query}“:"]
    if result["fallback"]:
        lines.append("")
        lines.append(
            "Kein Zettel erreicht den Orakel-Schwellenwert. "
            "Die Antwort stützt sich auf den nächsten Treffer."
        )
    for hit in result["hits"][:3]:
        lines.extend(["", f"**[[{hit['title']}]]** · {hit['probability'] * 100:.1f} %", "", hit["excerpt"]])
    return "\n".join(lines).strip()


def grover_prompt(query, result):
    context = "\n\n---\n\n".join(
        f"Zettel: [[{hit['title']}]]\n{hit['excerpt']}" for hit in result["hits"][:3]
    )
    return f"""Du antwortest aus einem Zettelkasten. Die Auszüge sind die vom Grover-Orakel verstärkten Treffer.

Beantworte die Frage auf Deutsch, knapp und direkt. Stütze jede Aussage auf die Auszüge und nenne die Zettel als [[Titel]]. Wenn die Auszüge die Frage nicht tragen, sage das in einem Satz.

FRAGE:
{query}

VERSTÄRKTE ZETTEL:
{context}
"""


def stream_ollama(prompt, model):
    request = urllib.request.Request(
        "http://localhost:11434/api/generate",
        data=json.dumps({"model": model, "prompt": prompt, "stream": True}).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        for raw_line in response:
            if not raw_line.strip():
                continue
            payload = json.loads(raw_line.decode("utf-8"))
            piece = payload.get("response", "")
            if piece:
                yield piece
            if payload.get("done"):
                break


def count_label(count, singular, plural):
    return f"{count} {singular if count == 1 else plural}"


def meta_line(result):
    if result["empty"]:
        return f"Grover-Antwort · {count_label(result['note_count'], 'Zettel', 'Zettel')} durchsucht · kein Treffer"
    return (
        f"Grover-Antwort · {count_label(result['iterations'], 'Iteration', 'Iterationen')} · "
        f"{count_label(result['target_count'], 'Orakel-Treffer', 'Orakel-Treffer')} "
        f"von {result['note_count']} Zetteln"
    )


def render_sources(result):
    if result["empty"]:
        return
    with st.expander(f"Verstärkte Zettel ({len(result['hits'])})"):
        for hit in result["hits"]:
            st.markdown(f"**{hit['title']}** · {hit['probability'] * 100:.1f} %")
            st.caption(hit["path"])
            if hit["tags"]:
                st.caption(" ".join(hit["tags"]))
            if hit["orphan"]:
                st.caption("Verwaiste Notiz")
            if hit["missing"]:
                linked = ", ".join(f"[[{title}]]" for title in hit["missing"])
                st.caption(f"Fehlende Querverbindungen: {linked}")
            st.markdown(hit["excerpt"])
        if result["chart"]:
            frame = pd.DataFrame(result["chart"])
            st.bar_chart(
                frame,
                x="Zettel",
                y="Wahrscheinlichkeit",
                color="Farbe",
                sort=False,
            )


def save_answer(vault_path, query, content, hits):
    stamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    safe_query = re.sub(r"[^\w\-]+", "_", query, flags=re.UNICODE).strip("_")[:30] or "frage"
    filename = f"{stamp}_Grover_{safe_query}.md"
    destination = os.path.join(vault_path, filename)
    parent_links = " ".join(f"[[{hit['title']}]]" for hit in hits[:3])
    with open(destination, "w", encoding="utf-8") as handle:
        handle.write("---\n")
        handle.write(f"title: Grover-Antwort zu {query}\n")
        handle.write(f"created: {datetime.datetime.now().isoformat()}\n")
        handle.write("tags: [ai-generated, zettelkasten, grover-studio]\n")
        handle.write("---\n\n")
        handle.write(f"# Grover-Antwort: {query}\n\n")
        if parent_links:
            handle.write(f"**Verknüpfte Quellnotizen:** {parent_links}\n\n")
        handle.write(content)
        handle.write("\n")
    return filename


st.set_page_config(
    page_title="Grover-Antworten",
    page_icon="🕸️",
    layout="wide",
)

st.title("Grover-Antworten")
st.caption("Stelle eine Frage. Grover verstärkt die passenden Zettel und antwortet daraus.")

st.session_state.setdefault("messages", [])
st.session_state.setdefault("pending_question", None)

st.sidebar.header("Obsidian Vault")
default_vault = os.path.expanduser("~/Documents/zettel")
vault_path = st.sidebar.text_input("Pfad zum Obsidian Vault", value=default_vault, key="vault_path")

if st.session_state.get("vault_path_seen") != vault_path:
    st.session_state.vault_path_seen = vault_path
    st.session_state.messages = []
    st.session_state.pending_question = None

st.sidebar.header("Suche")
search_mode = st.sidebar.segmented_control(
    "Suchmodus",
    ["Semantisch (TF-IDF)", "Stichwort"],
    default="Semantisch (TF-IDF)",
    required=True,
    key="search_mode",
)
if search_mode is None:
    search_mode = "Semantisch (TF-IDF)"
threshold = st.sidebar.slider("Grover-Schwellenwert", 0.1, 0.9, 0.30, step=0.05, key="threshold")
gap_sim_threshold = st.sidebar.slider(
    "Missing-Link-Schwellenwert", 0.15, 0.80, 0.25, step=0.05, key="gap_threshold"
)

st.sidebar.header("Antwort")
llm_provider = st.sidebar.selectbox(
    "Formulierung",
    ["Lokales Ollama", "Direkt aus den Zetteln"],
    key="llm_provider",
)
ollama_model = st.sidebar.text_input("Ollama-Modell", value="llama3", key="ollama_model")
use_ollama = llm_provider == "Lokales Ollama"

fingerprint = vault_fingerprint(vault_path)
notes = load_vault(vault_path, fingerprint)
gap_data = cached_gaps(vault_path, fingerprint, gap_sim_threshold) if notes else {"missing_links": [], "orphans": []}

st.sidebar.markdown("---")
st.sidebar.markdown(
    f"**{len(notes)}** Notizen · **{len(gap_data.get('missing_links', []))}** fehlende Links · "
    f"**{len(gap_data.get('orphans', []))}** verwaiste Zettel"
)
if not HAS_SKLEARN and search_mode.startswith("Semantisch"):
    st.sidebar.caption("TF-IDF ist nicht installiert. Die Suche läuft über Stichworte.")

if st.sidebar.button("Vault neu einlesen", icon=":material/refresh:"):
    load_vault.clear()
    cached_gaps.clear()
    st.rerun()

if notes and st.sidebar.button("Gap-Report schreiben", icon=":material/link:"):
    report_path = os.path.join(vault_path, "00_Knowledge_Gap_Report.md")
    try:
        with open(report_path, "w", encoding="utf-8") as handle:
            handle.write("# Knowledge Gap Report\n\n")
            handle.write(f"**Analysierte Notizen:** {len(notes)}\n\n")
            handle.write("## Empfohlene Querverbindungen\n")
            handle.write("| Ähnlichkeit | Notiz A | Notiz B |\n|---|---|---|\n")
            for link in gap_data.get("missing_links", [])[:15]:
                handle.write(
                    f"| {link['similarity'] * 100:.1f}% | [[{link['note_a']}]] | [[{link['note_b']}]] |\n"
                )
            handle.write("\n## Verwaiste Notizen\n")
            for orphan in gap_data.get("orphans", []):
                handle.write(f"- [[{orphan}]]\n")
        st.sidebar.success("Gespeichert als `00_Knowledge_Gap_Report.md`.")
    except OSError as error:
        st.sidebar.error(f"Bericht konnte nicht geschrieben werden: {error}")

if st.sidebar.button("Gespräch leeren", icon=":material/delete:"):
    st.session_state.messages = []
    st.session_state.pending_question = None
    st.rerun()

question = st.session_state.pop("pending_question", None)
typed = st.chat_input(
    "Frage an den Zettelkasten",
    submit_mode="disable",
    disabled=not notes,
)
if typed:
    question = typed.strip()

if not notes:
    st.info(f"Im Pfad `{vault_path}` liegen noch keine Notizen. Trag den Vault links ein.")
elif not st.session_state.messages and question is None:
    st.caption("Beispielfragen")
    selected = st.pills("Beispielfragen", list(SUGGESTIONS), label_visibility="collapsed")
    if selected:
        st.session_state.pending_question = SUGGESTIONS[selected]
        st.rerun()

for index, message in enumerate(st.session_state.messages):
    with st.chat_message(message["role"]):
        if message["role"] == "assistant":
            st.caption(message.get("meta", "Grover-Antwort"))
        st.markdown(message["content"])
        if message["role"] == "assistant" and message.get("result"):
            render_sources(message["result"])
            is_last = index == len(st.session_state.messages) - 1
            if is_last and st.button("Antwort als Zettel speichern", key=f"save_{index}", icon=":material/save:"):
                try:
                    filename = save_answer(
                        vault_path,
                        message.get("query", "Frage"),
                        message["content"],
                        message["result"].get("hits", []),
                    )
                    st.success(f"Gespeichert als `{filename}`.")
                except OSError as error:
                    st.error(f"Zettel konnte nicht gespeichert werden: {error}")

if question and notes:
    retrieval_query = question
    previous_questions = [item["content"] for item in st.session_state.messages if item["role"] == "user"]
    if previous_questions and len(query_terms(question)) < 3:
        retrieval_query = f"{previous_questions[-1]} {question}"

    st.session_state.messages.append({"role": "user", "content": question})
    with st.chat_message("user"):
        st.markdown(question)

    with st.chat_message("assistant"):
        with st.status("Zettel bewerten", type="step") as scoring:
            result = retrieve_grover(retrieval_query, notes, search_mode, threshold, gap_data)
            if result["empty"]:
                scoring.update(label="Kein passender Zettel", state="complete")
            elif result["fallback"]:
                scoring.update(label="Unter dem Schwellenwert, nächster Zettel", state="complete")
            else:
                scoring.update(
                    label=f"{result['target_count']} Zettel über dem Schwellenwert",
                    state="complete",
                )

        with st.status("Grover-Verstärkung", type="step") as amplification:
            if result["empty"]:
                amplification.update(label="Keine Amplitude verstärkt", state="complete")
            else:
                amplification.update(
                    label=(
                        f"{count_label(result['iterations'], 'Iteration', 'Iterationen')}, "
                        f"{count_label(len(result['hits']), 'Treffer', 'Treffer')}"
                    ),
                    state="complete",
                )

        st.caption(meta_line(result))
        if use_ollama and not result["empty"]:
            prompt = grover_prompt(question, result)
            try:
                answer = st.write_stream(stream_ollama(prompt, ollama_model))
                if not str(answer).strip():
                    answer = format_grover_answer(question, result)
                    st.markdown(answer)
            except Exception as error:
                answer = format_grover_answer(question, result)
                answer += f"\n\nOllama hat nicht geantwortet ({error}). Gezeigt sind die verstärkten Stellen."
                st.markdown(answer)
        else:
            answer = format_grover_answer(question, result)
            st.markdown(answer)

        render_sources(result)

    st.session_state.messages.append({
        "role": "assistant",
        "content": answer,
        "meta": meta_line(result),
        "query": question,
        "result": result,
    })
    st.rerun()
