import numpy as np
import re
from pathlib import Path

# ---------- small prime helper ----------
def first_n_primes(n):
    n = int(n)
    primes = []
    x = 2
    while len(primes) < n:
        is_p = True
        r = int(np.sqrt(x))
        for p in primes:
            if p > r:
                break
            if x % p == 0:
                is_p = False
                break
        if is_p:
            primes.append(x)
        x += 1 if x == 2 else 2
    return np.array(primes, dtype=np.int64)

def eabc_class_mod12(p):
    r = int(p % 12)
    return r if r in (1, 5, 7, 11) else None

# ---------- BM/EABC predicate (you can replace this!) ----------
def bm_predicate_factory(primes, mode="threads_v1"):
    """
    Returns a predicate(i)->bool based on primes[i], primes[i+1], ...
    mode:
      - "threads_v1": marks indices with small gaps and certain mod-12 patterns
    """
    primes = np.asarray(primes, dtype=np.int64)

    if mode == "threads_v1":
        def predicate(i):
            p  = int(primes[i])
            p1 = int(primes[i+1])
            gap = p1 - p

            c  = eabc_class_mod12(p)
            c1 = eabc_class_mod12(p1)

            # "Faden": very small gaps, plus optional class pattern
            is_thread = (gap <= 6)  # catches twins (2), cousin (4), sexy (6)
            class_pattern = (c in (1,11)) and (c1 in (1,11))  # example pattern

            # mark if either "thread" or special pattern triggers
            return bool(is_thread or class_pattern)

        return predicate

    raise ValueError("Unknown mode")


# ---------- Textverzeichnis: Dokumente laden und Ähnlichkeit ----------
def get_text_files(dir_path, extensions=(".txt", ".md", ".rtf")):
    """Sammelt Textdateien in einem Verzeichnis (rekursiv optional ausgeschaltet)."""
    dir_path = Path(dir_path).expanduser().resolve()
    if not dir_path.is_dir():
        raise FileNotFoundError(f"Verzeichnis nicht gefunden: {dir_path}")
    files = []
    for ext in extensions:
        files.extend(dir_path.glob(f"*{ext}"))
    return sorted([str(p) for p in files])


def read_text(file_path, encodings=("utf-8", "utf-8-sig", "latin-1")):
    """Liest Text aus einer Datei mit Fallback-Encodings (macOS-kompatibel)."""
    path = Path(file_path)
    if not path.is_file():
        return ""
    raw = path.read_bytes()
    for enc in encodings:
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def tokenize(text):
    """Einfache Tokenisierung: Kleinbuchstaben, Wörter aus Buchstaben/Ziffern."""
    text = (text or "").lower()
    return re.findall(r"[a-zäöüß0-9]+", text)


def word_counts(tokens):
    """Zählt Tokens -> dict (Wort -> Anzahl)."""
    out = {}
    for t in tokens:
        out[t] = out.get(t, 0) + 1
    return out


def cosine_similarity(vec_a, vec_b):
    """Ähnlichkeit zweier Bag-of-Words als Kosinus (0..1). Keys = alle Keys aus beiden."""
    if not vec_a or not vec_b:
        return 0.0
    keys = set(vec_a) | set(vec_b)
    a = np.array([vec_a.get(k, 0) for k in keys], dtype=np.float64)
    b = np.array([vec_b.get(k, 0) for k in keys], dtype=np.float64)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def build_doc_vectors(file_paths):
    """Liefert Liste von (path, word_count_dict) für jede Datei."""
    return [(p, word_counts(tokenize(read_text(p)))) for p in file_paths]


def marked_by_similarity(doc_vectors, query_vector, threshold=0.15):
    """
    Markiert Indizes, deren Dokument-Ähnlichkeit zum Query >= threshold ist.
    query_vector = word_counts dict (z. B. von einem Query-Dokument).
    """
    marked = []
    for path, vec in doc_vectors:
        sim = cosine_similarity(vec, query_vector)
        marked.append(sim >= threshold)
    return np.array(marked, dtype=bool)


# ---------- Grover core ----------
def grover_bm_search(N=1024, shots=2000, seed=1, oracle_mode="threads_v1"):
    """
    Grover-Suche auf BM/EABC-Prädikat (Primzahlen, Faden-Kriterien).
    Die Implementierung kann anhand der Phonizität (Korrektheit) geprüft werden:
    theoretische P(marked), Top-Indizes in marked, empirische Trefferquote.
    """
    if N & (N - 1) != 0:
        raise ValueError("N must be a power of 2 (e.g., 1024).")

    rng = np.random.default_rng(seed)

    # We need enough primes for predicate look-ahead i+1
    primes = first_n_primes(N + 10)
    predicate = bm_predicate_factory(primes, mode=oracle_mode)

    marked = np.array([predicate(i) for i in range(N)], dtype=bool)
    M = int(marked.sum())
    if M == 0:
        raise RuntimeError("Oracle marked set is empty; adjust predicate.")

    # init uniform superposition
    psi = np.ones(N, dtype=np.complex128) / np.sqrt(N)

    # Grover iteration count for M solutions
    k = int(np.floor((np.pi / 4.0) * np.sqrt(N / M)))
    k = max(1, k)

    print(f"Grover BM-search: N={N}, marked M={M}, iterations k={k}")

    def oracle(state):
        state[marked] *= -1
        return state

    def diffusion(state):
        mean_amp = np.mean(state)
        return 2.0 * mean_amp - state

    for _ in range(k):
        psi = oracle(psi)
        psi = diffusion(psi)

    probs = np.abs(psi)**2

    # ---------- Phonizität / Korrektheit ----------
    # Theoretische Erfolgswahrscheinlichkeit (Grover nach k Iterationen)
    prob_marked_theory = np.sum(probs[marked])
    print(f"Theoretische P(marked) nach k={k} Iterationen: {prob_marked_theory:.4f}")

    top = np.argsort(probs)[-12:][::-1]
    top_marked_count = int(np.sum(marked[top]))
    print(f"Phonizität: {top_marked_count}/12 der Top-Indizes sind markiert.")
    # Strenge Prüfung nur bei hoher Verstärkung (Grover-Regime M << N)
    if prob_marked_theory > 0.5 and top_marked_count < 12:
        bad = top[~marked[top]]
        raise AssertionError(f"Phonizität verletzt: P(marked)={prob_marked_theory:.4f}, aber Top-Indizes {bad.tolist()} nicht in marked.")

    # empirical sampling
    samples = rng.choice(N, size=shots, p=probs)
    hit_mask = marked[samples]
    hits = int(hit_mask.sum())

    print(f"Empirische Trefferquote (shots={shots}): {hits/shots:.4f}")
    print("Top outcomes (index, prob, marked, prime, gap):")
    for idx in top:
        p = int(primes[idx])
        gap = int(primes[idx+1] - primes[idx])
        print(f"  {idx:4d}  {probs[idx]:.6f}  marked={bool(marked[idx])}  p={p}  gap={gap}")

    return {
        "probs": probs,
        "samples": samples,
        "marked": marked,
        "primes": primes,
        "k": k,
        "M": M,
        "prob_marked_theory": prob_marked_theory,
    }


# ---------- Grover für Textverzeichnisse (ähnliche Dokumente) ----------
def grover_text_search(
    dir_path,
    query_path=None,
    extensions=(".txt", ".md", ".rtf"),
    similarity_threshold=0.15,
    N=1024,
    shots=2000,
    seed=1,
    top_display=12,
):
    """
    Grover-Suche in einem Textverzeichnis: findet Dokumente, die einem
    Query-Dokument ähnlich sind (Bag-of-Words, Kosinus-Ähnlichkeit).
    Für normale Textverzeichnisse auf dem Mac (z. B. ~/Documents, ~/Desktop).
    Nutzung: ähnliche Dokumente finden, um z. B. Vorlagen für neue, inhaltlich
    ähnliche Dokumente zu haben.

    - dir_path: Pfad zum Verzeichnis mit Textdateien
    - query_path: Pfad zu einem Referenzdokument (ähnliche werden gesucht).
      Wenn None: erstes Dokument im Verzeichnis dient als Query.
    - similarity_threshold: Mindest-Kosinus-Ähnlichkeit zum Query (0..1)
    - N: Anzahl Slots (muss Zweierpotenz sein); Verzeichnis wird ggf. aufgefüllt
    """
    if N & (N - 1) != 0:
        raise ValueError("N must be a power of 2 (e.g., 256, 1024).")

    rng = np.random.default_rng(seed)
    paths = get_text_files(dir_path, extensions=extensions)
    if not paths:
        raise FileNotFoundError(f"Keine Textdateien in {dir_path} (Erweiterungen: {extensions})")

    # Auf N Slots bringen: Auffüllen mit Dummy-Pfaden (leere Dokumente)
    while len(paths) < N:
        paths.append(None)  # Dummy = leeres Dokument
    paths = paths[:N]
    doc_vectors = []
    for p in paths:
        if p is not None:
            doc_vectors.append((p, word_counts(tokenize(read_text(p)))))
        else:
            doc_vectors.append((None, {}))

    # Query-Vektor: von query_path oder vom ersten echten Dokument
    if query_path:
        query_vector = word_counts(tokenize(read_text(query_path)))
    else:
        query_vector = doc_vectors[0][1] if doc_vectors[0][0] else {}

    # Ähnlichkeiten und Markierung
    similarities = np.array(
        [cosine_similarity(vec, query_vector) for _, vec in doc_vectors],
        dtype=np.float64,
    )
    marked = similarities >= similarity_threshold
    M = int(marked.sum())
    if M == 0:
        # Fallback: markiere die top-k nach Roh-Ähnlichkeit (z. B. top 10 %)
        k_mark = max(1, N // 10)
        top_idx = np.argsort(similarities)[-k_mark:]
        marked = np.zeros(N, dtype=bool)
        marked[top_idx] = True
        M = int(marked.sum())
        print(f"Hinweis: Kein Dokument >= {similarity_threshold}; markiere die {M} ähnlichsten.")

    # Uniform superposition
    psi = np.ones(N, dtype=np.complex128) / np.sqrt(N)
    k = max(1, int(np.floor((np.pi / 4.0) * np.sqrt(N / M))))

    def oracle(state):
        state[marked] *= -1
        return state

    def diffusion(state):
        mean_amp = np.mean(state)
        return 2.0 * mean_amp - state

    for _ in range(k):
        psi = oracle(psi)
        psi = diffusion(psi)

    probs = np.abs(psi) ** 2
    prob_marked_theory = np.sum(probs[marked])

    print(f"Grover Text-Suche: Verzeichnis={dir_path}, N={N}, ähnlich M={M}, k={k}")
    print(f"Theoretische P(ähnlich) nach k Iterationen: {prob_marked_theory:.4f}")

    top_idx = np.argsort(probs)[-top_display:][::-1]
    samples = rng.choice(N, size=shots, p=probs)
    hits = int(np.sum(marked[samples]))

    print(f"Empirische Trefferquote (shots={shots}): {hits/shots:.4f}")
    print("Top ähnliche Dokumente (Index, Wahrscheinlichkeit, Pfad, Ähnlichkeit):")
    for idx in top_idx:
        path = paths[idx] if idx < len(paths) else None
        path_str = (path or "(Dummy)")[:60] + ("..." if (path and len(path) > 60) else "")
        print(f"  {idx:4d}  {probs[idx]:.6f}  {path_str}  sim={similarities[idx]:.3f}")

    return {
        "probs": probs,
        "paths": paths,
        "marked": marked,
        "similarities": similarities,
        "k": k,
        "M": M,
        "prob_marked_theory": prob_marked_theory,
        "top_indices": top_idx,
    }


if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 2 and sys.argv[1].strip().lower() in ("--dir", "-d", "--text"):
        # Textverzeichnis: python "Grover EABC.py" --dir /pfad/zum/ordner [query.txt]
        dir_path = sys.argv[2] if len(sys.argv) > 2 else "."
        query_path = sys.argv[3] if len(sys.argv) > 3 else None
        grover_text_search(dir_path, query_path=query_path, N=256, similarity_threshold=0.1)
    else:
        grover_bm_search(N=1024, shots=3000, seed=42, oracle_mode="threads_v1")