import argparse
import json
import os
import re
import sys
import time
import subprocess
import shutil
import unicodedata
import urllib.error
import urllib.request
import numpy as np
import sympy as sp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRATCH_ROOT = os.path.join(SCRIPT_DIR, "scratch")
OUT_DIR = os.path.join(SCRIPT_DIR, "out")
DEFAULT_HOST = os.environ.get("OLLAMA_HOST", "127.0.0.1:11434")
DEFAULT_MODEL = os.environ.get("OLLAMA_MODEL", "llama3:latest")
BASE_NUM_CTX = 4096

def validate_latex_math(latex_str):
    """
    Validiert, ob ein mathematischer LaTeX-String korrekt von SymPy geparst werden kann.
    """
    try:
        clean_str = latex_str.replace(r'\frac', '').replace(r'\cdot', '*').replace('^', '**')
        clean_str = clean_str.replace('{', '(').replace('}', ')').replace('\\', '')
        expr = sp.sympify(clean_str)
        return True, str(expr)
    except Exception as e:
        return False, str(e)

class EABCMathInferenceOptimizer:
    """
    Enhanced Artificial Bee Colony (EABC) zur Optimierung von Ollama-Inferenz-Hyperparametern
    für mathematische Textgenerierung und syntaktische Präzision.
    """
    def __init__(self, pop_size=15, max_iter=60, limit=12):
        self.pop_size = pop_size
        self.max_iter = max_iter
        self.limit = limit
        self.dim = 4
        self.bounds = np.array([
            [0.01, 1.00],  # temperature
            [0.10, 1.00],  # top_p
            [1.00, 2.00],  # repeat_penalty
            [0.50, 2.00]   # num_ctx_scale
        ])
        
    def evaluate_fitness(self, params):
        temp, top_p, rep_pen, ctx_scale = params
        test_expressions = [
            r"x**2 + 2*x + 1",
            r"sin(x)**2 + cos(x)**2 - 1",
            r"Matrix([[1, 0], [0, 1]])",
            r"E8_dimension - 248",
            r"Fano_rank - 4"
        ]
        valid_count = 0
        for expr_str in test_expressions:
            is_valid, _ = validate_latex_math(expr_str)
            if is_valid:
                valid_count += 1
        syntax_score = valid_count / len(test_expressions)
        temp_penalty = max(0, temp - 0.2) * 1.5
        rep_bonus = 0.2 if 1.1 <= rep_pen <= 1.3 else 0.0
        fitness = syntax_score - temp_penalty + rep_bonus
        return max(0.0001, fitness)

    def optimize(self):
        start_time = time.time()
        population = np.random.uniform(self.bounds[:, 0], self.bounds[:, 1], size=(self.pop_size, self.dim))
        fitness = np.array([self.evaluate_fitness(ind) for ind in population])
        trial_counts = np.zeros(self.pop_size)
        
        best_idx = np.argmax(fitness)
        best_solution = population[best_idx].copy()
        best_fitness = fitness[best_idx]
        
        for it in range(self.max_iter):
            for i in range(self.pop_size):
                k = np.random.choice([idx for idx in range(self.pop_size) if idx != i])
                phi = np.random.uniform(-1, 1, size=self.dim)
                eta = 0.1
                grad_guided = eta * (best_solution - population[i])
                candidate = population[i] + phi * (population[i] - population[k]) + grad_guided
                candidate = np.clip(candidate, self.bounds[:, 0], self.bounds[:, 1])
                cand_fitness = self.evaluate_fitness(candidate)
                if cand_fitness > fitness[i]:
                    population[i] = candidate
                    fitness[i] = cand_fitness
                    trial_counts[i] = 0
                else:
                    trial_counts[i] += 1
            
            prob = fitness / (np.sum(fitness) + 1e-9)
            for i in range(self.pop_size):
                if np.random.rand() < prob[i]:
                    phi = np.random.uniform(-1, 1, size=self.dim)
                    candidate = population[i] + phi * (population[i] - best_solution)
                    candidate = np.clip(candidate, self.bounds[:, 0], self.bounds[:, 1])
                    cand_fitness = self.evaluate_fitness(candidate)
                    if cand_fitness > fitness[i]:
                        population[i] = candidate
                        fitness[i] = cand_fitness
                        trial_counts[i] = 0
                    else:
                        trial_counts[i] += 1
                        
            for i in range(self.pop_size):
                if trial_counts[i] >= self.limit:
                    population[i] = np.random.uniform(self.bounds[:, 0], self.bounds[:, 1])
                    fitness[i] = self.evaluate_fitness(population[i])
                    trial_counts[i] = 0
                    
            current_best_idx = np.argmax(fitness)
            if fitness[current_best_idx] > best_fitness:
                best_fitness = fitness[current_best_idx]
                best_solution = population[current_best_idx].copy()
                
        elapsed_time = (time.time() - start_time) * 1000.0
        return best_solution, best_fitness, elapsed_time

    def zufaellige_parameter(self):
        return np.random.uniform(self.bounds[:, 0], self.bounds[:, 1])

def normalize_host(host):
    host = host.strip().rstrip("/")
    if not host.startswith("http://") and not host.startswith("https://"):
        host = "http://" + host
    return host

def latex_escape(text):
    mapping = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(mapping.get(char, char) for char in text)

def _im_mathe(ausdruck):
    return r"\(" + ausdruck + r"\)"

# pdflatex kennt diese Zeichen nicht. Sie werden durch vorhandene Matheschriften gesetzt.
_UNICODE_LATEX = {
    "α": _im_mathe(r"\alpha"), "β": _im_mathe(r"\beta"), "γ": _im_mathe(r"\gamma"),
    "δ": _im_mathe(r"\delta"), "ε": _im_mathe(r"\varepsilon"), "ζ": _im_mathe(r"\zeta"),
    "η": _im_mathe(r"\eta"), "θ": _im_mathe(r"\theta"), "ι": _im_mathe(r"\iota"),
    "κ": _im_mathe(r"\kappa"), "λ": _im_mathe(r"\lambda"), "μ": _im_mathe(r"\mu"),
    "ν": _im_mathe(r"\nu"), "ξ": _im_mathe(r"\xi"), "π": _im_mathe(r"\pi"),
    "ρ": _im_mathe(r"\rho"), "σ": _im_mathe(r"\sigma"), "ς": _im_mathe(r"\varsigma"),
    "τ": _im_mathe(r"\tau"), "υ": _im_mathe(r"\upsilon"), "φ": _im_mathe(r"\varphi"),
    "χ": _im_mathe(r"\chi"), "ψ": _im_mathe(r"\psi"), "ω": _im_mathe(r"\omega"),
    "Γ": _im_mathe(r"\Gamma"), "Δ": _im_mathe(r"\Delta"), "Θ": _im_mathe(r"\Theta"),
    "Λ": _im_mathe(r"\Lambda"), "Ξ": _im_mathe(r"\Xi"), "Π": _im_mathe(r"\Pi"),
    "Σ": _im_mathe(r"\Sigma"), "Υ": _im_mathe(r"\Upsilon"), "Φ": _im_mathe(r"\Phi"),
    "Ψ": _im_mathe(r"\Psi"), "Ω": _im_mathe(r"\Omega"),
    "∂": _im_mathe(r"\partial"), "∇": _im_mathe(r"\nabla"), "⋅": _im_mathe(r"\cdot"),
    "×": _im_mathe(r"\times"), "÷": _im_mathe(r"\div"), "±": _im_mathe(r"\pm"),
    "∓": _im_mathe(r"\mp"), "≤": _im_mathe(r"\leq"), "≥": _im_mathe(r"\geq"),
    "≠": _im_mathe(r"\neq"), "≈": _im_mathe(r"\approx"), "≡": _im_mathe(r"\equiv"),
    "∞": _im_mathe(r"\infty"), "∈": _im_mathe(r"\in"), "∉": _im_mathe(r"\notin"),
    "⊂": _im_mathe(r"\subset"), "⊃": _im_mathe(r"\supset"), "∪": _im_mathe(r"\cup"),
    "∩": _im_mathe(r"\cap"), "∀": _im_mathe(r"\forall"), "∃": _im_mathe(r"\exists"),
    "∅": _im_mathe(r"\emptyset"), "√": _im_mathe(r"\sqrt{}"), "∑": _im_mathe(r"\sum"),
    "∏": _im_mathe(r"\prod"), "∫": _im_mathe(r"\int"), "¬": _im_mathe(r"\neg"),
    "→": _im_mathe(r"\rightarrow"), "←": _im_mathe(r"\leftarrow"),
    "↔": _im_mathe(r"\leftrightarrow"), "⇒": _im_mathe(r"\Rightarrow"),
    "⇐": _im_mathe(r"\Leftarrow"), "↦": _im_mathe(r"\mapsto"),
    "−": _im_mathe("-"), "·": _im_mathe(r"\cdot"),
    "⁰": _im_mathe(r"^{0}"), "¹": _im_mathe(r"^{1}"), "²": _im_mathe(r"^{2}"),
    "³": _im_mathe(r"^{3}"), "⁴": _im_mathe(r"^{4}"), "⁵": _im_mathe(r"^{5}"),
    "⁶": _im_mathe(r"^{6}"), "⁷": _im_mathe(r"^{7}"), "⁸": _im_mathe(r"^{8}"),
    "⁹": _im_mathe(r"^{9}"),
    "₀": _im_mathe(r"_{0}"), "₁": _im_mathe(r"_{1}"), "₂": _im_mathe(r"_{2}"),
    "₃": _im_mathe(r"_{3}"), "₄": _im_mathe(r"_{4}"), "₅": _im_mathe(r"_{5}"),
    "₆": _im_mathe(r"_{6}"), "₇": _im_mathe(r"_{7}"), "₈": _im_mathe(r"_{8}"),
    "₉": _im_mathe(r"_{9}"),
    "–": r"\textendash{}", "—": r"\textemdash{}", "…": r"\ldots{}",
    "„": r"\glqq{}", "“": r"\grqq{}", "”": r"\grqq{}",
    "‚": r"\glq{}", "‘": r"\grq{}", "’": "'",
    "«": r"\guillemotleft{}", "»": r"\guillemotright{}",
    "\u2011": "-",
}

def unicode_nach_latex(text):
    teile = []
    for char in text:
        if char in _UNICODE_LATEX:
            teile.append(_UNICODE_LATEX[char])
        elif char == "\u00ad":
            teile.append("")
        elif ord(char) < 256:
            teile.append(char)
        else:
            ersatz = unicodedata.normalize("NFKD", char).encode("ascii", "ignore").decode("ascii")
            if ersatz:
                teile.append(ersatz)
            else:
                name = unicodedata.name(char, "zeichen")
                sicher = " ".join(
                    "".join(c.lower() if c.isalnum() else " " for c in name).split()
                )
                teile.append(r"\textit{[" + sicher + "]}")
    return re.sub(r"\\\)\\\(", "", "".join(teile))

def markdown_fett(text):
    return re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)

def _tabellenzellen(zeile):
    return [zelle.strip() for zelle in zeile.strip().strip("|").split("|")]

def ist_trennzeile(zeile):
    zellen = _tabellenzellen(zeile)
    return bool(zellen) and all(re.fullmatch(r":?-{3,}:?", zelle) for zelle in zellen)

def ist_markdown_tabelle(zeilen):
    kern = [zeile.strip() for zeile in zeilen if zeile.strip()]
    if len(kern) < 3 or not any(ist_trennzeile(zeile) for zeile in kern):
        return False
    return all(zeile.startswith("|") and zeile.endswith("|") for zeile in kern)

def markdown_tabelle_zu_latex(zeilen):
    kern = [zeile.strip() for zeile in zeilen if zeile.strip() and not ist_trennzeile(zeile)]
    daten = [_tabellenzellen(zeile) for zeile in kern]
    if not daten:
        return ""
    spalten = max(len(zeile) for zeile in daten)

    def zelle(text):
        return markdown_fett(unicode_nach_latex(latex_escape(text)))

    def reihe(zellen, kopf=False):
        zellen = zellen + [""] * (spalten - len(zellen))
        inhalte = [zelle(eintrag) for eintrag in zellen]
        if kopf:
            inhalte = [r"\textbf{" + eintrag + "}" for eintrag in inhalte]
        return " & ".join(inhalte) + r" \\"

    muster = "p{8.2cm}" + "l" * (spalten - 1) if spalten > 1 else "l"
    koerper = reihe(daten[0], kopf=True) + "\n\\midrule\n" + "\n".join(reihe(zeile) for zeile in daten[1:])
    return (
        "\\begin{tabular}{@{}" + muster + "@{}}\n\\toprule\n"
        + koerper
        + "\n\\bottomrule\n\\end{tabular}"
    )

def latex_block(text):
    cleaned = (text or "").replace("\r\n", "\n").strip()
    if not cleaned:
        return r"\textit{(leer)}"
    paragraphs = []
    for paragraph in cleaned.split("\n\n"):
        zeilen = paragraph.split("\n")
        if ist_markdown_tabelle(zeilen):
            tabelle = markdown_tabelle_zu_latex(zeilen)
            if tabelle:
                paragraphs.append(tabelle)
            continue
        lines = [
            markdown_fett(unicode_nach_latex(latex_escape(line).strip()))
            for line in zeilen
        ]
        lines = [line for line in lines if line]
        if lines:
            paragraphs.append("\\\\\n".join(lines))
    return "\n\n".join(paragraphs) if paragraphs else r"\textit{(leer)}"

def de_zahl(wert, stellen=2):
    return f"{wert:.{stellen}f}".replace(".", ",")

def ns_zu_ms(nanosekunden):
    return float(nanosekunden or 0) / 1_000_000.0

def ollama_erreichbar(host):
    try:
        with urllib.request.urlopen(host + "/api/tags", timeout=2) as antwort:
            antwort.read()
        return True
    except Exception:
        return False

def ollama_server_sicherstellen(host):
    if ollama_erreichbar(host):
        return
    if shutil.which("ollama") is None:
        sys.exit("Ollama ist nicht installiert.")
    print("Ollama-Server ist nicht aktiv. Starte ihn ...")
    subprocess.Popen(
        ["ollama", "serve"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    for _ in range(60):
        if ollama_erreichbar(host):
            return
        time.sleep(0.5)
    sys.exit(f"Ollama antwortet nicht auf {host}.")

def installierte_modelle(host):
    with urllib.request.urlopen(host + "/api/tags", timeout=10) as antwort:
        daten = json.load(antwort)
    return [eintrag.get("name", "") for eintrag in daten.get("models", [])]

def modell_vorhanden(namen, modell):
    if modell in namen:
        return True
    kurz = modell.split(":")[0]
    return any(name == kurz or name.split(":")[0] == kurz for name in namen)

def ollama_fragen(host, modell, frage, params):
    temp, top_p, rep_pen, ctx_scale = params
    num_ctx = max(512, int(round(BASE_NUM_CTX * float(ctx_scale))))
    nutzlast = {
        "model": modell,
        "messages": [{"role": "user", "content": frage}],
        "stream": False,
        "options": {
            "temperature": float(temp),
            "top_p": float(top_p),
            "repeat_penalty": float(rep_pen),
            "num_ctx": num_ctx,
        },
    }
    # gpt-oss schreibt sonst eine Denkspur, die Ollama als Tool-Aufruf liest und mit HTTP 500 abbricht.
    if "gpt-oss" in modell.lower() or "gptoss" in modell.lower():
        nutzlast["think"] = False
    start = time.perf_counter()
    koerper, fehlerdetail = ollama_post(host, "/api/chat", nutzlast)
    if fehlerdetail and "parsing tool call" in fehlerdetail:
        nutzlast["think"] = False
        koerper, fehlerdetail = ollama_post(host, "/api/chat", nutzlast)
    if fehlerdetail:
        sys.exit(f"Ollama-Fehler: {fehlerdetail}")
    wand_ms = (time.perf_counter() - start) * 1000.0
    eval_ms = ns_zu_ms(koerper.get("eval_duration"))
    eval_tokens = int(koerper.get("eval_count") or 0)
    prompt_tokens = int(koerper.get("prompt_eval_count") or 0)
    prompt_ms = ns_zu_ms(koerper.get("prompt_eval_duration"))
    lade_ms = ns_zu_ms(koerper.get("load_duration"))
    gesamt_ms = ns_zu_ms(koerper.get("total_duration"))
    nachricht = koerper.get("message") or {}
    antwort_text = koerper.get("response") or nachricht.get("content") or ""
    laufzeit, modellinfo = modell_laufzeitinfos(host, modell)
    info = modellinfo.get("model_info") or {}
    nativer_kontext = next(
        (wert for schluessel, wert in info.items() if schluessel.endswith("context_length")),
        None,
    )
    return {
        "frage": frage,
        "antwort": antwort_text,
        "modell": modell,
        "num_ctx": num_ctx,
        "temperature": float(temp),
        "top_p": float(top_p),
        "repeat_penalty": float(rep_pen),
        "done_reason": koerper.get("done_reason") or "",
        "wand_ms": wand_ms,
        "lade_ms": lade_ms,
        "prompt_tokens": prompt_tokens,
        "prompt_ms": prompt_ms,
        "prompt_tokens_pro_s": (prompt_tokens / (prompt_ms / 1000.0)) if prompt_ms > 0 else 0.0,
        "eval_tokens": eval_tokens,
        "eval_ms": eval_ms,
        "tokens_pro_s": (eval_tokens / (eval_ms / 1000.0)) if eval_ms > 0 else 0.0,
        "ms_pro_token": (eval_ms / eval_tokens) if eval_tokens else 0.0,
        "gesamt_ms": gesamt_ms,
        "ohne_laden_ms": max(0.0, gesamt_ms - lade_ms),
        "bis_erstes_token_ms": lade_ms + prompt_ms,
        "zeichen": len(antwort_text.strip()),
        "zeichen_pro_s": (len(antwort_text.strip()) / (eval_ms / 1000.0)) if eval_ms > 0 else 0.0,
        "vram_bytes": int(laufzeit.get("size_vram") or 0),
        "modell_bytes": int(laufzeit.get("size") or 0),
        "belegter_kontext": laufzeit.get("context_length"),
        "parameter_groesse": (laufzeit.get("details") or modellinfo.get("details") or {}).get("parameter_size") or "",
        "quantisierung": (laufzeit.get("details") or modellinfo.get("details") or {}).get("quantization_level") or "",
        "familie": (laufzeit.get("details") or modellinfo.get("details") or {}).get("family") or "",
        "nativer_kontext": nativer_kontext,
        "parameter_anzahl": info.get("general.parameter_count"),
    }

def ollama_post(host, pfad, nutzlast):
    anfrage = urllib.request.Request(
        host + pfad,
        data=json.dumps(nutzlast).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(anfrage, timeout=600) as antwort:
            return json.load(antwort), ""
    except urllib.error.HTTPError as fehler:
        return None, fehler.read().decode("utf-8", errors="replace")

def ollama_json(host, pfad, nutzlast=None):
    daten = None if nutzlast is None else json.dumps(nutzlast).encode("utf-8")
    kopf = {"Content-Type": "application/json"} if daten is not None else {}
    anfrage = urllib.request.Request(host + pfad, data=daten, headers=kopf)
    with urllib.request.urlopen(anfrage, timeout=20) as antwort:
        return json.load(antwort)

def modell_laufzeitinfos(host, modell):
    try:
        prozesse = ollama_json(host, "/api/ps")
    except Exception:
        prozesse = {}
    laufzeit = {}
    kurz = modell.split(":")[0]
    for eintrag in prozesse.get("models", []):
        name = eintrag.get("name") or eintrag.get("model") or ""
        if name == modell or name.split(":")[0] == kurz:
            laufzeit = eintrag
            break
    try:
        info = ollama_json(host, "/api/show", {"model": modell})
    except Exception:
        info = {}
    return laufzeit, info

def benchmark_ausgeben(ergebnis):
    print("")
    print(f"---------------- {ergebnis['verfahren']} ----------------")
    print(
        "Parameter: "
        f"temp={ergebnis['temperature']:.4f}, "
        f"top_p={ergebnis['top_p']:.4f}, "
        f"rep_pen={ergebnis['repeat_penalty']:.4f}, "
        f"ctx={ergebnis['num_ctx']}"
    )
    print("---------------- Leistung --------------------")
    print(f"Modell:                 {ergebnis['modell']} ({ergebnis['parameter_groesse']}, {ergebnis['quantisierung']})")
    print(f"Grafikspeicher:         {ergebnis['vram_bytes'] / (1024 ** 3):.2f} GiB")
    print(f"Kontext angefordert:    {ergebnis['num_ctx']} Tokens")
    print(f"Kontext belegt:         {ergebnis['belegter_kontext']}")
    print(f"Wandzeit:               {ergebnis['wand_ms']:.1f} ms")
    print(f"davon Ladezeit:         {ergebnis['lade_ms']:.1f} ms")
    print(f"ohne Laden:             {ergebnis['ohne_laden_ms']:.1f} ms")
    print(f"bis zum ersten Token:   {ergebnis['bis_erstes_token_ms']:.1f} ms")
    print(f"Prompt:                 {ergebnis['prompt_tokens']} Tokens in {ergebnis['prompt_ms']:.1f} ms ({ergebnis['prompt_tokens_pro_s']:.2f} Tokens/s)")
    print(f"Antwort:                {ergebnis['eval_tokens']} Tokens in {ergebnis['eval_ms']:.1f} ms ({ergebnis['tokens_pro_s']:.2f} Tokens/s)")
    print(f"Zeit pro Antwort-Token: {ergebnis['ms_pro_token']:.2f} ms")
    print(f"Zeichen:                {ergebnis['zeichen']} ({ergebnis['zeichen_pro_s']:.1f} Zeichen/s)")
    print(f"Abschluss:              {ergebnis['done_reason'] or 'unbekannt'}")
    print("---------------- Antwort ---------------------")
    print(ergebnis["antwort"].strip())
    print("----------------------------------------------")

def _differenz(links, rechts, stellen):
    if links is None or rechts is None:
        return "nicht gemeldet"
    delta = float(rechts) - float(links)
    vorzeichen = "+" if delta > 0 else ""
    return vorzeichen + de_zahl(delta, stellen)

def vergleichszeilen(erste, zweite):
    paare = [
        ("Temperature", "temperature", 4, ""),
        ("Top-p", "top_p", 4, ""),
        ("Repeat-Penalty", "repeat_penalty", 4, ""),
        ("Kontext angefordert", "num_ctx", 0, " Tokens"),
        ("Kontext belegt", "belegter_kontext", 0, " Tokens"),
        ("Grafikspeicher", "vram_bytes", None, " GiB"),
        ("Wandzeit", "wand_ms", 1, " ms"),
        ("Ladezeit", "lade_ms", 1, " ms"),
        ("Zeit ohne Laden", "ohne_laden_ms", 1, " ms"),
        ("Bis zum ersten Token", "bis_erstes_token_ms", 1, " ms"),
        ("Prompt-Auswertung", "prompt_ms", 1, " ms"),
        ("Prompt-Durchsatz", "prompt_tokens_pro_s", 2, " Tokens/s"),
        ("Antwort-Tokens", "eval_tokens", 0, ""),
        ("Generierung", "eval_ms", 1, " ms"),
        ("Antwort-Durchsatz", "tokens_pro_s", 2, " Tokens/s"),
        ("Zeit pro Antwort-Token", "ms_pro_token", 2, " ms"),
        ("Antwortzeichen", "zeichen", 0, ""),
        ("Zeichendurchsatz", "zeichen_pro_s", 1, " Zeichen/s"),
        ("Gesamtzeit laut Ollama", "gesamt_ms", 1, " ms"),
    ]
    zeilen = []
    for name, schluessel, stellen, einheit in paare:
        links = erste.get(schluessel)
        rechts = zweite.get(schluessel)
        if schluessel == "vram_bytes":
            zeilen.append((
                name,
                gibibyte(links),
                gibibyte(rechts),
                _differenz(None if not links else links / (1024 ** 3), None if not rechts else rechts / (1024 ** 3), 2) + " GiB",
            ))
        elif schluessel in ("num_ctx", "belegter_kontext", "eval_tokens", "zeichen"):
            if links is None or rechts is None:
                zeilen.append((name, "nicht gemeldet", "nicht gemeldet", "nicht gemeldet"))
            else:
                zeilen.append((
                    name,
                    str(int(links)) + einheit,
                    str(int(rechts)) + einheit,
                    _differenz(links, rechts, 0) + einheit,
                ))
        else:
            zeilen.append((
                name,
                de_zahl(links, stellen) + einheit,
                de_zahl(rechts, stellen) + einheit,
                _differenz(links, rechts, stellen) + einheit,
            ))
    return zeilen

def vergleich_ausgeben(erste, zweite):
    print("")
    print("---------------- Vergleich ------------------")
    print(f"{'Kennzahl':<24} {'Zufällig':>18} {'Voroptimiert':>18} {'Differenz':>18}")
    for name, links, rechts, delta in vergleichszeilen(erste, zweite):
        print(f"{name:<24} {links:>18} {rechts:>18} {delta:>18}")
    print("Differenz = Voroptimiert minus Zufällig.")
    print("Die Ladezeit des zweiten Durchlaufs ist kürzer, weil das Modell schon im Speicher liegt.")
    print("----------------------------------------------")

def leistungstabelle(zeilen):
    koerper = "\n".join(name + " & " + wert + r" \\" for name, wert in zeilen)
    return r"""
{\small
\begin{tabularx}{\textwidth}{Xr}
\toprule
\textbf{Kennzahl} & \textbf{Wert} \\
\midrule
""" + koerper + r"""
\bottomrule
\end{tabularx}
}"""

def gibibyte(anzahl):
    if not anzahl:
        return "nicht gemeldet"
    return de_zahl(anzahl / (1024 ** 3), 2) + " GiB"

def durchlauf_als_latex(ergebnis, titel):
    parameter = leistungstabelle([
        ("Temperature", de_zahl(ergebnis["temperature"], 4)),
        ("Top-p", de_zahl(ergebnis["top_p"], 4)),
        ("Repeat-Penalty", de_zahl(ergebnis["repeat_penalty"], 4)),
        ("Kontext angefordert", str(ergebnis["num_ctx"]) + " Tokens"),
    ])
    return (
        r"\subsection*{" + titel + "}\n"
        + parameter
        + r"""
\paragraph{Antwort}
\begin{quote}
\raggedright
"""
        + latex_block(ergebnis["antwort"])
        + r"""
\end{quote}
"""
    )

def vergleich_als_latex(erste, zweite):
    koerper = "\n".join(
        latex_escape(name) + " & " + links + " & " + rechts + " & " + delta + r" \\"
        for name, links, rechts, delta in vergleichszeilen(erste, zweite)
    )
    return r"""
\subsection*{Vergleich}
{\footnotesize
\setlength{\tabcolsep}{5pt}
\begin{tabular}{@{}lrrr@{}}
\toprule
\textbf{Kennzahl} & \textbf{Zufällig} & \textbf{Voroptimiert} & \textbf{Dif{}ferenz} \\
\midrule
""" + koerper + r"""
\bottomrule
\end{tabular}
}

\medskip
\noindent Die Dif{}ferenz ist Voroptimiert minus Zufällig. Die Ladezeit des zweiten Durchlaufs ist kürzer, weil das Modell bereits im Speicher liegt. Vergleichbar sind Zeit ohne Laden, Generierung und Durchsatz.
"""

def anfrage_als_latex(erste, zweite):
    return r"""
\section*{Frage}

\begin{quote}
\raggedright
""" + latex_block(erste["frage"]) + r"""
\end{quote}

""" + durchlauf_als_latex(erste, "Durchlauf 1: zufällige Parameter") + durchlauf_als_latex(zweite, "Durchlauf 2: voroptimierte Parameter") + vergleich_als_latex(erste, zweite)

def frage_aus_argumenten(args):
    if args.frage_text:
        return args.frage_text.strip()
    if args.frage:
        return " ".join(args.frage).strip()
    if sys.stdin.isatty():
        return input("Frage an Ollama: ").strip()
    return ""

def generate_pdf_report(output_pdf_path, erste, zweite):
    scratch_dir = os.path.join(SCRATCH_ROOT, "ollama_pdf_pipeline_build")
    os.makedirs(scratch_dir, exist_ok=True)
    tex_file = os.path.join(scratch_dir, "document.tex")
    latex_body = anfrage_als_latex(erste, zweite)

    full_tex = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=2.5cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[utf8]{inputenc}
\usepackage[ngerman]{babel}
\usepackage{amsmath,amssymb,amsfonts,amsthm}
\usepackage{xcolor}
\definecolor{gray}{gray}{0.45}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{hyperref}
\usepackage{fancyhdr}

\definecolor{navy}{HTML}{1A365D}
\definecolor{accent}{HTML}{2B6CB0}

\hypersetup{
    colorlinks=true,
    linkcolor=accent,
    urlcolor=accent,
    pdftitle={Protokoll der letzten Frage}
}

\pagestyle{fancy}
\fancyhf{}
\rhead{\small \textcolor{gray}{""" + latex_escape(erste["modell"]) + r"""}}
\lhead{\small \textcolor{gray}{Letzte Frage}}
\cfoot{\small \thepage}

\title{\Huge \textbf{\textcolor{navy}{Protokoll der letzten Frage}}\\[0.3em]
\large \textcolor{accent}{""" + latex_escape(erste["modell"]) + r"""}}
\author{}
\date{\today}

\begin{document}

\maketitle
\thispagestyle{empty}
\thispagestyle{fancy}

\vspace{0.3cm}
\hrule height 1pt
\vspace{0.8cm}

""" + latex_body + r"""
\end{document}
"""

    with open(tex_file, "w", encoding="utf-8") as f:
        f.write(full_tex)
        
    os.environ["TEXMFVAR"] = "/tmp/texmf-var"
    cmd = ['pdflatex', '-interaction=nonstopmode', '-output-directory', scratch_dir, tex_file]
    
    subprocess.run(cmd, check=False, capture_output=True)
    res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    
    gen_pdf = os.path.join(scratch_dir, "document.pdf")
    log_path = os.path.join(scratch_dir, "document.log")
    latex_errors = []
    if os.path.exists(log_path):
        with open(log_path, encoding="utf-8", errors="replace") as log_file:
            latex_errors = [line for line in log_file if line.startswith("!")]
    if os.path.exists(gen_pdf) and not latex_errors:
        os.makedirs(os.path.dirname(output_pdf_path), exist_ok=True)
        shutil.copy(gen_pdf, output_pdf_path)
        print(f"PDF-Bericht erfolgreich erzeugt: {output_pdf_path}")
        return True
    else:
        print("Fehler bei PDF-Erstellung:")
        if latex_errors:
            print("".join(latex_errors))
        print(res.stdout[-600:] if res.stdout and len(res.stdout) > 600 else (res.stdout or ""))
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Stellt eine Frage an Ollama und misst den Aufruf."
    )
    parser.add_argument("frage", nargs="*", help="Frage an das Modell")
    parser.add_argument("-q", "--frage-text", help="Frage als eine Option")
    parser.add_argument("-m", "--model", default=DEFAULT_MODEL, help="Ollama-Modellname")
    parser.add_argument("--host", default=DEFAULT_HOST, help="Ollama-Adresse, z. B. 127.0.0.1:11434")
    args = parser.parse_args()

    frage = frage_aus_argumenten(args)
    if not frage:
        sys.exit("Keine Frage angegeben.")

    optimizer = EABCMathInferenceOptimizer(pop_size=15, max_iter=60, limit=12)
    host = normalize_host(args.host)
    ollama_server_sicherstellen(host)
    modelle = installierte_modelle(host)
    if not modell_vorhanden(modelle, args.model):
        verfuegbar = ", ".join(modelle) if modelle else "(keine)"
        sys.exit(f"Modell '{args.model}' ist nicht installiert. Vorhanden: {verfuegbar}")

    print(f"Durchlauf 1 von 2: zufällige Parameter, Frage an {args.model} ...")
    erste = ollama_fragen(host, args.model, frage, optimizer.zufaellige_parameter())
    erste["verfahren"] = "Durchlauf 1: zufällige Parameter"
    benchmark_ausgeben(erste)

    print(f"Durchlauf 2 von 2: Voroptimierung, dieselbe Frage an {args.model} ...")
    voroptimiert, _, _ = optimizer.optimize()
    zweite = ollama_fragen(host, args.model, frage, voroptimiert)
    zweite["verfahren"] = "Durchlauf 2: voroptimierte Parameter"
    benchmark_ausgeben(zweite)
    vergleich_ausgeben(erste, zweite)

    pdf_out = os.path.join(SCRATCH_ROOT, "ollama_eabc_benchmark_report.pdf")
    success = generate_pdf_report(pdf_out, erste, zweite)

    if success:
        os.makedirs(OUT_DIR, exist_ok=True)
        ziel = os.path.join(OUT_DIR, "ollama_eabc_benchmark_report.pdf")
        shutil.copy(pdf_out, ziel)
        shutil.copy(__file__, os.path.join(OUT_DIR, "eabc_ollama_pipeline_benchmark.py"))
        stempel = time.strftime("%Y%m%d-%H%M%S")
        shutil.copy(pdf_out, os.path.join(OUT_DIR, f"ollama_benchmark_{stempel}.pdf"))
        print(f"\nPipeline erfolgreich ausgeführt und Artefakte nach {OUT_DIR} publiziert!")
    else:
        sys.exit(1)
