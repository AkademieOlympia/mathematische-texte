"""
EABC-Topos: Ollama-Anbindung mit Kaluza-Klein-Sieb
=================================================
Lokale Ollama-Modelle, klassifiziert durch das Sieb Omega.
Textdateien des Verzeichnisses werden in Tokenabschnitte zerlegt.
Ein Grover-Orakel verstaerkt die passenden Abschnitte, bevor sie
als Kontext an das Modell gehen.

1_Omega gilt nur zugleich fuer drei Bedingungen: der Dienst hat eine
nichtleere Antwort geliefert, det(G) > 0, und die Galois-Norm des
Zeugen in Z[phi] ist exakt -1. Ein nicht erreichbarer Dienst bleibt
0_Omega; es wird kein Antworttext erzeugt.
"""

from __future__ import annotations

import json
import math
import os
import re
import urllib.error
import urllib.request
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional


@dataclass(frozen=True)
class Icosian:
    """Element a + b*phi von Z[phi], phi = (1 + sqrt(5)) / 2.

    Die Galois-Konjugation ist phi |-> 1 - phi. Die Norm
    N(a + b*phi) = a^2 + a*b - b^2 liegt in Z. Fuer phi = 0 + 1*phi
    gilt N(phi) = -1, also -N(phi) = 1.
    """

    a: int
    b: int

    def __post_init__(self) -> None:
        if isinstance(self.a, bool) or isinstance(self.b, bool):
            raise TypeError("Koeffizienten von Z[phi] sind ganzzahlig, keine Wahrheitswerte.")
        if not isinstance(self.a, int) or not isinstance(self.b, int):
            raise TypeError("Koeffizienten von Z[phi] sind ganzzahlig.")

    def conjugate(self) -> "Icosian":
        return Icosian(self.a + self.b, -self.b)

    def __mul__(self, other: "Icosian") -> "Icosian":
        # phi^2 = phi + 1
        return Icosian(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def norm(self) -> int:
        product = self * self.conjugate()
        if product.b != 0:
            raise RuntimeError("Norm verlaesst Z[phi] nicht nach Z.")
        return product.a


PHI = Icosian(0, 1)


@dataclass(frozen=True)
class SieveVerdict:
    label: str
    principal: bool
    reason: str

    def as_text(self) -> str:
        return f"{self.label} ({self.reason})"


def classify(norm: int, det_g: float) -> SieveVerdict:
    """Unterobjekt-Klassifikator Omega.

    det(G) <= 0 ist das leere Sieb. Andernfalls ist das Hauptsieb genau
    die Normschale N = -1.
    """
    if det_g <= 0.0:
        return SieveVerdict("0_Omega", False, "Singularitaet det(G) <= 0")
    if norm != -1:
        return SieveVerdict("0_Omega", False, "Norm-Verletzung N != -1")
    return SieveVerdict("1_Omega", True, "Hauptsieb")


class ToposGuardrail:
    """Klassifiziert einen geometrischen Zustand, nicht den Wortlaut."""

    def __init__(self, witness: Icosian = PHI) -> None:
        self.witness = witness
        if witness == PHI and witness.norm() != -1:
            raise RuntimeError("Die Identitaet -N(phi) = 1 ist gebrochen.")

    def evaluate(
        self,
        det_g: float,
        witness: Optional[Icosian] = None,
    ) -> SieveVerdict:
        element = self.witness if witness is None else witness
        return classify(element.norm(), det_g)


_TOKEN = re.compile(r"[0-9A-Za-zÄÖÜäöüß]+")
_CAMEL = re.compile(r"[A-ZÄÖÜ]+(?=[A-ZÄÖÜ][a-zäöüß])|[A-ZÄÖÜ]?[a-zäöüß]+|[A-ZÄÖÜ]+|\d+")
_STOP = frozenset({
    "aber", "als", "also", "and", "an", "am", "auf", "aus", "auch", "bei",
    "bin", "bis", "das", "dass", "dem", "den", "der", "des", "die", "dies",
    "diese", "dieser", "dieses", "ein", "eine", "einem", "einen", "einer",
    "eines", "er", "es", "for", "from", "fuer", "für", "hat", "im", "in",
    "ist", "man", "mit", "nach", "nicht", "nur", "oder", "of", "on", "sich",
    "sie", "sind", "so", "that", "the", "to", "und", "von", "vor", "was",
    "wie", "wir", "with", "zu", "zum", "zur", "über",
})
_SKIP_DIRS = frozenset({
    ".git", ".venv", ".claude", ".agents", ".cursor",
    "node_modules", "__pycache__", "scratch", "out", "external",
})
_TEXT_EXTENSIONS = frozenset({".tex", ".txt", ".md"})


def content_tokens(text: str) -> list[str]:
    """Wort- und Zahlentoken ohne Funktionswoerter.

    Binnenmajuskeln werden getrennt, damit ein Dateiname wie
    RiemannNullstellen die Token riemann und nullstellen ergibt.
    """
    tokens = []
    for match in _TOKEN.finditer(text):
        parts = _CAMEL.findall(match.group(0)) or [match.group(0)]
        for part in parts:
            token = (
                part.casefold()
                .replace("ä", "ae")
                .replace("ö", "oe")
                .replace("ü", "ue")
                .replace("ß", "ss")
            )
            if len(token) <= 1 or token in _STOP:
                continue
            if token.isdigit() and len(token) > 4:
                continue
            tokens.append(token)
    return tokens


@dataclass(frozen=True)
class TextSpan:
    source: str
    text: str
    tokens: tuple[str, ...]
    name_tokens: tuple[str, ...] = ()


@dataclass(frozen=True)
class GroverRetrieval:
    spans: list[TextSpan]
    iterations: int
    target_count: int
    span_count: int
    fallback: bool

    def as_dict(self) -> dict[str, Any]:
        return {
            "iterations": self.iterations,
            "target_count": self.target_count,
            "span_count": self.span_count,
            "fallback": self.fallback,
        }


def grover_amplify(
    scores: list[float],
    threshold: float,
) -> tuple[list[float], int, list[int]]:
    """Grover-Iteration ueber Tokenabschnitte.

    Das Orakel kippt das Vorzeichen der markierten Amplituden.
    Die Diffusion spiegelt an der mittleren Amplitude. Die optimale
    Schrittzahl ist (pi/4) * sqrt(N/M).
    """
    count = len(scores)
    if count == 0:
        return [], 0, []
    targets = [index for index, score in enumerate(scores) if score >= threshold]
    if not targets:
        targets = [max(range(count), key=lambda index: scores[index])]
    marked = len(targets)
    iterations = max(1, int(round((math.pi / 4.0) * math.sqrt(count / marked))))
    amplitude = 1.0 / math.sqrt(count)
    state = [amplitude] * count
    for _ in range(iterations):
        for index in targets:
            state[index] = -state[index]
        mean = sum(state) / count
        twice_mean = 2.0 * mean
        for index in range(count):
            state[index] = twice_mean - state[index]
    return [value * value for value in state], iterations, targets


class DirectoryTokenCorpus:
    """Tokenindex ueber die Textdateien eines Verzeichnisses."""

    def __init__(
        self,
        root: Path,
        chunk_words: int = 280,
        overlap_words: int = 40,
    ) -> None:
        self.root = root
        self.chunk_words = chunk_words
        self.overlap_words = overlap_words
        self.spans: list[TextSpan] = []
        self.file_count = 0
        self.token_count = 0
        self._idf: dict[str, float] = {}

    def build(self) -> "DirectoryTokenCorpus":
        documents = 0
        for path in self._files():
            text = self._read(path)
            if not text.strip():
                continue
            documents += 1
            source = path.relative_to(self.root).as_posix()
            name_tokens = tuple(dict.fromkeys(content_tokens(path.stem)))
            first = True
            for window in self._windows(text):
                body = content_tokens(window)
                self.token_count += len(body)
                tokens = tuple(dict.fromkeys(body))
                if not tokens and not (first and name_tokens):
                    first = False
                    continue
                self.spans.append(TextSpan(
                    source,
                    window,
                    tokens,
                    name_tokens if first else (),
                ))
                first = False
        self.file_count = documents
        df: Counter[str] = Counter()
        for span in self.spans:
            df.update(set(span.tokens).union(span.name_tokens))
        total = max(1, len(self.spans))
        self._idf = {
            token: math.log((total + 1) / (count + 1)) + 1.0
            for token, count in df.items()
        }
        return self

    def retrieve(
        self,
        question: str,
        limit: int = 4,
        threshold: float = 0.30,
    ) -> GroverRetrieval:
        query = content_tokens(question)
        empty = GroverRetrieval([], 0, 0, len(self.spans), False)
        if not query or not self.spans:
            return empty
        weights = Counter(query)
        scores = [self._score_span(span, weights) for span in self.spans]
        peak = max(scores)
        if peak <= 0.0:
            return empty
        normalized = [score / peak for score in scores]
        raw_targets = [index for index, score in enumerate(normalized) if score >= threshold]
        probabilities, iterations, targets = grover_amplify(normalized, threshold)
        best: dict[str, tuple[float, float, int]] = {}
        for index in targets:
            span = self.spans[index]
            score = normalized[index]
            probability = probabilities[index]
            current = best.get(span.source)
            if current is None or (score, probability) > (current[0], current[1]):
                best[span.source] = (score, probability, index)
        ranked = sorted(best.values(), key=lambda item: (-item[0], -item[1], item[2]))
        chosen = [self.spans[index] for _, _, index in ranked[:limit]]
        return GroverRetrieval(
            chosen,
            iterations,
            len(targets),
            len(self.spans),
            not raw_targets,
        )

    def _score_span(self, span: TextSpan, weights: Counter[str]) -> float:
        body = set(span.tokens)
        names = set(span.name_tokens)
        score = 0.0
        for token, count in weights.items():
            weight = self._idf.get(token, 1.0) * (1.0 + math.log(count))
            if token in body:
                score += weight
                if token in names:
                    score += 0.15 * weight
            elif token in names:
                score += 0.35 * weight
        return score

    def as_prompt(self, spans: list[TextSpan], question: str) -> str:
        if not spans:
            return question
        blocks = [
            f"[Quelle {number}: {span.source}]\n{span.text}"
            for number, span in enumerate(spans, start=1)
        ]
        context = "\n\n".join(blocks)
        return (
            "Die folgenden Tokenabschnitte wurden vom Grover-Orakel verstaerkt. "
            "Stuetze die Antwort darauf, nenne die Quelldatei, und trenne das von Modellannahmen. "
            "Wenn die Abschnitte die Frage nicht tragen, sage das.\n\n"
            f"{context}\n\nFrage: {question}"
        )

    def _files(self):
        for dirpath, dirnames, filenames in os.walk(self.root):
            dirnames[:] = [
                name for name in dirnames
                if name not in _SKIP_DIRS and not name.startswith(".")
            ]
            for name in filenames:
                if name.startswith("."):
                    continue
                path = Path(dirpath) / name
                if path.suffix.lower() in _TEXT_EXTENSIONS:
                    yield path

    def _read(self, path: Path) -> str:
        raw = path.read_bytes()
        if b"\x00" in raw[:4096]:
            return ""
        text = raw.decode("utf-8", errors="replace")
        if path.suffix.lower() != ".tex":
            return text
        kept = []
        for line in text.splitlines():
            if line.lstrip().startswith("%"):
                continue
            kept.append(line)
        return "\n".join(kept)

    def _windows(self, text: str) -> list[str]:
        words = text.split()
        if not words:
            return []
        if len(words) <= self.chunk_words:
            return [" ".join(words)]
        step = max(1, self.chunk_words - self.overlap_words)
        windows = []
        for start in range(0, len(words), step):
            piece = words[start:start + self.chunk_words]
            if len(piece) < self.overlap_words and windows:
                break
            windows.append(" ".join(piece))
            if start + self.chunk_words >= len(words):
                break
        return windows


class OllamaKaluzaPipeline:
    """HTTP-Schnittstelle zu einem lokalen Ollama-Daemon, danach Sieb Omega."""

    SYSTEM = (
        "Du antwortest als mathematischer Assistent im EABC-Rahmen. "
        "Trenne feststehende Identitaeten, etwa -N(phi) = 1 in Q(sqrt(5)), "
        "von offenen Modellannahmen. Erfinde keine Rechenergebnisse."
    )

    def __init__(
        self,
        host: str = "http://localhost:11434",
        default_model: str = "llama3",
        timeout: float = 120.0,
        witness: Icosian = PHI,
        corpus: Optional[DirectoryTokenCorpus] = None,
    ) -> None:
        self.host = host.rstrip("/")
        self.default_model = default_model
        self.timeout = timeout
        self.guardrail = ToposGuardrail(witness)
        self.corpus = corpus

    def generate(
        self,
        prompt: str,
        model: Optional[str] = None,
        det_g: float = 1.0,
        witness: Optional[Icosian] = None,
    ) -> dict[str, Any]:
        selected = model or self.default_model
        element = self.guardrail.witness if witness is None else witness
        if not prompt.strip():
            return self._record(
                selected,
                raw_text="",
                det_g=0.0,
                witness=element,
                transport="rejected",
                error="Leerer Prompt.",
            )
        retrieval = self.corpus.retrieve(prompt) if self.corpus is not None else None
        spans = retrieval.spans if retrieval is not None else []
        sources = list(dict.fromkeys(span.source for span in spans))
        grover = retrieval.as_dict() if retrieval is not None else None
        full_prompt = self.corpus.as_prompt(spans, prompt) if self.corpus is not None else prompt
        try:
            raw_text = self._request(full_prompt, selected)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
            return self._record(
                selected,
                raw_text="",
                det_g=0.0,
                witness=element,
                transport="error",
                error=str(exc),
                sources=sources,
                grover=grover,
            )
        if not raw_text.strip():
            return self._record(
                selected,
                raw_text="",
                det_g=0.0,
                witness=element,
                transport="empty",
                error="Leere Modellantwort.",
                sources=sources,
                grover=grover,
            )
        return self._record(
            selected,
            raw_text=raw_text,
            det_g=det_g,
            witness=element,
            transport="live",
            error=None,
            sources=sources,
            grover=grover,
        )

    def _request(self, prompt: str, model: str) -> str:
        payload = {
            "model": model,
            "prompt": prompt,
            "system": self.SYSTEM,
            "stream": False,
        }
        request = urllib.request.Request(
            f"{self.host}/api/generate",
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            body = json.loads(response.read().decode("utf-8"))
        text = body.get("response", "")
        if not isinstance(text, str):
            raise json.JSONDecodeError("Feld 'response' ist kein String.", "", 0)
        return text

    def _record(
        self,
        model: str,
        raw_text: str,
        det_g: float,
        witness: Icosian,
        transport: str,
        error: Optional[str],
        sources: Optional[list[str]] = None,
        grover: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        verdict = self.guardrail.evaluate(det_g, witness)
        # Eine lebende Antwort mit gebrochener Norm oder det(G) <= 0 bleibt Rohtext.
        # refined_text gibt es nur im Hauptsieb.
        accepted = verdict.principal and transport == "live"
        return {
            "sieve_classification": verdict.as_text(),
            "is_artifact_free": accepted,
            "refined_text": raw_text if accepted else None,
            "raw_text": raw_text,
            "model_used": model,
            "transport": transport,
            "error": error,
            "galois_norm": witness.norm(),
            "det_g": det_g,
            "sources": sources or [],
            "grover": grover,
        }


def frage_einlesen() -> Optional[str]:
    """Liest eine beliebige Frage. Leere Eingabe oder Ende-Befehl schliesst den Dialog."""
    try:
        text = input("Frage: ").strip()
    except (EOFError, KeyboardInterrupt):
        print()
        return None
    if not text or text.lower() in {"q", "quit", "exit", "ende"}:
        return None
    return text


def antwort_ausgeben(result: dict[str, Any]) -> None:
    print(f"Modell: {result['model_used']}")
    print(f"Transport: {result['transport']}")
    print(f"Status: {result['sieve_classification']}")
    if result["error"]:
        print(f"Fehler: {result['error']}")
    grover = result.get("grover") or {}
    if grover:
        print(
            "Grover: "
            f"{grover['iterations']} Iterationen, "
            f"{grover['target_count']} Orakel-Treffer von {grover['span_count']} Abschnitten"
        )
    sources = result.get("sources") or []
    print("Quellen: " + (", ".join(sources) if sources else "keine passende Textstelle"))
    text = result["refined_text"] if result["refined_text"] is not None else result["raw_text"]
    print(f"Antwort: {text}")


if __name__ == "__main__":
    phi = Icosian(0, 1)
    if phi.norm() != -1 or -phi.norm() != 1:
        raise SystemExit("Galois-Identitaet -N(phi) = 1 nicht erfuellt.")

    corpus = DirectoryTokenCorpus(Path(__file__).resolve().parent).build()
    pipeline = OllamaKaluzaPipeline(default_model="llama3", corpus=corpus)
    print(
        "Korpus: "
        f"{corpus.file_count} Dateien, {len(corpus.spans)} Abschnitte, {corpus.token_count} Token."
    )
    print("EABC-Ollama. Beliebige Frage eingeben. Leere Zeile oder 'ende' beendet den Dialog.")
    while True:
        frage = frage_einlesen()
        if frage is None:
            break
        antwort_ausgeben(pipeline.generate(frage))
        print()
