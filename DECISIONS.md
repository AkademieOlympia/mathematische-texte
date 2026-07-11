# DECISIONS — #Energiedoku / Bamberger Labor

> **Leitprinzip:** Nicht „das Modell ist bewiesen“, sondern: **Jede Behauptung ist einer reproduzierbaren Prüfspur zugeordnet.**

Stand: Juli 2026 (Synthesis 14 Seiten, Bridges 4 Seiten).

---

## Dokumentenarchitektur

| Dokument | Rolle | Inhalt |
|---|---|---|
| `Bridges2026_Submission.pdf` | **Schaufenster** | Visuelle/topologische Lesart (H12, Fano, Spatial Forcing als Motiv) |
| `Synthesis_2026.pdf` | **Tresor** | Invarianten, Audits, Methodik, explizite Schicht-Trennung |

Reviewer sollen nicht rhetorisch überzeugt, sondern in die **Auditorrolle** versetzt werden (Reproduktionsbefehle im Tresor).

---

## Entscheidung: Burn-in $p_{\min} = 101$

**Burn-in 101 is methodological, not theorem-level.**

- $101$ ist ein **spektral motivierter Stabilitätsfilter** (Einschwingen oberhalb transiente Kleinst-Primzahlen; erste dreistellige Primzahl).
- $255 = \mu(8) = 4\cdot 8^2 - 1$ ist das **modulare $n=8$-Scharnier** (Dual-Bit-Lemma, Schicht 1).
- Numerisch: $101 < 255$ — der Cutoff liegt **nicht** oberhalb der Mersenne-Grenze.
- Diese Ebenen dürfen **nicht vermischt** werden (kein „Cutoff wegen Scharnier“, kein „Scharnier wegen Cutoff“).

Implementierung: `eabc_mi_audit_v11_3.py`, `eabc_cluster_ellipse_h12.py`.

---

## Fünf epistemische Schichten

| Schicht | Status | Beispiel / Artefakt |
|---|---|---|
| **1 — Mathematische Identität** | Hart | $D_{\mathrm{KL}} \ge 0$ + `RuntimeError` (v11.3); $\mathcal{N}_{\text{crit}}=\{7,8\}$, $255=4\cdot 8^2-1$ (`eabc_dual_bit_inclusion_check.py`) |
| **2 — Numerischer Befund** | Reproduzierbar, empirisch | MI-Werte, $r=-0.7767$, spektrales Plateau (`synthesis_spectral_plateau_de.tex`, `zeros6.npy`) |
| **3 — Methodische Setzung** | Begründet, offenlegbar | $p_{\min}=101$, $B_{\mathrm{perm}}=10{,}000$ |
| **4 — Heuristische Deutung** | Markiert | Phase-Alternation $P(6)>P(0)$; Spatial-Forcing-Lesart |
| **5 — Forschungsprogramm** | Offen | Asymptotisches Forcing, $\varepsilon_{\mathrm{asymp}}>0$ (Dirichlet-$L$-Horizont) |

Was in Schicht 1 scheitert, darf nicht als struktureller Befund weitergeführt werden.

---

## Spektrales Plateau (Schicht 2)

Audit-Kette:

```
Terminal (Prim Konjungierte.py)
  ↔ Tabelle in synthesis_spectral_plateau_de.tex (§ spectral-plateau)
  ↔ zeros6.npy (2 001 052 Nullstellen)
  ↔ methodische Interpretation (Endlich-$x$-Rest, kein Trunkierungsbug)
```

Eingefrorene Referenzwerte (vollständige Nullstellenliste):

| $a$ | $\psi(a)$ | $\widehat{\psi}_{\mathrm{spektral}}(a)$ | $\Delta$ |
|---:|---:|---:|---:|
| 11 | 10.229909 | 9.030966 | −1.199 |
| 101 | 98.660432 | 96.352864 | −2.308 (~2.3 %) |

Fejér-/Lanczos-Glättung ändert das Plateau praktisch nicht → **kein** Gibbs-only-Problem.

---

## Reproduktion (Kern)

```bash
python3 eabc_mi_audit_v11_3.py --max-x 1000000 --B 10000
python3 eabc_dual_bit_inclusion_check.py
python3 eabc_cluster_ellipse_h12.py --plot
python3 "Hurrwitz Matrix.py"
python3 "Prim Konjungierte.py" --prime 101 --num-zeros 2001052 --compare-kernels
pdflatex Synthesis_2026.tex
pdflatex Bridges2026_Submission.tex
```

---

## Post-freeze program notes (Layer 4/5)

- Mannheim/eabc conformal synthesis is tracked as a Layer 4/5 program note in `notes/mannheim_eabc_conformal_synthesis_de.md`; it is explicitly outside the frozen audit dossier and does not extend the audited Layer 1/2 claims.
- Octonionic Collatz proof attempt lives on branch `post-freeze/octonionic-collatz-proof-attempt` under `post-freeze/octonionic-collatz-proof-attempt/`; O1/O2 are reproducible via `verify_o1_o2.py`; no Collatz closure claim; outside audit-freeze-2026.

---

## Was bewusst *nicht* behauptet wird

- Vollständiger analytischer HL-Beweis für Phase-Alternation mod 30
- Exakte Konvergenz $\widehat{\psi}_{\mathrm{spektral}}(x) \to \psi(x)$ bei endlichem $x$ ohne Riemann–Siegel-Rest
- „Beweis“ des gesamten Spatial-Forcing-Programms über alle Algebra-Körper

Stärke des Systems: **Prüfbarkeit**, nicht vermeintliche Unangreifbarkeit.
