# EABC: Diskrete Geometrie — kanonische Formulierung

**Status:** Kanonischer Geometrie-Layer (paper-ready Definitionen)  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59), `collatz/eabc-euklidische-hebung`  
**Tao-Labels:** Definition | Theorem | Analogie | Modellabbildung | Hypothese | Vermutung | Forschungsfrage

**Rolle:** **Primäres** EABC-Dokument. Alle Detailhypothesen (Zirkulation, Übergangsraum, Wigner-Legacy) verweisen hierher. Numerik: `collatz_eabc_hodge_eabc.py` (`Phi_E`, `flux_density_limit`, `synthesis_report`).

**Querverweise (Detail, nicht duplizieren):**
| Thema | Detaildokument | Code |
|-------|----------------|------|
| Zählgrößen $N_\pm$, Fehlerterm | `collatz_eabc_zirkulationshypothese.md` | `collatz_eabc_holonomie_fehlerterm.py` |
| Hodge, $L_{\mathrm{mag}}$ | `collatz_eabc_uebergangsraum.md` | `collatz_eabc_hodge_eabc.py` |
| Signierte Maßstruktur | `collatz_eabc_signierte_massstruktur.md` | `collatz_eabc_wigner_field.py` |
| Epistemik Physik vs. Modell | `collatz_eabc_epistemik_physik.md` | — |
| Wigner-Legacy (4-Pfad) | `collatz_eabc_wigner_analog.md` | `collatz_eabc_wigner_field.py` |
| Algebraischer Nebenzweig | `collatz_eabc_euklidische_hebung.md` | `collatz_eabc_euklid_hebung.py` |

---

## 0. These (boxed)

$$\boxed{\;\text{EABC = Theorie des priminduzierten Flusses auf dem Kreisgraphen } C_4\cong S^1.\;}$$

**Nicht** primär: „wie oft tritt ABCE vs. CEAB auf?" (Zählstatistik).  
**Sondern:** ob die Primstruktur einen **bevorzugten Windungssinn** auf dem elementaren Zyklus induziert.

---

## 1. Zentrale Frage

$$\boxed{\;\text{Erzeugt die Primstruktur eine bevorzugte Windungsrichtung auf } C_4\cong S^1\text{?}\;}$$

---

## 2. Paper-Definitionen (boxed)

**Kreisgraph:**
$$G_E = (V, E), \qquad V = \{E, A, B, C\}.$$

**Orientierte Kantenmengen** (Vorwärts- bzw. Rückwärtsorientierung auf $C_4$):
$$E^+ = \{EA, AB, BC, CE\}, \qquad E^- = \{EC, CB, BA, AE\}.$$

**Elementarer Vorwärtszyklus:**
$$E \xrightarrow{EA} A \xrightarrow{AB} B \xrightarrow{BC} C \xrightarrow{CE} E.$$

**5-Fenster-Zählgrößen** (Prim-Obergrenze $X$):
$$N_+(X) := \#\{n : p_{n+4}\le X,\; C_n^{(5)} = \mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n : p_{n+4}\le X,\; C_n^{(5)} = \mathrm{CEABC}\},$$
mit $C_n^{(5)} := (X_n, X_{n+1}, X_{n+2}, X_{n+3}, X_{n+4})$ und $X_k := \kappa(p_k)$.

**Zirkulation und normierte Observable:**
$$C_E(X) := N_+(X) - N_-(X), \qquad S_E(X) := N_+(X) + N_-(X),$$
$$W_E(X) := \frac{C_E(X)}{S_E(X)} = \frac{N_+(X) - N_-(X)}{N_+(X) + N_-(X)}.$$

**Normalisierter arithmetischer Magnetfluss:**
$$\Phi_E := \lim_{X\to\infty} \frac{N_+(X) - N_-(X)}{N_+(X) + N_-(X)} = \lim_{X\to\infty} \frac{C_E(X)}{S_E(X)} = \lim_{X\to\infty} W_E(X).$$

| $\Phi_E$ | Lesart |
|----------|--------|
| $\Phi_E = 0$ | keine bevorzugte Orientierung (asymptotische Symmetrie) |
| $\Phi_E \neq 0$ | stabile **arithmetische Orientierungsklasse** |

$$\boxed{\;\Phi_E \;\stackrel{?}{\neq}\; 0.\;}$$

**Label:** $G_E$, $E^\pm$, $N_\pm$, $C_E$, $S_E$, $W_E(X)$, $\Phi_E$ = **Definition**; $\Phi_E \neq 0$ = **Vermutung**.

**Legacy-Aliase:** $D_E := C_E$; $\mathrm{Hol}_E := \Phi_E$; `flux_density` = $W_E(X)$ = $S_E^{-1}C_E$ in Python.

**Numerik:** `Phi_E(X)` / `flux_density_limit(X)` in `collatz_eabc_hodge_eabc.py`.

---

## 3. Drei Ebenen

$$\boxed{\;\text{Zählung} \;\to\; \text{orientierter Fluss} \;\to\; \text{harmonische Kohomologieklasse.}\;}$$

| Stufe | Objekt | Lesart |
|-------|--------|--------|
| 1 — Zählung | $N_\pm(X)$ | wie oft welche Orientierung im Primfenster |
| 2 — orientierter Fluss | $\omega_E$ auf $E^+ \cup E^-$, $C_E = \oint_\gamma \omega_E$ | priminduzierte 1-Form auf Kanten |
| 3 — harmonische Klasse | $h \in H^1(C_4,\mathbb{Z})$, $\langle\omega_E, h\rangle$ | Kohomologie-Anteil des Flusses |

Topologie: $C_4 \cong S^1$, $H^1(S^1) \cong \mathbb{Z}$. Der Generator $h$ **existiert trivial**; die Frage ist, ob die Primfolge ihn **asymmetrisch** besetzt.

---

## 4. Vermutungen (boxed)

$$\boxed{\;\text{EABC-Vermutung: Der priminduzierte EABC-Fluss trägt eine nichttriviale harmonische Komponente.}\;}$$

$$\boxed{\;\langle\omega_E,\, h\rangle \neq 0,\quad h \text{ kanonische harmonische 1-Form auf dem Kreisgraphen.}\;}$$

- Existenz von $h$: **trivial** ($H^1 \cong \mathbb{Z}$).
- Nichttrivialität: Primfolge besetzt die Klasse **asymptotisch asymmetrisch** — nicht lokale Kantengeometrie allein.

Äquivalent zur Flussdichte-Vermutung:
$$\Phi_E = \lim_{X\to\infty} W_E(X) \;\stackrel{?}{\neq}\; 0 \quad\Longleftrightarrow\quad \lim_{X\to\infty} \langle\omega_E, h\rangle / \|\omega_E\| \;\stackrel{?}{\neq}\; 0.$$

**Numerik:** `harmonic_holonomy_component`, `inner_product_omega_h`.

---

## 5. Physikalische Analogie (nicht Wigner-primär)

**Primäre Analogien** (didaktisch, kein Physikanspruch):

| Analogie | EABC-Objekt | Eigenschaft |
|----------|-------------|-------------|
| Aharonov–Bohm | $\Phi_E$, $C_E = \oint \omega_E$ | lokal unsichtbar, global messbar |
| Sagnac | $\gamma^+$ vs. $\gamma^-$ auf $C_4$ | gegenläufige Umläufe, Schleifendefekt |
| magnetischer Laplace | $L_{\mathrm{mag}} = D - U$, $U_{ij} = A_{ij} e^{i\theta_{ij}}$ | chirale Phasen aus Orientierung |

$$\boxed{\;\text{Lokal unsichtbar, global messbarer Schleifendefekt — nicht „wer altert langsamer?".}\;}$$

**Sekundär (Legacy):** Wigner-Quasi-Wahrscheinlichkeit, 4-Pfad-Zählung $\#\mathrm{ABCE}-\#\mathrm{CEAB}$ — `collatz_eabc_wigner_analog.md`.

**Label:** AB / Sagnac / $L_{\mathrm{mag}}$ = **Analogie** / **Modellabbildung**; kein SRT, kein $c$ — `collatz_eabc_epistemik_physik.md`.

---

## 6. Harmonische Klasse und Hodge

Kanonischer Generator $h \in H^1(S^1)$: alle Vorwärtskanten $E^+$ tragen $+1$, normiert.

Priminduzierte Kanten-1-Form:
$$\omega_E = \sum_{e \in E^+ \cup E^-} w_e \cdot e, \qquad C_E = \sum_{e \in E^+} w_e - \sum_{e \in E^-} w_e.$$

**Hodge-Zerlegung** (Stub auf $C_4$):
$$\omega_E = d\phi + \delta\psi + h_{\mathrm{harm}}.$$

**Magnetischer Laplace:**
$$L_{\mathrm{mag}} = D - U, \qquad U_{ij} = A_{ij}\, e^{i\theta_{ij}},$$
$\theta_{ij}$ aus ABCEA/CEABC-Orientierungen.

**Numerik:** `harmonic_form_c4`, `discrete_hodge_decomposition`, `magnetic_laplacian`.

---

## 7. Epistemische Karte

| Inhalt | Label |
|--------|-------|
| $G_E$, $E^\pm$, $N_\pm$, $C_E$, $W_E(X)$, $\Phi_E$ | **Definition** |
| $C_4 \cong S^1$, $H^1 \cong \mathbb{Z}$, Existenz von $h$ | **Theorem** |
| Priminduktion von $\omega_E$ | **Experiment** |
| $\Phi_E \neq 0$, $\langle\omega_E,h\rangle \neq 0$ | **Vermutung** |
| AB / Sagnac / $L_{\mathrm{mag}}$ | **Analogie** / **Modellabbildung** |
| Wigner-4-Pfad-Legacy | **Analogie** (sekundär) |
| $\zeta$-Eigenmoden | **Forschungsfrage** (spekulativ) |

---

## 8. Roadmap: Code vs. offen

### Implementiert

| Funktion | Modul |
|----------|-------|
| $N_\pm$, $C_E$, $S_E$, $W_E(X)$ | `collatz_eabc_holonomie_fehlerterm.py` |
| `Phi_E`, `flux_density_limit`, `synthesis_report` | `collatz_eabc_hodge_eabc.py` |
| $E^+$, $E^-$ (`C4_EDGE_POSITIVE`, `C4_EDGE_NEGATIVE`) | `collatz_eabc_hodge_eabc.py` |
| $\langle\omega_E,h\rangle$, $L_{\mathrm{mag}}$ | `collatz_eabc_hodge_eabc.py` |

### Offen

- Asymptotischer Beweis $\Phi_E = 0$ oder $\neq 0$
- Lean: $\Phi_E$, harmonische Paarung, $L_{\mathrm{mag}}$
- Vollständige diskrete Hodge-Theorie auf erweitertem Übergangsraum

---

## 9. Experimente: Oktaeder-Umgebung

**Design:** `collatz_eabc_oktaeder_test.md` (**Modellabbildung**)

Regulärer Oktaeder $O_6 \subset \mathbb{R}^3$: C4-Äquator ($E,A,B,C$ auf $\pm e_1,\pm e_2$), polare Lift-Achse
$P^+$ (ABCEA) / $P^-$ (CEABC), Oktonion-Schalen-Gewicht $r_8(p)$ via $\theta_3^8$.

| Test | Observable | Erwartung |
|------|------------|-----------|
| Sanity | $\Phi_{\mathrm{oct,eq}}$ vs. $\Phi_E$ | exakte Übereinstimmung |
| Schale | $\Phi_{\mathrm{oct,shell}} = \sum \omega\, r_8(p) / \sum r_8(p)$ | gleiches Vorzeichen wie $\Phi_E$ |
| Harmonisch | $\langle\omega,h\rangle$ auf Äquator | $\neq 0$ bei endlichem $X$ |
| Pol | $P^+$ vs. $P^-$ Fluss | ABCEA-Präferenz |
| Spektral | $\mathrm{Spec}(L_{\mathrm{mag}})$ auf $O_6$ | near-zero-Moden |

```bash
python3 collatz_eabc_oktaeder_flux_test.py --max-p 1000000
pytest tests/test_eabc_oktaeder_flux.py -q
```

**JSON:** `collatz_eabc_oktaeder_flux_test.json`

---

## Artefakte

```bash
python3 collatz_eabc_hodge_eabc.py --max-p 1000000
python3 -c "from collatz_eabc_hodge_eabc import Phi_E; import json; print(json.dumps(Phi_E(10**6), indent=2))"
pytest tests/test_eabc_hodge_eabc.py -q
```

**JSON:** `collatz_eabc_diskrete_geometrie_synthesis.json` (via `synthesis_report`)
