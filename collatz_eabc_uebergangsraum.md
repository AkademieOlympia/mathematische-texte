# EABC: Übergangsraum-Geometrie ($C_4 \cong S^1$, Hodge, magnetischer Laplace)

**Status:** Kanonische Geometrie + Hypothese + Experiment  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59), `collatz/eabc-euklidische-hebung`  
**Tao-Labels:** Definition | Analogie | Hypothese | Modellabbildung | Vermutung | Forschungsfrage

**Querverweise:**
- `collatz_eabc_zirkulationshypothese.md` — kanonisch: $N_\pm$, $C_E$, $D_E$, $S_E$, **Zentralvermutung Flussdichte**
- `collatz_eabc_signierte_massstruktur.md` — signierte Maßstruktur, Wigner-Negativität
- `collatz_eabc_wigner_analog.md` — historische Wigner-Lesart (sekundär zu Fluss/Wilson)
- `collatz_eabc_zirkulation_spektral.md` — $\mathrm{Spec}(L_E)$, diskrete 1-Form $\alpha$
- `collatz_eabc_hodge_eabc.py` — Numerik: $h$, $\langle\omega_E,h\rangle$, $L_{\mathrm{mag}}$, `flux_density_limit`
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 0. Paradigmenwechsel

$$D_E,\; C_E,\; S_E,\; W_E \;\Rightarrow\; \text{Frage: welcher zugrundeliegende Raum?}$$

**Shift:** Geometrie der **Übergänge** (Kanten, 1-Formen, Zyklen auf $S^1$), nicht Punktfolgen $p_n$.

$$\boxed{\;\text{Fundamentale Objekte sind Kanten } EA, AB, BC, CE \text{ auf } C_4,\text{ nicht Knoten }\{E,A,B,C\}\text{ allein.}\;}$$

---

## 1. Übergangsraum: Kantenfundament

**Elementarer gerichteter Zyklus:**
$$E \xrightarrow{EA} A \xrightarrow{AB} B \xrightarrow{BC} C \xrightarrow{CE} E.$$

| Ebene | Objekt | Lesart |
|-------|--------|--------|
| 0-Formen | $\phi$ auf $V=\{E,A,B,C\}$ | Knotenpotentiale |
| **1-Formen** | $\omega_E = w_{EA}\,EA + w_{AB}\,AB + w_{BC}\,BC + w_{CE}\,CE$ | **primäre Observable** |
| 2-Formen (Stub) | $C_4$ als einzige Zelle | geschlossene Orientierung |

**Label:** Kantenraum = **Definition**; Primfolge induziert $\omega_E$ = **Experiment**.

---

## 2. Topologie: $C_4 \cong S^1$, $H^1(S^1)\cong\mathbb{Z}$

Der $C_4$-Zyklus ist homöomorph zu $S^1$. Es existiert **trivial** ein nichttrivialer harmonischer 1-Formen-Generator $h$.

$$\boxed{\;\text{Schlüssel: nicht } h\neq 0,\text{ sondern priminduziertes } \langle\omega_E, h\rangle \neq 0 \text{ asymptotisch.}\;}$$

**Lesart:** Bevorzugte Orientierung der harmonischen Klasse durch die Primfolge — nicht lokale Kantengeometrie allein.

**Numerik:** `harmonic_form_c4`, `inner_product_omega_h` in `collatz_eabc_hodge_eabc.py`.

---

## 3. AB-Fluss: $C_E = N_+ - N_- = \Phi_E$

$$C_E(X) = N_+(X) - N_-(X) = \oint_{\gamma} \omega_E = \Phi_E$$

**Analogie (didaktisch):** diskreter Aharonov–Bohm-Fluss — globale Orientierungsbilanz, **keine** lokale Krümmung der Kantenlängen.

| Objekt | Formel | Bedeutung |
|--------|--------|-----------|
| $C_E$ | $N_+-N_-$ | Gesamtfluss |
| $S_E$ | $C_E/(N_++N_-)$ | normierte Flussdichte |
| $\widetilde{D}_E$ | $C_E/\sqrt{N_++N_-}$ | skaliertes Rauschen |

**Label:** $\Phi_E$ = **Definition**; AB-Analogie = **Analogie**.

---

## 4. $W_E$-Neulesung: Wilson-Fluss, nicht Quasi-Wahrscheinlichkeit

**4-Block (Pfad):**
$$S_W(N) = \frac{W_E(N)}{N_+^{(4)}+N_-^{(4)}} = \frac{N_+^{(4)}-N_-^{(4)}}{N_+^{(4)}+N_-^{(4)}} \approx \tanh\Theta_E.$$

**5-Block (Zyklus):** $S_E = C_E/(N_++N_-)$.

$$\boxed{\;S_W,\; S_E \;\text{sind normierte Flussdichten / Wilson-Schleifenwerte — nicht Wigner-Quasi-Wahrscheinlichkeiten.}\;}$$

**Verknüpfung:** `collatz_eabc_signierte_massstruktur.md` §1, §7; `collatz_eabc_wigner_analog.md` (historische Analogie).

---

## 5. Reeller Laplace $L = D - W$ und Spektrum

Aus der $4\times 4$-Übergangsmatrix $W_{ij}=\sum_n \chi_i(n)\chi_j(n+1)$:
$$L = D - W,\qquad D_{ii}=\sum_j W_{ij},\qquad \lambda_1\le\cdots\le\lambda_4.$$

**Near-zero $\lambda_2\approx 0$:** fast erhaltene chirale Struktur (Dirac-/Index-/Hodge-Analog).

**Numerik:** `laplacian_from_W` in `collatz_eabc_hodge_eabc.py`; Vergleich `collatz_eabc_graph_laplacian.py`.

---

## 6. Magnetischer Laplace $L_{\mathrm{mag}} = D - U$

Für das chirale Programm eignet sich der **magnetische** Laplace besser als der reelle $L=D-W$:

$$L_{\mathrm{mag}} = D - U,\qquad U_{ij} = A_{ij}\,e^{i\theta_{ij}},$$

wobei $\theta_{ij}$ aus ABCEA/CEABC-Orientierungen bzw. $W_E(i,j;N)$ stammt.

**Label:** $L_{\mathrm{mag}}$ = **Modellabbildung**; physikalisches Magnetfeld = **nicht behauptet**.

**Numerik:** `magnetic_laplacian`, `magnetic_phase_matrix`.

---

## 7. Hodge-Zerlegung und harmonische Holonomie

$$\omega_E = d\phi + \delta\psi + h.$$

Auf $C_4$ trägt der **harmonische** Anteil $h$ geschlossene Schleifen → Holonomie.

**Offene Frage:** Existiert ein nichttrivialer harmonischer 1-Form-Anteil in $\omega_E$, der mit $C_E$ korreliert?

**Numerik:** `discrete_hodge_decomposition`, `harmonic_holonomy_component`.

---

## 8. Sagnac-Struktur

$$\Delta_E \propto \oint_{\gamma} \omega_E$$

dieselbe Struktur wie orientiertes Linienintegral auf geschlossenem $C_4$-Loop — **Zirkulation**, nicht Korrelation.

**Verknüpfung:** `collatz_eabc_sagnac.md` (Intuition only).

---

## 9. ZENTRALE VERMUTUNG (boxed)

$$\boxed{\;\omega_E = \text{priminduzierter orientierter Fluss auf } C_4\cong S^1.\;}$$

$$\boxed{\;\lim_{X\to\infty}\frac{C_E(X)}{\#\{\text{erkannte Zyklen}\le X\}} = \lim_{X\to\infty} S_E(X) \;\stackrel{?}{\neq}\; 0.\;}$$

| Grenzfall | Folgerung |
|-----------|-----------|
| $\lim S_E = 0$ | keine bevorzugte globale Orientierung (Hauptterm-Symmetrie) |
| $\lim S_E \neq 0$ | nichttriviale **arithmetische Orientierungsklasse** |

**Numerik:** `flux_density_limit`, `flux_density_series` — Schätzung bei $10^3\ldots 10^6$.

**Label:** Orientierungsklassen-Vermutung = **Vermutung** / **Forschungsfrage**.

---

## 10. Spekulativ: $\zeta$-Nullstellen als Eigenmoden

**Spekulation (nicht kanonisch):** Nullstellen der $\zeta$-Funktion als Eigenmoden der Übergangsgeometrie $\omega_E$ — nicht nur $\pi(x)$-Oszillationen.

**Label:** $\zeta$-Eigenmoden = **Forschungsfrage** (spekulativ).

---

## 11. Artefakte

```bash
python3 collatz_eabc_hodge_eabc.py --max-p 1000000
pytest tests/test_eabc_hodge_eabc.py -q
```

**JSON:** `collatz_eabc_hodge_eabc.json`
