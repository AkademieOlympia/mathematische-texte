# EABC: Diskrete Geometrie — kanonische Synthese

**Status:** Kanonischer Geometrie-Layer **vor** der nächsten konjekturellen Schicht  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59), `collatz/eabc-euklidische-hebung`  
**Tao-Labels:** Definition | Theorem | Analogie | Modellabbildung | Hypothese | Vermutung | Forschungsfrage

**Rolle:** Dieses Dokument **bündelt** alle geometrischen Fäden in **eine** kanonische Lesart. Detaildokumente bleiben Quellen; hier steht die **Programmarchitektur** der diskreten Geometrie auf dem EABC-Übergangsraum.

**Querverweise (integriert, nicht dupliziert):**
| Thema | Detaildokument | Code |
|-------|----------------|------|
| Zirkulation, $N_\pm$, $C_E$ | `collatz_eabc_zirkulationshypothese.md` | `collatz_eabc_holonomie_fehlerterm.py` |
| Signierte Maßstruktur | `collatz_eabc_signierte_massstruktur.md` | `collatz_eabc_hodge_eabc.py` |
| Übergangsraum / Hodge / $L_{\mathrm{mag}}$ | `collatz_eabc_uebergangsraum.md` | `collatz_eabc_hodge_eabc.py` |
| Holonomie-Stufen A/B/C | `collatz_eabc_holonomie_stufen.md` | `collatz_eabc_D_growth.py` |
| Epistemik Physik vs. Modell | `collatz_eabc_epistemik_physik.md` | — |
| Wigner / Wilson-Neulesung | `collatz_eabc_wigner_analog.md` | `collatz_eabc_wigner_field.py` |
| Chirale Polarisation | `collatz_eabc_chirale_polarisation.md` | `collatz_eabc_chirale_transport.py` |
| Potentialgeometrie Bohm/AB/Berry | `collatz_eabc_potential_geometrie.md` | `collatz_eabc_potential_phase.py` |
| Brachistochrone / Fermat | `collatz_eabc_brachistochrone.md` | `collatz_eabc_brachistochrone.py` |
| Kritische Abbildung | `collatz_eabc_kritische_abbildung.md` | `collatz_eabc_kritische_abbildung.py` |
| Spektralgeometrie | `collatz_eabc_zirkulation_spektral.md` | `collatz_eabc_graph_laplacian.py` |
| Algebraischer Nebenzweig | `collatz_eabc_euklidische_hebung.md` | `collatz_eabc_euklid_hebung.py` |

**Numerik (Synthese):** `synthesis_report(X)` in `collatz_eabc_hodge_eabc.py` → `collatz_eabc_diskrete_geometrie_synthesis.json`

---

## 1. Programmshift: Zählung → Geometrie auf gerichtetem $C_4 \cong S^1$

**Alt:** Primfolge $p_n$ als Punktstatistik; Korrelationen $E(a,b)$; Zählungen $\#(\cdots)$.

**Neu:** Geometrie der **gerichteten Übergänge** auf dem elementaren Zyklus
$$E \xrightarrow{EA} A \xrightarrow{AB} B \xrightarrow{BC} C \xrightarrow{CE} E \;\cong\; S^1.$$

Topologie: $H^1(S^1) \cong \mathbb{Z}$. Es existiert **immer** ein nichttrivialer harmonischer 1-Formen-Generator $h$.

$$\boxed{\;\text{Shift: } E(a,b) \;\to\; \oint_\gamma \omega_E \quad\text{(Paarstatistik } \to \text{ Zyklusgeometrie).}\;}$$

**Label:** $C_4 \cong S^1$, $H^1 \cong \mathbb{Z}$ = **Definition** / **Theorem** (Topologie); Priminduktion von $\omega_E$ = **Experiment**.

---

## 2. Fundamentalraum: Kanten $\{EA, AB, BC, CE\}$, nicht Knoten

$$\boxed{\;\text{Fundamentale Objekte sind Kanten auf } C_4,\text{ nicht Knoten }\{E,A,B,C\}\text{ allein.}\;}$$

| Ebene | Objekt | Rolle |
|-------|--------|-------|
| 0-Formen | $\phi$ auf $V=\{E,A,B,C\}$ | Knotenpotentiale (sekundär) |
| **1-Formen** | $\omega_E = w_{EA}\,EA + w_{AB}\,AB + w_{BC}\,BC + w_{CE}\,CE$ | **primäre Observable** |
| 2-Formen | $C_4$ als einzige Zelle | geschlossene Orientierung |

Der **Übergangsraum** (`collatz_eabc_uebergangsraum.md` §1) ist der Kantenraum des $C_4$-Komplexes. Der gerichtete Zustandsgraph $G_E=(V,E,w)$ aus der $W$-Matrix (`collatz_eabc_signierte_massstruktur.md` §2) lebt auf denselben Übergängen — Netzwerktheorie, nicht klassische Zahlentheorie.

**Numerik:** `edge_incidence_matrix_c4`, `omega_edge_from_holonomy`.

---

## 3. Observablen-Hierarchie: $D_E$, $C_E$, $S_E$, $W_E$ — Flussdichte, keine Wahrscheinlichkeit

| Symbol | Träger | Formel | Lesart |
|--------|--------|--------|--------|
| $D_E$ | 5-Zyklus | $N_+ - N_-$ | Zirkulationsdefekt = $C_E$ |
| $C_E$ | 5-Zyklus | $\oint_\gamma \omega_E = \Phi_E$ | diskreter AB-Fluss |
| $S_E$ | 5-Zyklus | $C_E/(N_++N_-)$ | normierte Flussdichte |
| $W_E$ | 4-Pfad | $\#\mathrm{ABCE}-\#\mathrm{CEAB}$ | Pfadorientierung |
| $S_W$ | 4-Pfad | $W_E/(N_+^{(4)}+N_-^{(4)})$ | normierte Pfad-Flussdichte |

**Neulesung (kanonisch):**
$$S_W,\; S_E \;\approx\; \tanh\Theta_E \quad\text{(Wilson-Schleifenwerte / normierte Flussdichten — \textbf{nicht} Quasi-Wahrscheinlichkeiten).}$$

4-Block ($W_E$, $S_W$) und 5-Block ($D_E=C_E$, $S_E$) strikt trennen (`collatz_eabc_wigner_analog.md` §0). Die historische Wigner-Analogie bleibt didaktisch; die **geometrische** Lesart ist Wilson-Fluss auf $S^1$.

**Label:** $C_E$, $S_E$, $S_W$ = **Definition**; $\tanh\Theta_E$-Analogie = **Analogie**; Wigner-Quasi-Wahrscheinlichkeit = **sekundäre Analogie**.

---

## 4. Harmonische Klasse $h$; Schlüssel: $\langle\omega_E, h\rangle \neq 0$ priminduziert

Auf $C_4 \cong S^1$ existiert der kanonische harmonische Generator $h \in H^1(S^1)$ (alle Vorwärtskanten $+1$, normiert).

$$\boxed{\;\text{Nicht } h \neq 0 \text{ (trivial), sondern priminduziertes } \langle\omega_E, h\rangle \neq 0 \text{ asymptotisch.}\;}$$

Die Primfolge induziert $\omega_E$ als Kanten-1-Form; die Paarung misst den **Fluss entlang der harmonischen Klasse** — bevorzugte Orientierung von ABCEA vs. CEABC.

**Hodge-Zerlegung** (Stub auf $C_4$):
$$\omega_E = d\phi + \delta\psi + h_{\mathrm{harm}}.$$

**Numerik:** `harmonic_form_c4`, `inner_product_omega_h`, `harmonic_holonomy_component`, `discrete_hodge_decomposition`.

**Label:** $h$ = **Definition**; $\langle\omega_E,h\rangle \neq 0$ bei endlichem $X$ = **Experiment**; asymptotisches $\neq 0$ = **Vermutung** (§7).

---

## 5. AB-Fluss und Sagnac: $C_E = \Phi_E$; $\oint \omega_E$

$$C_E(X) = N_+(X) - N_-(X) = \oint_{\gamma} \omega_E = \Phi_E.$$

**Analogie (didaktisch):** diskreter Aharonov–Bohm-Fluss — globale Orientierungsbilanz, **keine** lokale Kantenkrümmung (`collatz_eabc_potential_geometrie.md` §2).

**Sagnac-Struktur:**
$$\Delta_E \propto \oint_{\gamma} \omega_E$$
dieselbe Struktur wie orientiertes Linienintegral auf geschlossenem $C_4$-Loop (`collatz_eabc_sagnac.md` — Intuition only; Kern = Zirkulationshypothese).

Gegenläufige Orientierungen $\gamma^+$ (ABCEA) und $\gamma^-$ (CEABC) auf demselben $C_4$-Kreis; Chiralität als diskrete Helizität $\lambda=\pm 1$ (`collatz_eabc_chirale_polarisation.md`).

**Label:** $\Phi_E = C_E$ = **Definition**; AB-/Sagnac-Bild = **Analogie**.

---

## 6. Magnetischer Laplace $L_{\mathrm{mag}} = D - U$

Reeller Graph-Laplace (Knotenraum):
$$L = D - W, \qquad \lambda_1 \le \cdots \le \lambda_4.$$

Near-zero $\lambda_2 \approx 0$: fast erhaltene chirale Struktur (Dirac-/Index-/Hodge-Analog, `collatz_eabc_uebergangsraum.md` §5).

**Magnetischer Laplace** (besser für chirales Programm):
$$L_{\mathrm{mag}} = D - U, \qquad U_{ij} = A_{ij}\, e^{i\theta_{ij}},$$
mit $\theta_{ij}$ aus ABCEA/CEABC-Orientierungen bzw. $W_E(i,j;N)$.

**Numerik:** `laplacian_from_W`, `magnetic_laplacian`, `magnetic_laplacian_eigenvalues`.

**Label:** $L=D-W$ = **Definition**; $L_{\mathrm{mag}}$ = **Modellabbildung**; physikalisches Magnetfeld = **nicht behauptet**.

---

## 7. ZENTRALE VERMUTUNG (boxed)

$$\boxed{\;\omega_E = \text{priminduzierter orientierter Fluss auf } C_4 \cong S^1.\;}$$

$$\boxed{\;\lim_{X\to\infty} \frac{C_E(X)}{\#\{\text{erkannte Zyklen} \le X\}} = \lim_{X\to\infty} S_E(X) \;\stackrel{?}{\neq}\; 0.\;}$$

| Grenzfall | Folgerung |
|-----------|-----------|
| $\lim S_E = 0$ | keine bevorzugte globale Orientierung (Hauptterm-Symmetrie) |
| $\lim S_E \neq 0$ | nichttriviale **arithmetische Orientierungsklasse** |

**Konkurrenz:** Hauptterm-Vermutung $N_+ \sim N_-$ ⇒ $S_E \to 0$ (`collatz_eabc_zirkulationshypothese.md` §4). Die Orientierungsklassen-Vermutung ist die **geometrische** Umformulierung: trägt die Primfolge einen stabilen Fluss entlang $h$?

**Numerik:** `flux_density_limit`, `flux_density_series`, `synthesis_report`.

**Label:** Orientierungsklassen-Vermutung = **Vermutung** / **Forschungsfrage**.

---

## 8. Korkenzieher-Metapher (Modellabbildung)

Ein **Korkenzieher** ist ein orientierter Umlauf mit **aufgewickelter Phase**: jede Umdrehung akkumuliert einen festen Winkel, bevor ein neuer konzeptioneller Zustand (Herausziehen des Pfropfens) möglich wird. Auf $C_4 \cong S^1$ entspricht das der **Phasenakkumulation entlang des Zyklus** — $\oint \omega_E$, $\langle\omega_E, h\rangle$, $S_E \approx \tanh\Theta_E$ — **bevor** die nächste Schicht (arithmetische Orientierungsklasse, Lean-Beweis, $\zeta$-Eigenmoden) formal wird. Die Metapher bündelt: gerichteter Umlauf (ABCEA vs. CEABC), harmonische Klasse $h \in H^1(S^1)$, und den diskreten Fluss $\Phi_E$ als **gespeicherte Windung**, nicht als lokale Knotenstatistik.

**Label:** Korkenzieher = **Modellabbildung** (didaktisch); keine Physikbehauptung.

---

## 9. Epistemische Karte

| Inhalt | Label | Wo |
|--------|-------|-----|
| $C_4 \cong S^1$, $H^1 \cong \mathbb{Z}$ | Theorem / Definition | §1, `collatz_eabc_uebergangsraum.md` |
| Kantenfundament $\{EA,AB,BC,CE\}$ | Definition | §2 |
| $C_E = N_+ - N_-$, $S_E = C_E/N_{\mathrm{cycles}}$ | Definition | `collatz_eabc_zirkulationshypothese.md` |
| $S_W, S_E \approx \tanh\Theta_E$ (Wilson) | Analogie | §3, `collatz_eabc_signierte_massstruktur.md` |
| Wigner-Quasi-Wahrscheinlichkeit | Analogie (sekundär) | `collatz_eabc_wigner_analog.md` |
| $\langle\omega_E, h\rangle$ priminduziert | Experiment / Hypothese | §4 |
| AB-Fluss, Sagnac $\oint\omega$ | Analogie | §5, `collatz_eabc_potential_geometrie.md` |
| $L_{\mathrm{mag}} = D - U$ | Modellabbildung | §6 |
| **Zentralvermutung** $\lim S_E \neq 0$ | Vermutung | §7 |
| Korkenzieher | Modellabbildung | §8 |
| Bohm/AB/Berry, Brachistochrone | Analogie / Modell | `collatz_eabc_potential_geometrie.md`, `collatz_eabc_brachistochrone.md` |
| Kritische Abbildung $s_v(x)$ | Modellabbildung | `collatz_eabc_kritische_abbildung.md` |
| Chirale Helizität $\lambda=\pm 1$ | Modellabbildung | `collatz_eabc_chirale_polarisation.md` |
| Arithmetische Wigner-Negativität | Hypothese / offen | `collatz_eabc_signierte_massstruktur.md` §7 |
| Kein SRT, kein $c$, kein Zwillingsparadoxon | Abgrenzung | `collatz_eabc_epistemik_physik.md` |
| $\zeta$-Nullstellen als Eigenmoden | Forschungsfrage (spekulativ) | `collatz_eabc_uebergangsraum.md` §10 |

---

## 10. Roadmap: Lean/Code vs. offen

### Implementiert (Python)

| Funktion | Modul |
|----------|-------|
| $N_\pm$, $C_E$, $S_E$, $\widetilde{D}_E$ | `collatz_eabc_holonomie_fehlerterm.py` |
| $W_E(i,j;N)$, Informationsüberschuss | `collatz_eabc_wigner_field.py` |
| $h$, $\langle\omega_E,h\rangle$, Hodge-Stub | `collatz_eabc_hodge_eabc.py` |
| `flux_density_limit`, `synthesis_report` | `collatz_eabc_hodge_eabc.py` |
| $L_{\mathrm{mag}}$, `magnetic_laplacian_eigenvalues` | `collatz_eabc_hodge_eabc.py` |
| $\mathrm{Spec}(L_E)$ | `collatz_eabc_graph_laplacian.py` |
| $\phi_R/\phi_L$, chiraler Transport | `collatz_eabc_chirale_transport.py` |
| Bohm/AB/Berry-Stubs | `collatz_eabc_potential_phase.py` |
| Brachistochrone $T_R/T_L$ | `collatz_eabc_brachistochrone.py` |
| Kritische Abbildung, Holonomie-Sensor | `collatz_eabc_kritische_abbildung.py` |
| Wachstum Fall A/B/C | `collatz_eabc_D_growth.py` |

### Lean (PR #54 Kern)

| Datei | Inhalt |
|-------|--------|
| `CollatzEabc/HolonomieFehlerterm.lean` | Holonomie-Fehlerterm-Stub |
| `CollatzEabc/Core.lean`, `Kappa.lean`, `Mod12Matrix.lean` | $\kappa$, mod-$12$-Struktur |
| `CollatzEabc/Open.lean` | offene Ziele |

### Offen

- Asymptotischer Beweis $\lim S_E = 0$ oder $\neq 0$ (`collatz_eabc_holonomie_beweisversuch.md`)
- Vollständige diskrete Hodge-Theorie auf erweitertem Übergangsraum (nicht nur $C_4$-Stub)
- Lean: $\langle\omega_E,h\rangle$, Flussdichte, $L_{\mathrm{mag}}$
- Arithmetische Wigner-Negativität: Kantenfeld aus Marginalen rekonstruierbar?
- $\zeta$-Eigenmoden-Hypothese (spekulativ)

---

## Artefakte

```bash
python3 collatz_eabc_hodge_eabc.py --max-p 1000000
python3 -c "from collatz_eabc_hodge_eabc import synthesis_report; import json; print(json.dumps(synthesis_report(10**6), indent=2))"
pytest tests/test_eabc_hodge_eabc.py -q
```

**JSON:** `collatz_eabc_diskrete_geometrie_synthesis.json` (via `synthesis_report`)
