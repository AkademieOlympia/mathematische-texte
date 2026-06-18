# EABC-Fehlerterm-Hypothese (Teilhypothese)

**Status:** Teilhypothese — eingebettet in `collatz_eabc_zirkulationshypothese.md` §5  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment

**Primäre Hypothese:** `collatz_eabc_zirkulationshypothese.md` — kanonische EABC-Zirkulationshypothese. Dieses Dokument vertieft den **Fehlerterm** $D_E(X)$ und die normalisierte Observable $\widetilde{D}_E(X)$.

**Querverweise:**
- `collatz_eabc_epistemik_physik.md` — epistemische Drei Ebenen; zentrale Frage $\mathrm{Hol}_E$ und Prime-Race-Chiralität
- `collatz_eabc_zirkulationshypothese.md` — **kanonische Hypothese:** $N_\pm$, $C_E$, $D_E$, $S_E$, Hauptvermutung, Prime-Race-Box
- `collatz_eabc_zirkulation_spektral.md` — Spektralgeometrie, diskrete 1-Form $\alpha$, $\mathrm{Spec}(L_E)$
- `collatz_eabc_sagnac.md` — **Intuition only:** Sagnac-Bild für $\gamma^\pm$
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — Numerik $N_\pm$, $D_E$, $\widetilde{D}_E$
- `collatz_eabc_sagnac_circulation.py` — $C_E(X)$, $\omega(e)$, $\alpha$
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 0. Einbettung in die Zirkulationshypothese

Die Fehlerterm-Hypothese ist **Teil** der kanonischen EABC-Zirkulationshypothese (`collatz_eabc_zirkulationshypothese.md`):

$$\boxed{\;D_E(X)\;\text{ist ein Prime Race zwischen zwei Orientierungen desselben EABC-Zyklus.}\;}$$

Setup, $N_\pm$, $C_E$, $S_E$ und Hauptvermutung: dort §1–4. Hier: Vertiefung von $D_E$, $\widetilde{D}_E$ und experimentelle Fragen.

---

## 1. Definitionen (Fehlerterm)

**Zählgrößen** (für Prim-Obergrenze $X$):
$$N_+(X) := \#\{n:\,p_{n+4}\le X,\; C_n=\mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n:\,p_{n+4}\le X,\; C_n=\mathrm{CEABC}\}.$$

**Fehlerterm / Zirkulations-Differenz:**
$$D_E(X) := N_+(X) - N_-(X) = C_E(X) = \Delta_E(X).$$

**Normalisierte Observable:**
$$S_E(X) := \frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}.$$

**Normalisierter Fehlerterm:**
$$\widetilde{D}_E(X) := \frac{N_+(X)-N_-(X)}{\sqrt{N_+(X)+N_-(X)}} = \frac{D_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

**Label:** $D_E$, $\widetilde{D}_E$ = **Definition**.

---

## 2. Fehlerterm-Hypothese (stärker)

**Fehlerterm-Hypothese (stärker).** $D_E$ trägt einen **nichttrivialen Chebyshev-artigen Bias**, gesteuert durch die **Nullstellen der Dirichlet-$L$-Funktionen modulo $12$** (bzw. der zugehörigen Charaktere auf $(\mathbb{Z}/12\mathbb{Z})^\times$).

Qualitativ: wie beim klassischen Chebyshev-Bias mod $4$ kann die **absolute Differenz** vorzeichenbehaftet und strukturiert oszillieren, während der **normierte Hauptterm** $S_E(X)\to 0$.

**Label:** Fehlerterm-Hypothese = **Hypothese** (stärker als Hauptvermutung allein).

---

## 3. Zentrale Frage

**Zentrale Frage.** Verhält sich $\widetilde{D}_E(X)$ wie **reines Rauschen** (z. B. $|\widetilde{D}_E|=O(1)$ ohne stabile Vorzeichenstruktur), oder zeigt sie **stabile Vorzeichenasymmetrie / Oszillationsstruktur** gegenüber Isotropie- und Shuffle-Nullmodellen?

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py` — Zeitreihe $\widetilde{D}_E(X)$ über $X\in\{10^3,10^4,\ldots\}$; Vergleich mit Chebyshev mod $4$.

**Label:** Grenzverhalten von $\widetilde{D}_E$ = **Experiment** / offene **Forschungsfrage**.

---

## 4. Interpretation

**Interpretation.** Die empirische ABCEA/CEABC-Asymmetrie wird als **Fehlerterm-/Bias-Phänomen** reklassifiziert, nicht als permanente Hauptterm-Holonomie ($S_E\to 0$).

$$\boxed{\;\mathrm{Hol}_E = 0\;\text{im Hauptterm, aber}\;D_E(X)\;\text{kann nichttrivial sein.}\;}$$

---

## 5. Sekundäre Analogien

**Bell/CHSH** (`collatz_eabc_bell_holonomie.md`): Persistenter Fehlerterm $\Delta_E(X)\neq 0$ bei $S_E\to 0$ entspricht **nicht-faktorisierbaren** Holonomie-Resten auf $G_E$ — **sekundäre** Analogie, kein primärer Kern.

**Numerik:** `de_bell_combined_report(X)` in `collatz_eabc_bell_inequality_test.py` und `collatz_eabc_holonomie_fehlerterm.py`.

**Lean:** `CollatzEabc/HolonomieFehlerterm.lean` — kombinatorische $N_\pm$/`D_E` auf endlichen Listen (**bewiesen**); Primzählung und CHSH-LHV-Bound (**`sorry`**).

---

## 6. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $D_E$, $\widetilde{D}_E$ | **Definition** |
| $D_E$ mit $L$-Funktions-Nullstellen mod $12$ | **Hypothese** |
| Verhalten von $\widetilde{D}_E$ | **Experiment** |
| Bell/CHSH | **Analogie** (sekundär) |

---

*Teilhypothese: Fehlerterm $D_E$ innerhalb der kanonischen Zirkulationshypothese `collatz_eabc_zirkulationshypothese.md`. Spektralgeometrie: `collatz_eabc_zirkulation_spektral.md`.*
