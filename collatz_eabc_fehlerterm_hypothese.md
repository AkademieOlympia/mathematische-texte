# EABC-Holonomie und Fehlerterm-Hypothese (Endform)

**Status:** Kanonische Endform (Definition + Hauptvermutung + Fehlerterm-Hypothese)  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment

**Querverweise:**
- `collatz_eabc_sagnac.md` — **kanonische primäre Metapher:** $\gamma^\pm$, $\Delta_E$, $S_E$ (Sagnac-Observable)
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch (Lemma-Skizzen, HL-Argument, Satzskizze)
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH; $P_{\mathrm{same}}$, CHSH auf $G_E$ (nicht primäre Metapher)
- `collatz_eabc_bell_inequality_test.py` / `.json` — Numerik $P_{\mathrm{same}}^{\mathrm{win}}$, $\mathcal{B}_{\mathrm{marg}}$, $S_{\mathrm{EABC}}$, **`de_bell_combined_report`**
- `collatz_eabc_core/CollatzEabc/HolonomieFehlerterm.lean` — Lean: Lückenmuster, `N_\pm`/`D_E` auf Listen; Prime/CHSH-Skeleton
- `collatz_eabc_transport.md` — Übergangsgraph $G_E$, Transport $T_n$
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — Numerik $N_\pm$, $D_E$, $\widetilde{D}_E$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54

---

## 0. Setup

**Restklassen modulo $12$:**
$$E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11\pmod{12}.$$

**Primfolge-Labels:**
$$X_n := \kappa(p_n),\qquad C_n^{(5)} := (X_n,X_{n+1},X_{n+2},X_{n+3},X_{n+4}).$$

(Wir schreiben kurz $C_n$ für das geschlossene 5-Fenster.)

**Zählgrößen** (für Prim-Obergrenze $X$):
$$N_+(X) := \#\{n:\,p_{n+4}\le X,\; C_n=\mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n:\,p_{n+4}\le X,\; C_n=\mathrm{CEABC}\}.$$

**Legacy-Aliase:** $N_{\mathrm{ABCEA}}:=N_+$, $N_{\mathrm{CEABC}}:=N_-$.

**Sagnac-Lesart** (`collatz_eabc_sagnac.md`):
$$\Delta_E(X) := N_+(X) - N_-(X) = D_E(X),\qquad
S_E(X) := \frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}.$$
ABCEA = $\gamma^+$ ($A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$); CEABC = $\gamma^-$ ($C\!\to\!E\!\to\!A\!\to\!B\!\to\!C$). $S_E\to 0$: rotationsfrei; $S_E\neq 0$: bevorzugte Orientierung auf endlichem $X$.

---

## 1. Hauptvermutung (Hauptterm)

**Hauptvermutung.** Asymptotische Symmetrie der gegenläufigen Zyklusorientierungen:
$$N_+(X) \sim N_-(X)\qquad (X\to\infty),$$
und damit
$$\mathrm{Hol}_E := \lim_{X\to\infty}\frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)} = 0.$$

**Label:** Hauptvermutung = **Vermutung** (unter HL-Äquidistribution / mod-$12$-Symmetrie; siehe Beweisversuch).

---

## 2. Fehlerterm (EABC-Sagnac-Fehlerterm)

**Definition (Fehlerterm / Sagnac-Differenz).**
$$\Delta_E(X) := N_+(X) - N_-(X),\qquad D_E(X):=\Delta_E(X).$$

Für endliche $X$ gilt typischerweise $\Delta_E(X)\neq 0$, auch wenn $\mathrm{Hol}_E=0$. Interessant sind die **Fluktuationen** $\Delta_E$, $\widetilde{D}_E$, nicht der Hauptterm $\lim \Delta_E/(N_++N_-)=0$.

$$\boxed{\;\Delta_E(X) = D_E(X)\;\text{ist der EABC-Sagnac-Fehlerterm.}\;}$$

**Label:** $\Delta_E$, $D_E$ = **Definition**.

---

## 3. Fehlerterm-Hypothese (stärker)

**Fehlerterm-Hypothese (stärker).** $D_E$ trägt einen **nichttrivialen Chebyshev-artigen Bias**, gesteuert durch die **Nullstellen der Dirichlet-$L$-Funktionen modulo $12$** (bzw. der zugehörigen Charaktere auf $(\mathbb{Z}/12\mathbb{Z})^\times$).

Qualitativ: wie beim klassischen Chebyshev-Bias mod $4$ kann die **absolute Differenz** vorzeichenbehaftet und strukturiert oszillieren, während der **normierte Hauptterm** $\mathrm{Hol}_E$ verschwindet.

**Label:** Fehlerterm-Hypothese = **Hypothese** (stärker als Hauptvermutung allein).

---

## 4. Normalisierte Observable

**Definition (normalisierter Fehlerterm).**
$$\widetilde{D}_E(X) := \frac{N_+(X)-N_-(X)}{\sqrt{N_+(X)+N_-(X)}} = \frac{D_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

Äquivalent: $\chi_{\mathrm{Hol}}(X) = D_E(X)/(N_+(X)+N_-(X))$ und $\widetilde{D}_E(X) = \chi_{\mathrm{Hol}}(X)\sqrt{N_+(X)+N_-(X)}$.

**Label:** $\widetilde{D}_E$ = **Definition**.

---

## 5. Zentrale Frage

**Zentrale Frage.** Verhält sich $\widetilde{D}_E(X)$ wie **reines Rauschen** (z. B. $|\widetilde{D}_E|=O(1)$ ohne stabile Vorzeichenstruktur), oder zeigt sie **stabile Vorzeichenasymmetrie / Oszillationsstruktur** gegenüber Isotropie- und Shuffle-Nullmodellen?

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py` — Zeitreihe $\widetilde{D}_E(X)$ über $X\in\{10^3,10^4,\ldots\}$; Vergleich mit Chebyshev mod $4$.

**Label:** Grenzverhalten von $\widetilde{D}_E$ = **Experiment** / offene **Forschungsfrage**.

---

## 6. Interpretation

**Interpretation.** EABC-Holonomie ist **sekundärer Bias** im „Prime Race“ zwischen den **gegenläufigen Orientierungen** desselben geschlossenen 4-Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$ — **nicht** ein permanenter Hauptterm.

Die empirische ABCEA/CEABC-Asymmetrie wird damit als **Fehlerterm-/Bias-Phänomen** reklassifiziert, nicht als $\mathrm{Hol}_E\neq 0$.

---

## 7. Boxed Schluss

$$\boxed{\;\mathrm{Hol}_E = 0\;\text{im Hauptterm, aber}\;D_E(X)\;\text{kann nichttrivial sein.}\;}$$

$$\boxed{\;\text{ABCE/CEAB-Asymmetrie} \;\Rightarrow\; \text{Fehlerterm/Bias, nicht Hauptterm-Holonomie.}\;}$$

---

## 8. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $N_\pm$, $D_E$, $\widetilde{D}_E$ | **Definition** |
| $N_+\sim N_-$ $\Rightarrow$ $\mathrm{Hol}_E=0$ | **Vermutung** (Hauptvermutung) |
| $D_E$ mit $L$-Funktions-Nullstellen mod $12$ | **Hypothese** (Fehlerterm-Hypothese) |
| Verhalten von $\widetilde{D}_E$ | **Experiment** |
| $\mathrm{Hol}_E\neq 0$ (stark, legacy) | **Hypothese** (`collatz_eabc_zyklus_holonomie.md` §7) |

---

## 9. Primäre Metapher: Sagnac; sekundär: Bell/CHSH

**Kanonische Metapher (Sagnac):** ABCEA vs. CEABC sind gegenläufige geschlossene Wege $\gamma^\pm$; $\Delta_E=N_+-N_-$ ist die Sagnac-Differenz, $S_E=\Delta_E/(N_++N_-)$ die normalisierte Observable (analog $\Delta T/T$). Details: `collatz_eabc_sagnac.md`.

$$\boxed{\;\text{ABCEA gegen CEABC ist eine Sagnac-Observable.}\;}$$

**Sekundäre Analogie (Bell/CHSH, kein Theorem):** Persistenter Fehlerterm $\Delta_E(X)\neq 0$ bei $S_E\to 0$ entspricht **nicht-faktorisierbaren** Holonomie-Resten auf $G_E$ — parallel zu $|S_{\mathrm{EABC}}|>2$ als arithmetisches Korrelationssignal (`collatz_eabc_bell_holonomie.md` §12).

| Größe | Rolle | Label |
|-------|-------|-------|
| $\chi_{\mathrm{Hol}}(X)=D_E/(N_++N_-)$ | normierter Hauptterm → $0$ | **Vermutung** |
| $\widetilde{D}_E(X)$ | skaleninvariantes Bias-Signal | **Definition** / **Experiment** |
| $S_{\mathrm{EABC}}$ | CHSH-Analog auf gemeinsamem ABCE-Träger | **Definition** / **Experiment** |
| $P_{\mathrm{same}}^{\mathrm{hol}}$ | Pfad↔Holonomie-Kohärenz | **Definition** / **Hypothese** |

**Numerik:** `de_bell_combined_report(X)` in `collatz_eabc_bell_inequality_test.py` und `collatz_eabc_holonomie_fehlerterm.py` — $D_E$, $\widetilde{D}_E$, $S_{\mathrm{EABC}}$ am selben $X$.

**Lean:** `CollatzEabc/HolonomieFehlerterm.lean` — kombinatorische $N_\pm$/`D_E` auf endlichen Listen (**bewiesen**); Primzählung und CHSH-LHV-Bound (**`sorry`**).

---

*Kanonsiche Endform: Dieses Dokument fixiert die **operative Lesart** für PR #54. Primäre Metapher: `collatz_eabc_sagnac.md`. Details und Beweisskizzen: `collatz_eabc_holonomie_beweisversuch.md`. Sekundäre Bell-Brücke: `collatz_eabc_bell_holonomie.md` §12.*
