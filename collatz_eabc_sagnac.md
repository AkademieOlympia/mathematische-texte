# EABC-Sagnac-Observable (kanonische Metapher)

**Status:** Definition + Vermutung + Hypothese + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

**Epistemische Abgrenzung:** Dieses Dokument ist **keine** Relativitätstheorie-Behauptung. Es überträgt die **kombinatorische Sagnac-Logik** (gegenläufige geschlossene Wege, Laufzeitdifferenz) auf den **diskreten EABC-Transport** entlang der Primfolge. Die Observable misst **arithmetische Orientierungsasymmetrie** zwischen ABCEA und CEABC — nicht Rotation im physikalischen Sinn.

**Querverweise:**
- `collatz_eabc_fehlerterm_hypothese.md` — **kanonische Endform:** $N_\pm$, $\Delta_E=D_E$, $S_E$, Hauptvermutung, Fehlerterm-Hypothese
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_transport.md` — $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH (nicht primäre Metapher)
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — Numerik $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$; `sagnac_report()`
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch: $\mathrm{Hol}_E=0$, Fehlerterm $\Delta_E$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54

---

## 0. Sagnac vs. Bell

| Metapher | Objekt | Eignung für EABC |
|----------|--------|------------------|
| **Bell/CHSH** | $E(a,b)$-Korrelationen zwischen **Messachsen** | **nicht** natürlich — drei binäre Lesarten auf Fenstern, kein geschlossener Weg |
| **Sagnac** | $\gamma^+$ vs. $\gamma^-$ auf **geschlossenem Pfad**, $\Delta T = T(\gamma^+)-T(\gamma^-)$ | **natürlich** — ABCEA vs. CEABC als gegenläufige Orientierungen desselben 4-Zyklus |

$$\boxed{\;\text{Primäre Metapher: Sagnac (gegenläufige Zyklen), nicht Bell (Achsenkorrelationen).}\;}$$

**Sekundär:** Bell/CHSH bleibt als kombinatorische Konsistenzprüfung auf $G_E$ dokumentiert (`collatz_eabc_bell_holonomie.md`).

---

## 1. Definitionen

**Restklassen modulo $12$:**
$$E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11\pmod{12}.$$

**Primfolge-Labels:**
$$X_n := \kappa(p_n),\qquad C_n^{(5)} := (X_n,X_{n+1},X_{n+2},X_{n+3},X_{n+4}).$$

**Zählgrößen** (für Prim-Obergrenze $X$):
$$N_+(X) := N_{\mathrm{ABCEA}}(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{ABCEA}\},$$
$$N_-(X) := N_{\mathrm{CEABC}}(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{CEABC}\}.$$

**Orientierte Wege (Sagnac-Pfade):**
- **$\gamma^+$ (positiv):** ABCEA = $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$
- **$\gamma^-$ (negativ):** CEABC = $C\!\to\!E\!\to\!A\!\to\!B\!\to\!C$

**Sagnac-Differenz (EABC-Sagnac-Fehlerterm):**
$$\Delta_E(X) := N_+(X) - N_-(X) = D_E(X).$$

**Normalisierte Sagnac-Observable** (analog $\Delta T/T$):
$$S_E(X) := \frac{N_+(X) - N_-(X)}{N_+(X) + N_-(X)} = \frac{\Delta_E(X)}{N_+(X)+N_-(X)}.$$

**Legacy-Aliase:** $\chi_{\mathrm{Hol}}(X):=S_E(X)$; $\mathrm{Hol}_E:=\lim_{X\to\infty} S_E(X)$.

**Interpretation:**
- $S_E(X)\to 0$: **rotationsfrei** (asymptotische Symmetrie der gegenläufigen Orientierungen)
- $S_E(X)\neq 0$: **bevorzugte Orientierung** auf endlichem $X$ (Fehlerterm, nicht Hauptterm)

**Normalisierter Fehlerterm (Fluktuationen):**
$$\widetilde{D}_E(X) := \frac{\Delta_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

**Label:** $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$ = **Definition**.

---

## 2. Hauptvermutung und Fehlerterm

**Hauptvermutung (Hauptterm).**
$$N_+(X) \sim N_-(X)\qquad (X\to\infty),$$
und damit
$$\lim_{X\to\infty}\frac{\Delta_E(X)}{N_+(X)+N_-(X)} = \lim_{X\to\infty} S_E(X) = \mathrm{Hol}_E = 0.$$

**Label:** Hauptvermutung = **Vermutung**.

**Fehlerterm-Hypothese (stärker).** $\Delta_E$ trägt einen **nichttrivialen Chebyshev-artigen Bias**, gesteuert durch Nullstellen der Dirichlet-$L$-Funktionen modulo $12$. Für endliche $X$ gilt typischerweise $\Delta_E(X)\neq 0$, auch wenn $\mathrm{Hol}_E=0$.

**Interessant:** nicht der Hauptterm $\lim \Delta_E/(N_++N_-)=0$, sondern die **Fluktuationen** $\Delta_E$, $\widetilde{D}_E$.

**Label:** Fehlerterm-Hypothese = **Hypothese**.

---

## 3. Boxed Hierarchie

$$\boxed{\;\text{Primzahlen} \;\to\; \text{EABC-Klassen} \;\to\; \text{Transport} \;\to\; \text{Zyklen} \;\to\; \text{Sagnac-Observable} \;\to\; \Delta_E(X)\;}$$

| Stufe | Objekt | Symbol |
|------:|--------|--------|
| 1 | Primzahlen $p_n>3$ | $p_n$ |
| 2 | EABC-Klasse | $X_n=\kappa(p_n)$ |
| 3 | Transport / Kante | $\tau_n=(X_n,X_{n+1})$, $G_E$ |
| 4 | Geschlossener Zyklus | $C_n^{(5)}$, ABCEA / CEABC |
| 5 | Sagnac-Observable | $S_E(X)$, $\gamma^\pm$ |
| 6 | Fehlerterm | $\Delta_E(X)=D_E(X)$ |

---

## 4. Boxed Schluss

$$\boxed{\;\text{ABCEA gegen CEABC ist eine Sagnac-Observable.}\;}$$

$$\boxed{\;\Delta_E(X) = D_E(X)\;\text{ist der EABC-Sagnac-Fehlerterm.}\;}$$

$$\boxed{\;\mathrm{Hol}_E = 0\;\text{im Hauptterm, aber}\;\Delta_E(X)\;\text{kann nichttrivial sein.}\;}$$

---

## 5. Python-Symbolzuordnung

| LaTeX | Python (kanonisch) | Legacy-Alias |
|-------|-------------------|--------------|
| $N_+(X)$ | `N_plus` | `N_ABCEA` |
| $N_-(X)$ | `N_minus` | `N_CEABC` |
| $\Delta_E(X)$ | `Delta_E` | `D_E` |
| $S_E(X)$ | `S_E` | `chi_Hol` |
| $\widetilde{D}_E(X)$ | `D_tilde_E` | — |
| Sagnac-Report | `sagnac_report()` | — |

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py` — `sagnac_report(X)` liefert $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$ am selben $X$.

---

## 6. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$ | **Definition** |
| ABCEA vs. CEABC als $\gamma^\pm$ | **Definition** (Sagnac-Analogie) |
| $N_+\sim N_-$ $\Rightarrow$ $S_E\to 0$ | **Vermutung** (Hauptvermutung) |
| $\Delta_E$ mit $L$-Funktions-Nullstellen mod $12$ | **Hypothese** (Fehlerterm-Hypothese) |
| Verhalten von $\widetilde{D}_E$ | **Experiment** |
| Bell/CHSH als Korrelations-Analogie | **Analogie** (sekundär, `collatz_eabc_bell_holonomie.md`) |
| Physikalische Sagnac-Rotation | **nicht behauptet** |

---

*Kanonsiche Lesart: Der EABC-Transport liefert zwei gegenläufige geschlossene Wege $\gamma^\pm$ auf demselben 4-Zyklus. Ihre Zähl-Differenz $\Delta_E$ ist der **Sagnac-Fehlerterm**; die normalisierte Observable $S_E$ misst Orientierungsasymmetrie analog $\Delta T/T$. Bell/CHSH bleibt als sekundäre Konsistenz-Metapher dokumentiert.*
