# EABC-Sagnac-Observable (kanonische Metapher)

**Status:** Definition + Vermutung + Hypothese + Experiment  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

**Epistemische Abgrenzung:** Dieses Dokument ist **keine** Relativitätstheorie-Behauptung. Es überträgt die **kombinatorische Sagnac-Logik** (gegenläufige geschlossene Wege, Laufzeitdifferenz) auf den **diskreten EABC-Transport** entlang der Primfolge. Die Observable misst **arithmetische Orientierungsasymmetrie** zwischen ABCEA und CEABC — nicht Rotation im physikalischen Sinn.

**Querverweise:**
- `collatz_eabc_fehlerterm_hypothese.md` — **kanonische Endform:** $N_\pm$, $\Delta_E=D_E$, $S_E$, Hauptvermutung, Fehlerterm-Hypothese
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_transport.md` — $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH (nicht primäre Metapher)
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — Numerik $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$; `sagnac_report()`
- `collatz_eabc_sagnac_circulation.py` / `.json` — Zirkulation $C_E(X)$, Kantenorientierung $\omega$, diskrete 1-Form $A$
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch: $\mathrm{Hol}_E=0$, Fehlerterm $\Delta_E$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 0. Sagnac vs. Bell

| Metapher | Objekt | Eignung für EABC |
|----------|--------|------------------|
| **Bell/CHSH** | $E(a,b)$-Korrelationen zwischen **Messachsen** $(a,b)$, $(a,b')$, … | **nicht** natürlich — erzwingt **verschiedene Messkontexte**; für EABC künstlich |
| **Holonomie** | Paralleltransport entlang geschlossenem Pfad | **teilweise** — korrekte Stufe, aber ohne explizite Orientierungsgegenüberstellung |
| **Sagnac** | $\gamma^+$ vs. $\gamma^-$ auf **demselben geometrischen Kreis**, $\Delta T = T(\gamma^+)-T(\gamma^-)$ | **natürlich** — ABCEA vs. CEABC als $\gamma^\pm$ auf demselben 4-Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$ |

**Warum Sagnac $>$ Bell für EABC:**
- Bell braucht verschiedene Messkontexte $(a,b)$, $(a,b')$, … — für EABC **erzwungen**, nicht intrinsisch.
- EABC: **dieselbe** Primfolge, **dieselben** Restklassen — kein Kontextwechsel.
- Sagnac: $\gamma^+$ vs. $\gamma^-$ auf **demselben** geometrischen Kreis — passt exakt zu ABCEA vs. CEABC.

$$\boxed{\;\text{Primäre Metapher: Sagnac (gegenläufige Zyklen), nicht Bell (Achsenkorrelationen).}\;}$$

**Sekundär:** Bell/CHSH bleibt als kombinatorische Konsistenzprüfung auf $G_E$ dokumentiert (`collatz_eabc_bell_holonomie.md`).

### 0.1 Forschungsprogression

$$\kappa(p) \;\to\; Q\text{-Signatur} \;\to\; \chi_E \;\to\; D_E(X)=N_+-N_- \;\to\; C_E(X)$$

$D_E(X)=N_+(X)-N_-(X)$ ist die **erste geometrische Flussgröße** auf dem EABC-Transportgraphen; $C_E(X)$ formalisiert dieselbe Größe als **diskrete Zirkulation** (§8).

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

## 5. Graphformulierung $\Omega(C)$

**Gerichteter Transportgraph:**
$$G_E=(V,E),\qquad V=\{E,A,B,C\}.$$

**Kantenorientierung** auf kanonischen Zykluskanten des 4-Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$:
$$\omega(e)\in\{+1,-1,0\},\qquad
\omega(A\!\to\!B)=\omega(B\!\to\!C)=\omega(C\!\to\!E)=\omega(E\!\to\!A)=+1,$$
$$\omega(B\!\to\!A)=\omega(C\!\to\!B)=\omega(E\!\to\!C)=\omega(A\!\to\!E)=-1,$$
sonst $\omega=0$. Lückenmuster mod $12$: $(2,4,2,4)$.

**Geschlossener Zyklus** $C=(v_1,v_2,v_3,v_4,v_1)$:
$$\Omega(C) := \prod_{i=1}^{4}\omega(v_i,v_{i+1})\qquad (v_5:=v_1).$$

**Sagnac-Paare auf demselben 4-Zyklus** $(A,B,C,E)$:
- **ABCEA** ($\gamma^+$): $\omega(\gamma)=+1$
- **CEABC** ($\gamma^-$): $\omega(\gamma)=-1$

Beide Wörter tragen dasselbe Lückenmuster $(2,4,2,4)$; die Orientierung ist **zyklisch verschoben**, nicht nur rotiert im Kantenprodukt.

**Label:** $G_E$, $\omega(e)$, $\Omega(C)$ = **Definition**.

---

## 6. Diskrete 1-Form

**1-Form** auf den vier Zykluskanten:
$$A(E,A),\; A(A,B),\; A(B,C),\; A(C,E).$$

**Normierung** (diskrete Zirkulation mit Sagnac-Orientierung $\omega(\gamma)$):
$$A(i\!\to\!j) := \frac{\omega(i,j)}{4}\quad\text{auf Zykluskanten},\qquad
\oint_{\gamma} A := \omega(\gamma)\sum_{e\in\gamma}|A(e)|.$$

Dann:
$$\oint_{\mathrm{ABCEA}} A = +1,\qquad
\oint_{\mathrm{CEABC}} A = -\oint_{\mathrm{ABCEA}} A.$$

**Label:** $A$ = **Definition**.

---

## 7. Analytische Lesart

Unter **HL-Äquidistribution** und mod-$12$-Symmetrie:
$$N_+(X)\sim N_-(X)\quad\Rightarrow\quad S_E(X)\to 0$$
— **nichtrotierendes Interferometer** (asymptotische Symmetrie der gegenläufigen Orientierungen).

Die **Information** liegt in den **Fluktuationen** $\Delta_E$, $\widetilde{D}_E$ (Prime-Race-Fehlerterm), nicht im Hauptterm $\mathrm{Hol}_E=0$.

**Label:** $S_E\to 0$ = **Vermutung**; Fluktuationsstruktur = **Hypothese** / **Experiment**.

---

## 8. Zirkulation $C_E(X)$

**Definition (diskrete EABC-Zirkulation).**
$$C_E(X) := \sum_{\gamma\ \mathrm{erkannt}} \omega(\gamma),$$
wobei die Summe über alle erkannten 5-Fenster $C_n^{(5)}$ mit $\omega(\gamma)=+1$ (ABCEA) bzw. $\omega(\gamma)=-1$ (CEABC) läuft und $p_{n+4}\le X$.

**Identität:**
$$C_E(X) = N_+(X) - N_-(X) = \Delta_E(X) = D_E(X).$$

**Normalisierung:**
$$S_E(X) = \frac{C_E(X)}{N_+(X)+N_-(X)}.$$

**Nächster Schritt (Forschungsvorschlag):** Untersuchen, ob der Fehlerterm $C_E(X)$ eine **Chebyshev- / mod-$q$- / Dirichlet-$L$-Bias-Struktur** trägt (analog Prime Race mod $4$).

**Label:** $C_E$ = **Definition**; Bias-Struktur = **Hypothese**.

**Experiment:** `collatz_eabc_sagnac_circulation.py` — `circulation_C_E(X)`, `circulation_report(X)`.

---

## 9. Normschalen-Brücke

**Spektralgeometrische Kette** (`collatz_eabc_quaternion_mass_hypothese.md`):
$$\Sigma_n \;\to\; \mu_n \;\to\; G_n$$
mit EABC-Transport auf der Schale $\Sigma_n$.

**Brücke:** $D_E(n)$ bzw. $C_E(X)$ als **diskrete Zirkulation auf der Normschale** — die Sagnac-Differenz misst Orientierungsasymmetrie des EABC-Transports entlang der Primfolge, analog zu einem Defekt $D(n)=I(\mu_n)-I_{\mathrm{ref}}(n)$ auf $\Sigma_n$.

**Label:** Schalen-Brücke = **Analogie** / **Forschungsfrage** (kein Theorem).

---

## 10. Drei Analogien

| Analogie | Physik/Mathematik | EABC-Objekt | Eignung |
|----------|-------------------|-------------|---------|
| **Bell** | Korrelationen $E(a,b)$ zwischen Messachsen | $P_{\mathrm{same}}$, CHSH auf $G_E$ | sekundär — erzwingt Kontextwechsel |
| **Holonomie** | Paralleltransport, $\Omega_{\mathrm{Hol}}$ | $\chi_{\mathrm{Hol}}$, $\mathrm{Hol}_E$ | korrekte Stufe, ohne explizite $\gamma^\pm$-Gegenüberstellung |
| **Sagnac** | $\gamma^+$ vs. $\gamma^-$, $\Delta T$ | ABCEA vs. CEABC, $\Delta_E$, $C_E$ | **beste Passung** — gleicher Kreis, entgegengesetzte Orientierung |

$$\boxed{\;\text{Sagnac ist die beste Metapher: gleicher Zyklus, entgegengesetzte Orientierung, Fehlerterm in } \Delta_E.\;}$$

**Label:** Analogietabelle = **Analogie**.

---

## 11. Python-Symbolzuordnung

| LaTeX | Python (kanonisch) | Legacy-Alias |
|-------|-------------------|--------------|
| $N_+(X)$ | `N_plus` | `N_ABCEA` |
| $N_-(X)$ | `N_minus` | `N_CEABC` |
| $\Delta_E(X)$ | `Delta_E` | `D_E` |
| $S_E(X)$ | `S_E` | `chi_Hol` |
| $\widetilde{D}_E(X)$ | `D_tilde_E` | — |
| $C_E(X)$ | `C_E` | — |
| $\omega(e)$ | `edge_omega` | — |
| $A(i\!\to\!j)$ | `discrete_one_form` | — |
| $\Omega(C)$ | `cycle_omega_graph` | — |
| Sagnac-Report | `sagnac_report()` | — |
| Zirkulations-Report | `circulation_report()` | — |

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py` — `sagnac_report(X)`; `collatz_eabc_sagnac_circulation.py` — `circulation_C_E(X)`.

---

## 12. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$, $C_E$ | **Definition** |
| $G_E$, $\omega(e)$, $\Omega(C)$, diskrete 1-Form $A$ | **Definition** |
| ABCEA vs. CEABC als $\gamma^\pm$ | **Definition** (Sagnac-Analogie) |
| $N_+\sim N_-$ $\Rightarrow$ $S_E\to 0$ | **Vermutung** (Hauptvermutung) |
| $\Delta_E$, $C_E$ mit $L$-Funktions-Nullstellen mod $12$ | **Hypothese** (Fehlerterm-Hypothese) |
| Verhalten von $\widetilde{D}_E$ | **Experiment** |
| Bell / Holonomie / Sagnac | **Analogie** (Sagnac primär) |
| Normschalen-Brücke $\Sigma_n\to\mu_n\to G_n$ | **Analogie** / **Forschungsfrage** |
| Physikalische Sagnac-Rotation | **nicht behauptet** |

---

*Kanonsiche Lesart: Der EABC-Transport liefert zwei gegenläufige geschlossene Wege $\gamma^\pm$ auf demselben 4-Zyklus. Ihre Zähl-Differenz $\Delta_E$ ist der **Sagnac-Fehlerterm**; die normalisierte Observable $S_E$ misst Orientierungsasymmetrie analog $\Delta T/T$. Bell/CHSH bleibt als sekundäre Konsistenz-Metapher dokumentiert.*
