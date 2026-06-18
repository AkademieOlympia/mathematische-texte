# EABC-Sagnac-Observable (didaktischer Einstieg)

**Status:** Definition + Analogie + Verweis auf mathematischen Kern  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

**Mathematischer Kern:** `collatz_eabc_zirkulation_spektral.md` — **Korrelation → Zirkulation** (reine Mathematik: diskrete 1-Form, $H_1(C_4,\mathbb{Z})$, $\mathrm{Spec}(L_E)$). Dieses Dokument ist ein **didaktischer Einstieg** über gegenläufige geschlossene Wege; es behauptet **keine** physikalische Rotation oder Relativitätstheorie.

**Querverweise:**
- `collatz_eabc_zirkulation_spektral.md` — **kanonische reine Mathematik:** $E(a,b)$, $D_E$, $C_E$, $L_E$, Normschalen
- `collatz_eabc_fehlerterm_hypothese.md` — **kanonische Endform:** $N_\pm$, $\Delta_E=D_E$, $S_E$, Hauptvermutung, Fehlerterm-Hypothese
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_transport.md` — $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH (Kantenkorrelation $E(a,b)$)
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — Numerik $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$
- `collatz_eabc_sagnac_circulation.py` / `.json` — Zirkulation $C_E(X)$, Kantenorientierung $\omega$, diskrete 1-Form $\alpha$
- `collatz_eabc_graph_laplacian.py` — $\mathrm{Spec}(L_E)$, Spektrallücke
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch: $\mathrm{Hol}_E=0$, Fehlerterm $\Delta_E$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 0. Didaktik: gegenläufige Zyklen

Auf demselben $C_4$-Gerüst $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ betrachten wir zwei Orientierungen:
- **$\gamma^+$ (ABCEA):** eine Lesart des geschlossenen 5-Fensters
- **$\gamma^-$ (CEABC):** die zyklisch verschobene Gegenorientierung

Die **Sagnac-Lesart** benennt diese Gegenüberstellung bildhaft (gegenläufige Wege auf demselben Kreis). Mathematisch ist sie identisch mit der **Zirkulation** $C_E(X)$ in `collatz_eabc_zirkulation_spektral.md` §4.

$$\boxed{\;\text{Didaktik: Sagnac = gegenläufige } \gamma^\pm \text{; Kern: } C_E = D_E = \sum_\gamma \mathrm{sgn}(\gamma).\;}$$

| Lesart | Objekt | Dokument |
|--------|--------|----------|
| **Kern (rein mathematisch)** | Korrelation $E(a,b)$ → Zirkulation $C_E$ → $\mathrm{Spec}(L_E)$ | `collatz_eabc_zirkulation_spektral.md` |
| **Didaktik (dieses Dokument)** | $\gamma^+$ vs. $\gamma^-$, $\Delta_E$ | hier |
| **Sekundär** | Bell/CHSH auf $G_E$ | `collatz_eabc_bell_holonomie.md` |

---

## 1. Definitionen

**Restklassen modulo $12$:**
$$E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11\pmod{12}.$$

**Primfolge-Labels:**
$$X_n := \kappa(p_n),\qquad C_n^{(5)} := (X_n,X_{n+1},X_{n+2},X_{n+3},X_{n+4}).$$

**Zählgrößen** (für Prim-Obergrenze $X$):
$$N_+(X) := N_{\mathrm{ABCEA}}(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{ABCEA}\},$$
$$N_-(X) := N_{\mathrm{CEABC}}(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{CEABC}\}.$$

**Orientierte Wege:**
- **$\gamma^+$:** ABCEA
- **$\gamma^-$:** CEABC

**Zirkulations-Differenz (Fehlerterm):**
$$\Delta_E(X) := N_+(X) - N_-(X) = D_E(X) = C_E(X).$$

**Normalisierte Observable:**
$$S_E(X) := \frac{N_+(X) - N_-(X)}{N_+(X) + N_-(X)} = \frac{\Delta_E(X)}{N_+(X)+N_-(X)}.$$

**Legacy-Aliase:** $\chi_{\mathrm{Hol}}(X):=S_E(X)$; $\mathrm{Hol}_E:=\lim_{X\to\infty} S_E(X)$.

**Normalisierter Fehlerterm:**
$$\widetilde{D}_E(X) := \frac{\Delta_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

**Label:** $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$ = **Definition**.

---

## 2. Hauptvermutung und Fehlerterm

**Hauptvermutung (Hauptterm).**
$$N_+(X) \sim N_-(X)\qquad (X\to\infty),$$
und damit $\mathrm{Hol}_E = \lim_{X\to\infty} S_E(X) = 0$.

**Fehlerterm-Hypothese (stärker).** $\Delta_E$ trägt Chebyshev-artigen Bias (Nullstellen der Dirichlet-$L$-Funktionen mod $12$). Details: `collatz_eabc_fehlerterm_hypothese.md`.

**Label:** Hauptvermutung = **Vermutung**; Fehlerterm-Hypothese = **Hypothese**.

---

## 3. Hierarchie (Verweis auf Kern)

$$\boxed{\;\text{Primzahlen} \;\to\; \text{EABC-Klassen} \;\to\; \text{Transport} \;\to\; \text{Zirkulation } C_E \;\to\; D_E(X) \;\to\; \mathrm{Spec}(L_E).\;}$$

Vollständige Tabelle und Spektralgeometrie: `collatz_eabc_zirkulation_spektral.md` §9.

---

## 4. Graphformulierung $\Omega(C)$

**Gerichteter Transportgraph:** $G_E=(V,E)$, $V=\{E,A,B,C\}$.

**Kantenorientierung** auf dem kanonischen 4-Zyklus:
$$\omega(A\!\to\!B)=\omega(B\!\to\!C)=\omega(C\!\to\!E)=\omega(E\!\to\!A)=+1,$$
entgegengesetzte Kanten $-1$, sonst $0$. Lückenmuster mod $12$: $(2,4,2,4)$.

**Zyklus-Produkt:** $\Omega(C) := \prod_i \omega(v_i,v_{i+1})$.

**Label:** $G_E$, $\omega(e)$, $\Omega(C)$ = **Definition**.

---

## 5. Diskrete 1-Form und Zirkulation

$$\alpha(i\!\to\!j) := \frac{\omega(i,j)}{4},\qquad
C_E(X) := \sum_{\gamma\ \mathrm{erkannt}} \oint_\gamma \alpha = D_E(X).$$

**Experiment:** `collatz_eabc_sagnac_circulation.py` — `circulation_C_E(X)`.

**Label:** $\alpha$, $C_E$ = **Definition**.

---

## 6. Boxed Schluss (didaktisch)

$$\boxed{\;\text{ABCEA gegen CEABC: zwei Orientierungen in } H_1(C_4,\mathbb{Z}).\;}$$

$$\boxed{\;\Delta_E(X) = D_E(X) = C_E(X)\;\text{ist der EABC-Fehlerterm.}\;}$$

$$\boxed{\;\mathrm{Hol}_E = 0\;\text{im Hauptterm, aber}\;\Delta_E(X)\;\text{kann nichttrivial sein.}\;}$$

Forschungsfragen (Spektralgeometrie): `collatz_eabc_zirkulation_spektral.md` §8.

---

## 7. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $C_E(X)$, $D_E(X)$ | `C_E`, `D_E` | `collatz_eabc_sagnac_circulation` |
| $\omega(e)$, $\alpha$ | `edge_omega`, `discrete_one_form` | `collatz_eabc_sagnac_circulation` |
| $\mathrm{Spec}(L_E)$ | `eigenvalues_symmetrized` | `collatz_eabc_graph_laplacian` |

---

## 8. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $N_\pm$, $\Delta_E$, $S_E$, $\widetilde{D}_E$, $C_E$ | **Definition** |
| $G_E$, $\omega(e)$, diskrete 1-Form $\alpha$ | **Definition** |
| $\gamma^\pm$ als $H_1(C_4,\mathbb{Z})$-Orientierungen | **Definition** |
| $N_+\sim N_-$ $\Rightarrow$ $S_E\to 0$ | **Vermutung** |
| Fehlerterm mit $L$-Funktionen mod $12$ | **Hypothese** |
| Sagnac-Bild (gegenläufige Wege) | **Analogie** (didaktisch) |
| Mathematischer Kern | `collatz_eabc_zirkulation_spektral.md` |
| Physikalische Rotation | **nicht behauptet** |

---

*Didaktischer Einstieg: Gegenläufige Zyklen $\gamma^\pm$ motivieren $C_E=D_E$. Der mathematische Kern — Korrelation, Zirkulation, Fehlerterm, Spektralgeometrie — steht in `collatz_eabc_zirkulation_spektral.md`.*
