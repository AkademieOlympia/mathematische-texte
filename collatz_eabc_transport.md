# EABC-Transport und gerichteter Übergangsgraph

**Status:** Definition + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Theorem | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_zyklus_holonomie.md` — **kanonisch:** $\chi_E(N)$ auf Primfolge-Gleitfenstern, $\mathrm{Hol}_E$, Hypothese
- `collatz_eabc_holonomie.md` — Vierlings-Orientierung $\omega(Q)$, $\chi_E^{\mathrm{quad}}(N)$, projektive Holonomie $\mathcal H_E$
- `collatz_eabc_holonomie_test.py` — Vierlings-$\chi_E^{\mathrm{quad}}$, Lean-Chiralität
- `collatz_eabc_transition_graph.py` / `.json` — Numerik (Übergangsmatrix, $\chi_E$, Nullmodelle, $\mathrm{Hol}_E$-Schätzung)
- `eabc_from_lean.py` — $\kappa=\texttt{class\_of}$, Rotation $t\colon E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$
- `collatz_eabc_invarianzprogramm.md` — globale Observable $\chi(x)$

---

## 1. Zentrale Frage: natürlicher Transport $T$

**Forschungsfrage.** Welcher **natürliche Transport** verbindet EABC-Signaturen entlang der Primzahlfolge?

**Kandidat $T$.** Für Primzahlen $p_n$ in aufsteigender Reihenfolge ($p_n>3$):
$$\kappa(p_n)\in\{E,A,B,C\},\qquad
\tau(p_n):=\bigl(\kappa(p_n),\,\kappa(p_{n+1})\bigr).$$

Der **fundamentale Zustand** ist die **gerichtete Kante** $\tau$, nicht der isolierte Punkt $\kappa(p)$.

**Label:** Kandidat-Transport = **Definition**; Existenz einer kanonischen Dynamik = **Forschungsfrage**.

---

## 2. Hierarchie

| Stufe | Objekt | Symbol |
|------:|--------|--------|
| 1 | Primzahl | $p_n$ |
| 2 | EABC-Klasse | $\kappa(p_n)$ |
| 3 | Übergang | $\tau(p_n)=(\kappa(p_n),\kappa(p_{n+1}))$ |
| 4 | Zyklus (4-Schritt) | Wort $W\in\{\mathrm{ABCE},\mathrm{CEAB}\}$ oder $t$-Pfad |
| 5 | Zyklus-Holonomie | $\omega(W)\in\{+1,-1\}$ |
| 6 | Grenzwert | $\mathcal H_E=\lim_{N\to\infty}\chi_E(N)$ bzw. Transport-Analogon |

Stufen 5–6 sind in `collatz_eabc_holonomie.md` für **Prim-Vierlinge** kanonisch; dieses Dokument hebt die **Transportebene** (Stufen 3–4) auf dieselbe Hierarchie.

---

## 3. Gerichteter EABC-Übergangsgraph

**Definition (Zustandsraum).** $V_4=\{E,A,B,C\}$.

**Definition (gerichtete Kanten).** Für Primzahlen $p_n,p_{n+1}>3$:
$$i\to j \text{ gezählt, wenn }\kappa(p_n)=i,\;\kappa(p_{n+1})=j.$$

**Definition (Übergangsmatrix).**
$$T_{ij}(N):=\#\{n:\,\kappa(p_n)=i,\,\kappa(p_{n+1})=j,\,p_{n+1}\le N\}.$$

**Definition (Kantenfrequenz).**
$$f_{ij}(N):=\frac{T_{ij}(N)}{\sum_{k,\ell}T_{k\ell}(N)}.$$

**Experiment:** `collatz_eabc_transition_graph.py::transition_matrix`.

**Label:** $T_{ij}$, $f_{ij}$ = **Definition**; empirische Werte = **Experiment**.

---

## 4. Kanonische Rotation $t$ und Kantendynamik

Aus `eabc_from_lean.py`:
$$t\colon E\mapsto A,\; A\mapsto B,\; B\mapsto C,\; C\mapsto E,\qquad t^4=\mathrm{id}.$$

**Definition (vorwärts-kanonische Kante).** $i\to j$ ist **$t$-aligniert**, wenn $j=t(i)$.

**Dynamik auf Kanten:** $E\!\to\!A$, $A\!\to\!B$, $B\!\to\!C$, $C\!\to\!E$.

**Label:** $t$ = **Definition** (Lean-konsistent).

---

## 5. Zyklus-Holonomie auf 4-Schritten

### 5.1 Wort-Zyklen (Vierlings-Chiralität)

**Definition (Wortfenster).** Für vier aufeinanderfolgende Prim-Klassen $(c_0,c_1,c_2,c_3)$:
$$W(c_0,c_1,c_2,c_3):=c_0c_1c_2c_3\in\{E,A,B,C\}^4.$$

**Definition (Orientierung).**
$$\omega(W)=+1\Leftrightarrow W=\mathrm{ABCE},\qquad
\omega(W)=-1\Leftrightarrow W=\mathrm{CEAB}.$$

### 5.2 $t$-Zyklus vs. $t^{-1}$-Zyklus

**Vorwärtszyklus (ABCE-Phase):** vier $t$-alignierte Schritte:
$$c_{k+1}=t(c_k)\quad (k=0,1,2,3),\qquad \Omega=+1.$$

**Rückwärtszyklus (CEAB-Phase):** $c_{k+1}=t^{-1}(c_k)$:
$$\Omega=-1.$$

**Label:** $\omega$, $\Omega$ = **Definition**.

---

## 6. Transport-Chiralität $\chi_E(N)$ (Gleitfenster)

**Kanonsiche Definition:** `collatz_eabc_zyklus_holonomie.md` §4.

$$\chi_E(N):=
\frac{\displaystyle\sum_{n:\,p_{n+3}\le N}\Omega(Q_n)}
{\displaystyle\#\{n:\,p_{n+3}\le N,\;\Omega(Q_n)\neq 0\}}
\in[-1,1],$$

mit $Q_n=(X_n,X_{n+1},X_{n+2},X_{n+3})$ auf der Primfolge und $\Omega=+1$ (ABCE), $-1$ (CEAB), $0$ sonst.

**Legacy-Notation:** $\chi_{\mathrm{trans}}(N):=\chi_E(N)$ (dieses Dokument, frühere Version).

**Vergleichsträger (Vierlinge):**
$$\chi_E^{\mathrm{quad}}(N)=\frac{\#\mathrm{ABCE}-\#\mathrm{CEAB}}{\#\mathrm{ABCE}+\#\mathrm{CEAB}}$$

auf **Prim-Vierlingen** $Q(p)=(p,p{+}2,p{+}6,p{+}8)$ (`collatz_eabc_holonomie.md` §4).

**Experiment:** `collatz_eabc_transition_graph.py::chi_E_sliding`, `chi_sliding_vs_quadruplet`.

---

## 7. Boxed Forschungsfrage

$$\boxed{\;\text{Besitzt der gerichtete EABC-Übergangsgraph eine nichttriviale Zyklus-Holonomie?}\;}$$

**Präzisierung (ehrlich).**

| Observable | Träger | Was zählt? |
|------------|--------|------------|
| $\chi_E(N)$ | Prim-Vierlinge $Q(p)$ | ABCE/CEAB auf vier **arithmetisch gekoppelten** Beinen |
| $\chi_{\mathrm{trans}}(N)$ | aufeinanderfolgende Prim-Klassen | ABCE/CEAB auf **Primfolgen-Fenstern** |
| $\chi_{\mathrm{t\text{-}cycle}}(N)$ | $t$-bzw. $t^{-1}$-Pfade | kanonische Rotation vs. Inversion |

**Theorem ($V_4$).** Algebraische Produkte auf $V_4$ sind trivial assoziativ (`collatz_eabc_holonomie.md` §1). Holonomie ist **projektiv / transportiert**, nicht $V_4$-Klammertheorie.

**Heuristik.** Signifikante Abweichung von $\chi_{\mathrm{trans}}$ gegenüber einem **Marginal-Nullmodell** (Permutation der Klassenfolge bei erhaltener Häufigkeitsverteilung) wäre ein Hinweis auf **nichttriviale Transport-Holonomie** — kein Beweis.

**Experiment:** Shuffle-Null in `collatz_eabc_transition_graph.py`.

---

## 8. Stationäre Verteilung und Ergodizität

**Definition.** Eine Verteilung $\pi$ auf $V_4$ ist **stationär**, wenn $\pi T=\pi$ (Zeilenvektor, $T$ zeilenstochastisch normalisiert).

**Definition (ergodisch).** Die Kette ist **ergodisch**, wenn sie irreduzibel und aperiodisch ist; dann ist $\pi$ eindeutig.

**Experiment:** `stationary_distribution`, Irreduzibilitäts-Check in `collatz_eabc_transition_graph.py`.

**Label:** stationäre Verteilung = **Definition**; empirische Ergodizität = **Experiment**.

---

## 9. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $\tau(p_n)=(\kappa(p_n),\kappa(p_{n+1}))$ | **Definition** |
| $T_{ij}$, $f_{ij}$ | **Definition** |
| $t$-Rotation $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ | **Definition** |
| $\chi_{\mathrm{trans}}$ vs. $\chi_E$ — verwandt, **nicht identisch** | **Experiment** |
| ABCE-Dominanz auf Vierlingen = ABCE-Dominanz auf Primübergängen | **Forschungsfrage** |
| Shuffle-Null Signifikanz | **Experiment** |
| $\mathcal H_E=\lim\chi_E(N)$ | **Definition** / Grenzwert offen |

---

*Kanonsiche Notiz: Der Übergangsgraph operationalisiert **Transport** entlang der Primzahlordnung. Vierlings-$\chi_E$ misst **lokale arithmetische Pakete** — die ehrliche Gegenüberstellung ist Kern von PR #54 Stufe 5–6.*
