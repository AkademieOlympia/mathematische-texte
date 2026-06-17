# EABC-Transport und gerichteter Übergangsgraph

**Status:** Definition + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Theorem | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_zyklus_holonomie.md` — **kanonisch:** $\chi_{\mathrm{path}}$ (4-Block), $\chi_{\mathrm{hol}}^{(5)}$ (5-Block), $\mathrm{Hol}_E$
- `collatz_eabc_holonomie.md` — Vierlings-Orientierung $\omega(Q)$, $\chi_E^{\mathrm{quad}}(N)$, projektive Holonomie $\mathcal H_E$
- `collatz_eabc_holonomie_test.py` — Vierlings-$\chi_E^{\mathrm{quad}}$, Lean-Chiralität
- `collatz_eabc_transition_graph.py` / `.json` — Numerik (Übergangsmatrix, $\chi_{\mathrm{path}}$, $\chi_{\mathrm{hol}}^{(5)}$, Nullmodelle)
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
| 4 | 4-Block-Pfad | $Q_n^{(4)}$, Wort ABCE/CEAB — **orientierter Transportpfad** |
| 5 | Pfad-Chiralität | $\chi_{\mathrm{path}}(N)$, $\Omega_{\mathrm{path}}\in\{+1,-1,0\}$ |
| 6 | 5-Block-Zyklus | $Q_n^{(5)}$, ABCEA/CEABC — **geschlossene Holonomie** |
| 7 | Zyklus-Chiralität / Grenzwert | $\chi_{\mathrm{hol}}^{(5)}(N)$, $\mathrm{Hol}_E=\lim\chi_{\mathrm{hol}}^{(5)}$ |

Stufen 4–5 (Pfad) und 6–7 (Holonomie) sind in `collatz_eabc_zyklus_holonomie.md` kanonisch; dieses Dokument hebt die **Transportebene** (Stufen 3–4) auf dieselbe Hierarchie.

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

## 5. Pfadorientierung auf 4-Schritten (kein geschlossener Zyklus)

$$\boxed{\;\text{Viererblock} = \text{orientierter Transportpfad.}\;}$$

### 5.1 Wort-Pfade (ABCE / CEAB)

**Definition (4-Fenster).** Für vier aufeinanderfolgende Prim-Klassen $(c_0,c_1,c_2,c_3)$:
$$W^{(4)}(c_0,c_1,c_2,c_3):=c_0c_1c_2c_3\in\{E,A,B,C\}^4.$$

**Definition (Pfadorientierung).**
$$\chi_{\mathrm{path}}(W^{(4)})=+1\Leftrightarrow W^{(4)}=\mathrm{ABCE},\qquad
\chi_{\mathrm{path}}(W^{(4)})=-1\Leftrightarrow W^{(4)}=\mathrm{CEAB}.$$

Dies ist **EABC-Pfadorientierung**, **nicht** Holonomie im strengen Sinn — der Pfad schließt nicht.

### 5.2 $t$-Zyklus vs. $t^{-1}$-Zyklus (kanonische Rotation)

**Vorwärts-$t$-Pfad:** $c_{k+1}=t(c_k)$ für $k=0,1,2$ (drei Schritte, vier Knoten).

**Rückwärts-$t^{-1}$-Pfad:** $c_{k+1}=t^{-1}(c_k)$.

**Label:** $\chi_{\mathrm{path}}$, $t$-Pfade = **Definition**.

---

## 6. Zyklus-Holonomie auf 5-Schritten

$$\boxed{\;\text{Fünferblock} = \text{geschlossener Zyklus / Holonomie.}\;}$$

Geschlossener Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$:
$$\Omega^{(5)}(W^{(5)})=+1\Leftrightarrow W^{(5)}=\mathrm{ABCEA},\qquad
\Omega^{(5)}(W^{(5)})=-1\Leftrightarrow W^{(5)}=\mathrm{CEABC}.$$

**Label:** $\Omega^{(5)}$, $\chi_{\mathrm{hol}}^{(5)}$ = **Definition** (`collatz_eabc_zyklus_holonomie.md` §5).

---

## 7. Pfad-Chiralität $\chi_{\mathrm{path}}(N)$ und Holonomie $\chi_{\mathrm{hol}}^{(5)}(N)$

**Kanonsiche Definition:** `collatz_eabc_zyklus_holonomie.md` §4–5.

$$\chi_{\mathrm{path}}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+3}\le N}\chi_{\mathrm{path}}(Q_n^{(4)})}
{\displaystyle\#\{n:\,p_{n+3}\le N,\;\chi_{\mathrm{path}}(Q_n^{(4)})\neq 0\}}
\in[-1,1],$$

$$\chi_{\mathrm{hol}}^{(5)}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+4}\le N}\Omega^{(5)}(Q_n^{(5)})}
{\displaystyle\#\{n:\,p_{n+4}\le N,\;\Omega^{(5)}(Q_n^{(5)})\neq 0\}}
\in[-1,1].$$

**Legacy-Notation:** $\chi_{\mathrm{trans}}(N):=\chi_{\mathrm{path}}(N)$; $\chi_E(N):=\chi_{\mathrm{path}}(N)$ (ältere Bezeichnung).

**Vergleichsträger (Vierlinge):**
$$\chi_E^{\mathrm{quad}}(N)=\frac{\#\mathrm{ABCE}-\#\mathrm{CEAB}}{\#\mathrm{ABCE}+\#\mathrm{CEAB}}$$

auf **Prim-Vierlingen** $Q(p)=(p,p{+}2,p{+}6,p{+}8)$ (`collatz_eabc_holonomie.md` §4).

**Experiment:** `collatz_eabc_transition_graph.py::chi_path_sliding`, `chi_hol_sliding`, `chi_path_vs_hol`.

---

## 8. Boxed Forschungsfrage

$$\boxed{\;\text{Besitzt der gerichtete EABC-Übergangsgraph eine nichttriviale 5-Block-Zyklus-Holonomie }\mathrm{Hol}_E\neq 0\text{?}\;}$$

**Präzisierung (ehrlich).**

| Observable | Block | Träger | Was zählt? |
|------------|------:|--------|------------|
| $\chi_{\mathrm{path}}(N)$ | 4 | aufeinanderfolgende Prim-Klassen | ABCE/CEAB — **Pfadorientierung** |
| $\chi_{\mathrm{hol}}^{(5)}(N)$ | 5 | aufeinanderfolgende Prim-Klassen | ABCEA/CEABC — **geschlossene Holonomie** |
| $\chi_E^{\mathrm{quad}}(N)$ | 4 | Prim-Vierlinge $Q(p)$ | ABCE/CEAB auf **arithmetisch gekoppelten** Beinen |
| $\chi_{\mathrm{t\text{-}cycle}}(N)$ | 4 | $t$-bzw. $t^{-1}$-Pfade | kanonische Rotation vs. Inversion |

**Theorem ($V_4$).** Algebraische Produkte auf $V_4$ sind trivial assoziativ (`collatz_eabc_holonomie.md` §1). Holonomie ist **projektiv / transportiert**, nicht $V_4$-Klammertheorie.

**Heuristik.** Signifikante Abweichung von $\chi_{\mathrm{hol}}^{(5)}$ gegenüber einem **Marginal-Nullmodell** wäre ein Hinweis auf **nichttriviale Zyklus-Holonomie** — kein Beweis. $\chi_{\mathrm{path}}$ dient als verwandte Pfad-Observable, misst aber **keine** geschlossene Holonomie.

**Experiment:** Shuffle-Null in `collatz_eabc_transition_graph.py`.

---

## 9. Stationäre Verteilung und Ergodizität

**Definition.** Eine Verteilung $\pi$ auf $V_4$ ist **stationär**, wenn $\pi T=\pi$ (Zeilenvektor, $T$ zeilenstochastisch normalisiert).

**Definition (ergodisch).** Die Kette ist **ergodisch**, wenn sie irreduzibel und aperiodisch ist; dann ist $\pi$ eindeutig.

**Experiment:** `stationary_distribution`, Irreduzibilitäts-Check in `collatz_eabc_transition_graph.py`.

**Label:** stationäre Verteilung = **Definition**; empirische Ergodizität = **Experiment**.

---

## 10. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $\tau(p_n)=(\kappa(p_n),\kappa(p_{n+1}))$ | **Definition** |
| $T_{ij}$, $f_{ij}$ | **Definition** |
| $t$-Rotation $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ | **Definition** |
| $\chi_{\mathrm{path}}$ — 4-Block-Pfad, **keine** Holonomie | **Definition** |
| $\chi_{\mathrm{hol}}^{(5)}$, $\mathrm{Hol}_E$ — 5-Block-Zyklus | **Definition** / **Hypothese** |
| $\chi_{\mathrm{path}}$ vs.\ $\chi_{\mathrm{hol}}^{(5)}$ vs.\ $\chi_E^{\mathrm{quad}}$ | **Experiment** |
| ABCE-Dominanz auf Vierlingen = ABCEA-Dominanz auf Primübergängen | **Forschungsfrage** |
| Shuffle-/Isotropie-Null Signifikanz | **Experiment** |

---

*Kanonsiche Notiz: Der Übergangsgraph operationalisiert **Transport** entlang der Primzahlordnung. Viererblöcke messen **Pfade**; Fünferblöcke schließen den Zyklus. Vierlings-$\chi_E^{\mathrm{quad}}$ misst **lokale arithmetische Pakete** — die ehrliche Gegenüberstellung ist Kern von PR #54.*
