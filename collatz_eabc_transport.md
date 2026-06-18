# EABC-Transport und gerichteter Übergangsgraph

**Status:** Definition + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Theorem | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_zyklus_holonomie.md` — **kanonisch:** Klasse→Kante→Pfad→Zyklus→Holonomie; $\chi_{\mathrm{Pfad}}$, $\chi_{\mathrm{Hol}}$, $\mathrm{Hol}_E$
- `collatz_eabc_holonomie.md` — Vierlings-Orientierung $\omega(Q)$, $\chi_E^{\mathrm{quad}}(N)$, projektive Holonomie $\mathcal H_E$
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch: $\mathrm{Hol}_E=0$, Fehlerterm $D_E$
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — Numerik $D_E$, $\widetilde{D}_E$, Chebyshev-Vergleich
- `collatz_eabc_holonomie_test.py` — Vierlings-$\chi_E^{\mathrm{quad}}$, Lean-Chiralität
- `collatz_eabc_transition_graph.py` / `.json` — Numerik (Übergangsmatrix, $\chi_{\mathrm{Pfad}}$, $\chi_{\mathrm{Hol}}$, Nullmodelle)
- `eabc_from_lean.py` — $\kappa=\texttt{class\_of}$, Rotation $t\colon E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$
- `collatz_eabc_invarianzprogramm.md` — globale Observable $\chi(x)$

---

## 1. Zentrale Frage: natürlicher Transport $T$

**Forschungsfrage.** Welcher **natürliche Transport** verbindet EABC-Signaturen entlang der Primzahlfolge?

**Kandidat $T$.** Für Primzahlen $p_n$ in aufsteigender Reihenfolge ($p_n>3$):
$$\kappa(p_n)\in\{E,A,B,C\},\qquad
\tau_n:=\bigl(X_n,\,X_{n+1}\bigr)=\bigl(\kappa(p_n),\,\kappa(p_{n+1})\bigr),$$
mit $X_n:=\kappa(p_n)$.

Der **fundamentale Zustand** ist die **gerichtete Kante** $\tau_n$, nicht der isolierte Punkt $X_n$.

**Label:** Kandidat-Transport = **Definition**; Existenz einer kanonischen Dynamik = **Forschungsfrage**.

---

## 2. Hierarchie

$$\boxed{\;\text{Klasse} \;\to\; \text{Kante} \;\to\; \text{Pfad} \;\to\; \text{Zyklus} \;\to\; \text{Holonomie}\;}$$

| Stufe | Objekt | Symbol |
|------:|--------|--------|
| 1 | **Klasse** | $X_n=\kappa(p_n)$ |
| 2 | **Kante** | $\tau_n=(X_n,X_{n+1})$ |
| 3 | **Pfad** | $P_n^{(4)}=(X_n,X_{n+1},X_{n+2},X_{n+3})$; $\Omega_{\mathrm{Pfad}}$ (ABCE/CEAB) |
| 4 | Pfad-Chiralität | $\chi_{\mathrm{Pfad}}(N)$ |
| 5 | **Zyklus** | $C_n^{(5)}$; ABCEA/CEABC |
| 6 | **Holonomie** | $\Omega_{\mathrm{Hol}}(C_n^{(5)})$, $\chi_{\mathrm{Hol}}(N)$, $\mathrm{Hol}_E=\lim\chi_{\mathrm{Hol}}$ |

Stufen 3–4 (Pfad) und 5–6 (Zyklus/Holonomie) sind in `collatz_eabc_zyklus_holonomie.md` kanonisch; dieses Dokument hebt die **Transportebene** (Stufen 2–3) auf dieselbe Hierarchie.

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

$$\boxed{\;\text{Pfad} = \text{orientierter Transportpfad.}\;}$$

### 5.1 Wort-Pfade (ABCE / CEAB)

**Definition (4-Fenster / Pfad).**
$$P_n^{(4)}=(X_n,X_{n+1},X_{n+2},X_{n+3}),\qquad
W^{(4)}:=X_nX_{n+1}X_{n+2}X_{n+3}\in\{E,A,B,C\}^4.$$

**Definition (Pfadorientierung).**
$$\Omega_{\mathrm{Pfad}}(P_n^{(4)})=+1\Leftrightarrow W^{(4)}=\mathrm{ABCE},\qquad
\Omega_{\mathrm{Pfad}}(P_n^{(4)})=-1\Leftrightarrow W^{(4)}=\mathrm{CEAB}.$$

Dies ist **EABC-Pfadorientierung**, **nicht** Holonomie im strengen Sinn — der Pfad schließt nicht.

### 5.2 $t$-Zyklus vs. $t^{-1}$-Zyklus (kanonische Rotation)

**Vorwärts-$t$-Pfad:** $X_{k+1}=t(X_k)$ für $k=0,1,2$ (drei Schritte, vier Knoten).

**Rückwärts-$t^{-1}$-Pfad:** $X_{k+1}=t^{-1}(X_k)$.

**Label:** $\Omega_{\mathrm{Pfad}}$, $t$-Pfade = **Definition**.

---

## 6. Zyklus-Holonomie auf 5-Schritten

$$\boxed{\;\text{Zyklus} = \text{geschlossener Zyklus; Holonomie beginnt hier.}\;}$$

Geschlossener Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$:
$$\Omega_{\mathrm{Hol}}(C_n^{(5)})=+1\Leftrightarrow C_n^{(5)}=\mathrm{ABCEA},\qquad
\Omega_{\mathrm{Hol}}(C_n^{(5)})=-1\Leftrightarrow C_n^{(5)}=\mathrm{CEABC}.$$

**Label:** $\Omega_{\mathrm{Hol}}$, $\chi_{\mathrm{Hol}}$ = **Definition** (`collatz_eabc_zyklus_holonomie.md` §5).

---

## 7. Pfad-Chiralität $\chi_{\mathrm{Pfad}}(N)$ und Holonomie $\chi_{\mathrm{Hol}}(N)$

**Kanonsiche Definition:** `collatz_eabc_zyklus_holonomie.md` §4–5.

$$\chi_{\mathrm{Pfad}}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+3}\le N}\Omega_{\mathrm{Pfad}}(P_n^{(4)})}
{\displaystyle\#\{n:\,p_{n+3}\le N,\;\Omega_{\mathrm{Pfad}}(P_n^{(4)})\neq 0\}}
\in[-1,1],$$

$$\chi_{\mathrm{Hol}}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+4}\le N}\Omega_{\mathrm{Hol}}(C_n^{(5)})}
{\displaystyle\#\{n:\,p_{n+4}\le N,\;\Omega_{\mathrm{Hol}}(C_n^{(5)})\neq 0\}}
\in[-1,1].$$

**Legacy-Notation:** $\chi_{\mathrm{path}}:=\chi_{\mathrm{Pfad}}$; $\chi_{\mathrm{hol}}^{(5)}:=\chi_{\mathrm{Hol}}$; $\chi_{\mathrm{trans}}:=\chi_{\mathrm{Pfad}}$.

**Vergleichsträger (Vierlinge):**
$$\chi_E^{\mathrm{quad}}(N)=\frac{\#\mathrm{ABCE}-\#\mathrm{CEAB}}{\#\mathrm{ABCE}+\#\mathrm{CEAB}}$$

auf **Prim-Vierlingen** $Q(p)=(p,p{+}2,p{+}6,p{+}8)$ (`collatz_eabc_holonomie.md` §4).

**Experiment:** `collatz_eabc_transition_graph.py::chi_pfad_sliding`, `chi_hol_sliding`, `chi_pfad_vs_hol`.

---

## 8. Boxed Forschungsfrage

$$\boxed{\;\text{Besitzt der gerichtete EABC-Übergangsgraph eine nichttriviale Zyklus-Holonomie }\mathrm{Hol}_E\neq 0\text{?}\;}$$

**Präzisierung (ehrlich).**

| Observable | Block | Träger | Was zählt? |
|------------|------:|--------|------------|
| $\chi_{\mathrm{Pfad}}(N)$ | 4 | aufeinanderfolgende Prim-Klassen $P_n^{(4)}$ | ABCE/CEAB — **Pfadorientierung** |
| $\chi_{\mathrm{Hol}}(N)$ | 5 | aufeinanderfolgende Prim-Klassen $C_n^{(5)}$ | ABCEA/CEABC — **geschlossene Holonomie** |
| $\chi_E^{\mathrm{quad}}(N)$ | 4 | Prim-Vierlinge $Q(p)$ | ABCE/CEAB auf **arithmetisch gekoppelten** Beinen |
| $\chi_{\mathrm{t\text{-}cycle}}(N)$ | 4 | $t$-bzw. $t^{-1}$-Pfade | kanonische Rotation vs. Inversion |

**Theorem ($V_4$).** Algebraische Produkte auf $V_4$ sind trivial assoziativ (`collatz_eabc_holonomie.md` §1). Holonomie ist **projektiv / transportiert**, nicht $V_4$-Klammertheorie.

**Heuristik.** Signifikante Abweichung von $\chi_{\mathrm{Hol}}$ gegenüber einem **Marginal-Nullmodell** wäre ein Hinweis auf **nichttriviale Zyklus-Holonomie** — kein Beweis. $\chi_{\mathrm{Pfad}}$ dient als verwandte Pfad-Observable, misst aber **keine** geschlossene Holonomie.

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
| $\tau_n=(X_n,X_{n+1})$ — **Kante** | **Definition** |
| $T_{ij}$, $f_{ij}$ | **Definition** |
| $t$-Rotation $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ | **Definition** |
| $P_n^{(4)}$, $\Omega_{\mathrm{Pfad}}$ — Pfad, **keine** Holonomie | **Definition** |
| $C_n^{(5)}$, $\Omega_{\mathrm{Hol}}$, $\chi_{\mathrm{Hol}}$, $\mathrm{Hol}_E$ | **Definition** / **Hypothese** |
| $\chi_{\mathrm{Pfad}}$ vs.\ $\chi_{\mathrm{Hol}}$ vs.\ $\chi_E^{\mathrm{quad}}$ | **Experiment** |
| ABCE-Dominanz auf Vierlingen = ABCEA-Dominanz auf Primübergängen | **Forschungsfrage** |
| Shuffle-/Isotropie-Null Signifikanz | **Experiment** |

---

*Kanonsiche Notiz: Der Übergangsgraph operationalisiert **Transport** entlang der Primzahlordnung. Pfade ($P_n^{(4)}$) messen **gerichtete 4-Ketten**; Zyklen ($C_n^{(5)}$) schließen den Zyklus und tragen Holonomie. Vierlings-$\chi_E^{\mathrm{quad}}$ misst **lokale arithmetische Pakete** — die ehrliche Gegenüberstellung ist Kern von PR #54.*
