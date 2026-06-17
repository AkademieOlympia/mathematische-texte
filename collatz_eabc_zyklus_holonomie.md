# EABC-Zyklus-Holonomie-Hypothese (kanonisch)

**Status:** Definition (Stufen 1–6) + Hypothese (Stufe 7) + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Hypothese | Experiment | Theorem

**Querverweise:**
- `collatz_eabc_transport.md` — gerichteter Übergangsgraph $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_holonomie.md` — Vierlings-Orientierung $\omega(Q)$, $\chi_E^{\mathrm{quad}}(N)$ (arithmetische Vierlinge)
- `collatz_eabc_transition_graph.py` / `.json` — Numerik: $\chi_{\mathrm{path}}$, $\chi_{\mathrm{hol}}^{(5)}$, Nullmodelle
- `collatz_eabc_holonomie_test.py` — $\chi_E^{\mathrm{quad}}(N)$ auf Prim-Vierlingen (Vergleichsträger)
- `eabc_from_lean.py` — $\kappa=\texttt{class\_of}$, Rotation $t$

---

## 1. EABC-Klassifikation $\kappa(p)$

**Definition.** Für Primzahlen $p>3$:
$$\kappa(p)\in\{E,A,B,C\},\qquad X_n := \kappa(p_n),$$
wobei $(p_n)_{n\ge 1}$ die aufsteigende Folge aller Primzahlen $>3$ ist.

**Label:** $\kappa$, $X_n$ = **Definition** (Stufe 1).

**Implementierung:** `eabc_from_lean.py::class_of`, `EABC.lean`.

---

## 2. Transport $T_n$ und gerichteter Graph $G_E$

**Definition (Transportoperator).** Entlang der Primzahlordnung:
$$T_n\colon X_n \longrightarrow X_{n+1},\qquad T_n := \tau(p_n) = \bigl(\kappa(p_n),\,\kappa(p_{n+1})\bigr).$$

**Definition (Zustandsraum).** $V_4=\{E,A,B,C\}$.

**Definition (gerichteter Übergangsgraph $G_E$).** Kante $i\to j$ wird gezählt, wenn $\kappa(p_n)=i$ und $\kappa(p_{n+1})=j$.

Der **fundamentale Zustand** ist die **gerichtete Kante** $\tau(p_n)$, nicht der isolierte Punkt $\kappa(p_n)$.

**Label:** $T_n$, $G_E$, $\tau$ = **Definition** (Stufe 2–3).

**Experiment:** `collatz_eabc_transition_graph.py::transition_matrix_report`.

---

## 3. Viererblock $Q_n^{(4)}$ — orientierter Transportpfad

$$\boxed{\;\text{Viererblock} = \text{orientierter Transportpfad (kein geschlossener Zyklus).}\;}$$

**Definition (4-Fenster auf der Primfolge).**
$$Q_n^{(4)} := (X_n,\,X_{n+1},\,X_{n+2},\,X_{n+3})
=\bigl(\kappa(p_n),\,\kappa(p_{n+1}),\,\kappa(p_{n+2}),\,\kappa(p_{n+3})\bigr).$$

**Definition (Pfadorientierung).** ABCE und CEAB sind **EABC-Pfadorientierungen** entlang vier aufeinanderfolgender Transportkanten — **nicht** Holonomie im strengen DG-Sinn (der Pfad schließt nicht):
$$\chi_{\mathrm{path}}(Q_n^{(4)}):=
\begin{cases}
+1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}=\mathrm{ABCE},\\[4pt]
-1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}=\mathrm{CEAB},\\[4pt]
0 & \text{sonst.}
\end{cases}$$

**Legacy-Symbol:** $\Omega_{\mathrm{path}}(Q_n^{(4)}) := \chi_{\mathrm{path}}(Q_n^{(4)})$ (frühere Notation $\Omega$ auf 4-Fenstern).

**Wichtig:** $Q_n^{(4)}$ läuft über **aufeinanderfolgende Prim-Klassen** (Gleitfenster), **nicht** über arithmetische Prim-Vierlinge $Q(p)=(p,p{+}2,p{+}6,p{+}8)$.

**Label:** $Q_n^{(4)}$, $\chi_{\mathrm{path}}$ = **Definition** (Stufe 4).

---

## 4. Pfad-Chiralität $\chi_{\mathrm{path}}(N)$

Sei $N$ eine Primzahl-Obergrenze. Fenster $Q_n^{(4)}$ mit $p_{n+3}\le N$ werden berücksichtigt.

**Definition (4-Block-Pfad-Chiralität).**
$$\chi_{\mathrm{path}}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+3}\le N}\chi_{\mathrm{path}}(Q_n^{(4)})}
{\displaystyle\#\{n:\,p_{n+3}\le N,\;\chi_{\mathrm{path}}(Q_n^{(4)})\neq 0\}}
\in[-1,1],$$
mit der Konvention $\chi_{\mathrm{path}}(N)=0$, falls der Nenner $0$ ist.

**Legacy-Alias:** $\chi_E(N):=\chi_{\mathrm{path}}(N)$ (ältere Bezeichnung vor PR #54-Korrektur).

**Experiment:** `collatz_eabc_transition_graph.py::chi_path_sliding`.

**Label:** $\chi_{\mathrm{path}}(N)$ = **Definition** (Stufe 5).

---

## 5. Fünferblock $Q_n^{(5)}$ — geschlossener Zyklus / Holonomie

$$\boxed{\;\text{Fünferblock} = \text{geschlossener Zyklus / Holonomie.}\;}$$

**Definition (5-Fenster auf der Primfolge).**
$$Q_n^{(5)} := (X_n,\,X_{n+1},\,X_{n+2},\,X_{n+3},\,X_{n+4}).$$

Der **geschlossene EABC-Zyklus** ist $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$. Entlang der Primfolge erscheinen die Wörter **ABCEA** (positiv) und **CEABC** (negativ, umgekehrte Zyklusorientierung).

**Definition (Zyklus-Holonomie $\Omega^{(5)}$).**
$$\Omega^{(5)}(Q_n^{(5)})=
\begin{cases}
+1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}X_{n+4}=\mathrm{ABCEA},\\[4pt]
-1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}X_{n+4}=\mathrm{CEABC},\\[4pt]
0 & \text{sonst.}
\end{cases}$$

**Definition (5-Block-Zyklus-Chiralität).**
$$\chi_{\mathrm{hol}}^{(5)}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+4}\le N}\Omega^{(5)}(Q_n^{(5)})}
{\displaystyle\#\{n:\,p_{n+4}\le N,\;\Omega^{(5)}(Q_n^{(5)})\neq 0\}}
\in[-1,1],$$
mit $\chi_{\mathrm{hol}}^{(5)}(N)=0$, falls der Nenner $0$ ist.

**Experiment:** `collatz_eabc_transition_graph.py::chi_hol_sliding`.

**Label:** $Q_n^{(5)}$, $\Omega^{(5)}$, $\chi_{\mathrm{hol}}^{(5)}$ = **Definition** (Stufe 6).

---

## 6. Projektive Zyklus-Holonomie $\mathrm{Hol}_E$

**Definition (Grenzwert).**
$$\boxed{\;\mathrm{Hol}_E := \lim_{N\to\infty}\chi_{\mathrm{hol}}^{(5)}(N)\;}$$
sofern der Grenzwert existiert und $\neq 0$ ist — dann als **projektive Zyklus-Holonomie** des geschlossenen 5-Schritt-Zyklus interpretiert.

**Nicht** $\lim\chi_{\mathrm{path}}(N)$: der 4-Block misst Pfadorientierung, keine geschlossene Holonomie.

Falls $\mathrm{Hol}_E=0$, bleibt die Frage offen, ob $\chi_{\mathrm{hol}}^{(5)}(N)$ **stabile, nicht-zufällige Fluktuationen** trägt (schwächere Lesart).

**Experiment:** `collatz_eabc_transition_graph.py::hol_E_estimates` (Stichproben bei $N=10^5$, $10^6$).

**Label:** $\mathrm{Hol}_E$ = **Definition** (Grenzwert auf 5-Blöcken); Existenz und $\neq 0$ = **Hypothese**.

---

## 7. Boxed Hypothese (stark)

$$\boxed{\;\mathrm{Hol}_E\neq 0\;}$$

**Präzisierung.** Eine nichttriviale asymptotische Orientierungsasymmetrie der geschlossenen ABCEA/CEABC-Zyklen entlang der Primzahl-Transportkette — Hinweis auf **projektive Holonomie** des gerichteten EABC-Übergangsgraphen, kein Beweis im DG-Sinn.

**Schwächere Lesart (falls $\mathrm{Hol}_E=0$):**
$$\boxed{\;\chi_{\mathrm{hol}}^{(5)}(N)\text{ zeigt stabile, nicht-zufällige Fluktuationen gegenüber Isotropie- und Shuffle-Nullmodellen.}\;}$$

**Nullmodelle (Experiment):**
- **Marginal-Shuffle** und **Isotropie-Null** jeweils für $\chi_{\mathrm{path}}$ und $\chi_{\mathrm{hol}}^{(5)}$ (`shuffle_null_chi_path`, `isotropy_null_chi_hol`).

**Label:** $\mathrm{Hol}_E\neq 0$ = **Hypothese**; Nulltests = **Experiment**.

---

## 8. Abgrenzung der Träger

| Symbol | Block | Träger | Was misst es? |
|--------|------:|--------|---------------|
| $\chi_{\mathrm{path}}(N)$ | 4 | Primfolge-Gleitfenster $Q_n^{(4)}$ | **Pfadorientierung** ABCE/CEAB |
| $\chi_{\mathrm{hol}}^{(5)}(N)$ | 5 | Primfolge-Gleitfenster $Q_n^{(5)}$ | **Zyklus-Holonomie** ABCEA/CEABC |
| $\chi_E^{\mathrm{quad}}(N)$ | 4 | Prim-Vierlinge $Q(p)$ | arithmetische Vierlings-Orientierung |

**Label:** Trägervergleich = **Experiment**.

---

## 9. Boxed Grundsatz: Geometrie beginnt bei Übergängen

$$\boxed{\;\text{EABC-Geometrie beginnt bei den Übergängen }T_n,\text{ nicht bei isolierten Klassen }X_n.\;}$$

**Begründung.** Der Übergangsgraph $G_E$ operationalisiert **Transport** entlang der Primzahlordnung. Viererblöcke messen **gerichtete Pfade**; Fünferblöcke schließen den kanonischen Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$ und tragen die holonomische Lesart.

**Label:** **Definition** (epistemischer Grundsatz).

---

## 10. Hierarchie und Dokumentenverknüpfung

```
collatz_eabc_holonomie.md          Stufe 1–4: κ, σ, ω, χ_E^quad (Vierlinge)
        │
        ▼
collatz_eabc_transport.md          Stufe 3–5: G_E, T_n, χ_path
        │
        ▼
collatz_eabc_zyklus_holonomie.md   Stufe 1–7: Q^(4), χ_path; Q^(5), Ω^(5), χ_hol^(5), Hol_E
```

| Stufe | Objekt | Dokument |
|------:|--------|----------|
| 1 | $\kappa(p)$, $X_n$ | hier §1; `collatz_eabc_holonomie.md` |
| 2–3 | $T_n$, $G_E$ | hier §2; `collatz_eabc_transport.md` §3 |
| 4 | $Q_n^{(4)}$, $\chi_{\mathrm{path}}$ | hier §3 |
| 5 | $\chi_{\mathrm{path}}(N)$ | hier §4 |
| 6 | $Q_n^{(5)}$, $\Omega^{(5)}$, $\chi_{\mathrm{hol}}^{(5)}$ | hier §5 |
| 7 | $\mathrm{Hol}_E=\lim\chi_{\mathrm{hol}}^{(5)}$ | hier §6–7 |
| — | $\chi_E^{\mathrm{quad}}(N)$ | `collatz_eabc_holonomie.md` §4 |

---

## 11. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $\kappa(p)$, $X_n=\kappa(p_n)$ | **Definition** |
| $T_n\colon X_n\to X_{n+1}$, $G_E$ | **Definition** |
| $Q_n^{(4)}$, $\chi_{\mathrm{path}}$ — Pfadorientierung, **keine** Holonomie | **Definition** |
| $\chi_{\mathrm{path}}(N)$ auf Gleitfenstern | **Definition** |
| $Q_n^{(5)}$, $\Omega^{(5)}$ — geschlossener Zyklus | **Definition** |
| $\chi_{\mathrm{hol}}^{(5)}(N)$ | **Definition** |
| $\mathrm{Hol}_E=\lim\chi_{\mathrm{hol}}^{(5)}(N)$ | **Definition** |
| $\mathrm{Hol}_E\neq 0$ | **Hypothese** |
| $\chi_{\mathrm{path}}$ vs.\ $\chi_{\mathrm{hol}}^{(5)}$ vs.\ $\chi_E^{\mathrm{quad}}$ | **Experiment** |
| EABC-Geometrie bei Übergängen | **Definition** (Grundsatz) |
| $V_4$ assoziativ, keine $V_4$-Holonomie | **Theorem** (`collatz_eabc_holonomie.md`) |

---

*Kanonsiche Notiz: Viererblöcke (ABCE/CEAB) sind **orientierte Transportpfade**; echte Zyklus-Holonomie erfordert den **Fünferblock** (ABCEA/CEABC) mit geschlossenem Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$. Vierlings-$\chi_E^{\mathrm{quad}}$ bleibt der arithmetische Vergleichsträger auf Prim-Vierlingen.*
