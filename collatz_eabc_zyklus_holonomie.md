# EABC-Zyklus-Holonomie-Hypothese (kanonisch)

**Status:** Definition (Stufen 1–5) + Hypothese (Stufe 6) + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Hypothese | Experiment | Theorem

**Querverweise:**
- `collatz_eabc_transport.md` — gerichteter Übergangsgraph $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_holonomie.md` — Vierlings-Orientierung $\omega(Q)$, $\chi_E^{\mathrm{quad}}(N)$ (arithmetische Vierlinge)
- `collatz_eabc_transition_graph.py` / `.json` — Numerik: $\chi_E(N)$ auf Primfolge-Gleitfenstern, $\mathrm{Hol}_E$-Schätzung, Nullmodelle
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

## 3. Gleitfenster-Quartett $Q_n$ und Zyklus-Orientierung $\Omega$

**Definition (4-Fenster auf der Primfolge).**
$$Q_n := (X_n,\,X_{n+1},\,X_{n+2},\,X_{n+3})
=\bigl(\kappa(p_n),\,\kappa(p_{n+1}),\,\kappa(p_{n+2}),\,\kappa(p_{n+3})\bigr).$$

**Definition (Zyklus-Orientierung).**
$$\Omega(Q_n)=
\begin{cases}
+1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}=\mathrm{ABCE},\\[4pt]
-1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}=\mathrm{CEAB},\\[4pt]
0 & \text{sonst.}
\end{cases}$$

**Wichtig:** $Q_n$ läuft über **aufeinanderfolgende Prim-Klassen** (Gleitfenster), **nicht** über arithmetische Prim-Vierlinge $Q(p)=(p,p{+}2,p{+}6,p{+}8)$.

**Label:** $Q_n$, $\Omega$ = **Definition** (Stufe 4).

---

## 4. Zyklus-Chiralität $\chi_E(N)$

Sei $N$ eine Primzahl-Obergrenze. Fenster $Q_n$ mit $p_{n+3}\le N$ werden berücksichtigt.

**Definition (EABC-Zyklus-Chiralität).**
$$\boxed{\;
\chi_E(N):=
\frac{\displaystyle\sum_{n:\,p_{n+3}\le N}\Omega(Q_n)}
{\displaystyle\#\{n:\,p_{n+3}\le N,\;\Omega(Q_n)\neq 0\}}
\in[-1,1],\;}
$$
mit der Konvention $\chi_E(N)=0$, falls der Nenner $0$ ist.

**Äquivalente Darstellung.** Mit $N_+(N):=\#\{n:\,p_{n+3}\le N,\,\Omega(Q_n)=+1\}$ und $N_-(N):=\#\{n:\,p_{n+3}\le N,\,\Omega(Q_n)=-1\}$:
$$\chi_E(N)=\frac{N_+(N)-N_- (N)}{N_+(N)+N_-(N)}\quad\text{(falls }N_++N_->0\text{).}$$

**Abgrenzung (Vierlings-Träger).** Auf Prim-Vierlingen $Q(p)$ definiert `collatz_eabc_holonomie.md` eine **andere** Observable $\chi_E^{\mathrm{quad}}(N)$ mit gleicher ABCE/CEAB-Formel, aber **arithmetisch gekoppelten** Beinen — verwandt, **nicht identisch**.

| Symbol | Träger | Indizierung |
|--------|--------|-------------|
| $\chi_E(N)$ | Primfolge-Gleitfenster $Q_n$ | $p_{n+3}\le N$ |
| $\chi_E^{\mathrm{quad}}(N)$ | Prim-Vierlinge $Q(p)$ | Vierlingsstart $p\le N$ |

**Experiment:** `collatz_eabc_transition_graph.py::chi_E_sliding`.

**Label:** $\chi_E(N)$ (Gleitfenster) = **Definition** (Stufe 5).

---

## 5. Projektive Zyklus-Holonomie $\mathrm{Hol}_E$

**Definition (Grenzwert).**
$$\boxed{\;\mathrm{Hol}_E := \lim_{N\to\infty}\chi_E(N)\;}$$
sofern der Grenzwert existiert und $\neq 0$ ist — dann als **projektive Zyklus-Holonomie** interpretiert.

Falls $\mathrm{Hol}_E=0$, bleibt die Frage offen, ob $\chi_E(N)$ **stabile, nicht-zufällige Fluktuationen** trägt (schwächere Lesart).

**Experiment:** `collatz_eabc_transition_graph.py::hol_E_estimates` (Stichproben bei $N=10^5$, $10^6$).

**Label:** $\mathrm{Hol}_E$ = **Definition** (Grenzwert); Existenz und $\neq 0$ = **Hypothese**.

---

## 6. Boxed Hypothese (stark)

$$\boxed{\;\mathrm{Hol}_E\neq 0\;}$$

**Präzisierung.** Eine nichttriviale asymptotische Orientierungsasymmetrie der ABCE/CEAB-Zyklen entlang der Primzahl-Transportkette — Hinweis auf **projektive Holonomie** des gerichteten EABC-Übergangsgraphen, kein Beweis im DG-Sinn.

**Schwächere Lesart (falls $\mathrm{Hol}_E=0$):**
$$\boxed{\;\chi_E(N)\text{ zeigt stabile, nicht-zufällige Fluktuationen gegenüber Isotropie- und Shuffle-Nullmodellen.}\;}$$

**Nullmodelle (Experiment):**
- **Marginal-Shuffle:** Permutation der Klassenfolge bei erhaltener Häufigkeitsverteilung (`shuffle_null_chi_E`).
- **Isotropie-Null:** Zufällige Permutation $\sigma\in S_4$ auf $\{E,A,B,C\}$ vor Fensterzählung (`isotropy_null_chi_E`).

**Label:** $\mathrm{Hol}_E\neq 0$ = **Hypothese**; Nulltests = **Experiment**.

---

## 7. Boxed Grundsatz: Geometrie beginnt bei Übergängen

$$\boxed{\;\text{EABC-Geometrie beginnt bei den Übergängen }T_n,\text{ nicht bei isolierten Klassen }X_n.\;}$$

**Begründung.** Der Übergangsgraph $G_E$ und die Gleitfenster $Q_n$ operationalisieren **Transport** entlang der Primzahlordnung. Isolierte Punkte $\kappa(p)$ tragen keine Orientierungsinformation; erst Kanten $\tau(p_n)$ und 4-Schritt-Zyklen $\Omega(Q_n)$ eröffnen die holonomische Lesart.

**Label:** **Definition** (epistemischer Grundsatz).

---

## 8. Hierarchie und Dokumentenverknüpfung

```
collatz_eabc_holonomie.md          Stufe 1–4: κ, σ, ω, χ_E^quad (Vierlinge)
        │
        ▼
collatz_eabc_transport.md          Stufe 3–5: G_E, T_n, χ_trans
        │
        ▼
collatz_eabc_zyklus_holonomie.md   Stufe 1–7: Q_n, Ω, χ_E (Gleitfenster), Hol_E
```

| Stufe | Objekt | Dokument |
|------:|--------|----------|
| 1 | $\kappa(p)$, $X_n$ | hier §1; `collatz_eabc_holonomie.md` |
| 2–3 | $T_n$, $G_E$ | hier §2; `collatz_eabc_transport.md` §3 |
| 4 | $Q_n$, $\Omega$ | hier §3 |
| 5 | $\chi_E(N)$ Gleitfenster | hier §4 |
| 6 | $\mathrm{Hol}_E$ | hier §5–6 |
| — | $\chi_E^{\mathrm{quad}}(N)$ | `collatz_eabc_holonomie.md` §4 |

---

## 9. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $\kappa(p)$, $X_n=\kappa(p_n)$ | **Definition** |
| $T_n\colon X_n\to X_{n+1}$, $G_E$ | **Definition** |
| $Q_n$, $\Omega(Q_n)\in\{+1,-1,0\}$ | **Definition** |
| $\chi_E(N)=\sum\Omega/\#\{\Omega\neq 0\}$ auf Gleitfenstern | **Definition** |
| $\mathrm{Hol}_E=\lim\chi_E(N)$ | **Definition** |
| $\mathrm{Hol}_E\neq 0$ | **Hypothese** |
| $\chi_E$ vs.\ $\chi_E^{\mathrm{quad}}$ — verwandt, nicht identisch | **Experiment** |
| EABC-Geometrie bei Übergängen | **Definition** (Grundsatz) |
| $V_4$ assoziativ, keine $V_4$-Holonomie | **Theorem** (`collatz_eabc_holonomie.md`) |

---

*Kanonsiche Notiz: Dieses Dokument fixiert $\chi_E(N)$ auf **Primfolge-Gleitfenstern** als Zyklus-Holonomie-Kandidat. Vierlings-$\chi_E^{\mathrm{quad}}$ bleibt in `collatz_eabc_holonomie.md` der arithmetische Vergleichsträger — beide messen ABCE/CEAB-Orientierung in verschiedenen Räumen.*
