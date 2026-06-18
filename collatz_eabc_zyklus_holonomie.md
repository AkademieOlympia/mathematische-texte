# EABC-Zyklus-Holonomie-Hypothese (kanonisch)

**Status:** Definition (Stufen 1–6) + Hypothese (Stufe 7) + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Hypothese | Experiment | Theorem

**Querverweise:**
- `collatz_eabc_fehlerterm_hypothese.md` — **kanonische Endform:** $N_\pm$, Hauptvermutung, Fehlerterm-Hypothese, $\widetilde{D}_E$
- `collatz_eabc_bell_holonomie.md` — Bell-Ungleichung als Holonomie-/Zykluskonsistenz; $P_{\mathrm{same}}$, CHSH-Analog (kein QM)
- `collatz_eabc_bell_inequality_test.py` / `.json` — Numerik Bell-Summen, $G_E$-Vergleich
- `collatz_eabc_transport.md` — gerichteter Übergangsgraph $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_holonomie.md` — Vierlings-Orientierung $\omega(Q)$, $\chi_E^{\mathrm{quad}}(N)$ (arithmetische Vierlinge)
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch: mod-$12$-Symmetrie, $\mathrm{Hol}_E=0$, Fehlerterm $D_E$
- `collatz_eabc_transition_graph.py` / `.json` — Numerik: $\chi_{\mathrm{Pfad}}$, $\chi_{\mathrm{Hol}}$, Nullmodelle
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — $N_\pm$, $D_E$, $\widetilde{D}_E$, Lückenmuster $(2,4,2,4)$
- `collatz_eabc_holonomie_test.py` — $\chi_E^{\mathrm{quad}}(N)$ auf Prim-Vierlingen (Vergleichsträger)
- `eabc_from_lean.py` — $\kappa=\texttt{class\_of}$, Rotation $t$

---

## 0. Boxed Hierarchie

$$\boxed{\;\text{Klasse} \;\to\; \text{Kante} \;\to\; \text{Pfad} \;\to\; \text{Zyklus} \;\to\; \text{Holonomie}\;}$$

| Stufe | Objekt | Symbol |
|------:|--------|--------|
| 1 | **Klasse** | $X_n=\kappa(p_n)$ |
| 2 | **Kante** | $\tau_n=(X_n,X_{n+1})$ |
| 3 | **Pfad** | $P_n^{(4)}=(X_n,X_{n+1},X_{n+2},X_{n+3})$ |
| 4 | **Zyklus** | $C_n^{(5)}=(X_n,X_{n+1},X_{n+2},X_{n+3},X_{n+4})$ |
| 5 | **Holonomie** | $\Omega_{\mathrm{Hol}}(C_n^{(5)})$, $\chi_{\mathrm{Hol}}(N)$, $\mathrm{Hol}_E$ |

**Kernbotschaft:** ABCE/CEAB-Befunde bleiben gültig als **Pfadorientierung** ($\Omega_{\mathrm{Pfad}}$); echte **Holonomie** beginnt erst beim geschlossenen Zyklus ($C_n^{(5)}$, ABCEA/CEABC).

---

## 1. EABC-Klassifikation $\kappa(p)$ — Klasse

**Definition.** Für Primzahlen $p>3$:
$$\kappa(p)\in\{E,A,B,C\},\qquad X_n := \kappa(p_n),$$
wobei $(p_n)_{n\ge 1}$ die aufsteigende Folge aller Primzahlen $>3$ ist.

**Label:** $\kappa$, $X_n$ = **Definition** (Stufe 1 — Klasse).

**Implementierung:** `eabc_from_lean.py::class_of`, `EABC.lean`.

---

## 2. Transport $T_n$ und Kante $\tau_n$

**Definition (Transportoperator).** Entlang der Primzahlordnung:
$$T_n\colon X_n \longrightarrow X_{n+1},\qquad
\tau_n := \bigl(X_n,\,X_{n+1}\bigr)
=\bigl(\kappa(p_n),\,\kappa(p_{n+1})\bigr).$$

**Definition (Zustandsraum).** $V_4=\{E,A,B,C\}$.

**Definition (gerichteter Übergangsgraph $G_E$).** Kante $i\to j$ wird gezählt, wenn $\kappa(p_n)=i$ und $\kappa(p_{n+1})=j$.

Der **fundamentale Zustand** ist die **gerichtete Kante** $\tau_n$, nicht der isolierte Punkt $X_n$.

**Label:** $T_n$, $G_E$, $\tau_n$ = **Definition** (Stufe 2 — Kante).

**Experiment:** `collatz_eabc_transition_graph.py::transition_matrix_report`.

---

## 3. Pfad $P_n^{(4)}$ — Pfadorientierung

$$\boxed{\;\text{Pfad} = \text{orientierter Transportpfad (kein geschlossener Zyklus).}\;}$$

**Definition (4-Fenster auf der Primfolge).**
$$P_n^{(4)} := (X_n,\,X_{n+1},\,X_{n+2},\,X_{n+3})
=\bigl(\kappa(p_n),\,\kappa(p_{n+1}),\,\kappa(p_{n+2}),\,\kappa(p_{n+3})\bigr).$$

**Definition (Pfadorientierung).** ABCE und CEAB sind **EABC-Pfadorientierungen** entlang vier aufeinanderfolgender Transportkanten — **nicht** Holonomie im strengen DG-Sinn (der Pfad schließt nicht):
$$\Omega_{\mathrm{Pfad}}(P_n^{(4)}):=
\begin{cases}
+1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}=\mathrm{ABCE},\\[4pt]
-1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}=\mathrm{CEAB},\\[4pt]
0 & \text{sonst.}
\end{cases}$$

**Wichtig:** $P_n^{(4)}$ läuft über **aufeinanderfolgende Prim-Klassen** (Gleitfenster), **nicht** über arithmetische Prim-Vierlinge $Q(p)=(p,p{+}2,p{+}6,p{+}8)$.

**Label:** $P_n^{(4)}$, $\Omega_{\mathrm{Pfad}}$ = **Definition** (Stufe 3 — Pfad).

**Legacy-Aliase:** $Q_n^{(4)}:=P_n^{(4)}$; $\chi_{\mathrm{path}}:=\Omega_{\mathrm{Pfad}}$ (Python: `omega_path`, `chi_path`).

---

## 4. Pfad-Chiralität $\chi_{\mathrm{Pfad}}(N)$

Sei $N$ eine Primzahl-Obergrenze. Pfade $P_n^{(4)}$ mit $p_{n+3}\le N$ werden berücksichtigt.

**Definition (gerichtete 4-Pfade).**
$$\chi_{\mathrm{Pfad}}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+3}\le N}\Omega_{\mathrm{Pfad}}(P_n^{(4)})}
{\displaystyle\#\{n:\,p_{n+3}\le N,\;\Omega_{\mathrm{Pfad}}(P_n^{(4)})\neq 0\}}
\in[-1,1],$$
mit der Konvention $\chi_{\mathrm{Pfad}}(N)=0$, falls der Nenner $0$ ist.

**Experiment:** `collatz_eabc_transition_graph.py::chi_pfad_sliding`.

**Label:** $\chi_{\mathrm{Pfad}}(N)$ = **Definition** (Observable für gerichtete 4-Pfade).

**Legacy-Alias:** $\chi_{\mathrm{path}}(N):=\chi_{\mathrm{Pfad}}(N)$.

---

## 5. Zyklus $C_n^{(5)}$ — geschlossener Zyklus / Holonomie

$$\boxed{\;\text{Zyklus} = \text{geschlossener Transportzyklus; Holonomie beginnt hier.}\;}$$

**Definition (5-Fenster auf der Primfolge).**
$$C_n^{(5)} := (X_n,\,X_{n+1},\,X_{n+2},\,X_{n+3},\,X_{n+4}).$$

Der **geschlossene EABC-Zyklus** ist $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$. Entlang der Primfolge erscheinen die Wörter **ABCEA** (positiv) und **CEABC** (negativ, umgekehrte Zyklusorientierung).

**Definition (Zyklus-Holonomie $\Omega_{\mathrm{Hol}}$).**
$$\Omega_{\mathrm{Hol}}(C_n^{(5)})=
\begin{cases}
+1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}X_{n+4}=\mathrm{ABCEA},\\[4pt]
-1 & \text{wenn } X_nX_{n+1}X_{n+2}X_{n+3}X_{n+4}=\mathrm{CEABC},\\[4pt]
0 & \text{sonst.}
\end{cases}$$

**Definition (Holonomie-Observable auf echten 5-Zyklen).**
$$\chi_{\mathrm{Hol}}(N):=
\frac{\displaystyle\sum_{n:\,p_{n+4}\le N}\Omega_{\mathrm{Hol}}(C_n^{(5)})}
{\displaystyle\#\{n:\,p_{n+4}\le N,\;\Omega_{\mathrm{Hol}}(C_n^{(5)})\neq 0\}}
\in[-1,1],$$
mit $\chi_{\mathrm{Hol}}(N)=0$, falls der Nenner $0$ ist.

**Experiment:** `collatz_eabc_transition_graph.py::chi_hol_sliding`.

**Label:** $C_n^{(5)}$, $\Omega_{\mathrm{Hol}}$, $\chi_{\mathrm{Hol}}$ = **Definition** (Stufe 4–5 — Zyklus / Holonomie).

**Legacy-Aliase:** $Q_n^{(5)}:=C_n^{(5)}$; $\Omega^{(5)}:=\Omega_{\mathrm{Hol}}$; $\chi_{\mathrm{hol}}^{(5)}:=\chi_{\mathrm{Hol}}$.

---

## 6. Zählgrößen $N_\pm$ und projektive Holonomie $\mathrm{Hol}_E$

**Definition ($N_+$, $N_-$).** Für Prim-Obergrenze $X$:
$$N_+(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{CEABC}\}.$$

**Definition (Fehlerterm und Normalisierung).**
$$D_E(X) := N_+(X)-N_-(X),\qquad
\widetilde{D}_E(X) := \frac{D_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

**Definition (Grenzwert / Hauptterm).**
$$\chi_{\mathrm{Hol}}(X)=\frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)},\qquad
\boxed{\;\mathrm{Hol}_E := \lim_{X\to\infty}\chi_{\mathrm{Hol}}(X)\;}$$
sofern der Grenzwert existiert.

**Nicht** $\lim\chi_{\mathrm{Pfad}}(N)$: der 4-Pfad misst Pfadorientierung, keine geschlossene Holonomie.

### 6.1 Hauptvermutung (konservativ)

$$\boxed{\;N_+(X)\sim N_-(X)\;\Rightarrow\;\mathrm{Hol}_E=0\;}$$

**Begründung (Skizze):** `collatz_eabc_holonomie_beweisversuch.md` — mod-$12$-Symmetrie, gemeinsames Lückenmuster $(2,4,2,4)$, HL-Äquidistribution.

**Label:** Hauptvermutung = **Vermutung**.

### 6.2 Fehlerterm-Hypothese (stärker)

Selbst bei $\mathrm{Hol}_E=0$ bleibt $D_E(X)$ nichttrivial. **Fehlerterm-Hypothese:** Chebyshev-artiger Bias in $D_E$, gesteuert durch Nullstellen der Dirichlet-$L$-Funktionen mod $12$.

**Zentrale Frage:** Verhält sich $\widetilde{D}_E$ wie reines Rauschen oder zeigt stabile Vorzeichenasymmetrie / Oszillation?

**Analytischer Beweisversuch:** `collatz_eabc_holonomie_beweisversuch.md`; **Endform:** `collatz_eabc_fehlerterm_hypothese.md`.

**Experiment:** `collatz_eabc_transition_graph.py::hol_E_estimates`; `collatz_eabc_holonomie_fehlerterm.py` ($N_\pm$, $D_E$, $\widetilde{D}_E$-Zeitreihe).

**Label:** Fehlerterm-Hypothese = **Hypothese**; $\mathrm{Hol}_E$ als Grenzwert = **Definition**.

---

## 7. Boxed Hypothese (stark, legacy)

$$\boxed{\;\mathrm{Hol}_E\neq 0\;}$$

**Präzisierung.** Eine nichttriviale asymptotische Orientierungsasymmetrie der geschlossenen ABCEA/CEABC-Zyklen entlang der Primzahl-Transportkette — Hinweis auf **projektive Holonomie** des gerichteten EABC-Übergangsgraphen, kein Beweis im DG-Sinn.

**Schwächere Lesart (falls $\mathrm{Hol}_E=0$, operative Endform):**
$$\boxed{\;\widetilde{D}_E(X)\text{ zeigt strukturierte Fluktuationen (Bias/Oszillation), nicht permanente Hauptterm-Holonomie.}\;}$$
Siehe `collatz_eabc_fehlerterm_hypothese.md` §5–7.

**Nullmodelle (Experiment):**
- **Marginal-Shuffle** und **Isotropie-Null** jeweils für $\chi_{\mathrm{Pfad}}$ und $\chi_{\mathrm{Hol}}$ (`shuffle_null_chi_pfad`, `isotropy_null_chi_hol`).

**Label:** $\mathrm{Hol}_E\neq 0$ = **Hypothese**; Nulltests = **Experiment**.

---

## 8. Abgrenzung der Träger

| Symbol | Block | Träger | Was misst es? |
|--------|------:|--------|---------------|
| $\chi_{\mathrm{Pfad}}(N)$ | 4 | Primfolge-Gleitfenster $P_n^{(4)}$ | **Pfadorientierung** ABCE/CEAB |
| $\chi_{\mathrm{Hol}}(N)$ | 5 | Primfolge-Gleitfenster $C_n^{(5)}$ | **Zyklus-Holonomie** ABCEA/CEABC |
| $\chi_E^{\mathrm{quad}}(N)$ | 4 | Prim-Vierlinge $Q(p)$ | arithmetische Vierlings-Orientierung |

**Label:** Trägervergleich = **Experiment**.

---

## 9. Boxed Grundsatz: Geometrie beginnt bei Übergängen

$$\boxed{\;\text{EABC-Geometrie beginnt bei den Übergängen }\tau_n,\text{ nicht bei isolierten Klassen }X_n.\;}$$

**Begründung.** Der Übergangsgraph $G_E$ operationalisiert **Transport** entlang der Primzahlordnung. Pfade messen **gerichtete 4-Ketten**; Zyklen schließen den kanonischen Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$ und tragen die holonomische Lesart.

**Label:** **Definition** (epistemischer Grundsatz).

---

## 10. Hierarchie und Dokumentenverknüpfung

```
collatz_eabc_holonomie.md          Stufe 1–4: κ, σ, ω, χ_E^quad (Vierlinge)
        │
        ▼
collatz_eabc_transport.md          Stufe 2–4: G_E, τ_n, P_n^(4), χ_Pfad
        │
        ▼
collatz_eabc_zyklus_holonomie.md   Klasse→Kante→Pfad→Zyklus→Holonomie
```

| Stufe | Objekt | Dokument |
|------:|--------|----------|
| 1 | $X_n=\kappa(p_n)$ — **Klasse** | hier §1; `collatz_eabc_holonomie.md` |
| 2 | $\tau_n=(X_n,X_{n+1})$ — **Kante** | hier §2; `collatz_eabc_transport.md` §3 |
| 3 | $P_n^{(4)}$, $\Omega_{\mathrm{Pfad}}$ — **Pfad** | hier §3 |
| 4 | $\chi_{\mathrm{Pfad}}(N)$ | hier §4 |
| 5 | $C_n^{(5)}$, $\Omega_{\mathrm{Hol}}$, $\chi_{\mathrm{Hol}}(N)$ — **Zyklus / Holonomie** | hier §5 |
| 6 | $N_\pm$, $D_E$, $\widetilde{D}_E$, $\mathrm{Hol}_E$ | hier §6; `collatz_eabc_fehlerterm_hypothese.md` |
| — | $\chi_E^{\mathrm{quad}}(N)$ | `collatz_eabc_holonomie.md` §4 |

---

## 11. Python-Symbolzuordnung

| LaTeX | Python (kanonisch) | Legacy-Alias |
|-------|-------------------|--------------|
| $\Omega_{\mathrm{Pfad}}(P_n^{(4)})$ | `omega_pfad` | `omega_path`, `omega_window` |
| $\Omega_{\mathrm{Hol}}(C_n^{(5)})$ | `omega_hol` | `omega_5` |
| $\chi_{\mathrm{Pfad}}(N)$ | `chi_pfad_sliding` | `chi_path_sliding`, `chi_E_sliding` |
| $\chi_{\mathrm{Hol}}(N)$ | `chi_hol_sliding` | — |
| $N_+(X)$, $N_-(X)$ | `N_plus`, `N_minus` | `N_ABCEA`, `N_CEABC` |
| $D_E$, $\widetilde{D}_E$ | `D_E`, `D_tilde_E` | — |
| Vergleich Pfad vs. Holonomie | `chi_pfad_vs_hol` | `chi_path_vs_hol` |

---

## 12. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $\kappa(p)$, $X_n=\kappa(p_n)$ — **Klasse** | **Definition** |
| $\tau_n=(X_n,X_{n+1})$, $G_E$ — **Kante** | **Definition** |
| $P_n^{(4)}$, $\Omega_{\mathrm{Pfad}}$ — Pfadorientierung, **keine** Holonomie | **Definition** |
| $\chi_{\mathrm{Pfad}}(N)$ auf Gleitfenstern | **Definition** |
| $C_n^{(5)}$, $\Omega_{\mathrm{Hol}}$ — geschlossener Zyklus | **Definition** |
| $\chi_{\mathrm{Hol}}(N)$ | **Definition** |
| $N_\pm$, $D_E$, $\widetilde{D}_E$ | **Definition** |
| $N_+\sim N_-$ $\Rightarrow$ $\mathrm{Hol}_E=0$ | **Vermutung** (Hauptvermutung) |
| Fehlerterm-Hypothese ($L$-Funktionen mod $12$) | **Hypothese** |
| $\mathrm{Hol}_E=\lim\chi_{\mathrm{Hol}}(N)$ | **Definition** |
| $\mathrm{Hol}_E\neq 0$ | **Hypothese** |
| $\chi_{\mathrm{Pfad}}$ vs.\ $\chi_{\mathrm{Hol}}$ vs.\ $\chi_E^{\mathrm{quad}}$ | **Experiment** |
| EABC-Geometrie bei Übergängen | **Definition** (Grundsatz) |
| $V_4$ assoziativ, keine $V_4$-Holonomie | **Theorem** (`collatz_eabc_holonomie.md`) |

---

*Kanonsiche Notiz: ABCE/CEAB auf $P_n^{(4)}$ sind **orientierte Transportpfade** ($\Omega_{\mathrm{Pfad}}$); echte Zyklus-Holonomie erfordert $C_n^{(5)}$ (ABCEA/CEABC) mit geschlossenem Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$. Vierlings-$\chi_E^{\mathrm{quad}}$ bleibt der arithmetische Vergleichsträger auf Prim-Vierlingen.*
