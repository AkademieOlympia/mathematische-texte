# EABC-Bell-Ungleichung als Holonomie- und Zykluskonsistenz

**Status:** Analogie + Definition + Hypothese + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Analogie | Definition | Hypothese | Theorem | Experiment

**Epistemische Abgrenzung:** Dieses Dokument ist **keine** Quantenphysik-Behauptung. Es überträgt die **kombinatorische Taubenloch-Logik** der Bell-Tripel-Ungleichung auf den **diskreten EABC-Transportgraphen** $G_E$ und die Holonomie-Hierarchie. Verletzungen sind **arithmetische Korrelationsphänomene** (Bias, Fehlerterm $D_E$), nicht „nichtlokalität“ im QM-Sinn.

**Querverweise:**
- `collatz_eabc_zyklus_holonomie.md` — Klasse→Kante→Pfad→Zyklus→Holonomie; $\Omega_{\mathrm{Pfad}}$, $\Omega_{\mathrm{Hol}}$
- `collatz_eabc_fehlerterm_hypothese.md` — $N_\pm$, $D_E$, $\widetilde{D}_E$, Hauptvermutung $\mathrm{Hol}_E=0$
- `collatz_eabc_transport.md` — $G_E$, $T_{ij}$, $t$-Rotation
- `collatz_eabc_bell_inequality_test.py` / `.json` — Numerik $P_{\mathrm{same}}$, CHSH-Analog, Vergleich $G_E$, **gemeinsamer Report mit $D_E$**
- `collatz_eabc_holonomie_fehlerterm.py` — $D_E$, $\widetilde{D}_E$, Chebyshev-Vergleich, **$S_{\mathrm{EABC}}$-Brücke**
- `collatz_eabc_core/CollatzEabc/HolonomieFehlerterm.lean` — Lean Phase 1–2 (Lückenmuster bewiesen; Prime/CHSH-Skeleton)
- `collatz_mathlib_eabc_kandidaten.md` — Mathlib `PrimeCounting` / `DirichletCharacter` (noch nicht für Bell/$D_E$)
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54

---

## 0. Boxed Kernidee

$$\boxed{\;\text{Bell-Tripel} \;\leadsto\; \text{drei binäre EABC-Lesarten auf demselben Transportfenster} \;\Rightarrow\; \text{Taubenloch-Konsistenz } \ge 1.\;}$$

$$\boxed{\;\text{„Verletzung''} = \text{marginalisierte oder kontextgemischte Lesarten} \;\not\equiv\; \text{gemeinsame versteckte Zustände pro Fenster.}\;}$$

---

## 1. Bell im Mengen-/Wahrscheinlichkeitskalkül (Referenz)

**Setting (klassisch, lokal, deterministisch pro Zustand).** Drei binäre Observablen $A,B,C\in\{0,1\}$ werden am **selben** versteckten Zustand $\lambda$ gelesen. Für festes $\lambda$ sind $a(\lambda),b(\lambda),c(\lambda)$ drei Bits.

**Taubenloch.** In $\{0,1\}$ gibt es unter drei Bits mindestens **ein gleiches Paar**:
$$\mathbb{1}[a=b]+\mathbb{1}[a=c]+\mathbb{1}[b=c]\;\ge\;1.$$

**Ungleichung.** Mit $P_{\mathrm{same}}(X,Y):=\Pr[X=Y]$ über die Verteilung von $\lambda$:
$$\boxed{\;P_{\mathrm{same}}(A,B)+P_{\mathrm{same}}(A,C)+P_{\mathrm{same}}(B,C)\;\ge\;1.\;}$$

**QM-Kontext (nur Referenz, nicht EABC-These).** Verschränkte Spins verletzen CHSH ($|S|\le 2$ klassisch, $2\sqrt2$ quantenmechanisch). EABC übernimmt **nur** die kombinatorische Struktur „drei binäre Lesarten + Konsistenzsumme“, nicht Hilbert-Raum oder Tensorprodukt.

**Label:** Bell-Tripel-Ungleichung = **Analogie** (Referenz aus klassischer Wahrscheinlichkeit).

---

## 2. EABC-Transport: Zustand und Kanten

**Definition (Zustandsraum).** $V_4=\{E,A,B,C\}$, $\kappa(p)\in V_4$ für Primzahlen $p>3$.

**Definition (Transportkette).** $X_n:=\kappa(p_n)$, Kante $\tau_n=(X_n,X_{n+1})$.

**Definition ($t$-Rotation).** $t\colon E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$, $t^4=\mathrm{id}$ (`eabc_from_lean.py`).

**Definition (binäre Kantenlesart).**
$$\sigma(n):=\mathbb{1}\bigl[X_{n+1}=t(X_n)\bigr]\in\{0,1\}.$$
$\sigma=1$ = **$t$-alignierte** (vorwärts-kanonische) Kante; $\sigma=0$ = nicht $t$-aligniert.

**Label:** $\sigma(n)$ = **Definition** (Stufe 2 — Kante).

---

## 3. Mapping Bell $(A,B,C)\to$ EABC-Dreieck

### 3.1 Kanonische Zuordnung

| Bell | EABC-Knoten | Rolle im 4-Zyklus $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ |
|------|-------------|----------------------------------------------------------|
| $A$ | $E$ | Start/Erzeuger-Knoten |
| $B$ | $A$ | zweiter Schritt |
| $C$ | $C$ | vierter Schritt (vor Rückkehr zu $E$) |

Alternativ **zyklische Lesart** $(A,B,C)\mapsto(A,B,C)$ auf dem Dreieck $A$–$B$–$C$ im Vierergraph; die **operative** Numerik in `collatz_eabc_bell_inequality_test.py` verwendet primär $(E,A,C)$ auf **ABCE-Fenstern** (siehe §4).

**Label:** Knoten-Mapping = **Definition** (Analogie-Brücke).

### 3.2 Gemeinsamer Träger (Lokalität analog)

**Definition (gemeinsames Fenster).** Ein **versteckter Zustand** im EABC-Sinn ist ein Gleitfenster
$$W_n^{(4)}=P_n^{(4)}=(X_n,X_{n+1},X_{n+2},X_{n+3})$$
der Länge $4$ auf der Primfolge.

**Definition (drei gleichzeitige Lesarten auf $W_n^{(4)}=\mathrm{ABCE}$).** Wort $\mathrm{ABCE}$ bedeutet $X_n=A$, $X_{n+1}=B$, $X_{n+2}=C$, $X_{n+3}=E$ (`collatz_eabc_zyklus_holonomie.md` §3):
$$\begin{aligned}
O_E(W_n) &:= \sigma(n+3),\\
O_A(W_n) &:= \sigma(n),\\
O_C(W_n) &:= \sigma(n+2),
\end{aligned}$$
wobei $\sigma(n+k)$ die $t$-Alignierung der Kante ab Position $n+k$ ist (für $O_E$ ist $p_{n+4}$ nötig).

Analog für $W_n^{(4)}=\mathrm{CEAB}$ mit umgekehrter Zyklusorientierung (Rollen vertauscht).

**Label:** gemeinsames Fenster = **Definition** (EABC-„Lokalität'' = ein Transportfenster).

---

## 4. $P_{\mathrm{same}}$-Analoge

### 4.1 Gleichfenster-Koinzidenz (streng)

**Definition.** Für Paare $(i,j)\in\{(E,A),(E,C),(A,C)\}$ auf ABCE-Fenstern:
$$P_{\mathrm{same}}^{\mathrm{win}}(i,j;N):=
\frac{\#\{n:\,P_n^{(4)}=\mathrm{ABCE},\;O_i=O_j\}}
{\#\{n:\,P_n^{(4)}=\mathrm{ABCE}\}}.$$

**Definition (Bell-Summe auf gemeinsamem Träger).**
$$\mathcal{B}_{\mathrm{win}}(N):=
P_{\mathrm{same}}^{\mathrm{win}}(E,A;N)+
P_{\mathrm{same}}^{\mathrm{win}}(E,C;N)+
P_{\mathrm{same}}^{\mathrm{win}}(A,C;N).$$

**Theorem (Taubenloch auf ABCE-Fenstern).** Für jedes einzelne Fenster $P_n^{(4)}=\mathrm{ABCE}$ gilt
$$\mathbb{1}[O_E=O_A]+\mathbb{1}[O_E=O_C]+\mathbb{1}[O_A=O_C]\ge 1,$$
also $\mathcal{B}_{\mathrm{win}}(N)\ge 1$ **exakt** (sofern ABCE-Fenster existieren).

**Beweis.** Drei Bits $(O_E,O_A,O_C)\in\{0,1\}^3$ — Taubenloch, wie §1.

**Label:** $\mathcal{B}_{\mathrm{win}}\ge 1$ = **Theorem** (kombinatorisch, kein Grenzwert).

### 4.2 Marginalisierte Koinzidenz (schwächer, kann „verletzen'')

**Definition.** Unabhängig von gemeinsamen Fenstern:
$$P_{\mathrm{same}}^{\mathrm{marg}}(i,j;N):=
\sum_{a,b\in\{0,1\}}\min\bigl(\Pr[O_i=a\mid X_n=i],\,\Pr[O_j=b\mid X_n=j]\bigr)\;\mathbb{1}[a=b],$$
wobei die Ränder über **alle** Kanten mit Startklasse $i$ bzw. $j$ gebildet werden (nicht nur ABCE).

**Hypothese (Bias-Lesart).** $\mathcal{B}_{\mathrm{marg}}:=\sum_{(i,j)}P_{\mathrm{same}}^{\mathrm{marg}}(i,j)$ kann **$<1$** sein, weil die Ränder **keinen gemeinsamen versteckten Zustand** pro Fenster repräsentieren — analog zu kontextabhängigen „Messungen'' ohne gemeinsame Realität.

**Label:** $P_{\mathrm{same}}^{\mathrm{marg}}$ = **Definition**; $\mathcal{B}_{\mathrm{marg}}<1$ = **Hypothese** / **Experiment**.

### 4.3 Holonomie-/Pfad-Lesart (Orientierungskoinzidenz)

**Definition (Pfad-Observable).** Auf $P_n^{(4)}$:
$$\widetilde{O}_{\mathrm{Pfad}}(n):=\frac{1+\Omega_{\mathrm{Pfad}}(P_n^{(4)})}{2}\in\{0,1\}\quad(\Omega_{\mathrm{Pfad}}\neq 0),$$
sonst undefiniert.

**Definition (Holonomie-Observable).** Auf $C_n^{(5)}$:
$$\widetilde{O}_{\mathrm{Hol}}(n):=\frac{1+\Omega_{\mathrm{Hol}}(C_n^{(5)})}{2}\in\{0,1\}\quad(\Omega_{\mathrm{Hol}}\neq 0).$$

**Definition (Orientierungs-$P_{\mathrm{same}}$).** Koinzidenz von $\widetilde{O}_{\mathrm{Pfad}}$ und $\widetilde{O}_{\mathrm{Hol}}$ auf **überlappenden** Indizes $n$ (4-Fenster in 5-Fenster eingebettet):
$$P_{\mathrm{same}}^{\mathrm{hol}}(N):=
\frac{\#\{n:\,\widetilde{O}_{\mathrm{Pfad}}(n)=\widetilde{O}_{\mathrm{Hol}}(n),\,\text{beide definiert}\}}
{\#\{n:\,\text{beide definiert}\}}.$$

**Interpretation.** Wenn $\chi_{\mathrm{Hol}}(N)\neq 0$ oder $D_E(X)\neq 0$, driftet die Pfad-Holonomie-Kohärenz — **Fehlerterm** statt globaler $V_4$-Assoziativität (`collatz_eabc_fehlerterm_hypothese.md`).

**Label:** $P_{\mathrm{same}}^{\mathrm{hol}}$ = **Definition**; Bezug zu $D_E$ = **Hypothese**.

---

## 5. Taubenloch $\Rightarrow$ Bell-Typ auf $G_E$

**Definition (Dreieck $E$–$A$–$C$ im Transportgraphen).** Kanten $\tau_{E\to A}$, $\tau_{A\to B}$, $\tau_{B\to C}$, $\tau_{C\to E}$ sind die **kanonischen** $t$-Kanten des 4-Zyklus (nicht jede Prim-Transition folgt ihnen).

**Satz (Zykluskonsistenz, kombinatorisch).** Auf jedem ABCE-Fenster ist die Tripel-Lesart $(\sigma(n),\sigma(n+1),\sigma(n+3))$ ein Bell-Zustand; die Ungleichung aus §1 gilt **punktweise**.

**Heuristik (global vs. lokal).** Die **Übergangsmatrix** $T_{ij}(N)$ von $G_E$ (`collatz_eabc_transport.md` §3) beschreibt **marginalisierte** Sprünge. Taubenloch gilt auf **gemeinsamen Fenstern**; $T_{ij}$ allein erzwingt **keine** Bell-Summe $\ge 1$ für marginalisierte $P_{\mathrm{same}}$ — das ist der EABC-Analogon zu „klassische Korrelationen ohne gemeinsame Realität pro Kontext''.

**Label:** Zykluskonsistenz = **Theorem** (Fenster); $G_E$-Marginalen = **Heuristik**.

---

## 6. Wann „verletzt'' EABC die klassische Schranke?

| Phänomen | EABC-Analog | Klassische Bell-Schranke |
|----------|-------------|--------------------------|
| $\mathcal{B}_{\mathrm{win}}$ auf ABCE | **nie** verletzt (Theorem) | immer $\ge 1$ |
| $\mathcal{B}_{\mathrm{marg}}$ | kann $<1$ sein | Verletzung = falsches Modell (kein gemeinsamer $\lambda$-Träger) |
| $D_E(X)\neq 0$, $\chi_{\mathrm{Hol}}(N)\neq 0$ | Pfad/Holonomie-Drift, Fehlerterm | **kein** QM; strukturierter Bias |
| $\widetilde{D}_E$ Oszillation | „Chebyshev-artig'' mod $12$ | Korrelations-„Anomalie'' ohne Hauptterm |

**Wichtig:** Eine empirische $\mathcal{B}_{\mathrm{marg}}<1$ **widerlegt nicht** das Taubenloch-Theorem auf Fenstern; sie zeigt, dass **marginalisierte** EABC-Observablen **nicht** als gemeinsame versteckte Variable pro Fenster interpretiert werden dürfen.

**Label:** Interpretationstabelle = **Hypothese** + **Experiment**.

---

## 7. CHSH-Analog auf vier EABC-Observablen (Zyklus $E$–$A$–$B$–$C$)

### 7.1 Referenz: klassisches CHSH

Bei zwei Messparteien mit je zwei Einstellungen ($a,a'$ bzw. $b,b'$) und Korrelationen
$$E(\alpha,\beta):=\mathbb{E}[A_\alpha B_\beta]\in[-1,1]$$
gilt für **lokale versteckte Variablen** (LHV):
$$\boxed{\;|S|\le 2,\qquad S:=E(a,b)-E(a,b')+E(a',b)+E(a',b').\;}$$
Quantenmechanik (nur Referenz): $|S|_{\max}=2\sqrt2$ bei verschränkten Spins.

**Label:** CHSH-LHV-Schranke = **Analogie** (Referenz, kein EABC-Theorem).

### 7.2 Mapping $(a,a',b,b')\to$ EABC-Lesarten

**Gemeinsamer Träger.** Index $n$ mit $P_n^{(4)}=\mathrm{ABCE}$, $p_{n+4}$ definiert, $\Omega_{\mathrm{Pfad}}(P_n^{(4)})\neq 0$, $\Omega_{\mathrm{Hol}}(C_n^{(5)})\neq 0$.

| CHSH | EABC-Observable | Stufe | Definition auf gemeinsamem $n$ (Wort $\mathrm{ABCE}$) |
|------|-----------------|-------|--------------------------------------------------------|
| $a$ | **Kante-$E$** | Kante | $\sigma(n+3)$ bei $X_{n+3}=E$ |
| $a'$ | **Pfad-$A$** | Pfad | $\widetilde{O}_{\mathrm{Pfad}}(n)=\frac{1+\Omega_{\mathrm{Pfad}}(P_n^{(4)})}{2}$ |
| $b$ | **Kante-$C$** | Kante | $\sigma(n+2)$ bei $X_{n+2}=C$ |
| $b'$ | **Holonomie** | Zyklus | $\widetilde{O}_{\mathrm{Hol}}(n)=\frac{1+\Omega_{\mathrm{Hol}}(C_n^{(5)})}{2}$ |

**Interpretation.** Alice liest am **$E$-Knoten** ($\sigma$ vs. Pfadorientierung am Fensterstart); Bob am **$C$-Knoten** ($\sigma$ vs. geschlossene Holonomie). Alle vier Werte stammen vom **selben** $n$ mit $P_n^{(4)}=\mathrm{ABCE}$ und definierter 5-Holonomie.

**Label:** CHSH$\to$EABC-Mapping = **Definition**.

### 7.3 Korrelationsfunktion $E^{\mathrm{EABC}}$

**Definition (empirische Korrelation auf gemeinsamem Träger).**
$$E^{\mathrm{EABC}}(\alpha,\beta):=2\,\Pr[O_\alpha=O_\beta]-1\in[-1,1],$$
wobei $O_a,O_{a'},O_b,O_{b'}$ die vier Bits aus §7.2 sind und $\Pr$ über alle gültigen Indizes $n$ gebildet wird.

**Definition (CHSH-Summe, Standardvorzeichen).**
$$\boxed{\;S_{\mathrm{EABC}}:=E^{\mathrm{EABC}}(a,b)-E^{\mathrm{EABC}}(a,b')+E^{\mathrm{EABC}}(a',b)+E^{\mathrm{EABC}}(a',b').\;}$$

**Definition (Betrag).** $|S_{\mathrm{EABC}}|$ wird mit Schranken verglichen:
- **LHV-Analog:** $|S_{\mathrm{EABC}}|\le 2$ (gilt, wenn vier Bits pro $n$ aus **einer** binären LHV pro Fenster faktorisieren);
- **QM-Referenz:** $2\sqrt2\approx 2{,}828$ (nur Vergleichszahl, **keine** EABC-Vorhersage).

**Label:** $E^{\mathrm{EABC}}$, $S_{\mathrm{EABC}}$ = **Definition**.

### 7.4 Wann wäre $|S_{\mathrm{EABC}}|>2$ „quantum-like''?

**Hypothese (arithmetisches Korrelationsüberschuss-Signal).** $|S_{\mathrm{EABC}}|>2$ tritt auf, wenn Pfad- und Holonomie-Lesarten **nicht faktorisierbar** über eine einzige lokale Variable pro Fenster sind — z. B. durch:

1. **Nichtassoziative Transportmischung:** $G_E$-Übergänge verletzen gleichzeitige $t$-Konsistenz auf $E$- und $C$-Seite;
2. **Holonomie-Fehlerterm:** $D_E(X)\neq 0$, $\chi_{\mathrm{Hol}}(N)\neq 0$ erzeugen **Pfad↔Holonomie-Drift** (`collatz_eabc_fehlerterm_hypothese.md`);
3. **Kontextvermischung:** marginale CHSH-Konstruktion über **verschiedene** Träger (nicht §7.2) — dann ist $|S|>2$ **kein** physikalisches Signal, sondern Modellfehler.

$$\boxed{\;\text{„Quantum-like'' in EABC} = \text{nicht-faktorisierbare Zyklus-Holonomie auf Prim-Transport, nicht QM.}\;}$$

**Label:** Interpretation $|S|>2$ = **Hypothese**; empirischer Wert = **Experiment**.

### 7.5 LHV-Schranke auf gemeinsamem Fenster (Skizze)

**Theorem-Skizze (4-Bit-LHV pro Fenster).** Fixiere $n$. Vier Bits $(o_a,o_{a'},o_b,o_{b'})\in\{0,1\}^4$ sind Funktionen **eines** $\lambda_n$. Die CHSH-Polynom-Bilanz für deterministische $\pm1$-Werte liefert $|S|\le 4$; mit $E=2P_{\mathrm{same}}-1$ folgt $|S|\le 2$ unter LHV-Faktorisierung Alice$\times$Bob.

**Abweichung in Daten:** Wenn empirisch $|S_{\mathrm{EABC}}|>2$ auf **gemeinsamem** Träger §7.2, interpretieren wir **projektive Holonomie** / strukturierten Bias — analog zu $D_E$, nicht zu Bell-Nichtlokalität.

**Label:** LHV-Skizze = **Analogie**; empirische Verletzung = **Experiment**.

### 7.6 Experiment

`collatz_eabc_bell_inequality_test.py::chsh_eabc_cycle_report` — Primfolge bis $10^6$; Ausgabe $E^{\mathrm{EABC}}(a,b)$, …, $S_{\mathrm{EABC}}$, $|S_{\mathrm{EABC}}|$ vs. $2$ und $2\sqrt2$.

**Label:** Numerik = **Experiment**.

---

## 8. Holonomie: Dreieck vs. Viereck

### 8.1 Dreieck $E$–$A$–$C$ (offen)

Auf ABCE-Fenster: drei Kantenbits $(\sigma_E,\sigma_A,\sigma_C)$ — **kein geschlossener Zyklus**, Taubenloch-Konsistenz (§4.1).

### 8.2 Viereck $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$ (geschlossen)

**Definition.** $\Omega_{\mathrm{Hol}}(C_n^{(5)})=\pm 1$ auf ABCEA/CEABC — **globale Orientierung** des geschlossenen Zyklus.

**Analogie zur DG-Holonomie.** Paralleles Transportieren entlang $t$ und Zurückkehren über Prim-Abweichungen erzeugt **nicht-triviale Fehlerphase** $D_E$ auch bei $\mathrm{Hol}_E=0$ im Hauptterm (`collatz_eabc_fehlerterm_hypothese.md` §7).

**Konsistenzbedingung (Hypothese).** Wenn Pfadbits auf dem Rand und $\Omega_{\mathrm{Hol}}$ im Inneren **gemeinsam** aus einem assoziativen $V_4$-Wert stammen, erwarten wir $P_{\mathrm{same}}^{\mathrm{hol}}\approx 1$. Abweichung misst **projektive Holonomie** / Fehlerterm.

$$\boxed{\;\text{Dreieck (4-Pfad): Bell-Taubenloch exakt} \;\neq\; \text{Viereck (5-Zyklus): globale Holonomie } \Omega_{\mathrm{Hol}}.\;}$$

**Label:** Dreieck/Viereck-Trennung = **Definition**; Konsistenz Pfad↔Holonomie = **Hypothese**.

---

## 9. Vergleich mit $G_E$-Übergangswahrscheinlichkeiten

**Definition.** Zeilenstochastische Matrix $P_{ij}=\Pr[X_{n+1}=j\mid X_n=i]$ aus `transition_counts`.

**Experiment.** `collatz_eabc_bell_inequality_test.py` vergleicht:
- empirische $P_{\mathrm{same}}^{\mathrm{win}}$ auf ABCE-Fenstern;
- $\mathcal{B}_{\mathrm{marg}}$ aus Kantenrändern;
- $P_{ij}$-Vorhersagen für $t$-alignierte Kante $i\to t(i)$.

**Label:** $G_E$-Vergleich = **Experiment**.

---

## 12. Bell/CHSH-Verletzungsanalog ↔ $D_E$-Bias (Brücke)

**Epistemische Abgrenzung:** Dieser Abschnitt ist eine **strukturierte Analogie** zwischen kombinatorischer Bell/CHSH-Logik auf $G_E$ und dem Holonomie-Fehlerterm $D_E$. Es wird **keine** physikalische Nichtlokalität und **kein** Theorem $D_E \Leftrightarrow |S_{\mathrm{EABC}}|$ behauptet.

### 12.1 Faktorisierbarkeit auf $G_E$

**Definition (lokaler Transport, faktorisierbar).** Auf einem gemeinsamen Fenster $W_n^{(4)}$ seien Observablen Funktionen **eines** versteckten Zustands $\lambda_n$ (das Fenster selbst). Dann gilt:
- Bell-Taubenloch: $\mathcal{B}_{\mathrm{win}}\ge 1$ (**Theorem**, §4.1);
- CHSH auf gemeinsamem Träger §7.2: $|S_{\mathrm{EABC}}|\le 2$ unter **LHV-Faktorisierung** (**Analogie** / Skizze §7.5).

**Interpretation.** Wenn die klassische Bell-/CHSH-Schranke auf **gemeinsamem** Träger erfüllt ist, ist der Transport auf $G_E$ in diesem Fenster **faktorisierbar** — alle binären Lesarten stammen aus derselben lokalen Realität $\lambda_n$.

**Label:** Faktorisierbarkeit = **Definition**; Schranken auf gemeinsamem Träger = **Theorem** (Taubenloch) bzw. **Analogie** (CHSH).

### 12.2 Persistenter $D_E$-Bias = nicht-faktorisierbare Holonomie-Reste

**Definition.** $D_E(X)=N_+(X)-N_-(X)$ misst die **absolute** Asymmetrie der geschlossenen Orientierungen ABCEA vs. CEABC (`collatz_eabc_fehlerterm_hypothese.md`).

**Hauptvermutung.** $\mathrm{Hol}_E=\lim_{X\to\infty} D_E/(N_++N_-)=0$ — der **normierte** Hauptterm verschwindet.

**Analogie zu $|S|>2$.** Endliches $D_E(X)\neq 0$ bei kleinem $|\chi_{\mathrm{Hol}}(X)|$ ist ein **arithmetisches Korrelations-Restsignal**: wie eine CHSH-„Verletzung'' entsteht es nicht aus punktweiser Taubenloch-Verletzung, sondern aus **nicht-faktorisierbarer** Mischung von Pfad-, Kanten- und Zyklus-Lesarten über die Primfolge — ein **Holonomie-Fehlerterm**, kein Beweis für $\mathrm{Hol}_E\neq 0$.

$$\boxed{\;D_E(X)\neq 0\;\text{ persistent}\;\leadsto\;\text{„non-factorizable'' Holonomie-Residuum (Analog zu }|S|>2\text{), nicht QM.}\;}$$

**Label:** $D_E\neq 0$ als Restsignal = **Hypothese** / **Experiment**; Gleichsetzung mit CHSH-Verletzung = **Analogie**.

### 12.3 Explizite Zuordnungstabelle

| Bell/CHSH-Größe | $D_E$-Seite | Verhalten (empirisch / theoretisch) | Label |
|-----------------|-------------|-------------------------------------|-------|
| $\mathcal{B}_{\mathrm{win}}\ge 1$ auf ABCE | Taubenloch unabhängig von $D_E$ | immer erfüllt | **Theorem** |
| $\mathcal{B}_{\mathrm{marg}}<1$ | marginale Ränder vs. Fenster-Holonomie | Kontextvermischung | **Experiment** |
| $P_{\mathrm{same}}^{\mathrm{hol}}$ drift | $\chi_{\mathrm{Hol}}(X)=D_E/(N_++N_-)$ **beschränkt** | Hauptterm $\to 0$ | **Hypothese** |
| $|S_{\mathrm{EABC}}|$ vs. $2$ | $|\widetilde{D}_E(X)|=|D_E|/\sqrt{N_++N_-}$ | Oszillation / Bias-Struktur | **Experiment** |
| $|S_{\mathrm{EABC}}|>2$ auf gemeinsamem Träger | Pfad↔Holonomie nicht faktorisierbar | strukturierter Bias | **Hypothese** |
| $|S|>2\sqrt2$ (QM) | — | **nicht** EABC-These | **nicht behauptet** |

**Lesart der Normalisierung.** $\widetilde{D}_E$ und $|S_{\mathrm{EABC}}|$ sind beides **skaleninvariante** Korrelationssignale auf wachsendem $X$: $\chi_{\mathrm{Hol}}$ bleibt klein (Hauptterm), während $D_E$ wachsen kann (Chebyshev-artiger Fehlerterm) — parallel zu $|S|>2$ bei erfülltem Taubenloch auf Einzelfenstern.

### 12.4 Lean- und Python-Schnittstelle

| Artefakt | Inhalt | Status |
|----------|--------|--------|
| `HolonomieFehlerterm.lean` | `wordABCEA`/`wordCEABC`, Lücken $(2,4,2,4)$, `N_plus`/`N_minus`/`D_E` auf Listen, Taubenloch | kombinatorisch **bewiesen** |
| `HolonomieFehlerterm.lean` | `N_plus_up_to`, `Hol_E_zero`, `chsh_lhv_bound_skel` | **`sorry`** / Skeleton |
| `collatz_eabc_bell_inequality_test.py` | `de_bell_combined_report(X)` — $S_{\mathrm{EABC}}$, $D_E$, $\widetilde{D}_E$ am selben $X$ | **Experiment** |
| `collatz_eabc_holonomie_fehlerterm.py` | gleiche kombinierte Zeile in `run_series` | **Experiment** |

**Label:** Lean/Python-Schnittstelle = **Definition** (operativ) + **Experiment** (Numerik).

---

## 13. Dokumentenverknüpfung (aktualisiert)

```
collatz_eabc_transport.md
        │
        ▼
collatz_eabc_zyklus_holonomie.md ──► collatz_eabc_fehlerterm_hypothese.md
        │                                      │
        │                                      │ D_E, widetilde{D}_E
        ▼                                      ▼
collatz_eabc_bell_holonomie.md ◄──── §12 Brücke Bell↔D_E
        │
        ├── collatz_eabc_bell_inequality_test.py  (S_EABC + D_E combined)
        ├── collatz_eabc_holonomie_fehlerterm.py
        └── CollatzEabc/HolonomieFehlerterm.lean
```

---

## 14. Epistemische Tabelle (erweitert)

| Aussage | Label |
|---------|-------|
| Bell-Tripel-Ungleichung (Referenz) | **Analogie** |
| $\sigma(n)$, $O_E,O_A,O_C$ auf ABCE-Fenster | **Definition** |
| $P_{\mathrm{same}}^{\mathrm{win}}$, $\mathcal{B}_{\mathrm{win}}$ | **Definition** |
| $\mathcal{B}_{\mathrm{win}}\ge 1$ | **Theorem** |
| $P_{\mathrm{same}}^{\mathrm{marg}}$, $\mathcal{B}_{\mathrm{marg}}$ | **Definition** / **Experiment** |
| $P_{\mathrm{same}}^{\mathrm{hol}}$, Pfad↔Holonomie | **Definition** / **Hypothese** |
| CHSH-Analog $S_{\mathrm{EABC}}$, $E^{\mathrm{EABC}}$ | **Definition** / **Analogie** |
| $|S_{\mathrm{EABC}}|>2$ auf gemeinsamem ABCE-Träger | **Experiment** / **Hypothese** |
| $|S|=2\sqrt2$ (QM) | **Referenz**, nicht EABC-These |
| $D_E$-Bias ↔ nicht-faktorisierbare Holonomie | **Analogie** / **Hypothese** |
| $\widetilde{D}_E$ ↔ $|S_{\mathrm{EABC}}|$ Skala | **Experiment** |
| $\chi_{\mathrm{Hol}}$ beschränkt bei wachsendem $D_E$ | **Hypothese** / **Experiment** |
| QM-Verschränkung / Nichtlokalität | **nicht behauptet** |
| Bezug $D_E$, $\chi_{\mathrm{Hol}}$ | **Hypothese** |
| `HolonomieFehlerterm.lean` Lückenmuster | **Theorem** |
| `HolonomieFehlerterm.lean` Prime/CHSH | **Skeleton** (`sorry`) |

---

*Kanonsiche Lesart: Die Bell-Ungleichung liefert auf EABC ein **Konsistenzkriterium für drei binäre Lesarten auf demselben Transportfenster** (Taubenloch, immer erfüllt). „Verletzungen'' entstehen nur bei **marginalisierten oder kontextvermischten** Observablen und werden mit **Holonomie-Fehlerterm** $D_E$ / $\chi_{\mathrm{Hol}}$ in Verbindung gebracht — kombinatorisch-statistisch, nicht quantenphysikalisch.*
