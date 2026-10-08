# v60 Anyon Braid Compiler — korrigierter Freeze

## Status

Dieses Dokument fixiert den **v60-Standard** für den Ising-Anyon-/Braid-Sektor des EABC-Resonatorprogramms. Es trennt strikt zwischen:

- exakten algebraischen Aussagen im logischen Fusionsraum,
- topologisch geschützten Pfadklassen,
- physikalischen Fehlerkanälen realer Hardware,
- offenen EABC-Brücken zu GaAs-, \(\nu=5/2\)- und Resonatorphysik.

Es wird **keine absolute Fehlerfreiheit** behauptet. Der Schutz gilt gegen homotopisch triviale Pfaddeformationen innerhalb eines lückenerhaltenden topologischen Sektors; Quasiteilchen-Poisoning, Gap-Schließung, endliche Anyon-Überlappung und Nichtadiabatik bleiben reale Fehlerquellen.

---

## 1. Logischer Ising-Anyon-Sektor

Für vier Ising-Anyonen mit festem Gesamtfusionskanal wird ein zweidimensionaler logischer Raum

\[
\mathcal H_L=\operatorname{span}\{|0_L\rangle,|1_L\rangle\}
\]

verwendet.

Bis auf eine globale Phase können die elementaren Braidgeneratoren in einer Standardkonvention als

\[
B_1=
\begin{pmatrix}
1&0\\
0&i
\end{pmatrix},
\]

und

\[
B_2=F B_1 F,
\qquad
F=\frac1{\sqrt2}
\begin{pmatrix}
1&1\\
1&-1
\end{pmatrix}
\]

geschrieben werden.

Globale Phasen sind für die logische Gatterwirkung irrelevant; daher wird im Folgenden die Äquivalenz

\[
U\sim e^{i\alpha}U
\]

verwendet.

---

## 2. Korrigierte Braid-Identitäten

Aus den Generatoren folgt

\[
B_1^2\sim Z,
\]

\[
B_2^2\sim X.
\]

Damit ergibt sich

\[
B_1^2B_2^2B_1^2
\sim
ZXZ
=
-X
\sim X.
\]

Die vorherige Identifikation

\[
\tau_1^2\tau_2^2\tau_1^2\sim H
\]

ist daher falsch.

Der v60-Freeze lautet:

\[
\boxed{
\tau_1^2\tau_2^2\tau_1^2\sim X.
}
\]

Für das Hadamard-Gatter gilt in derselben Konvention

\[
\boxed{
H\sim\tau_1\tau_2\tau_1,
}
\]

beziehungsweise bei vertauschter Braid-Konvention äquivalent

\[
H\sim\tau_2\tau_1\tau_2.
\]

Damit reduziert sich die minimale Hadamard-Braidfolge in diesem Modell von sechs elementaren Vertauschungen auf drei.

Weitere Clifford-Identitäten können entsprechend als Produkte der Braidgeneratoren kompiliert werden. Die exakte Sequenz hängt von der gewählten Basis- und Orientierungs-Konvention ab; jede Implementierung muss daher die verwendeten \(F\)- und \(R\)-Matrizen explizit angeben.

---

## 3. Clifford-Sektor und Universalitätsgrenze

Ising-Anyon-Braiding erzeugt einen Clifford-Sektor. Damit sind insbesondere Pauli-Operationen, Hadamard-Operationen und geeignete Clifford-Verschränkungen zugänglich.

Dies impliziert jedoch keine universelle Quantenberechnung.

Für Universalität ist eine zusätzliche Nicht-Clifford-Ressource erforderlich, beispielsweise:

- ein \(T\)-Gate,
- Magic-State-Injektion,
- eine nichttopologische kontrollierte Phase,
- oder eine andere explizit spezifizierte Zusatzressource.

Daher lautet der epistemisch korrekte Satz:

\[
\boxed{
\text{Ising-Braids liefern einen topologisch strukturierten Clifford-Sektor, aber keine vollständige Universalität.}
}
\]

Für metrologische Protokolle kann dieser Clifford-Sektor dennoch ausreichend sein, sofern die benötigte Zustandspräparation, Interferometrie und Auslese vollständig innerhalb dieses Sektors realisierbar ist.

---

## 4. Topologischer Schutz: korrekte Aussage

Die Schutzbehauptung des v60-Standards lautet:

\[
\boxed{
\text{Braiding schützt die logische Operation gegen homotopisch triviale Pfadstörungen,}
}
\]

solange

\[
\Delta_{\mathrm{top}}>0
\]

entlang der gesamten Deformation offen bleibt und kein unerwünschtes Anyon den logischen Fusionssektor verändert.

Nicht geschützt sind insbesondere:

\[
\boxed{
\text{Poisoning},\quad
\text{Gap-Schließung},\quad
\text{Nichtadiabatik},\quad
\text{endliche Anyon-Überlappung},\quad
\text{Mess- und Initialisierungsfehler}.
}
\]

Topologische Invarianz bedeutet daher nicht

\[
F_{\mathrm{gate}}=1
\]

für reale Hardware.

---

## 5. Realistisches Fehlerbudget

Für ein einzelnes logisches Gatter wird die Infidelity in erster Ordnung durch

\[
\boxed{
1-F_{\mathrm{gate}}
\lesssim
\epsilon_{\mathrm{topo}}
+
\epsilon_{\mathrm{poison}}
+
\epsilon_{\mathrm{dyn}}
+
\epsilon_{\mathrm{meas}}
}
\]

modelliert.

### 5.1 Endliche topologische Überlappung

Für räumlich getrennte Anyonen mit charakteristischer Distanz \(L\) und Korrelationslänge \(\xi\) wird als asymptotische Zielstruktur

\[
\boxed{
\epsilon_{\mathrm{topo}}
\le
C_{\mathrm{topo}} e^{-L/\xi}
}
\]

verwendet.

Dies beschreibt die Tatsache, dass die logische Ununterscheidbarkeit bei endlicher Separation nur exponentiell genau ist.

### 5.2 Quasiteilchen-Poisoning

Sei \(\Delta_{\mathrm{pair}}\) die minimale Energie zur Erzeugung eines für den logischen Sektor schädlichen Quasiteilchenpaares.

Als einfaches thermisch aktiviertes Modell wird

\[
\boxed{
\Gamma_{\mathrm{poison}}(T)
=
\Gamma_0
\exp\!\left(
-\frac{\Delta_{\mathrm{pair}}}{k_B T}
\right)
}
\]

angesetzt.

Für eine effektive Gatterdauer \(t_g\) lautet die Wahrscheinlichkeit mindestens eines Poisoning-Ereignisses

\[
\boxed{
P_{\mathrm{poison}}(T,t_g)
=
1-
\exp\!\left[
-\Gamma_{\mathrm{poison}}(T)t_g
\right].
}
\]

Im schwachen Fehlerregime

\[
\Gamma_{\mathrm{poison}}t_g\ll1
\]

folgt

\[
\boxed{
P_{\mathrm{poison}}
\simeq
\Gamma_0 t_g
e^{-\Delta_{\mathrm{pair}}/(k_B T)}.
}
\]

Damit kann

\[
\epsilon_{\mathrm{poison}}
\approx
P_{\mathrm{poison}}
\]

als erste Näherung verwendet werden.

### 5.3 Dynamische / nichtadiabatische Fehler

Für eine minimale dynamische Gap-Skala \(\Delta_{\min}\) und eine zeitabhängige Kontrollbahn \(\lambda(t)\) muss ein Adiabatik-Kriterium separat nachgewiesen werden.

Schematisch:

\[
\boxed{
\epsilon_{\mathrm{dyn}}
\lesssim
C_{\mathrm{ad}}
\max_t
\frac{
|\langle m(t)|\dot H(t)|n(t)\rangle|^2
}{
\Delta_{\min}(t)^4
}
}
\]

in einer geeigneten, modellabhängigen Form.

Die konkrete Schranke ist aus dem tatsächlichen Hamiltonoperator herzuleiten.

---

## 6. Fidelity-Untergrenze

Unter Addition der führenden Fehlerbeiträge erhält man die operative v60-Untergrenze

\[
\boxed{
F_{\mathrm{gate}}(T,L,t_g)
\gtrsim
1
-
C_{\mathrm{topo}}e^{-L/\xi}
-
\Gamma_0 t_g
e^{-\Delta_{\mathrm{pair}}/(k_BT)}
-
\epsilon_{\mathrm{dyn}}
-
\epsilon_{\mathrm{meas}}.
}
\]

Diese Gleichung ist als **physikalische Zielabschätzung** zu verstehen, nicht als universeller Satz ohne konkrete mikroskopische Herleitung.

---

## 7. Arrhenius-Grenztemperatur

Fordere für den thermischen Poisoning-Beitrag

\[
P_{\mathrm{poison}}\le\varepsilon_{\mathrm{th}}.
\]

Im linearen schwachen-Fehler-Regime genügt

\[
\Gamma_0t_g
e^{-\Delta_{\mathrm{pair}}/(k_BT)}
\le
\varepsilon_{\mathrm{th}}.
\]

Falls

\[
\Gamma_0t_g>\varepsilon_{\mathrm{th}},
\]

folgt

\[
\boxed{
T
\le
T_{\max}
=
\frac{
\Delta_{\mathrm{pair}}
}{
k_B
\ln\!\left(
\frac{\Gamma_0t_g}{\varepsilon_{\mathrm{th}}}
\right)
}.
}
\]

Dies ist die falsifizierbare obere Temperaturgrenze des einfachen Arrhenius-Modells.

Die Größe \(T_{\max}\) hängt explizit ab von:

\[
\Delta_{\mathrm{pair}},
\qquad
\Gamma_0,
\qquad
t_g,
\qquad
\varepsilon_{\mathrm{th}}.
\]

Ohne experimentell oder mikroskopisch begründete Werte für diese Größen darf keine konkrete Millikelvin-Zahl als Vorhersage eingefroren werden.

---

## 8. Einfluss der korrigierten Hadamard-Sequenz auf das Fehlerbudget

Sei \(t_{\mathrm{braid}}\) die Dauer eines elementaren Braid-Schritts.

Die alte, fehlerhafte Hadamard-Kompilierung hätte formal

\[
t_H^{(6)}=6t_{\mathrm{braid}}
\]

beansprucht.

Die korrigierte Sequenz benötigt

\[
\boxed{
t_H^{(3)}=3t_{\mathrm{braid}}.
}
\]

Im Arrhenius-Poisoning-Modell halbiert sich damit im kleinen Fehlerregime der führende thermische Beitrag:

\[
P_{\mathrm{poison}}^{(3)}
\simeq
\frac12
P_{\mathrm{poison}}^{(6)}.
\]

Dies ist ein realer Vorteil der kürzeren Compiler-Sequenz, sofern jeder elementare Braid-Schritt vergleichbare Dauer und vergleichbares thermisches Risiko besitzt.

---

## 9. Mehrere Poisoning-Kanäle

Für mehrere unabhängige thermische Kanäle \(a\) mit

\[
\Gamma_a(T)
=
\Gamma_{0,a}
e^{-\Delta_a/(k_BT)}
\]

gilt

\[
\boxed{
\Gamma_{\mathrm{tot}}(T)
=
\sum_a \Gamma_a(T).
}
\]

Dann

\[
P_{\mathrm{poison}}
=
1-e^{-\Gamma_{\mathrm{tot}}t_g}.
\]

Im schwachen Fehlerregime:

\[
P_{\mathrm{poison}}
\simeq
t_g
\sum_a
\Gamma_{0,a}
e^{-\Delta_a/(k_BT)}.
\]

Die niedrigste relevante Aktivierungsenergie dominiert asymptotisch bei tiefen Temperaturen, sofern die Präfaktoren nicht anomal groß oder klein sind.

---

## 10. v60-Gates für eine experimentelle Zertifizierung

Der Anyonen-Sektor gilt erst dann als physikalisch zertifiziert, wenn die folgenden Gates separat bestanden sind:

\[
\boxed{
G_{\mathrm{braid}}:
\text{korrekte Ising-Braid-Darstellung und Compiler-Identitäten}
}
\]

\[
\boxed{
G_{\mathrm{gap}}:
\Delta_{\mathrm{top}}>0
\text{ entlang aller physikalischen Braidpfade}
}
\]

\[
\boxed{
G_{\mathrm{sep}}:
C_{\mathrm{topo}}e^{-L/\xi}
<
\varepsilon_{\mathrm{sep}}
}
\]

\[
\boxed{
G_{\mathrm{poison}}:
P_{\mathrm{poison}}(T,t_g)
<
\varepsilon_{\mathrm{poison}}
}
\]

\[
\boxed{
G_{\mathrm{ad}}:
\epsilon_{\mathrm{dyn}}
<
\varepsilon_{\mathrm{dyn}}
}
\]

\[
\boxed{
G_{\mathrm{read}}:
\epsilon_{\mathrm{meas}}
<
\varepsilon_{\mathrm{meas,max}}
}
\]

und, sofern das Protokoll quantenmetrologische Überlegenheit beansprucht,

\[
\boxed{
G_{\mathrm{QFI}}:
F_Q^{\mathrm{full}}
>
F_Q^{\mathrm{reference}}
}
\]

bei gleicher Ressourcenbilanz.

---

## 11. Epistemischer Status

### Exakt im algebraischen Modell

- Ising-Fusionsraum ist zweidimensional nach Festlegung des Gesamtfusionskanals.
- Die Braidgeneratoren erzeugen einen Clifford-Sektor.
- In der verwendeten Konvention gilt
  \[
  \tau_1^2\sim Z,\qquad \tau_2^2\sim X.
  \]
- Damit
  \[
  \tau_1^2\tau_2^2\tau_1^2\sim X.
  \]
- Das Hadamard-Gatter kann als
  \[
  H\sim\tau_1\tau_2\tau_1
  \]
  dargestellt werden, bis auf globale Phase und Konventionswahl.

### Physikalisch bedingt

- Exponentielle Unterdrückung endlicher Anyon-Überlappung.
- Arrhenius-Form des Poisoning-Kanals.
- konkrete Gap-Größe \(\Delta_{\mathrm{pair}}\).
- konkrete Korrelationslänge \(\xi\).
- konkrete Attempt-Rate \(\Gamma_0\).
- konkrete Braiding-Zeit \(t_g\).
- tatsächliche Fehlerbudgets im GaAs-\(\nu=5/2\)-System.

### EABC-spezifisch offen

- Identifikation des v44/v60-Resonators mit einem physikalisch realisierten Ising-/Moore–Read-Codesektor.
- Kopplung an den Barnes-Wall-/de-Branges-Sektor.
- Ableitung der benötigten Gaps aus der EABC-Gittergeometrie.
- metrologische QFI-Überlegenheit nach vollständiger Verlust- und Noise-Bilanz.

---

## 12. Freeze-Satz

Der v60-Standard lautet:

\[
\boxed{
\text{Topologisches Braiding schützt die logische Gatterklasse gegen homotopisch triviale Pfadstörungen,}
}
\]

aber

\[
\boxed{
\text{reale Gatter-Fidelity bleibt durch endliche Überlappung, Poisoning, Nichtadiabatik und Auslesefehler begrenzt.}
}
\]

Die experimentell relevante Zielgröße ist daher nicht absolute Fehlerfreiheit, sondern eine nachprüfbare Schranke

\[
\boxed{
1-F_{\mathrm{gate}}
\le
\varepsilon_{\mathrm{budget}}
}
\]

mit separat zertifizierten Beiträgen.

Die obere Betriebstemperatur des einfachen thermischen Poisoning-Modells wird durch

\[
\boxed{
T_{\max}
=
\frac{
\Delta_{\mathrm{pair}}
}{
k_B\ln(\Gamma_0t_g/\varepsilon_{\mathrm{th}})
}
}
\]

bestimmt und ist erst nach Einsetzen experimentell begründeter Parameter numerisch zu fixieren.
