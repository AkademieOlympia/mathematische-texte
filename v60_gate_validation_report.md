# v60 Gate Validation Report

## Abschließendes Kapitel: v60-Zertifizierbarkeitssatz
### *Theorema de Coherentia Geometriae* \(\mathfrak G_{60}\)

> **Epistemischer Status.** Der folgende „Satz“ ist ein **Zertifizierbarkeitssatz innerhalb der v60-Spezifikation**: Er definiert notwendige und gemeinsam hinreichende PASS/FAIL-Bedingungen für die *v60-Zertifizierung*. Er ist **kein Beweis**, dass ein konkreter EABC-Laboraufbau diese Bedingungen bereits erfüllt, und kein Beweis der offenen Brückensätze \(WT_1,WT_2,G6\).

Ein physikalischer EABC-v60-Demonstrator gilt genau dann **im Sinne der v60-Spezifikation** als topologisch, thermisch und metrologisch zertifiziert, wenn

\[
\boxed{
\mathfrak G_{60}
=
\prod_{j\in\{\mathrm{mode},\eta,T,\mathrm{coh},\mathrm{EP},\mathrm{QFI}\}}G_j
=1,
\qquad G_j\in\{0,1\}.
}
\]

Äquivalent:

\[
\boxed{
\mathfrak G_{60}=1
\iff
G_{\rm mode}\land G_\eta\land G_T\land
G_{\rm coh}\land G_{\rm EP}\land G_{\rm QFI}.
}
\]

Die Produktdarstellung kodiert eine **nicht-kompensatorische Konjunktion**: Ein einziges FAIL erzwingt \(\mathfrak G_{60}=0\).

---

## I. Prinzip der Nicht-Kompensation

Für jedes \(j\) gilt

\[
G_j=0\Longrightarrow\mathfrak G_{60}=0.
\]

Es existiert innerhalb des Zertifikators keine gewichtete Ersatzmetrik, durch die ein FAIL eines Gates durch Übererfüllung eines anderen Gates kompensiert werden könnte.

### 1. Verlust-Gate

Für beobachtetes Squeezing mit Gesamtwirkungsgrad \(\eta_{\rm tot}\) gilt im idealisierten Verlustmodell

\[
V_{\rm obs}=(1-\eta_{\rm tot})+\eta_{\rm tot}e^{-2r}.
\]

Daraus folgt als notwendige Bedingung für ein Ziel \(S_{\rm target}\):

\[
\boxed{
1-\eta_{\rm tot}<10^{-S_{\rm target}/10}.
}
\]

Also

\[
G_\eta=
\mathbf 1\!\left[
1-\eta_{\rm tot}<10^{-S_{\rm target}/10}
\right].
\]

Für

\[
S_{\rm target}=19.9\,{\rm dB}
\]

ist

\[
10^{-1.99}\simeq0.01023,
\]

also

\[
\boxed{\eta_{\rm tot}>0.98977}
\]

bzw. ein Gesamtverlust von weniger als ungefähr \(1.023\%\).

Ein PASS von \(G_{\rm mode}\) ersetzt dieses Gesamtverlustkriterium nicht.

### 2. Moden- und Phasenfehler sind getrennte Flanken

Das Moden-Gate kontrolliert räumliche Überlappung, beispielsweise über

\[
\eta_{\rm mode}
=
|\langle u_{\rm LO},u_{\rm sq}\rangle|^2.
\]

Phasenjitter ist davon logisch getrennt. Für einen kleinen Quadraturfehler \(\delta\phi\) gilt schematisch

\[
V_{\rm meas}
=
V_X\cos^2\delta\phi
+
V_P\sin^2\delta\phi,
\]

so dass bei starkem Squeezing wegen \(V_P\gg V_X\) bereits kleine Phasenfehler relevant werden können.

**Wichtig:** Ein FAIL von \(G_T\) ist nicht ohne zusätzlichen dynamischen Kopplungssatz automatisch identisch mit Phasenmischung des antisqueezed Quadraturrauschens. Poisoning und Quadratur-Phasenjitter sind zunächst verschiedene Fehlerkanäle. Eine EABC-Identifikation beider Mechanismen benötigt einen expliziten Hamilton-/Liouvillian-Kopplungssatz.

---

## II. Verallgemeinertes thermisches Kontrolltor \(G_T\)

Unter der Modellannahme eines Poisson-Prozesses mit thermisch aktivierter Rate

\[
\Gamma_{\rm poison}(T)
=
\Gamma_0
\exp\!\left(
-\frac{\Delta_{\rm eff}}{k_BT}
\right)
\]

ist die Wahrscheinlichkeit mindestens eines Poisoning-Ereignisses während einer Gatterzeit \(t_g\)

\[
\boxed{
P_{\rm poison}(T,t_g)
=
1-
\exp\!\left[
-\Gamma_0t_g
e^{-\Delta_{\rm eff}/(k_BT)}
\right].
}
\]

Das Gate wird definiert durch

\[
\boxed{
G_T(T_{\rm op})
=
\mathbf1\!\left[
P_{\rm poison}(T_{\rm op},t_g)
<
\varepsilon_{\rm th}
\right].
}
\]

### Exakte Inversion innerhalb des Arrhenius-Poisson-Modells

Aus

\[
1-e^{-x}<\varepsilon_{\rm th},
\qquad
x=
\Gamma_0t_g e^{-\Delta_{\rm eff}/(k_BT)}
\]

folgt

\[
x<-\ln(1-\varepsilon_{\rm th}).
\]

Damit erhält man, sofern

\[
\Gamma_0t_g>
-\ln(1-\varepsilon_{\rm th}),
\]

die Temperaturbedingung

\[
\boxed{
T_{\rm op}
<
T_{\max}
=
\frac{\Delta_{\rm eff}}
{k_B
\ln\!\left[
\dfrac{\Gamma_0t_g}
{-\ln(1-\varepsilon_{\rm th})}
\right]}.
}
\]

Dies korrigiert die Notation der Grenzformel: Im Nenner steht
\(-\ln(1-\varepsilon_{\rm th})\), äquivalent zu
\(\ln[1/(1-\varepsilon_{\rm th})]\).

### Kleinfehlergrenze

Für

\[
x=\Gamma_{\rm poison}t_g\ll1
\]

gilt

\[
1-e^{-x}=x+O(x^2),
\]

und damit

\[
\boxed{
\Gamma_0t_g
e^{-\Delta_{\rm eff}/(k_BT_{\rm op})}
<
\varepsilon_{\rm th}.
}
\]

Für \(\varepsilon_{\rm th}\ll1\) gilt außerdem

\[
-\ln(1-\varepsilon_{\rm th})
=
\varepsilon_{\rm th}+O(\varepsilon_{\rm th}^2),
\]

wodurch die frühere Näherungsform von \(T_{\max}\) zurückgewonnen wird.

### Aktivierungs-Brücke

Die Größe \(\Delta_{\rm eff}\) ist die **für den tatsächlichen Poisoning-Kanal relevante effektive Aktivierungsenergie**. Eine gemessene Transport-Aktivierungslücke eines \(\nu=5/2\)-Zustands darf nicht ohne zusätzlichen experimentellen oder mikroskopischen Nachweis mit \(\Delta_{\rm eff}\) gleichgesetzt werden.

Insbesondere ist die Brücke

\[
\boxed{
\Delta_{\rm transport}
\stackrel{?}{=}
\Delta_{\rm pair}
\stackrel{?}{=}
\Delta_{\rm eff}
}
\]

offen und experimentell zu bestimmen.

---

## III. Die sechs v60-Gates

### \(G_{\rm mode}\): räumliche Modenkontrolle

\[
G_{\rm mode}
=
\mathbf1[L_{\rm spatial}\le L_{\rm spatial}^{\max}].
\]

Ein Wert wie \(0.4\%\) ist nur dann ein zertifizierter Grenzwert, wenn er für den konkreten Aufbau gemessen und als Spezifikation festgelegt wurde.

### \(G_\eta\): Gesamtverlust

\[
\boxed{
G_\eta
=
\mathbf1\!\left[
\eta_{\rm tot}>
1-10^{-S_{\rm target}/10}
\right].
}
\]

### \(G_T\): thermisches Poisoning

\[
\boxed{
G_T
=
\mathbf1[
P_{\rm poison}(T_{\rm op},t_g)<\varepsilon_{\rm th}
].
}
\]

### \(G_{\rm coh}\): offene Quantenkohärenz

\[
\boxed{
G_{\rm coh}
=
\mathbf1[C_q>1]
}
\]

für die im konkreten Liouvillian-Modell definierte effektive Cooperativity.

### \(G_{\rm EP}\): echte EP2-Zertifizierung

Ein EP2-PASS erfordert mindestens

\[
\lambda_1=\lambda_2=\lambda_{\rm EP},
\]

\[
\operatorname{algmult}(\lambda_{\rm EP})=2,
\qquad
\operatorname{geomult}(\lambda_{\rm EP})=1,
\]

sowie den separat definierten topologischen Umlaufindex

\[
\nu_{\rm EP}=1.
\]

Daher

\[
\boxed{
G_{\rm EP}
=
\mathbf1[
\text{EP2-Defektivität und geforderte Monodromie nachgewiesen}
].
}
\]

### \(G_{\rm QFI}\): ressourcenbereinigter Informationsgewinn

\[
\boxed{
G_{\rm QFI}
=
\mathbf1[
F_Q^{\rm full}>F_Q^{\rm ref}
]
}
\]

wobei Verluste, Postselektion, Detektion, Präparation, Messzeit und sonstige relevante Ressourcen in beiden Protokollen konsistent bilanziert werden müssen.

---

## IV. Epistemische Schichtung

\[
\boxed{
\begin{array}{ccl}
\mathbf{Level\ A}&:&
\text{exakte algebraische, gittertheoretische und funktionalanalytische Kerne},\\[2mm]
\mathbf{Level\ B}&:&
\text{numerisch oder experimentell bestimmte Gate-Parameter},\\[2mm]
\mathbf{Level\ C}&:&
\mathfrak G_{60}=1
\text{ für einen konkreten Laboraufbau},\\[2mm]
\mathbf{Level\ D}&:&
WT_1\land WT_2\land G6
\text{ als separat zu beweisende analytisch-spektrale Brücken},\\[2mm]
\mathbf{Level\ E}&:&
\text{übergeordnete EABC-/Hilbert--Pólya-Interpretation}.
\end{array}
}
\]

Für Level A ist entscheidend, nur tatsächlich bewiesene mathematische Aussagen einzuschließen. QRM-, Barnes-Wall-, Yang-Baxter-, Viazovska- oder de-Branges-Strukturen gehören nur mit den jeweils explizit bewiesenen Relationen in diesen Level; ihre physikalische Identifikation mit dem v60-Demonstrator ist eine zusätzliche Brücke.

Die zentrale Nichtimplikation lautet

\[
\boxed{
\mathfrak G_{60}=1
\;\not\Rightarrow\;
WT_1\land WT_2\land G6.
}
\]

Ebenso gilt ohne zusätzlichen Brückensatz nicht automatisch

\[
WT_1\land WT_2\land G6
\Rightarrow
\mathfrak G_{60}=1.
\]

Die mathematische und die experimentelle Beweislast bleiben getrennt.

---

## V. Bedeutung eines PASS

Ein experimenteller Nachweis

\[
\mathfrak G_{60}=1
\]

belegt **per Definition der v60-Zertifizierung**, dass der konkrete Aufbau sämtliche spezifizierten Moden-, Verlust-, Temperatur-, Kohärenz-, EP- und QFI-Grenzen gleichzeitig erfüllt.

Er beweist damit nicht automatisch:

- universelles topologisches Quantencomputing,
- eine konkrete Moore–Read-Identifikation des Trägerzustands,
- die Barnes-Wall- oder Viazovska-Brücke,
- \(WT_1\), \(WT_2\) oder \(G6\),
- eine Hilbert–Pólya-Realisierung,
- oder die Riemannsche Vermutung.

Ein PASS ist daher ein starkes physikalisches Zertifikat **innerhalb des definierten v60-Protokolls**, kein Beweis des gesamten EABC-Überbaus.

---

## VI. Bedeutung eines FAIL

Für

\[
\mathfrak G_{60}=0
\]

wird der Fehlervektor

\[
\boxed{
\mathbf f_{60}
=
(1-G_{\rm mode},
1-G_\eta,
1-G_T,
1-G_{\rm coh},
1-G_{\rm EP},
1-G_{\rm QFI})
}
\]

eingeführt.

Er lokalisiert unmittelbar, welche Zertifizierungsbedingungen nicht erfüllt wurden.

Dabei muss zwischen einem **FAIL des Demonstrators** und einer **Falsifikation einer mathematischen Hypothese** unterschieden werden.

Beispielsweise:

\[
G_{\rm EP}=0
\]

zeigt zunächst, dass die geforderte EP2-Struktur im untersuchten Aufbau und Parameterbereich nicht nachgewiesen wurde. Erst wenn eine mathematische Hypothese exakt behauptet, dass gerade dieser wohldefinierte Operator bei genau diesen Parametern zwingend eine EP2 besitzen muss, kann das Ergebnis diese Hypothese falsifizieren.

Analog bedeutet

\[
G_T=0
\]

zunächst, dass das thermische Fehlerbudget des konkreten Aufbaus überschritten wird. Daraus folgt nicht ohne zusätzliche Ableitung unmittelbar das Scheitern eines „16D-Freeze-Outs“.

---

## VII. Nicht-Kompensationssatz

**Satz (v60-Nicht-Kompensation).**  
Seien \(G_1,\ldots,G_6\in\{0,1\}\) und

\[
\mathfrak G_{60}=\prod_{j=1}^{6}G_j.
\]

Dann gilt

\[
\boxed{
\mathfrak G_{60}=1
\iff
G_j=1\quad\forall j.
}
\]

**Beweis.**  
Da jedes \(G_j\) nur die Werte \(0\) oder \(1\) annimmt, ist das Produkt genau dann \(1\), wenn jeder Faktor \(1\) ist. Enthält das Produkt mindestens einen Faktor \(0\), ist es \(0\). \(\square\)

Dieser algebraisch triviale Satz ist methodisch wesentlich: Er verhindert eine nachträgliche gewichtete Kompensation physikalisch unabhängiger Ausschlusskriterien.

---

## VIII. Epistemische Statushierarchie der topologischen Qubit-Architektur

### 1. Geltungsbereich: Clifford-Sektor statt Universalitätsbehauptung

Die v60-Architektur wird für den Ising-/Moore–Read-Zweig ausdrücklich als **vorgeschlagene topologische Qubit-Architektur im Clifford-Sektor** verstanden. Die logische Braid-Darstellung darf daher nicht als Nachweis universeller Quantenberechnung gelesen werden.

Für den algebraischen Ising-Sektor gilt innerhalb der verwendeten Konvention:

\[
\boxed{
\text{Ising-Braids}\longrightarrow\text{Clifford-Sektor},
\qquad
\text{Ising-Braids}\not\Rightarrow\text{Universalität}.
}
\]

Insbesondere können Pauli-, Hadamard- und geeignete Clifford-Verschränkungsoperationen durch Braidfolgen repräsentiert werden. Eine universelle Architektur benötigt darüber hinaus eine explizit spezifizierte Nicht-Clifford-Ressource, etwa ein \(T\)-Gate bzw. eine entsprechende Zustandsinjektion und Fehlerreduktion. Diese Zusatzschicht ist **nicht** Bestandteil des rein topologisch geschützten Ising-Braid-Sektors.

Ebenso wird kein absoluter Hardware-Schutz behauptet. Topologische Pfadinvarianz schützt nur innerhalb eines lückenerhaltenden Fusionssektors gegen die dafür zulässigen homotopisch trivialen Pfadstörungen. Poisoning, Gap-Schließung, Nichtadiabatik, endliche Anyon-Überlappung sowie Mess- und Initialisierungsfehler bleiben eigenständige physikalische Fehlerkanäle.

### 2. Status der fünf vertikalen Übergänge

Die schematische v60-Kette

\[
\text{GaAs}/\nu=5/2
\dashrightarrow
\text{Moore--Read/Ising}
\dashrightarrow
\Lambda_{\rm BW16}
\longrightarrow
\mathrm{QRM}(1,4)
\dashrightarrow
\text{Braid-Compiler}
\longrightarrow
\text{klassischer Audit}
\]

ist **kein einheitlicher Beweispfeil**. Jeder Übergang besitzt einen eigenen epistemischen Status:

1. **GaAs / \(\nu=5/2\) \(\dashrightarrow\) Moore–Read-/Ising-Sektor — Physical Bridge Conjecture.**  
   Die v60-Spezifikation setzt nicht voraus, dass die Moore–Read-/Ising-Identifikation eines konkreten GaAs-Trägerzustands bereits bewiesen ist. Sie ist experimentell separat zu zertifizieren. Ein PASS der übrigen v60-Gates ersetzt diesen Nachweis nicht.

2. **Moore–Read-/Ising-Sektor \(\dashrightarrow\Lambda_{\rm BW16}\) — EABC Bridge Hypothesis.**  
   Die Identifikation bzw. kontrollierte Repräsentation des physikalischen Fusionsraums durch einen Barnes–Wall-/16D-Gittersektor ist eine zusätzliche EABC-spezifische Brücke. Die Existenz des mathematischen Gitters allein beweist diese physikalische Kopplung nicht.

3. **\(\Lambda_{\rm BW16}\longrightarrow\mathrm{QRM}(1,4)\) — algebraischer Korrespondenzsektor, nur soweit explizit konstruiert.**  
   Barnes–Wall- und Reed–Muller-Strukturen dürfen auf Level A ausschließlich mit den im Modell tatsächlich hergeleiteten Code-/Gitterrelationen verwendet werden. Eine solche endliche algebraische Korrespondenz impliziert für sich genommen weder eine Moore–Read-Realisierung noch einen physikalischen Schutzmechanismus.

4. **\(\mathrm{QRM}(1,4)\dashrightarrow\) Braid-Compiler / Clifford-Sektor — Proposed Architectural Mapping.**  
   Die Compiler-Schicht spezifiziert, wie logische Operationen auf Braidfolgen abgebildet werden sollen. Für die verwendete Ising-Konvention müssen die \(F\)- und \(R\)-Matrizen sowie die resultierenden Braid-Identitäten explizit fixiert werden. Die Abbildung ist eine Architekturvorgabe und kein automatisches Korollar der Codegeometrie.

5. **Braid-/Hardware-Sektor \(\longrightarrow\) klassische Validierung — Numerical / Algorithmic Audit.**  
   Simulatoren, QFI-Auswertungen und der v60-Gate-Zertifikator prüfen die spezifizierten Modellbedingungen. Ein algorithmisches PASS ist ein Konsistenz- bzw. Zertifizierungsergebnis innerhalb der definierten Spezifikation; es ist weder ein physikalischer Existenzbeweis noch ein analytischer Beweis offener EABC-Brückensätze.

### 3. Nichtimplikationen

Die Statushierarchie wird durch folgende Nichtimplikationen eingefroren:

\[
\boxed{
\begin{aligned}
\text{Clifford-Braid-Darstellung}
&\not\Rightarrow \text{universelles Quantencomputing},\\
\Lambda_{\rm BW16}\text{-Algebra}
&\not\Rightarrow \text{physikalische Moore--Read-Kopplung},\\
\mathrm{QRM}(1,4)\text{-Codekorrespondenz}
&\not\Rightarrow \text{fehlerfreie Braid-Hardware},\\
\mathfrak G_{60}=1
&\not\Rightarrow WT_1\land WT_2\land G6.
\end{aligned}
}
\]

Umgekehrt falsifiziert ein experimentelles FAIL eines konkreten Demonstrators nicht automatisch die zugrunde liegende endliche Algebra. Es falsifiziert zunächst die betreffende operationalisierte Realisierungs- oder Zertifizierungsbedingung.

### 4. Epistemischer Freeze der Architekturabbildung

Für Abbildungen und Übersichten lautet der bevorzugte Titel daher:

\[
\boxed{\text{v60 Proposed Topological-Qubit Architecture (Clifford Sector)}}.
\]

Die graphische Darstellung ist als **Architektur- und Prüfdiagramm**, nicht als Darstellung einer bereits realisierten Hardwarekette zu lesen. Durchgezogene Pfeile dürfen nur für im jeweiligen Abschnitt explizit etablierte Konstruktionen oder Audit-Abbildungen verwendet werden; gestrichelte Pfeile markieren offene physikalische oder EABC-spezifische Brücken.

Der methodische Kernsatz lautet:

\[
\boxed{
\text{mathematische Konstruktion}
\neq
\text{physikalische Identifikation}
\neq
\text{experimentelle Zertifizierung}.
}
\]

---

## IX. v60-Freeze

Der experimentelle v60-Sektor wird abschließend durch folgende Aussage eingefroren:

\[
\boxed{
\begin{gathered}
\mathfrak G_{60}
=
G_{\rm mode}G_\eta G_TG_{\rm coh}G_{\rm EP}G_{\rm QFI},\\[1mm]
\mathfrak G_{60}=1
\iff
\text{alle sechs spezifizierten physikalischen Zertifizierungstore PASS}.
\end{gathered}
}
\]

Der Zertifikator ist nicht-kompensatorisch und falsifizierbar. Seine einzelnen Gates müssen mit vorab festgelegten Messgrößen, Schwellenwerten, Unsicherheiten und Ressourcenbilanzen operationalisiert werden.

\[
\boxed{
\text{Spezifikation}
\neq
\text{experimentelle Verifikation}
\neq
\text{analytischer Brückenbeweis}.
}
\]

Insbesondere bleibt

\[
\boxed{
\mathfrak G_{60}=1
\not\Rightarrow
WT_1\land WT_2\land G6.
}
\]

Damit ist v60 als **Zertifizierungsarchitektur** formal geschlossen. Die tatsächliche Zertifizierung eines Demonstrators beginnt erst mit der experimentellen Bestimmung aller sechs Gates.
