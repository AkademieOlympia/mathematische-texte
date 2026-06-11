\documentclass{article}

% --- Encoding: Wichtig für pdflatex ---
\usepackage{inputenc}
\usepackage{fontenc}
\usepackage{lmodern}

% --- Standard-Pakete ---
\usepackage{amsmath,amssymb,amsthm}
\usepackage{geometry}
\usepackage{booktabs}
\usepackage{listings}
\usepackage{hyperref}

\geometry{margin=2.5cm}

\lstset{
language=Python,
basicstyle=\ttfamily\small,
keywordstyle=\bfseries,
breaklines=true,
tabsize=4
}

% --- Theorem-Umgebungen ---
\theoremstyle{plain}
\newtheorem{satz}{Satz}
\newtheorem{lemma}{Lemma}
\newtheorem{korollar}{Korollar}
\newtheorem{vermutung}{Vermutung}

\theoremstyle{definition}
\newtheorem{definition}{Definition}
\newtheorem{beispiel}{Beispiel}

\theoremstyle{remark}
\newtheorem{bemerkung}{Bemerkung}

\newenvironment{beweis}{\begin{proof}[Beweis]}{\end{proof}}

\title{Quaternionische Zeitkristalle auf EABC-Tetraedern mit alternierender Parit\"atsphase}
\author{Thomas Hoffbauer}
\date{\today}

\begin{document}
\maketitle

\begin{abstract}
Wir untersuchen eine diskrete Dynamik auf einem tetraedrischen Zustandsraum mit Quaternionenstruktur (EABC-Modell). Numerische Experimente zeigen die Existenz diskreter Zeitkristall-Orbits mit Perioden
\[
T = 3,4,6,8,10,12,14,16,\dots
\]
Die Einf\"uhrung einer alternierenden Parit\"atsphase zwischen geraden und ungeraden Tetraedern f\"uhrt zu einer systematischen Transformation der beobachteten Perioden:
\[
T_{\mathrm{obs}}=
\begin{cases}
T,& T \text{ gerade}\\
2T,& T \text{ ungerade}
\end{cases}
\]
Dies deutet auf eine bipartite Tetraederstruktur hin, in der ontische Quaternionenorbits durch eine zus\"atzliche $\mathbb{Z}_2$-Symmetrie \"uberlagert werden.
\end{abstract}

\section{Geometrischer Ausgangspunkt}
Wir betrachten einen tetraedrischen Zustandsraum mit vier Komponenten
\[
(E,A,B,C)
\]
Die Darstellung kann als baryzentrischer Simplex interpretiert werden:
\[
\lambda_E+\lambda_A+\lambda_B+\lambda_C = 1
\]

\section{Quaternionische Darstellung}
Der Zustand wird als Quaternion geschrieben:
\[
Q = E + A i + B j + C k
\]
mit Norm
\[
\|Q\|^2 = E^2 + A^2 + B^2 + C^2
\]
Normierte Zust\"ande liegen auf der 3-Sph\"are $S^3$.

\section{Zeitentwicklung}
Die Dynamik wird durch Konjugation mit einem normierten Quaternion definiert:
\[
Q_{n+1} = U Q_n U^{-1}
\]
Da Quaternionen normerhaltend sind, gilt:
\[
\|Q_{n+1}\| = \|Q_n\|
\]

\section{Zeitkristall-Orbits}
Numerische Experimente zeigen diskrete Perioden:
\[
T = 3,4,6,8,10,12,14,16
\]
F\"ur Rotationen
\[
\theta = \frac{\pi}{k}
\]
ergibt sich empirisch:
\[
T = 2k
\]

\begin{satz}[Achsiale Periodenregel]
F\"ur einen Rotationsoperator
\[
U(\theta), \quad \theta=\frac{\pi}{k}
\]
ergibt sich generisch
\[
T=2k
\]
\end{satz}

\section{Triadischer Spezialorbit}
Der symmetrische Operator
\[
U_{\text{tri}} \propto 1+i+j+k
\]
erzeugt numerisch eine Periode
\[
T=3
\]

\section{Parit\"atsphase}
Wir erweitern den Zustand um eine bin\"are Phase:
\[
\epsilon_n \in \{+1,-1\}, \quad \epsilon_n = (-1)^n
\]

\begin{definition}
Der beobachtete Zustand sei
\[
\widetilde{Q}_n = \epsilon_n Q_n
\]
\end{definition}

\section{Beobachtete Periodizit\"at}
Angenommen:
\[
Q_{n+T} = Q_n
\]
Dann folgt:
\[
\widetilde{Q}_{n+T} = (-1)^T \widetilde{Q}_n
\]

\begin{satz}[Parit\"atsregel]
Die beobachtete Periode lautet:
\[
T_{\text{obs}} =
\begin{cases}
T & T \text{ gerade}\\
2T & T \text{ ungerade}
\end{cases}
\]
\end{satz}

\begin{beweis}
Einsetzen liefert:
\[
\widetilde{Q}_{n+T} = (-1)^{n+T}Q_{n+T} = (-1)^T \widetilde{Q}_n
\]
Ist $T$ gerade, folgt sofort Wiederkehr.
Ist $T$ ungerade, ist erst nach $2T$ Wiederkehr erreicht.
\end{beweis}

\section{Numerische Resultate}
Tests bis 100000 Iterationen ergeben:

\begin{center}
\begin{tabular}{ccc}
\toprule
Operator & ontische Periode & beobachtete Periode \\
\midrule
triadisch & 3 & 6 \\
$\pi/2$ & 4 & 4 \\
$\pi/3$ & 6 & 6 \\
$\pi/4$ & 8 & 8 \\
$\pi/5$ & 10 & 10 \\
$\pi/6$ & 12 & 12 \\
$\pi/7$ & 14 & 14 \\
$\pi/8$ & 16 & 16 \\
\bottomrule
\end{tabular}
\end{center}

\section{Erweiterter Zustandsraum}
\begin{definition}
Der parit\"atserweiterte Zustandsraum sei:
\[
\mathcal{X} = \mathbb{H} \times \mathbb{Z}_2
\]
\end{definition}

\section{Vermutung}
\begin{vermutung}
Alle ungeraden ontischen Zeitkristallperioden verdoppeln sich unter alternierender Tetraederparit\"at.
\end{vermutung}

\appendix
\section{Python-Code}
\begin{lstlisting}
import numpy as np

def quat_mult(q1,q2):
a1,b1,c1,d1 = q1
a2,b2,c2,d2 = q2
return np.array([
a1*a2-b1*b2-c1*c2-d1*d2,
a1*b2+b1*a2+c1*d2-d1*c2,
a1*c2-b1*d2+c1*a2+d1*b2,
a1*d2+b1*c2-c1*b2+d1*a2
])

def evolve(q0,u,steps):
u_inv = np.array([u[0],-u[1],-u[2],-u[3]])
q = q0
orbit = 
for _ in range(steps):
q = quat_mult(u,quat_mult(q,u_inv))
orbit.append(q)
return np.array(orbit)

def observed_orbit(base):
obs = []
for n,q in enumerate(base):
eps = -1 if n%2 else 1
obs.append(eps*q)
return np.array(obs)
\end{lstlisting}

\end{document}