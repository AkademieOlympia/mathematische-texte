import subprocess

latex_code = r"""\documentclass[11pt,a4paper,reqno]{amsart}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[german]{babel}
\usepackage{amsmath,amssymb,amsthm,amsfonts,mathtools}
\usepackage{hyperref}
\usepackage{geometry}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{tikz-cd}

\geometry{verbose,tmargin=2.5cm,bmargin=2.5cm,lmargin=2.5cm,rmargin=2.5cm}

\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{axiom}[theorem]{Axiom}
\newtheorem{example}[theorem]{Example}
\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}

\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.97,0.97,0.97}

\lstdefinestyle{sagestyle}{
 backgroundcolor=\color{backcolour}, 
 commentstyle=\color{codegreen},
 keywordstyle=\color{magenta},
 numberstyle=\tiny\color{codegray},
 stringstyle=\color{codepurple},
 basicstyle=\ttfamily\footnotesize,
 breakatwhitespace=false, 
 breaklines=true, 
 captionpos=b, 
 keepspaces=true, 
 numbers=left, 
 numbersep=5pt, 
 showspaces=false, 
 showstringspaces=false,
 showtabs=false, 
 tabsize=4
}
\lstset{style=sagestyle}

\newcommand{\M}{\mathcal{M}}
\newcommand{\N}{\mathbb{N}}
\newcommand{\Pset}{\mathbb{P}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\Fact}{\operatorname{Fact}}
\newcommand{\status}{\operatorname{status}}
\newcommand{\supp}{\operatorname{supp}}
\newcommand{\PABC}{P_{ABC}}

\title[Skalare Projektionen, Faserquotienten und Paritätsebenen]{Skalare Projektionen, Faserquotienten und EABC-Paritätsebenen in Modellräumen}
\author{Thomas Hoffbauer}
\address{Computer Institut Bamberg, Deutschland}
\email{research@cib.de}
\date{5. August 2026}

\begin{document}

\begin{abstract}
The purpose of this paper is not to extend the Fundamental Theorem of Arithmetic, but to formalize the compatibility between scalar factorization in $\mathbb{N}$ and multiple state representations in a structured model space $\mathcal{M}$. We analyze the algebraic interface between a model space $\mathcal{M}$ and the natural numbers $\mathbb{N}$ via scalar projections $\pi: \mathcal{M} \to \mathbb{N}$. We prove the canonical fiber quotient theorem $\mathcal{M}/{\equiv_\pi} \cong \operatorname{Im}(\pi)$, establish the universal property of fiber projections, and derive the scalar factorization invariant $\Phi_{\mathrm{fact}} = \Fact_\mathbb{N} \circ \pi$. Furthermore, we refine the multi-level decomposition of $\mathbb{N}_{>0}$ into carrier, parity, and residue levels under the $EABC$ prime classification framework, providing a fully audited SageMath 10.5 reference implementation.
\end{abstract}

\maketitle

\section{Einleitung und Motivation}
In strukturellen und algebraischen Modellräumen $\mathcal{M}$ treten häufig Zustandsrepräsentationen auf, deren skalare Invariante in den natürlichen Zahlen $\mathbb{N}$ liegt, während die Modellzustände selbst geometrische, nichtkommutative oder erweiterte algebraische Eigenschaften aufweisen.

Ziel dieser Arbeit ist es keineswegs, den Hauptsatz der Arithmetik zu erweitern oder zu modifizieren. Vielmehr wird der mathematische Kompatibilitätsrahmen formalisiert, der sicherstellt, dass Mehrfachdarstellungen von Zuständen im Modellraum $\mathcal{M}$ bei identischer Skalarprojektion $\pi(x) = \pi(y) = n \in \mathbb{N}$ in vollständiger Harmonie mit der eindeutigen Primfaktorzerlegung von $n$ stehen.

\section{Strukturierte Modellräume und Überkategorien}
Um den Begriff des Modellraums kategorial sauber einzuordnen, betrachten wir ein Modell $(M, \pi)$ nicht als isoliertes Objekt, sondern als Objekt einer Schnitt- bzw. Überkategorie (Slice Category).

\begin{definition}[Modellobjekt in der Überkategorie]
Es sei $\mathcal{C}$ eine algebraische Kategorie (z.\,B. die Kategorie der Monoide $\mathbf{Mon}$ oder Magmen $\mathbf{Mag}$) und $\mathbb{N}$ ein ausgezeichnetes Objekt in $\mathcal{C}$. Ein \emph{Modellraum mit skalaren Projektionen} ist ein Objekt $(M, \pi)$ der Überkategorie $\mathcal{C}/\mathbb{N}$, wobei $\pi: M \to \mathbb{N}$ ein $\mathcal{C}$-Morphismus ist. Ein Morphismus $f: (M, \pi_M) \to (M', \pi_{M'})$ in $\mathcal{C}/\mathbb{N}$ ist ein $\mathcal{C}$-Morphismus $f: M \to M'$, der das Diagramm
\begin{equation}
\begin{tikzcd}
M \arrow[rr, "f"] \arrow[rd, "\pi_M"'] & & M' \arrow[ld, "\pi_{M'}"] \\
& \mathbb{N} &
\end{tikzcd}
\end{equation}
kommutativ macht ($\pi_{M'} \circ f = \pi_M$).
\end{definition}

\section{Kanonische Quotientenstruktur und Universelle Eigenschaft}
Wir untersuchen die durch eine Projektion $\pi: M \to \mathbb{N}$ induzierte algebraische Zerlegung des Modellraums $M$.

\begin{definition}[Projektionsäquivalenz und Fasern]
Auf dem Modellraum $M$ wird die Äquivalenzrelation der \emph{Projektionsäquivalenz} $\equiv_\pi$ definiert durch:
\begin{equation}
x \equiv_\pi y \quad :\Longleftrightarrow \quad \pi(x) = \pi(y).
\end{equation}
Für ein $n \in \operatorname{Im}(\pi) \subseteq \mathbb{N}$ heißt die Äquivalenzklasse $M_n := \pi^{-1}(n) = \{ x \in M \mid \pi(x) = n \}$ die \emph{Projektionsfaser} von $n$.
\end{definition}

\begin{theorem}[Theorem A: Kanonischer Quotientensatz]
\label{thm:quotient}
Es sei $q: M \to M/{\equiv_\pi}$ die kanonische quotientenseitige Projektion $x \mapsto [x]_\pi$. Dann existiert eine kanonische Bijektion $\overline{\pi}: M/{\equiv_\pi} \xrightarrow{\,\cong\,} \operatorname{Im}(\pi)$ mit $\pi = \overline{\pi} \circ q$, mithin:
\begin{equation}
\boxed{ M/{\equiv_\pi} \;\cong\; \operatorname{Im}(\pi). }
\end{equation}
Insbesondere gilt für die Faserzerlegung von $M$:
\begin{equation}
M = \bigsqcup_{n \in \operatorname{Im}(\pi)} M_n.
\end{equation}
\end{theorem}

\begin{proof}
Die Wohldefiniertheit und Injektivität von $\overline{\pi}([x]_\pi) := \pi(x)$ folgt unmittelbar aus der Definition von $\equiv_\pi$. Die Surjektivität auf $\operatorname{Im}(\pi)$ ist trivial.
\end{proof}

\begin{theorem}[Theorem B: Universelle Eigenschaft des Faserquotienten]
\label{thm:universal}
Es sei $f: M \to X$ eine beliebige Abbildung in eine Menge $X$, die auf den $\pi$-Fasern konstant ist, d.\,h.:
\begin{equation}
\forall x,y \in M: \quad \pi(x) = \pi(y) \implies f(x) = f(y).
\end{equation}
Dann existiert genau eine Abbildung $\widetilde{f}: M/{\equiv_\pi} \to X$, sodass das folgende Diagramm kommutiert:
\begin{equation}
\begin{tikzcd}
M \arrow[r, "f"] \arrow[d, "q"'] & X \\
M/{\equiv_\pi} \arrow[ru, "\exists! \, \widetilde{f}"'] & 
\end{tikzcd}
\end{equation}
\end{theorem}

\begin{proof}
Wir definieren $\widetilde{f}([x]_\pi) := f(x)$. Da $f$ auf den Fasern von $\pi$ konstant ist, ist $\widetilde{f}$ wohldefiniert und eindeutig bestimmt.
\end{proof}

\section{Die Skalare Faktorisierungsinvariante}
Mithilfe der universellen Eigenschaft konstruieren wir die systematische Verbindung zum Hauptsatz der Arithmetik. Es sei $\mathbb{P} \subset \mathbb{N}$ die Menge der natürlichen Primzahlen. Jede Zahl $n \in \mathbb{N}$ besitzt eine eindeutige Bewertungsdarstellung:
\begin{equation}
\Fact_\mathbb{N}(n) := (v_p(n))_{p \in \mathbb{P}} \;\in\; \bigoplus_{p \in \mathbb{P}} \mathbb{N}_0.
\end{equation}

\begin{theorem}[Theorem C: Skalare Faktorisierungsinvariante]
Die Komposition $\Phi_{\mathrm{fact}} := \Fact_\mathbb{N} \circ \, \pi$ definiert eine Abbildung
\begin{equation}
\Phi_{\mathrm{fact}}: M \longrightarrow \bigoplus_{p \in \mathbb{P}} \mathbb{N}_0.
\end{equation}
Sie ist auf jeder Projektionsfaser $M_n$ konstant. Nach Theorem~\ref{thm:universal} induziert sie eine eindeutige Invariante $\widetilde{\Phi}_{\mathrm{fact}}: M/{\equiv_\pi} \to \bigoplus_{p \in \mathbb{P}} \mathbb{N}_0$.
\end{theorem}

\begin{corollary}[Schutzformel]
\label{cor:schutz}
Für beliebige Elemente $x,y \in M$ gilt:
\begin{equation}
\boxed{ \pi(x) = \pi(y) \;\;\not\Rightarrow\;\; x = y, \quad x \sim y \quad \text{oder} \quad \Fact_M(x) = \Fact_M(y). }
\end{equation}
Gleicher Skalarwert impliziert weder Identität noch Assoziiertheit oder Übereinstimmung der Faktorenstrukturen im Modellraum.
\end{corollary}

\section{Typentrennung und Exakte Algebraische Beispiele}

\begin{definition}[Semiprimes Modell-Tupel]
Ein semiprimes Modell-Tupel ist ein Element $(a,b \,; c)_M \;\in\; \mathbb{P}^2 \times M$. Das Semikolon trennt die arithmetischen Primfaktoren $(a,b) \in \mathbb{P}^2$ vom Modellzustand $c \in M$.
\end{definition}

\begin{proposition}[Korrektes Gegenbeispiel in $\mathbb{Z}[i]$]
\label{prop:gauss}
In den Gaußschen Zahlen $M = \mathbb{Z}[i]$ mit der quadratischen Norm $\pi(z) = N(z) = z\bar{z}$ betrachten wir $x = 5 \in \mathbb{Z}[i]$ und $y = 3+4i \in \mathbb{Z}[i]$:
\begin{enumerate}
 \item $N(x) = N(5) = 25$ und $N(y) = N(3+4i) = 3^2 + 4^2 = 25$, mithin $x \equiv_\pi y$.
 \item $x \neq y$ und $x \not\sim y$ (da $5 / (3+4i) = \frac{15-20i}{25} \notin \mathbb{Z}[i]$).
 \item Das Element $x = 5$ ist in $\mathbb{Z}[i]$ \textbf{reduzibel}: $5 = (2+i)(2-i)$.
 \item Die Faktorisierungsstrukturen unterscheiden sich: $5 = (2+i)(2-i)$ besitzt zwei verschiedene irreducible Faktoren der Norm 5, während $3+4i = (2+i)^2$ das Quadrat eines einzelnen irreduciblen Elements ist.
\end{enumerate}
\end{proposition}

\begin{example}[Zweites Normbeispiel: Eisenstein-Zahlen $\mathbb{Z}[\omega]$]
Es sei $\omega = e^{2\pi i / 3} = -\frac{1}{2} + \frac{\sqrt{3}}{2}i$. Im Ring der Eisenstein-Zahlen $M = \mathbb{Z}[\omega]$ mit der Norm $N(a+b\omega) = a^2 - ab + b^2$ zerfällt $7 = (3+\omega)(2-\omega)$. Für $x = 7$ gilt $N(x) = 49$. Für den Zustand $u = 3+8\omega$ gilt $N(u) = 3^2 - 3\cdot 8 + 8^2 = 9 - 24 + 64 = 49$. Somit sind $x = 7$ und $u = 3+8\omega$ projektionsäquivalent ($\pi(x) = \pi(u) = 49$), während ihre Zerlegungen $7 = (3+\omega)(2-\omega)$ und $3+8\omega = (1+\omega)(3+\omega)^2$ völlig verschieden sind.
\end{example}

\section{Das Geschärfte EABC-Drei-Ebenen-System}
Auf der Menge $\mathbb{N}_{>0}$ definieren wir für jede Zahl $n = 2^{\alpha(n)} 3^{\beta(n)} n^\sharp$ mit $\gcd(n^\sharp, 6) = 1$ die Einteilung der Primzahlen $p > 3$ in vier Familien:
\begin{align*}
\mathbb{P}_E &= \{p > 3: p \equiv 1 \pmod{12}\}, & \mathbb{P}_A &= \{p > 3: p \equiv 5 \pmod{12}\}, \\
\mathbb{P}_B &= \{p > 3: p \equiv 7 \pmod{12}\}, & \mathbb{P}_C &= \{p > 3: p \equiv 11 \pmod{12}\}.
\end{align*}

Daraus ergeben sich drei hierarchisch aufeinander aufbauende Ebene-Konzepte:

\subsection{1. Trägerebene (Familienträger)}
Der Familienträger $\supp_{EABC}(n) := \{ X \in \{E,A,B,C\} : \exists p \in \mathbb{P}_X \text{ mit } p \mid n \}$ induziert die disjunkte 16-Teilung $\mathbb{N}_{>0} = \bigsqcup_{S \subseteq \{E,A,B,C\}} \mathcal{N}_S$. Es gilt das Vereinigungsgesetz:
\begin{equation}
\supp_{EABC}(mn) = \supp_{EABC}(m) \cup \supp_{EABC}(n).
\end{equation}

\subsection{2. Paritätsebene (Exponentenraum Modulo 2)}
Für $X \in \{A,B,C\}$ sei $\Omega_X(n) := \sum_{p \in \mathbb{P}_X} v_p(n)$. Der Paritätsvektor $\varepsilon(n) := (\Omega_A(n) \bmod 2, \Omega_B(n) \bmod 2, \Omega_C(n) \bmod 2)^T \in \mathbb{F}_2^3$ induziert einen Epimorphismus $\varepsilon: (\mathbb{N}_{>0}, \cdot) \to (\mathbb{F}_2^3, \oplus)$ mit
\begin{equation}
\varepsilon(mn) = \varepsilon(m) \oplus \varepsilon(n).
\end{equation}

\subsection{3. Restklassenebene (Klein'sche Vierergruppe)}
Die Zuordnung $\varphi: \mathbb{F}_2^3 \to (\mathbb{Z}/12\mathbb{Z})^\times \cong V_4$ gegeben durch $\varphi(a,b,c) = 5^a 7^b 11^c \pmod{12}$ besitzt den 1-dimensionalen Kern:
\begin{equation}
\operatorname{ker}(\varphi) = \left\{ \begin{pmatrix} 0\\0\\0 \end{pmatrix}, \begin{pmatrix} 1\\1\\1 \end{pmatrix} \right\} \;\subset\; \mathbb{F}_2^3.
\end{equation}
Das Element $(1,1,1)^T$ entspricht dem vollen $ABC$-Tripel $p_A p_B p_C \equiv 5 \cdot 7 \cdot 11 = 385 \equiv 1 \pmod{12}$, welches im EABC-System als neutraler Baustein direkt in den $E$-Kern $K_E(n)$ absorbiert wird.

\section{Algorithmetische Referenzimplementierung in SageMath 10.5}

\begin{lstlisting}[language=Python, caption={SageMath 10.5 Skript zur auditierbaren Faser- und Äquivalenzprüfung}]
from sage.all import *

class SemiprimeModelSpace:
    def __init__(self, R):
        self.R = R
        self.P = Primes()

    def projection(self, element):
        return ZZ(element.norm())

    def compute_fiber(self, n, search_bound=20):
        fiber = []
        for x_val in range(-search_bound, search_bound + 1):
            for y_val in range(-search_bound, search_bound + 1):
                elem = self.R(x_val + y_val*I)
                if elem != 0 and self.projection(elem) == n:
                    fiber.append(elem)
        return fiber

    def verify_equivalence(self, a, b, x_elem, y_elem):
        assert a in self.P and b in self.P, "a und b muessen Primzahlen sein!"
        n_scalar = ZZ(a * b)
        
        pi_x = self.projection(x_elem)
        pi_y = self.projection(y_elem)
        
        is_equiv = (pi_x == pi_y == n_scalar)
        is_identical = (x_elem == y_elem)
        
        is_associated = False
        if y_elem != 0:
            quot = x_elem / y_elem
            is_associated = quot in self.R and (1/quot) in self.R

        return {
            "n_scalar": n_scalar,
            "Fact_N": factor(n_scalar),
            "pi(x)": pi_x,
            "pi(y)": pi_y,
            "equiv_pi": is_equiv,
            "identical": is_identical,
            "associated": is_associated,
            "x_factors": x_elem.factor(),
            "y_factors": y_elem.factor()
        }

if __name__ == "__main__":
    M_Gauss = SemiprimeModelSpace(ZZ[I])
    a_prime, b_prime = 5, 5
    x_state = ZZ[I](5)
    y_state = ZZ[I](3 + 4*I)
    
    res = M_Gauss.verify_equivalence(a_prime, b_prime, x_state, y_state)
    print("="*65)
    print(f"VERIFIKATION DER PROJEKTIONSAEQUIVALENZ FUER n = {res['n_scalar']}")
    print("="*65)
    print(f"1. Skalarwert n in N       : {res['n_scalar']}")
    print(f"2. Faktorisierung Fact_N(n): {res['Fact_N']}")
    print(f"3. Projektion pi(x)        : {res['pi(x)']}")
    print(f"4. Projektion pi(y)        : {res['pi(y)']}")
    print(f"5. Projektionsaequivalent?  : {res['equiv_pi']}")
    print(f"6. Identisch (x == y)?     : {res['identical']}")
    print(f"7. Assoziiert (x ~ y)?     : {res['associated']}")
    print(f"8. Faktorisierung von x    : {res['x_factors']}")
    print(f"9. Faktorisierung von y    : {res['y_factors']}")
    print("="*65)
\end{lstlisting}

\section{Fazit}
Die Differenzierung zwischen dem skalaren Invariantenwert $n \in \mathbb{N}$ und den Elementen der Modellfaser $M_n = \pi^{-1}(n)$ löst sämtliche Scheinwidersprüche bezüglich des Hauptsatzes der Arithmetik vollständig auf. Die Synthese aus kategorialem Faserquotienten, geschärfter Paritätsebene und auditiertem SageMath-Code bildet eine transparente und rigorose Grundlage für strukturierte Modellrepräsentationen.

\begin{thebibliography}{99}
\bibitem{Hoffbauer2026} T. Hoffbauer, \emph{On the Compatibility of Model-Space Representations with the Fundamental Theorem of Arithmetic}, arXiv preprint, August 2026.
\bibitem{Niven} I. Niven, H. S. Zuckerman, H. L. Montgomery, \emph{An Introduction to the Theory of Numbers}, 5th ed., John Wiley \& Sons, 1991.
\bibitem{Conway} J. H. Conway, D. A. Smith, \emph{On Quaternions and Octonions: Their Geometry, Arithmetic, and Symmetry}, A K Peters/CRC Press, 2003.
\end{thebibliography}

\end{document}
"""

with open("paper.tex", "w") as f:
    f.write(latex_code)

res = subprocess.run(["pdflatex", "-interaction=nonstopmode", "paper.tex"], capture_output=True, text=True)
res2 = subprocess.run(["pdflatex", "-interaction=nonstopmode", "paper.tex"], capture_output=True, text=True)

print("Return code:", res2.returncode)
if res2.returncode != 0:
    print("Output:", res2.stdout[-1500:])
else:
    print("PDF generated successfully.")