## Ikosier-Zahlen und Minkowski–Galois-Einbettung

Die Grundidee lautet:

> Eine Ikosier-Zahl ist zunächst ein vierdimensionales Quaternion mit Koeffizienten im Goldenen-Zahlkörper. Durch die gleichzeitige Betrachtung dieser Koeffizienten und ihrer Galois-Konjugierten wird daraus ein achtdimensionaler reeller Vektor.

### 1. Der Goldene Zahlkörper

Man arbeitet über

[
K=\mathbb Q(\sqrt5)
]

mit dem Goldenen Schnitt

[
\varphi=\frac{1+\sqrt5}{2}.
]

Der Zahlkörper besitzt zwei reelle Einbettungen:

[
\operatorname{id}:\sqrt5\longmapsto \sqrt5,
\qquad
\sigma:\sqrt5\longmapsto-\sqrt5.
]

Insbesondere gilt

[
\sigma(\varphi)
===============

# \bar\varphi

# \frac{1-\sqrt5}{2}

# 1-\varphi

-\frac1\varphi.
]

Ein Element

[
x=a+b\varphi
]

besitzt daher zwei reelle Erscheinungsformen:

[
x=a+b\varphi,
\qquad
\sigma(x)=a+b\bar\varphi.
]

Beispielsweise wird

[
\varphi
\quad\longmapsto\quad
(\varphi,\bar\varphi)
\approx
(1{,}618,-0{,}618).
]

Die beiden Zahlen sind keine voneinander unabhängigen Koordinaten. Sie sind die beiden Galois-Bilder desselben algebraischen Elements.

---

### 2. Was ist eine Ikosier-Zahl?

Ein Quaternion über (K) hat die Form

[
q=a_0+a_1i+a_2j+a_3k,
\qquad a_\nu\in K,
]

mit

[
i^2=j^2=k^2=ijk=-1.
]

Der **Ikosierring** (\mathbb I) ist eine besondere ganzzahlige Ordnung innerhalb der Quaternionenalgebra

[
\mathbb H(K)=\mathbb H\otimes_{\mathbb Q}\mathbb Q(\sqrt5).
]

Die PDF schreibt vereinfachend

[
q=a_0+a_1i+a_2j+a_3k,
\qquad a_\nu\in\mathbb Z[\varphi].
]

Streng genommen ist der Ikosierring nicht einfach die Menge aller Quaternionen mit beliebigen Koeffizienten aus (\mathbb Z[\varphi]); zusätzlich gelten bestimmte Gitter- beziehungsweise Kongruenzbedingungen. Für das geometrische Verständnis genügt zunächst aber die Vorstellung eines Quaternionengitters über dem Goldenen Ring.

Die 120 geeigneten Einheiten beziehungsweise Einheitswurzeln bilden die binäre Ikosaederstruktur und können als die 120 Wurzeln des nichtkristallographischen Wurzelsystems (H_4), also als Eckpunkte des 600-Zells, aufgefasst werden. Genau diesen Übergang beschreibt auch die PDF.  Die Verbindung zwischen Ikosiern, (H_4) und dem (E_8)-Gitter ist auch in der mathematischen Literatur etabliert. ([arXiv][1])

---

### 3. Zwei verschiedene Konjugationen

Hier müssen zwei Operationen streng getrennt werden.

#### Quaternionische Konjugation

Für

[
q=a_0+a_1i+a_2j+a_3k
]

ist

[
q^*=a_0-a_1i-a_2j-a_3k.
]

Damit erhält man die Quaternionennorm

[
n(q)=qq^*
=========

a_0^2+a_1^2+a_2^2+a_3^2
\in K.
]

#### Galois-Konjugation

Die Galois-Abbildung (\sigma) wirkt dagegen nur auf die Koeffizienten:

[
\sigma(q)
=========

\sigma(a_0)+\sigma(a_1)i+\sigma(a_2)j+\sigma(a_3)k.
]

Sie ändert also

[
\sqrt5\mapsto-\sqrt5,
\qquad
\varphi\mapsto\bar\varphi,
]

aber nicht die Quaternioneneinheiten (i,j,k).

Das ist entscheidend:

[
q^*
\neq
\sigma(q).
]

Die erste Operation gehört zur Quaternionenalgebra, die zweite zur algebraischen Zahlentheorie.

---

### 4. Die Minkowski–Galois-Einbettung

Da (K=\mathbb Q(\sqrt5)) zwei reelle Einbettungen besitzt, wird jedes Quaternion (q) gleichzeitig in beiden Einbettungen betrachtet:

[
\iota(q)
========

\frac1{\sqrt2}
\bigl(q,\sigma(q)\bigr).
]

Weil jedes Quaternion vier reelle Koordinaten besitzt, liegt das Paar in

[
\mathbb H\oplus\mathbb H
\cong
\mathbb R^4\oplus\mathbb R^4
\cong
\mathbb R^8.
]

Aus

[
q=(a_0,a_1,a_2,a_3)
]

entsteht konkret

[
\iota(q)
========

\frac1{\sqrt2}
\left(
a_0,a_1,a_2,a_3,,
\sigma(a_0),\sigma(a_1),\sigma(a_2),\sigma(a_3)
\right).
]

Man kann dies als zwei gekoppelte vierdimensionale Schatten desselben algebraischen Objekts verstehen:

[
\boxed{
q
\quad\rightsquigarrow\quad
\text{kanonischer (H_4)-Sektor}
\oplus
\text{Galois-konjugierter (H_4)-Sektor}
}
]

Die Dimension (8) wird also nicht künstlich hinzugefügt. Sie entsteht aus

[
4\text{ Quaternionenkoordinaten}
\times
2\text{ reelle Einbettungen}
============================

8.

]

---

### 5. Die Norm in acht Dimensionen

Für die skalierte Einbettung gilt

[
|\iota(q)|^2
============

\frac12
\left(
n(q)+\sigma(n(q))
\right).
]

Die rechte Seite ist die halbe Körpertrace der Quaternionennorm:

[
|\iota(q)|^2
============

\frac12
\operatorname{Tr}_{K/\mathbb Q}\bigl(n(q)\bigr).
]

Falls

[
n(q)=1,
]

dann ist auch

[
\sigma(n(q))=1
]

und damit

[
|\iota(q)|^2=1.
]

Hier liegt eine **Normierungsfrage der PDF**: Standardmäßig werden (E_8)-Wurzeln meist mit

[
|\alpha|^2=2
]

normiert. Die dort verwendete Abbildung mit dem Faktor (1/\sqrt2) liefert für eine Ikosier-Einheit zunächst Normquadrat (1). Für die übliche (E_8)-Normierung muss man daher entweder anschließend mit (\sqrt2) skalieren oder die Einbettung ohne diesen Vorfaktor verwenden.

Das ist kein struktureller Widerspruch, sondern eine nicht ausdrücklich getrennte Skalierungskonvention.

---

### 6. Von (H_4) zu (E_8)

Das Wurzelsystem (H_4) besitzt 120 Wurzeln. Das (E_8)-Wurzelsystem besitzt 240 Wurzeln. In einer geeigneten ikosianischen Darstellung kann man die 240 Wurzeln als zwei miteinander gekoppelte 120er-Sätze organisieren:

[
\Phi(E_8)
\sim
H_4\cup\varphi H_4.
]

Der zweite Satz ist dabei nicht einfach eine beliebige vergrößerte Kopie des ersten. Die Multiplikation mit (\varphi) wird unter Galois-Konjugation durch

[
\sigma(\varphi)=\bar\varphi=-\frac1\varphi
]

ergänzt. Gerade das Zusammenspiel von (\varphi) und (\bar\varphi) ermöglicht es, die nichtkristallographische fünfzählige (H_4)-Geometrie in einer achtdimensionalen kristallographischen Struktur unterzubringen.

Die PDF stellt dies auf Seite 5 als Wirkungskette dar:

[
\text{Ikosierring}
\longrightarrow
\text{Minkowski–Galois-Einbettung}
\longrightarrow
H_4\cup\varphi H_4
\longrightarrow
E_8.
]

Sie betont dabei zu Recht, dass es sich um eine **algebraisch-arithmetische Zerlegung** und nicht um eine topologische Zerlegung des acht-dimensionalen Raumes handelt. 

---

## Der anschauliche Kern

Man kann die Konstruktion auf einen Satz verdichten:

[
\boxed{
\text{Die Minkowski–Galois-Einbettung macht aus einem goldenen 4D-Quaternion zwei gekoppelte 4D-Bilder und damit einen 8D-Vektor.}
}
]

Oder geometrisch:

* (H_4) trägt die ikosaedrische Fünfergeometrie.
* (\varphi) erzeugt die goldene Skalierung.
* (\sigma) liefert die konjugierte Skalierung (\bar\varphi).
* Beide Sektoren gemeinsam können sich zu einer kristallographischen (E_8)-Struktur schließen.

Die Ikosier-Konstruktion ist damit keine bloße Analogie zu (E_8), sondern eine konkrete arithmetische Realisierung des (E_8)-Gitters aus Quaternionen über (\mathbb Q(\sqrt5)).

[1]: https://arxiv.org/pdf/2303.12219?utm_source=chatgpt.com "arXiv:2303.12219v2 [math.RA] 17 Apr 2023"
