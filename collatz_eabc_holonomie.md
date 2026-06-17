# EABC-Holonomie — Korrektur des diskreten Assoziators

**Status:** Theorem + Heuristik + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Theorem | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_discrete_associator.md` / `.py` — naive $V_4$-Assoziator $\mathfrak a\equiv 0$ (**Theorem**)
- `collatz_eabc_kommutator_assoziator.md` — Kommutator ($\mathbb{H}$) vs. Assoziator ($\mathbb{O}$)
- `collatz_eabc_oktonion_assoziator_spektralhypothese.md` — $\Gamma_E$-Projektion, $M_n^E(t)$
- `collatz_eabc_invarianzprogramm.md` — globale Observable $\chi(x)$
- `collatz_eabc_holonomie_test.py` / `.json` — Numerik (Orientierung $\omega$, $\chi$, Oktanion-$\Gamma_E$)
- `eabc_from_lean.py` / `EABC.lean` — Chiralität ABCE / CEAB

---

## 1. Kritische Korrektur: $V_4$ ist die Klein-Vierergruppe

**Definition ($\Phi$ auf $V_4$).** Mit $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11 \pmod{12}$:
$$\Phi(X,Y) := \mathrm{classOf}\bigl(\mathrm{residue}(X)\cdot\mathrm{residue}(Y)\bmod 12\bigr).$$

**Theorem ($V_4$-Struktur).** $\{1,5,7,11\}\subset(\mathbb Z/12\mathbb Z)^\times$ ist unter Multiplikation
isomorph zur **Klein-Vierergruppe** $V_4\cong\mathbb Z_2\times\mathbb Z_2$; $E$ ist neutral.

**Theorem (Assoziativität).** Jede Gruppe ist assoziativ. Daher für alle $X,Y,Z\in\{E,A,B,C\}$:
$$\Phi(\Phi(X,Y),Z)=\Phi(X,\Phi(Y,Z)),\qquad
\mathfrak a_{\mathrm{naiv}}(X,Y,Z):=\mathrm{sgn}\bigl(\Phi(\Phi(X,Y),Z),\Phi(X,\Phi(Y,Z))\bigr)=0.$$

**Experiment:** `collatz_eabc_discrete_associator.py::prove_v4_klein_associativity` — $4^3$ Tripel, alle $\mathfrak a=0$.

> **Ehrliche Antwort:** Der **naive** EABC-Assoziator auf $V_4$ allein misst **keine** Nichtassoziativität.
> Er ist **immer null** — kein Bug, sondern Algebra.

**Label:** $V_4$-Isomorphie und $\mathfrak a_{\mathrm{naiv}}\equiv 0$ = **Theorem**.

---

## 2. Echter EABC-„Assoziator“ = projektionsbasierte Holonomie

**Definition (Projektions-Holonomie).** Sei $\Gamma$ eine Projektion von einem **Trägerobjekt**
(Primzahl, Hurwitz-Punkt, Oktanion, Schalenblock) auf EABC-Koordinaten (glatt oder diskret).
Für ein Tripel $(x,y,z)$ im Träger:
$$\mathfrak a_E(x,y,z) := \Gamma\bigl((xy)z\bigr) - \Gamma\bigl(x(yz)\bigr).$$
Der Defekt entsteht **nach** der Multiplikation im höheren Ring, **nicht** durch $\Phi$ auf $V_4$ allein.

| Ebene | Träger | $\Gamma$ | Typischer Defekt |
|-------|--------|----------|------------------|
| Diskret sichtbar | $V_4$ | Klassen-ID | $\equiv 0$ (**Theorem**) |
| Quaternion | $\Sigma_p\subset\mathbb H_{\mathrm H}$ | $\Gamma_4$, Chiralität | Kommutator / $\chi_p$ |
| Oktanion | $\Sigma_n\subset\mathbb Z^8$ (Stub) | glatt-$\Gamma_E$ (16 Koordinaten) | $\neq 0$ typisch (**Experiment**) |
| Prim-Vierling | $Q=(p,p{+}2,p{+}6,p{+}8)$ | Orientierung $\omega(Q)$ | diskrete $\pm 1$ (**Definition**) |

**Label:** $\mathfrak a_E$ auf Trägern = **Definition**; $V_4$-Trivialität = **Theorem**.

> **Boxed (Heuristik, kein DG-Theorem).**
> $$\boxed{\;\text{Der EABC-Assoziator ist vermutlich keine Nichtassoziativität auf }V_4,
> \text{ sondern diskrete Holonomie: }\Gamma(\text{Weg}_1)\neq\Gamma(\text{Weg}_2).\;}$$
> Vergleiche zwei **Klammerwege** $(xy)z$ vs.\ $x(yz)$ **nach** Projektion — analog zu
> $\mathrm{Hol}(\gamma)=\exp(i\oint A)$ in der Eichtheorie.

---

## 3. Physik-Analogie: Projektionsverlust

In einem Eichtheorie-Bild:

| Algebraisches Objekt | Physikalische Lesart |
|----------------------|----------------------|
| Multiplikation in $\mathbb O$ | lokale Zusammensetzung in voller 8D-Struktur |
| Assoziator $[x,y,z]=(xy)z-x(yz)$ | Krümmung der Multiplikation (Oktanionen) |
| Projektion $\Gamma$ auf EABC | **Messkanal** — wie ein Detektor, der nur vier Familien sieht |
| $\Gamma((xy)z)-\Gamma(x(yz))$ | **Holonomie** — derselbe algebraische Endpunkt, verschiedene Klammerwege liefern verschiedene $\Gamma$-Spuren |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | Messung **nur** im Bild der Projektion auf eine **assoziative** Klein-Gruppe → immer 0 |

**Heuristik:** Der „Verlust“ beim Übergang $\mathbb H\to\mathbb O$ ist nicht nur Nichtassoziativität in $\mathbb O$,
sondern auch **Informationsverlust durch $\Gamma$**: zwei Wege können im Träger divergieren, während ihr
$V_4$-Schatten identisch bleibt.

**Label:** Physik-Analogie = **Heuristik** (kein bewiesenes DG-Modell).

---

## 4. ABCE / CEAB = Orientierung $\omega$, nicht Produkt

Aus `EABC.lean` / `eabc_from_lean.py`:

| Chiralität | Klassenfolge | Start $p\bmod 12$ | **Orientierung** $\omega(Q)$ |
|------------|--------------|-------------------|------------------------------|
| **ABCE** | $A,B,C,E$ | $5$ | $+1$ |
| **CEAB** | $C,E,A,B$ | $11$ | $-1$ |

**Definition.**
$$\omega(Q)\in\{+1,-1\},\qquad
\omega(Q)=+1\Leftrightarrow\sigma(Q)=\text{ABCE},\quad
\omega(Q)=-1\Leftrightarrow\sigma(Q)=\text{CEAB}.$$

**Wichtig:** ABCE und CEAB sind **zyklische Orientierungen** auf dem Vierling — **keine** Produkte
$A\cdot B\cdot C\cdot E$ in $V_4$. Die Heuristik „ABCE $\leftrightarrow$ linke Klammer, CEAB $\leftrightarrow$ rechte Klammer“
(`collatz_eabc_discrete_associator.md` §4) beschreibt eine **Lesart**, keinen algebraischen Defekt auf $V_4$.

**Label:** $\omega$ = **Definition**; Klammer-Link = **Heuristik**.

---

## 5. Verbindung zur Invariante $\chi$

### Globale $\chi$ (Invarianzprogramm)

$$\chi(x)=\frac{(E(x)+C(x))-(A(x)+B(x))}{\pi_{>3}(x)}
=\mathcal I_\chi(S(x)),\qquad
\mathcal I_\chi(e,a,b,c)=(e+c)-(a+b).$$

Misst die **chirale Asymmetrie** der Primverteilung in den vier EABC-Restklassen (Modus $\Phi_1$ auf $V_4$,
vgl.\ `collatz_eabc_invarianzprogramm.md` §8).

### Quadruplet-Holonomie

Für Prim-Vierlinge $Q$ mit $p\le N$:
$$\Phi_{\mathrm{quad}}(N):=\frac{1}{|Q_{\le N}|}\sum_{Q:\,p\le N}\omega(Q)
=\frac{\#\mathrm{ABCE}-\#\mathrm{CEAB}}{|Q_{\le N}|}\in[-1,1].$$

**Theorem (balancierte Beine).** Für jeden Vierling $Q$ gilt bei der **Bein-Chiralität**
$\chi_{\mathrm{leg}}(Q):=\#\{E,C\}-\#\{A,B\}$ auf den vier Primbeinen:
$$\chi_{\mathrm{leg}}(Q)=0\quad\text{für ABCE und CEAB}.$$
Die Orientierung $\omega$ trägt **keine** Bein-Asymmetrie — sie codiert die **zyklische Phase** des Vierlings.

**Heuristik (Holonomie $=$ $\chi$ oder Erweiterung).**

| Observable | Raum | Was misst sie? |
|------------|------|----------------|
| $\chi(x)$ | alle Primzahlen $>3$ | globale EC-vs-AB-Bilanz ($\Phi_1$-Richtung) |
| $\Phi_{\mathrm{quad}}(N)$ | Vierlingsstarts | ABCE-vs-CEAB-Orientierung (diskrete Holonomie) |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | Klassenprodukt | **nichts** (immer 0) |
| $\mathfrak a_E$ auf $\mathbb O$ | $\Gamma_E((xy)z)-\Gamma_E(x(yz))$ | echte Klammer-Holonomie im 8D-Stub |

**Ehrliches Fazit:** Quadruplet-Holonomie $\omega$ ist **nicht** identisch mit $\chi(x)$ — sie ist eine
**Erweiterung** auf der Vierlings-Untermenge: $\chi$ auf Vierlings**beinen** ist trivial ($0$),
während $\omega$ die **Orientierung** der Signatur unterscheidet. Beide sind „chiral“, aber in
**verschiedenen Koordinaten** (globale Primzählung vs.\ lokale Vierlingsphase).

**Experiment:** `collatz_eabc_holonomie_test.py` — `chi_global`, `chi_quad_legs`, `holonomy_flux_phi_quad`.

---

## 6. Vergleich zum Oktanion-Assoziator

| Größe | Domäne | Typ | Referenz |
|-------|--------|-----|----------|
| $[x,y,z]=(xy)z-x(yz)$ | $\mathbb O$ | algebraisch, typ.\ $\neq 0$ | `collatz_eabc_oktonion_associator.py` |
| $\|\Gamma_E((xy)z)-\Gamma_E(x(yz))\|$ | $\mathbb Z^8$-Stub | projektionsbasierte Holonomie | `eabc_associator_vector` |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | $\{E,A,B,C\}$ | **immer 0** | `collatz_eabc_discrete_associator.py` |
| $\omega(Q)$ | Prim-Vierlinge | diskrete $\pm 1$-Phase | dieses Dokument |

Der Oktanion-Assoziator misst **algebraische** Nichtassoziativität; die EABC-Holonomie misst den
**Projektionsdefekt** zwischen zwei Klammerwegen. Auf $V_4$ kollabiert beides zum Null-Assoziator;
in $\mathbb O$ und auf Vierlingen bleibt messbare Struktur.

---

## 7. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $V_4$ Klein-Gruppe, $\Phi$ assoziativ | **Theorem** |
| $\mathfrak a_{\mathrm{naiv}}\equiv 0$ auf $V_4$ | **Theorem** |
| $\mathfrak a_E=\Gamma((xy)z)-\Gamma(x(yz))$ auf Trägern | **Definition** |
| EABC-Assoziator = diskrete Holonomie | **Heuristik** |
| $\omega(\mathrm{ABCE})=+1$, $\omega(\mathrm{CEAB})=-1$ | **Definition** |
| $\chi_{\mathrm{leg}}(Q)=0$ für kanonische Vierlinge | **Theorem** |
| $\Phi_{\mathrm{quad}}$ vs.\ globale $\chi$ — verwandt, nicht identisch | **Experiment** |
| Physik-Analogie (Projektionsverlust) | **Heuristik** |

---

*Kanonsiche Notiz: Nach der Korrektur ist der diskrete $V_4$-Assoziator ein **Negativkontroll-Experiment** —
er bestätigt, dass Klammerdefekte **erst nach Hebung/Projektion** sichtbar werden. Die operative Forschung
verlagert sich auf $\Gamma_E$-Holonomie ($\mathbb O$) und $\omega$-Holonomie (Vierlinge), parallel zur globalen $\chi$.*
