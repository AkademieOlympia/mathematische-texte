# EABC-Holonomie (kanonisch)

**Status:** Theorem + Definition + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Theorem | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_discrete_associator.md` / `.py` — naive $V_4$-Klammerung $\mathfrak a_{\mathrm{naiv}}\equiv 0$ (**Negativkontrolle**, kein „EABC-Assoziator“)
- `collatz_eabc_kommutator_assoziator.md` — Kommutator ($\mathbb{H}$) vs. algebraischer Assoziator $[x,y,z]$ ($\mathbb{O}$)
- `collatz_eabc_oktonion_assoziator_spektralhypothese.md` — $\Gamma_E$-Projektion, $M_n^E(t)$ (Oktanion-Ebene)
- `collatz_eabc_invarianzprogramm.md` — globale Observable $\chi(x)$
- `collatz_eabc_holonomie_test.py` / `.json` — Numerik ($\chi_E$, $\omega$, $\mathcal H_E$, Vergleich zu $\chi$)
- `eabc_from_lean.py` / `EABC.lean` — Chiralität ABCE / CEAB

---

## 0. Kanonische Terminologie (PR #54)

> **Nicht** „EABC-Assoziator“ auf $V_4$-Ebene: $V_4$ ist assoziativ (Klein-Vierergruppe).

| Begriff | Ebene | Status |
|---------|-------|--------|
| **EABC-Holonomie** / **projektive EABC-Holonomie** | Träger + $\Gamma$ | **kanonisch** |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | mod-$12$-Klassenprodukt | **verworfen** als Chiralitätsmaß |
| $[x,y,z]=(xy)z-x(yz)$ | $\mathbb{O}$ | **algebraischer Assoziator** (bleibt gültig) |
| $\Delta\Gamma_E=\Gamma_E((xy)z)-\Gamma_E(x(yz))$ | $\mathbb{O}$-Stub | **projektive Holonomie** $\mathcal H_E$ |

> $$\boxed{\;\text{EABC-Chiralität ist projektive Holonomie, nicht algebraische Nichtassoziativität.}\;}$$

---

## 1. $V_4$ ist assoziativ — Negativkontrolle

**Definition ($\Phi$ auf $V_4$).** Mit $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11 \pmod{12}$:
$$\Phi(X,Y) := \mathrm{classOf}\bigl(\mathrm{residue}(X)\cdot\mathrm{residue}(Y)\bmod 12\bigr).$$

**Theorem ($V_4$-Struktur).** $\{1,5,7,11\}\subset(\mathbb Z/12\mathbb Z)^\times$ ist isomorph zur **Klein-Vierergruppe** $V_4$; $E$ ist neutral.

**Theorem (Assoziativität).** Für alle $X,Y,Z\in\{E,A,B,C\}$:
$$\Phi(\Phi(X,Y),Z)=\Phi(X,\Phi(Y,Z)),\qquad
\mathfrak a_{\mathrm{naiv}}(X,Y,Z):=\mathrm{sgn}\bigl(\Phi(\Phi(X,Y),Z),\Phi(X,\Phi(Y,Z))\bigr)=0.$$

**Experiment:** `collatz_eabc_discrete_associator.py::prove_v4_klein_associativity` — $4^3$ Tripel, alle $\mathfrak a_{\mathrm{naiv}}=0$.

**Label:** $V_4$-Isomorphie und $\mathfrak a_{\mathrm{naiv}}\equiv 0$ = **Theorem**.

---

## 2. Projektive EABC-Holonomie $\mathcal H_E$

**Definition (Rekonstruktionspfade).** Sei $\Gamma$ eine Projektion von einem **Trägerobjekt**
(Primzahl, Hurwitz-Punkt, Oktanion, Schalenblock, Vierling) auf EABC-Koordinaten.
Zwei **Rekonstruktionspfade** $\gamma_1,\gamma_2$ durch dieselbe Konfiguration sind z.\,B.:
- Klammerwege $(xy)z$ vs.\ $x(yz)$ mit gemeinsamem algebraischem Endpunkt in $\mathbb O$;
- zyklische Orientierungen **ABCE** vs.\ **CEAB** auf demselben Prim-Vierling $Q(p)$.

**Definition ($d_E$ auf $\Gamma$-Raum).** Für $\gamma,\eta$ im Bild von $\Gamma$ (z.\,B. $\mathbb{Z}^{16}$ bei glatt-$\Gamma_E$):
$$d_E(\gamma,\eta) := \|\gamma-\eta\|_2.$$

**Definition (projektive EABC-Holonomie).**
$$\boxed{\;\mathcal H_E(\gamma_1,\gamma_2) := d_E\bigl(\Gamma(\gamma_1),\,\Gamma(\gamma_2)\bigr).\;}$$

Spezialfall Oktanion-Stub ($\Gamma=\Gamma_E$):
$$\mathcal H_E\bigl((xy)z,\,x(yz)\bigr) = d_E\bigl(\Gamma_E((xy)z),\,\Gamma_E(x(yz))\bigr).$$

| Ebene | Träger | $\Gamma$ | $\mathcal H_E$ |
|-------|--------|----------|----------------|
| $V_4$ | Klassen-ID | $\mathrm{classOf}$ | $\equiv 0$ (**Theorem**) |
| Vierling | $Q=(p,p{+}2,p{+}6,p{+}8)$ | Orientierung $\omega(Q)$ | diskret $\in\{0,2\}$ in $\{\pm1\}$-Raum |
| $\mathbb{O}$ | $\Sigma_n\subset\mathbb{Z}^8$ | glatt-$\Gamma_E$ | typ.\ $\neq 0$ (**Experiment**) |

**Label:** $\mathcal H_E$, $d_E$ = **Definition**; $V_4$-Trivialität = **Theorem**.

---

## 3. Orientierung $\omega$ auf Prim-Vierlingen

Aus `EABC.lean` / `eabc_from_lean.py`:

| Chiralität | Klassenfolge | Start $p\bmod 12$ | **Orientierung** $\omega(Q)$ |
|------------|--------------|-------------------|------------------------------|
| **ABCE** | $A,B,C,E$ | $5$ | $+1$ |
| **CEAB** | $C,E,A,B$ | $11$ | $-1$ |

**Definition.**
$$\omega(Q)\in\{+1,-1\},\qquad
\omega(Q)=+1\Leftrightarrow\sigma(Q)=\text{ABCE},\quad
\omega(Q)=-1\Leftrightarrow\sigma(Q)=\text{CEAB}.$$

ABCE und CEAB sind **zyklische Orientierungen** — **keine** Produkte $A\cdot B\cdot C\cdot E$ in $V_4$.

**Label:** $\omega$ = **Definition**.

---

## 4. Vierlings-Chiralität $\chi_E(N)$

Für Prim-Vierlinge mit Startprimzahl $p\le N$:

**Definition.**
$$\boxed{\;\chi_E(N)=\frac{\#\mathrm{ABCE}_{\le N}-\#\mathrm{CEAB}_{\le N}}
{\#\mathrm{ABCE}_{\le N}+\#\mathrm{CEAB}_{\le N}}\in[-1,1].\;}$$

Äquivalent: $\chi_E(N)=\frac{1}{|Q_{\le N}|}\sum_{Q:\,p\le N}\omega(Q)$ (mittlere Orientierung).

**Experiment:** `collatz_eabc_holonomie_test.py::chi_E`.

---

## 5. Verbindung zur globalen Invariante $\chi(x)$

### Globale $\chi$ (Invarianzprogramm)

$$\chi(x)=\frac{(E(x)+C(x))-(A(x)+B(x))}{\pi_{>3}(x)}
=\mathcal I_\chi(S(x)),\qquad
\mathcal I_\chi(e,a,b,c)=(e+c)-(a+b).$$

Misst die **chirale Asymmetrie** aller Primzahlen $>3$ in den vier EABC-Restklassen.

### Vergleich $\chi_E$ vs.\ $\chi$

| Observable | Raum | Formel | Was misst sie? |
|------------|------|--------|----------------|
| $\chi(x)$ | alle Primzahlen $>3$ | $((E+C)-(A+B))/\pi_{>3}$ | globale EC-vs-AB-Bilanz |
| $\chi_E(N)$ | Vierlingsstarts $p\le N$ | $(\#\mathrm{ABCE}-\#\mathrm{CEAB})/(\#\mathrm{ABCE}+\#\mathrm{CEAB})$ | **projektive Holonomie** auf Vierlingen |
| $\chi_{\mathrm{leg}}(Q)$ | vier Beine eines Vierlings | $(\#E+\#C)-(\#A+\#B)$ | **$0$** für kanonische Vierlinge (**Theorem**) |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | Klassenprodukt | — | **nichts** (immer $0$) |

**Theorem (balancierte Beine).** Für jeden kanonischen Prim-Vierling $Q$:
$$\chi_{\mathrm{leg}}(Q)=0\quad\text{für ABCE und CEAB}.$$

**Experiment (ehrlicher Vergleich).** $\chi_E(N)$ und $\chi(x)$ sind **verwandte chirale Observablen in verschiedenen Räumen** — nicht identisch. $\chi_E$ codiert die ABCE/CEAB-**Phase**; $\chi$ trackt die globale Primzählung. Lean- und Python-Implementierung stimmen für $\omega$ überein (`tests/test_eabc_from_lean.py`, `tests/test_eabc_holonomie.py`).

---

## 6. Oktanion-Ebene: algebraischer Assoziator vs. projektive Holonomie

| Größe | Domäne | Typ | Referenz |
|-------|--------|-----|----------|
| $[x,y,z]=(xy)z-x(yz)$ | $\mathbb O$ | **algebraische** Nichtassoziativität | `collatz_eabc_oktonion_associator.py` |
| $\mathcal H_E((xy)z,x(yz))=d_E(\Gamma_E((xy)z),\Gamma_E(x(yz)))$ | $\mathbb{Z}^8$-Stub | **projektive Holonomie** | `eabc_associator_vector` |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | $\{E,A,B,C\}$ | **immer $0$** | `collatz_eabc_discrete_associator.py` |
| $\omega(Q)$, $\chi_E(N)$ | Prim-Vierlinge | diskrete Holonomie | dieses Dokument |

Der **algebraische** Oktanion-Assoziator $[x,y,z]$ und die **projektive** Holonomie $\mathcal H_E$ sind **verschiedene Ebenen**. Auf $V_4$ kollabiert $\mathfrak a_{\mathrm{naiv}}$; in $\mathbb O$ und auf Vierlingen bleibt messbare Struktur.

---

## 7. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $V_4$ Klein-Gruppe, $\Phi$ assoziativ | **Theorem** |
| $\mathfrak a_{\mathrm{naiv}}\equiv 0$ auf $V_4$ | **Theorem** |
| $\mathcal H_E(\gamma_1,\gamma_2)=d_E(\Gamma(\gamma_1),\Gamma(\gamma_2))$ | **Definition** |
| EABC-Chiralität = projektive Holonomie, nicht $V_4$-Nichtassoziativität | **Theorem** + **Definition** |
| $\omega(\mathrm{ABCE})=+1$, $\omega(\mathrm{CEAB})=-1$ | **Definition** |
| $\chi_E(N)$ Vierlings-Chiralität | **Definition** |
| $\chi_{\mathrm{leg}}(Q)=0$ für kanonische Vierlinge | **Theorem** |
| $\chi_E$ vs.\ globale $\chi$ — verwandt, nicht identisch | **Experiment** |
| $[x,y,z]$ auf $\mathbb O$ | **Theorem** (algebraisch, eigene Ebene) |

---

*Kanonsiche Notiz: Der diskrete $V_4$-„Assoziator“ ist eine **Negativkontrolle** — er bestätigt, dass Klammerdefekte **erst nach Hebung/Projektion** sichtbar werden. Operative Forschung: $\mathcal H_E$ auf $\mathbb O$ ($M_n^E$) und $\chi_E$, $\omega$ auf Vierlingen, parallel zur globalen $\chi$.*
