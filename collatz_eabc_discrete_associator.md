# Diskreter EABC-Assoziator auf $V_4$ — Negativkontrolle

> **Redirect (kanonisch, PR #54):** Der naive diskrete „EABC-Assoziator“ auf $V_4$ **verschwindet**
> ($\mathfrak a_{\mathrm{naiv}}\equiv 0$, Klein-Vierergruppe). Für Orientierung, Chiralität und die
> sechsstufige Begriffshierarchie siehe **`collatz_eabc_holonomie.md`** — $\omega$ (Stufe 3),
> $\chi_E(N)$ (Stufe 4); projektive EABC-Holonomie $\mathcal H_E$ = Zielobjekt (Stufe 6).
> Der Begriff „EABC-Assoziator“ ist auf $V_4$-Ebene **deprecated**; der **algebraische** Assoziator
> $[x,y,z]$ auf $\mathbb{O}$ bleibt in `collatz_eabc_oktonion_associator.py` gültig.

**Status:** Theorem (Trivialität) + Negativkontrolle  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Theorem | Heuristik | Experiment

**Querverweise:**
- `EABC.lean` / `eabc_from_lean.py` — Restklassen, Chiralität, Vierling $Q(p)$
- `collatz_eabc_holonomie.md` — Begriffshierarchie Stufe 1–6; $\omega$, $\chi_E$ (implementiert); $\mathcal H_E$ (Zielobjekt)
- `collatz_eabc_kommutator_assoziator.md` — Kommutator ($\mathbb{H}$) vs. Assoziator ($\mathbb{O}$)
- `collatz_eabc_discrete_associator.py` / `.json` — Rechnung und Tabellen
- `Global Lokal.py` — explizite $V_4$-Multiplikationstabelle (Referenz)

---

## 1. Motivation (Negativkontrolle)

Der volle **algebraische** Oktanion-Assoziator $[x,y,z]$ (`collatz_eabc_oktonion_associator.py`) ist auf $\mathbb{O}$
**typisch $\neq 0$**. Für die **sichtbare** EABC-Ebene — vier Familien $\{E,A,B,C\}$ modulo
$12$ — wurde ein diskreter Klammertest auf $V_4$ formuliert.

**Ergebnis (Theorem):** $\Phi$ ist assoziativ; $\mathfrak a_{\mathrm{naiv}}\equiv 0$.
Messbare Orientierung und Chiralität liegen in Stufe 3–4 (`collatz_eabc_holonomie.md`: $\omega$, $\chi_E$),
nicht in $V_4$-Nichtassoziativität. Wahre Holonomie (Stufe 6) ist Zielobjekt, noch nicht erreicht.

---

## 2. Definition von $\Phi$ (aus Lean/Python)

**Definition ($\Phi$ auf $V_4$).** Aus `EABC.lean` / `eabc_from_lean.py`:
$$E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11 \pmod{12}.$$
Setze
$$\Phi(X,Y) := \mathrm{classOf}\bigl(\mathrm{residue}(X)\cdot \mathrm{residue}(Y) \bmod 12\bigr).$$

**Theorem ($V_4$-Struktur).** Die Menge $\{1,5,7,11\}\subset (\mathbb Z/12\mathbb Z)^\times$
ist unter Multiplikation eine **Klein-Vierergruppe** $V_4$; $E$ ist neutral.
Die Tabelle stimmt mit `Global Lokal.py` und der quaternionischen Hebung
$E\leftrightarrow 1$, $A\leftrightarrow i$, $B\leftrightarrow j$, $C\leftrightarrow k$
als **kommutative** $V_4$-Schattenstruktur überein.

**Label:** $\Phi$ = **Definition**; $V_4$-Isomorphie = **Theorem**.

---

## 3. Diskreter Assoziator

**Definition.**
$$\mathfrak a(X,Y,Z) := \mathrm{sgn}\bigl(\Phi(\Phi(X,Y),Z),\;\Phi(X,\Phi(Y,Z))\bigr)\in\{-1,0,+1\},$$
mit $\mathrm{sgn}=0$ bei Gleichheit.

**Theorem (Ehrlichkeit auf $V_4$).** $\Phi$ ist **assoziativ** auf ganz $V_4$. Daher
$$\mathfrak a(X,Y,Z)=0 \quad \text{für alle } X,Y,Z\in\{E,A,B,C\},$$
einschließlich des reinen $\{A,B,C\}^3$-Subtripels.

**Experiment:** `collatz_eabc_discrete_associator.py::check_associativity` — $4^3$ Tripel,
alle $\mathfrak a=0$.

> **Ehrliche Antwort:** Auf der **vollständigen** $V_4$-Algebra ist der diskrete Assoziator
> **trivial**. Nichtassoziativität tritt erst in der Hebung $\mathbb{H}\to\mathbb{O}$ auf
> (`collatz_eabc_kommutator_assoziator.md` §2–§3).

---

## 4. Chiralität ABCE / CEAB — Klammerpräferenz (Heuristik)

Aus `EABC.lean` / `eabc_from_lean.py`:
- **ABCE**-Reihenfolge: $[A,B,C,E]$
- **CEAB**-Reihenfolge: $[C,E,A,B]$

**Heuristik (kein Theorem).**
| Chiralität | Klammerpräferenz | Lesart auf $(A,B,C)$ |
|------------|------------------|----------------------|
| **ABCE** | **links** | $\Phi(\Phi(A,B),C)$ |
| **CEAB** | **rechts** | $\Phi(A,\Phi(B,C))$ |

Primzahlvierlinge $Q(p)$ mit $p\equiv 5\pmod{12}$ tragen **ABCE**; $p\equiv 11\pmod{12}$
tragen **CEAB** (`collatz_phi_pref_test.py`, Lean-Tests in `tests/test_eabc_from_lean.py`).

**Experiment:** Für alle Vierlinge bis einer Grenze liefern beide Klammern **dasselbe**
$V_4$-Ergebnis — die Heuristik unterscheidet **Orientierung der Lesart**, nicht einen
algebraischen Defekt auf $V_4$.

**Label:** ABCE/CEAB-Klammerlink = **Heuristik**.

---

## 5. Primfolge-Observable $\mathcal A(N)$

**Definition (Experiment).** Für die ersten $N$ Primzahlen $p_1<\cdots<p_N$ mit
$X_k=\mathrm{classOf}(p_k)$:
$$\mathcal A(N) := \frac{1}{N-2}\sum_{k=2}^{N-1} \mathfrak a(X_{k-1},X_k,X_{k+1}).$$

**Experiment:** $\mathcal A(N)\equiv 0$ für getestete $N$ — konsistent mit $V_4$-Assoziativität.

**Label:** $\mathcal A(N)$ = **Experiment**; Stabilität $\neq 0$ = **nicht beobachtet** auf $V_4$.

---

## 6. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $\Phi$ via $\mathrm{residue}(X)\cdot\mathrm{residue}(Y)$ | **Definition** |
| $V_4$ assoziativ, $\mathfrak a\equiv 0$ | **Theorem** |
| ABCE ↔ links, CEAB ↔ rechts | **Heuristik** |
| $\mathcal A(N)$ auf Prim-EABC-Folge | **Experiment** |
| Nichttrivialer diskreter Assoziator auf vollem $V_4$ | **verworfen** (Theorem) |

---

*Kanonsiche Notiz: Der diskrete EABC-Assoziator ist der ehrliche $V_4$-Schatten des
Oktanion-Programms — er bestätigt Assoziativität auf der sichtbaren mod-$12$-Ebene und
reserviert Klammerdefekte für das **Holonomie-Zielobjekt** Stufe 6 (`collatz_eabc_holonomie.md`).*
