# EABC-Holonomie (kanonisch)

**Status:** Definition (Stufe 1–4) + Zielobjekt (Stufe 5–6) + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition (1–4) | Conjecture (5–6) | Forschungsfrage | Theorem | Experiment

**Querverweise:**
- `collatz_eabc_discrete_associator.md` / `.py` — naive $V_4$-Klammerung $\mathfrak a_{\mathrm{naiv}}\equiv 0$ (**Negativkontrolle**, kein „EABC-Assoziator“)
- `collatz_eabc_kommutator_assoziator.md` — Kommutator ($\mathbb{H}$) vs. algebraischer Assoziator $[x,y,z]$ ($\mathbb{O}$)
- `collatz_eabc_oktonion_assoziator_spektralhypothese.md` — $\Gamma_E$-Projektion, $M_n^E(t)$ (Oktanion-Ebene)
- `collatz_eabc_invarianzprogramm.md` — globale Observable $\chi(x)$; Stufe 1–2 ($\kappa$, $\sigma(Q)$)
- `collatz_eabc_holonomie_test.py` / `.json` — Numerik ($\chi_E$, $\omega$, $\mathcal H_E$-Stub, Vergleich zu $\chi$)
- `eabc_from_lean.py` / `EABC.lean` — Chiralität ABCE / CEAB

---

## 0. Epistemische Einordnung

### Wahre Holonomie (Differentialgeometrie)

**Holonomie** im klassischen Sinn: $\mathrm{Hol}(\gamma)=P_\gamma$ — Paralleltransport entlang einer **geschlossenen** Schleife $\gamma$ in einem Zustandsraum mit Verbindung. Voraussetzungen: Zustandsraum, Transportoperator, geschlossene Pfade, Vergleich von Transporten.

### Aktueller EABC-Status — noch **keine** wahre Holonomie

Was heute implementiert und gemessen ist, sind **Orientierung** und **Chiralität** auf Signaturen und Vierlingen — nicht yet Paralleltransport geschlossener Schleifen.

| Begriff (aktuell) | Status |
|-------------------|--------|
| **EABC-Orientierungsfunktion** $\omega(Q)$ | **Definition** (Stufe 3) |
| **EABC-Chiralitätsfunktion** $\chi_E(N)$ | **Definition** (Stufe 4) |
| ABCE/CEAB-Clustering auf Vierlingen | **Vorläufer** holonomer Struktur, **noch keine** Holonomie |
| **projektive EABC-Holonomie** $\mathcal H_E$ | **Zielobjekt** (langfristige Vision, Stufe 6) |

> $$\boxed{\;\text{EABC-Chiralität als beobachtete Orientierungsasymmetrie — nicht (noch) wahre Holonomie.}\;}$$

**Langfristige Vision:** projektive EABC-Holonomie — Vergleich geschlossener Transportpfade $\gamma_1,\gamma_2$ nach Definition eines natürlichen Transportoperators $T$ (Stufe 5–6).

> $$\boxed{\;\textbf{Forschungsfrage:}\;\text{Welcher natürliche Transportoperator verbindet zwei EABC-Signaturen?}\;}$$

---

## 1. Begriffshierarchie (Stufe 1–6)

Die sechs Stufen bauen aufeinander auf. Stufe 1–4 sind **definiert und implementiert**; Stufe 5–6 sind **offen** (Conjecture / Forschungsfrage).

| Stufe | Objekt | Formel / Inhalt | Tao-Label | Status |
|-------|--------|-----------------|-----------|--------|
| **1** | **EABC-Klassifikation** | $\kappa(p)\in\{E,A,B,C\}$ | **Definition** | implementiert (`EABC.lean`, `eabc_from_lean.py`) |
| **2** | **EABC-Signaturen** | $Q=(E,A,B,C)^4$, $\sigma(Q)\in\Sigma_4$ | **Definition** | implementiert (`collatz_eabc_invarianzprogramm.md` Def. 5) |
| **3** | **EABC-Orientierung** | $\omega(Q)\in\{-1,+1\}$; $\omega(\mathrm{ABCE})=+1$, $\omega(\mathrm{CEAB})=-1$ | **Definition** | implementiert (`EABC.lean`) |
| **4** | **EABC-Chiralität** | $\chi_E(N)=\dfrac{1}{N}\sum_{Q:\,p\le N}\omega(Q)$ | **Definition** | implementiert (`collatz_eabc_holonomie_test.py`) |
| **5** | **EABC-Transport** | $T:\,Q_k\mapsto Q_{k+1}$ | **Conjecture** | **NICHT DEFINIERT** — Forschungsfrage |
| **6** | **EABC-Holonomie** | Vergleich geschlossener Transportpfade $\gamma_1,\gamma_2$ via $T$ | **Conjecture** | Zukunft — benötigt Stufe 5 |

**Verbindung ABCE/CEAB:** Das beobachtete Clustering der Vierlingssignaturen in die zyklischen Orientierungen ABCE und CEAB ist ein **Vorläufer** holonomer Struktur — es misst Orientierungsasymmetrie (Stufe 3–4), noch keinen Paralleltransport (Stufe 5–6).

**Numerik:** `collatz_eabc_holonomie.json` → Feld `research_stage` (Stufe 1–4: `implemented`; Stufe 5–6: `open`).

---

## 2. $V_4$ ist assoziativ — Negativkontrolle

**Definition ($\Phi$ auf $V_4$).** Mit $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11 \pmod{12}$:
$$\Phi(X,Y) := \mathrm{classOf}\bigl(\mathrm{residue}(X)\cdot\mathrm{residue}(Y)\bmod 12\bigr).$$

**Theorem ($V_4$-Struktur).** $\{1,5,7,11\}\subset(\mathbb Z/12\mathbb Z)^\times$ ist isomorph zur **Klein-Vierergruppe** $V_4$; $E$ ist neutral.

**Theorem (Assoziativität).** Für alle $X,Y,Z\in\{E,A,B,C\}$:
$$\Phi(\Phi(X,Y),Z)=\Phi(X,\Phi(Y,Z)),\qquad
\mathfrak a_{\mathrm{naiv}}(X,Y,Z):=\mathrm{sgn}\bigl(\Phi(\Phi(X,Y),Z),\Phi(X,\Phi(Y,Z))\bigr)=0.$$

**Experiment:** `collatz_eabc_discrete_associator.py::prove_v4_klein_associativity` — $4^3$ Tripel, alle $\mathfrak a_{\mathrm{naiv}}=0$.

**Label:** $V_4$-Isomorphie und $\mathfrak a_{\mathrm{naiv}}\equiv 0$ = **Theorem**.

> **Nicht** „EABC-Assoziator“ auf $V_4$-Ebene: $V_4$ ist assoziativ (Klein-Vierergruppe). Der algebraische Assoziator $[x,y,z]$ auf $\mathbb{O}$ bleibt auf seiner Ebene gültig.

---

## 3. Orientierung $\omega$ auf Prim-Vierlingen (Stufe 3)

Aus `EABC.lean` / `eabc_from_lean.py`:

| Chiralität | Klassenfolge | Start $p\bmod 12$ | **Orientierung** $\omega(Q)$ |
|------------|--------------|-------------------|------------------------------|
| **ABCE** | $A,B,C,E$ | $5$ | $+1$ |
| **CEAB** | $C,E,A,B$ | $11$ | $-1$ |

**Definition (EABC-Orientierungsfunktion).**
$$\omega(Q)\in\{+1,-1\},\qquad
\omega(Q)=+1\Leftrightarrow\sigma(Q)=\text{ABCE},\quad
\omega(Q)=-1\Leftrightarrow\sigma(Q)=\text{CEAB}.$$

ABCE und CEAB sind **zyklische Orientierungen** — **keine** Produkte $A\cdot B\cdot C\cdot E$ in $V_4$.

**Label:** $\omega$ = **Definition** (Stufe 3).

---

## 4. Vierlings-Chiralität $\chi_E(N)$ (Stufe 4)

Für Prim-Vierlinge mit Startprimzahl $p\le N$:

**Definition (EABC-Chiralitätsfunktion).**
$$\boxed{\;\chi_E(N)=\frac{\#\mathrm{ABCE}_{\le N}-\#\mathrm{CEAB}_{\le N}}
{\#\mathrm{ABCE}_{\le N}+\#\mathrm{CEAB}_{\le N}}\in[-1,1].\;}$$

Äquivalent: $\chi_E(N)=\frac{1}{|Q_{\le N}|}\sum_{Q:\,p\le N}\omega(Q)$ (mittlere Orientierung).

$\chi_E(N)$ misst die **beobachtete Orientierungsasymmetrie** über Vierlinge — **nicht** (noch) Holonomie im DG-Sinn.

**Experiment:** `collatz_eabc_holonomie_test.py::chi_E`.

**Label:** $\chi_E$ = **Definition** (Stufe 4).

---

## 5. Zielobjekt: projektive EABC-Holonomie $\mathcal H_E$ (Stufe 6, Conjecture)

**Definition (Rekonstruktionspfade).** Sei $\Gamma$ eine Projektion von einem **Trägerobjekt**
(Primzahl, Hurwitz-Punkt, Oktanion, Schalenblock, Vierling) auf EABC-Koordinaten.
Zwei **Rekonstruktionspfade** $\gamma_1,\gamma_2$ durch dieselbe Konfiguration sind z.\,B.:
- Klammerwege $(xy)z$ vs.\ $x(yz)$ mit gemeinsamem algebraischem Endpunkt in $\mathbb O$;
- zyklische Orientierungen **ABCE** vs.\ **CEAB** auf demselben Prim-Vierling $Q(p)$.

**Definition ($d_E$ auf $\Gamma$-Raum).** Für $\gamma,\eta$ im Bild von $\Gamma$ (z.\,B. $\mathbb{Z}^{16}$ bei glatt-$\Gamma_E$):
$$d_E(\gamma,\eta) := \|\gamma-\eta\|_2.$$

**Definition (projektive EABC-Holonomie — Zielobjekt).**
$$\boxed{\;\mathcal H_E(\gamma_1,\gamma_2) := d_E\bigl(\Gamma(\gamma_1),\,\Gamma(\gamma_2)\bigr).\;}$$

Spezialfall Oktanion-Stub ($\Gamma=\Gamma_E$):
$$\mathcal H_E\bigl((xy)z,\,x(yz)\bigr) = d_E\bigl(\Gamma_E((xy)z),\,\Gamma_E(x(yz))\bigr).$$

| Ebene | Träger | $\Gamma$ | $\mathcal H_E$ | Status |
|-------|--------|----------|----------------|--------|
| $V_4$ | Klassen-ID | $\mathrm{classOf}$ | $\equiv 0$ | **Theorem** |
| Vierling | $Q=(p,p{+}2,p{+}6,p{+}8)$ | Orientierung $\omega(Q)$ | diskret $\in\{0,2\}$ in $\{\pm1\}$-Raum | Orientierung (Stufe 3), **nicht** Holonomie |
| $\mathbb{O}$ | $\Sigma_n\subset\mathbb{Z}^8$ | glatt-$\Gamma_E$ | typ.\ $\neq 0$ | **Experiment** (Stub) |

**Label:** $\mathcal H_E$, $d_E$ = **Definition** (Zielrahmen); vollständige Holonomie = **Conjecture** (Stufe 6, benötigt Transportoperator Stufe 5).

### Stufe 5 — EABC-Transport (offen)

**Conjecture / Forschungsfrage:** Existiert ein natürlicher Operator
$$T:\,Q_k \longrightarrow Q_{k+1}$$
der aufeinanderfolgende EABC-Signaturen entlang einer Collatz- oder Vierlingskette verbindet, sodass geschlossene Pfade $\gamma$ sinnvoll definierbar sind?

> $$\boxed{\;\textbf{Forschungsfrage:}\;\text{Welcher natürliche Transportoperator verbindet zwei EABC-Signaturen?}\;}$$

**Label:** Stufe 5 = **Conjecture**; Stufe 6 = **Conjecture** (Zielobjekt).

---

## 6. Verbindung zur globalen Invariante $\chi(x)$

### Globale $\chi$ (Invarianzprogramm)

$$\chi(x)=\frac{(E(x)+C(x))-(A(x)+B(x))}{\pi_{>3}(x)}
=\mathcal I_\chi(S(x)),\qquad
\mathcal I_\chi(e,a,b,c)=(e+c)-(a+b).$$

Misst die **chirale Asymmetrie** aller Primzahlen $>3$ in den vier EABC-Restklassen.

### Vergleich $\chi_E$ vs.\ $\chi$

| Observable | Raum | Formel | Was misst sie? |
|------------|------|--------|----------------|
| $\chi(x)$ | alle Primzahlen $>3$ | $((E+C)-(A+B))/\pi_{>3}$ | globale EC-vs-AB-Bilanz |
| $\chi_E(N)$ | Vierlingsstarts $p\le N$ | $(\#\mathrm{ABCE}-\#\mathrm{CEAB})/(\#\mathrm{ABCE}+\#\mathrm{CEAB})$ | **Orientierungsasymmetrie** (Stufe 4) |
| $\chi_{\mathrm{leg}}(Q)$ | vier Beine eines Vierlings | $(\#E+\#C)-(\#A+\#B)$ | **$0$** für kanonische Vierlinge (**Theorem**) |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | Klassenprodukt | — | **nichts** (immer $0$) |

**Theorem (balancierte Beine).** Für jeden kanonischen Prim-Vierling $Q$:
$$\chi_{\mathrm{leg}}(Q)=0\quad\text{für ABCE und CEAB}.$$

**Experiment (ehrlicher Vergleich).** $\chi_E(N)$ und $\chi(x)$ sind **verwandte chirale Observablen in verschiedenen Räumen** — nicht identisch. $\chi_E$ codiert die ABCE/CEAB-**Orientierung** (Stufe 3–4); $\chi$ trackt die globale Primzählung. Lean- und Python-Implementierung stimmen für $\omega$ überein (`tests/test_eabc_from_lean.py`, `tests/test_eabc_holonomie.py`).

---

## 7. Oktanion-Ebene: algebraischer Assoziator vs. projektive Holonomie (Ziel)

| Größe | Domäne | Typ | Referenz |
|-------|--------|-----|----------|
| $[x,y,z]=(xy)z-x(yz)$ | $\mathbb O$ | **algebraische** Nichtassoziativität | `collatz_eabc_oktonion_associator.py` |
| $\mathcal H_E((xy)z,x(yz))=d_E(\Gamma_E((xy)z),\Gamma_E(x(yz)))$ | $\mathbb{Z}^8$-Stub | **projektive Holonomie** (Zielobjekt) | `eabc_associator_vector` |
| $\mathfrak a_{\mathrm{naiv}}$ auf $V_4$ | $\{E,A,B,C\}$ | **immer $0$** | `collatz_eabc_discrete_associator.py` |
| $\omega(Q)$, $\chi_E(N)$ | Prim-Vierlinge | Orientierung / Chiralität (Stufe 3–4) | dieses Dokument |

Der **algebraische** Oktanion-Assoziator $[x,y,z]$ und die **projektive** Holonomie $\mathcal H_E$ sind **verschiedene Ebenen**. Auf $V_4$ kollabiert $\mathfrak a_{\mathrm{naiv}}$; in $\mathbb O$ und auf Vierlingen bleibt messbare Struktur — aktuell als Orientierung/Chiralität (Stufe 3–4), langfristig als Holonomie (Stufe 6).

---

## 8. Epistemische Tabelle

| Aussage | Stufe | Label |
|---------|-------|-------|
| $\kappa(p)\in\{E,A,B,C\}$ | 1 | **Definition** |
| $\sigma(Q)\in\Sigma_4$ | 2 | **Definition** |
| $\omega(\mathrm{ABCE})=+1$, $\omega(\mathrm{CEAB})=-1$ | 3 | **Definition** |
| $\chi_E(N)$ Vierlings-Chiralität | 4 | **Definition** |
| $T:\,Q_k\to Q_{k+1}$ | 5 | **Conjecture** (offen) |
| $\mathcal H_E(\gamma_1,\gamma_2)$ geschlossene Transportpfade | 6 | **Conjecture** (Zielobjekt) |
| Welcher natürliche Transportoperator? | 5 | **Forschungsfrage** |
| $V_4$ Klein-Gruppe, $\Phi$ assoziativ | — | **Theorem** |
| $\mathfrak a_{\mathrm{naiv}}\equiv 0$ auf $V_4$ | — | **Theorem** |
| $\mathcal H_E(\gamma_1,\gamma_2)=d_E(\Gamma(\gamma_1),\Gamma(\gamma_2))$ | 6 | **Definition** (Zielrahmen) |
| $\chi_{\mathrm{leg}}(Q)=0$ für kanonische Vierlinge | — | **Theorem** |
| $\chi_E$ vs.\ globale $\chi$ — verwandt, nicht identisch | 4 | **Experiment** |
| $[x,y,z]$ auf $\mathbb O$ | — | **Theorem** (algebraisch, eigene Ebene) |
| ABCE/CEAB-Clustering als Holonomie | — | **verworfen** — Vorläufer, nicht Holonomie |

---

*Kanonsiche Notiz: Der diskrete $V_4$-„Assoziator“ ist eine **Negativkontrolle** — er bestätigt, dass Klammerdefekte **erst nach Hebung/Projektion** sichtbar werden. Operative Forschung heute: $\omega$, $\chi_E$ (Stufe 3–4); langfristig: Transportoperator $T$ (Stufe 5) und $\mathcal H_E$ auf geschlossenen Pfaden (Stufe 6), parallel zur globalen $\chi$.*
