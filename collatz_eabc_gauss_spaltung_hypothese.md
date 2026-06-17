# Gauß–EABC-Spaltungshypothese (glatt-EABC)

**Status:** Forschungshypothese + Experiment  
**Kanonsiche Erweiterung von:** `collatz_eabc_normabstieg_hypothese.md` (§8–§9)  
**Experiment:** `collatz_eabc_gauss_spaltung_test.py` → `collatz_eabc_gauss_spaltung.json`  
**Vorgänger (korrigiert):** `collatz_eabc_gauss_faktor_eabc_test.py` — testete rohe mod-$12$-Legs
(höchstens eine EABC-sichtbare Leg pro Paar); **dieser Text ersetzt ihn als kanonischer Test.**

---

## 1. Definition (Gauß-Spaltung in $\mathbb{Z}[i]$)

**Definition 1 (split-Primzahl).** Für eine rationale Primzahl $p>3$ mit $p\equiv 1\pmod 4$
existiert eine **eindeutige** Darstellung
$$p = a^2 + b^2,\qquad 0 < a \le b,$$
und $p$ **zerlegt** in den Gaußschen Ganzzahlen:
$$p = (a+bi)(a-bi)\quad\text{in }\mathbb{Z}[i].$$

**Label:** **Theorem** (klassische algebraische Zahlentheorie).

---

## 2. Definition (glatt-EABC-Zerlegung)

**Definition 2 (glatter Kern).** Für $n\in\mathbb{N}_{>0}$ schreibe eindeutig
$$n = 2^{\alpha}\,3^{\beta}\, n',\qquad \gcd(n',6)=1.$$
Schreibe $(\alpha,\beta,n') = \mathrm{strip\_smooth}(n)$.

**Definition 3 ($\kappa$-Klassifikation).** Für $n'\equiv 1,5,7,11\pmod{12}$:
$$\kappa(n') \in \{E,A,B,C\},\quad E\equiv 1,\; A\equiv 5,\; B\equiv 7,\; C\equiv 11 \pmod{12}$$
(via `eabc_from_lean.class_of`).

**Lemma.** Ist $\gcd(n',6)=1$, so ist $\kappa(n')$ stets definiert
($n'\bmod 12 \in \{1,5,7,11\}$).

**Label:** **Definition** (EABC-Programm); Lemma = elementare Kongruenzrechnung.

---

## 3. Definition ($\Gamma$-Paar und Observablenraum)

Für split $p=a^2+b^2$ setze
$$\Gamma(p) := \bigl(\kappa(a'),\,\kappa(b')\bigr) \in \{E,A,B,C\}^2,$$
wobei $a=2^{\alpha_a}3^{\beta_a}a'$, $b=2^{\alpha_b}3^{\beta_b}b'$.

**Prinzip $N=(N_{\mathrm{glatt}},N_{\mathrm{EABC}})$:** Die glatten Faktoren $2^\alpha 3^\beta$
werden **vor** der EABC-Orientierung abgetrennt; nur die Kerne $a',b'$ tragen $\kappa$.

**Observablenraum:** $\Sigma := \{E,A,B,C\}^2$, $|\Sigma|=16$.

---

## 4. Forschungsfrage (Experiment)

Für $X\to\infty$, split-Primzahlen $p\le X$:
$$\mu_X(\gamma) := \frac{\#\{p\le X \text{ split} : \Gamma(p)=\gamma\}}{\#\{p\le X \text{ split}\}}, \quad \gamma\in\Sigma.$$

**Frage:** Ist $\mu_X$ **asymptotisch uniform** auf $\Sigma$ (je $1/16$), oder zeigen sich
**stabile Biases** / Orientierungspräferenzen?

---

## 5. Conjecture (Spaltungs-Orientierung)

$$\boxed{
\text{Die Gauß-Spaltung } p=(a+bi)(a-bi) \text{ trägt eine nichttriviale EABC-Orientierung }
\Gamma(p) \text{, die aus } p\equiv 1\pmod 4 \text{ allein nicht folgt.}
}$$

**Lesart:** Wenn $\mu_X(\gamma)\not\equiv 1/16$ stabil abweicht, ist $\mathbb{Z}[i]\to$EABC
eine **reelle** Brücke jenseits der trivialen bipartiten Zuordnung (§8 in
`collatz_eabc_normabstieg_hypothese.md`).

**Label:** **Conjecture** — empirisch zu prüfen; Uniformität **schwächt** die Konjektur.

---

## 6. Experiment (Tao-Stil)

**Setup.** Für $X\in\{10^4,10^5,10^6\}$ (erweiterbar):

1. Alle $p\le X$, $p\equiv 1\pmod 4$, $p>3$.
2. Kanonisches $(a,b)$; `strip_smooth` auf beiden Legs.
3. Zähle $\Gamma(p)$ über alle $16$ Klassen.
4. Schätze $\mu_X(\gamma)$; $\chi^2$ gegen Uniform $1/16$ (15 Freiheitsgrade).
5. Maximalabweichung $|\mu_X(\gamma)-1/16|$; Top-$k$-Abweichungen.
6. **Shuffle-Null:** marginals-erhaltende Permutation von $\kappa(a')$ und $\kappa(b')$
   unabhängig; Vergleich der beobachteten $\chi^2$ mit der Nullverteilung.

**Ehrliche Auswertung:**

| Ergebnis | Epistemik |
|----------|-----------|
| $\mu_X\approx 1/16$, $\chi^2$ klein | Konjektur **nicht gestützt** |
| stabile Abweichung + signifikante $\chi^2$ vs. Shuffle-Null | **Anker** für $\mathbb{Z}[i]\to$EABC |

---

## 7. Verhältnis zu früheren Tests

| Test | Objekt | Problem |
|------|--------|---------|
| `collatz_eabc_gauss_defekt_test.py` | $\kappa(p)$ vs. split/inert | **arithmetisch trivial** (mod $4\times 12$) |
| `collatz_eabc_gauss_faktor_eabc_test.py` | rohe $(a,b)\bmod 12$ | höchstens **eine** EABC-Leg sichtbar |
| **dieser Test** | $\Gamma(p)=(\kappa(a'),\kappa(b'))$ | beide Legs sichtbar; voller $16$-Klassen-Raum |

---

## 8. Zusammenfassung

Die **Gauß–EABC-Spaltungshypothese** fragt, ob die Idealzerlegung in $\mathbb{Z}[i]$
nach glatt-EABC-Projektion eine messbare Orientierungsverteilung $\Gamma(p)$ auf
$\{E,A,B,C\}^2$ erzeugt. Der kanonische Test ist `collatz_eabc_gauss_spaltung_test.py`.

---

**Eisenstein-Analogon:** `collatz_eabc_eisenstein_spaltung.md` — $\mathbb{Z}[\omega]$, mod-$3$-Bipartition
$E\cup B$ vs. $A\cup C$, Experiment `collatz_eabc_eisenstein_spaltung_test.py`.

---

*Epistemische Einordnung: Gauß-Spaltung = Theorem; $\Gamma$-Bias = offene Conjecture.*
