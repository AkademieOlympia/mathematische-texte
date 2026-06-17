# Eisenstein–EABC-Spaltungshypothese (glatt-EABC)

**Status:** Forschungshypothese + Experiment  
**Kanonsiche Erweiterung von:** `collatz_eabc_normabstieg_hypothese.md` (§8–§9b)  
**Gauß-Vergleich:** `collatz_eabc_gauss_spaltung_hypothese.md`  
**Experiment:** `collatz_eabc_eisenstein_spaltung_test.py` → `collatz_eabc_eisenstein_spaltung.json`

---

## 1. Definition (Eisenstein-Spaltung in $\mathbb{Z}[\omega]$)

Sei $\omega = e^{2\pi i/3}$ (primitive dritte Einheitswurzel). Für eine rationale Primzahl $p>3$:

| Bedingung | Verhalten in $\mathbb{Z}[\omega]$ |
|-----------|-----------------------------------|
| $p=3$ | **ramifiziert** ($3 = -\omega^2(1-\omega)^2$) |
| $p\equiv 1\pmod 3$ | **zerlegt** (split): $p = \pi\bar\pi$ |
| $p\equiv 2\pmod 3$ | **träge** (inert) |

**Theorem (Darstellung).** Für $p\equiv 1\pmod 3$, $p\neq 3$, existiert eine eindeutige
(kanonische) Darstellung
$$p = a^2 - ab + b^2,\qquad 0 < a \le b,$$
entsprechend $p = N(a+b\omega) = (a+b\omega)(a+b\bar\omega)$ in $\mathbb{Z}[\omega]$.

**Label:** **Theorem** (klassische algebraische Zahlentheorie).

---

## 2. Vergleich Gauß vs. Eisenstein

| Aspekt | Gauß $\mathbb{Z}[i]$ | Eisenstein $\mathbb{Z}[\omega]$ |
|--------|----------------------|----------------------------------|
| Norm | $a^2+b^2$ | $a^2-ab+b^2$ |
| Split-Bedingung | $p\equiv 1\pmod 4$ | $p\equiv 1\pmod 3$ |
| Inert | $p\equiv 3\pmod 4$ | $p\equiv 2\pmod 3$ |
| Ramifikation | $p=2$ | $p=3$ |
| Natürliche Peano-Schicht | **tetrahedral** $4n\pm 1$ | **hexagonal** $6n\pm 1$ |
| Grobe EABC-Bipartition | $E\cup A$ vs. $B\cup C$ | $E\cup B$ vs. $A\cup C$ |
| Modulo-Grundlage | $\equiv 1$ vs. $\equiv 3 \pmod 4$ | $\equiv 1$ vs. $\equiv 2 \pmod 3$ |

### 2.1 Mod-$3$-Bipartition (exakt, analog §8 Gauß)

EABC-Klassen mod $12$ projizieren auf mod $3$:

| Klasse | Rest mod $12$ | Rest mod $3$ |
|--------|---------------|--------------|
| $E$ | $1$ | $1$ |
| $A$ | $5$ | $2$ |
| $B$ | $7$ | $1$ |
| $C$ | $11$ | $2$ |

Damit gilt für $p>3$ mit definiertem $\kappa(p)$:

$$\boxed{
p\equiv 1\pmod 3 \;\Leftrightarrow\; \kappa(p)\in\{E,B\}
\quad\text{(hexagonal $6n\pm 1$, ungerade $\equiv 1$ mod $3$)},
}$$
$$\boxed{
p\equiv 2\pmod 3 \;\Leftrightarrow\; \kappa(p)\in\{A,C\}
\quad\text{(hexagonal $6n\pm 1$, ungerade $\equiv 2$ mod $3$)}.
}$$

**Wichtig:** Die Paarung $E\cup C$ vs. $A\cup B$ ist **nicht** die mod-$3$-Bipartition
($E\equiv 1$, $C\equiv 2$ gemischt). Die korrekte grobe Zuordnung ist
**$E\cup B$ vs. $A\cup C$** — strukturell parallel zu $E\cup A$ vs. $B\cup C$ bei Gauß.

**Label:** **arithmetisch exakt** (mod $3$ × mod $12$); als geometrische *Bedeutungs*-Brücke
**Heuristik** — dieselbe epistemische Warnung wie in §8 der Normabstiegshypothese.

### 2.2 Bezug zu §18 (hexagonal vs. tetrahedral)

- **Tetrahedral** ($\gcd(n,4)=1$): Reste $1,3,5,7,9,11\bmod 12$ — enthält $E,A,B,C$.
- **Hexagonal** ($\gcd(n,6)=1$): genau $\{1,5,7,11\}\bmod 12 = \{E,A,B,C\}$.

Eisenstein-Spaltung lebt auf der **hexagonalen** Schicht ($6n\pm 1$); Gauß-Spaltung filtert
daraus die **tetrahedralen** $4n+1$-Primzahlen ($E\cup A$). Die Eisenstein-Variante nutzt
damit die **volle** EABC-Vierergruppe in der Split/Inert-Grobklassifikation:
split sieht $E$ **und** $B$, nicht nur $E$ **oder** $A$.

---

## 3. Glatt-EABC-Zerlegung (Eisenstein)

**Definition (glatter Kern).** Wie beim Gauß-Test:
$$n = 2^{\alpha}\,3^{\beta}\, n',\qquad \gcd(n',6)=1,$$
$\kappa(n')\in\{E,A,B,C\}$ via `eabc_from_lean.class_of`.

**Diskussion $2^\alpha 3^\beta$ vs. nur $2^\alpha$:**

| Variante | Begründung | Kanonischer Test |
|----------|------------|------------------|
| $2^\alpha 3^\beta$ abziehen | Gleiche $N=(N_{\mathrm{glatt}},N_{\mathrm{EABC}})$-Trennung wie Gauß; vergleichbare $\Gamma$-Statistik | **ja** (dieser Test) |
| nur $2^\alpha$ abziehen | $3$ ramifiziert in $\mathbb{Z}[\omega]$; $3^\beta$ ist „Eisenstein-glatt“ | optional, nicht kanonisch |

Für $a',b'$ nach `strip_smooth`: stets $\gcd(a',6)=\gcd(b',6)=1$, also $\kappa$ definiert.

---

## 4. Definition ($\Gamma_E$-Paar)

Für split $p = a^2-ab+b^2$ setze
$$\Gamma_E(p) := \bigl(\kappa(a'),\,\kappa(b')\bigr) \in \{E,A,B,C\}^2.$$

**Observablenraum:** $\Sigma := \{E,A,B,C\}^2$, $|\Sigma|=16$.

---

## 5. Forschungsfrage

$$\mu_X(\gamma) := \frac{\#\{p\le X \text{ split in } \mathbb{Z}[\omega] : \Gamma_E(p)=\gamma\}}
{\#\{p\le X \text{ split}\}}, \quad \gamma\in\Sigma.$$

Ist $\mu_X$ asymptotisch uniform ($1/16$), oder zeigen sich stabile Biases jenseits der
trivialen mod-$3$-Grobzuordnung?

**Conjecture (Eisenstein-Spaltungs-Orientierung):** $\Gamma_E(p)$ trägt nichttriviale
EABC-Orientierung, nicht folgend aus $p\equiv 1\pmod 3$ allein.

---

## 6. Experiment

Setup analog `collatz_eabc_gauss_spaltung_test.py`:

1. $p\le X$, $p\equiv 1\pmod 3$, $p>3$.
2. Kanonisches $(a,b)$ via $4p=(2a-b)^2+3b^2$.
3. `strip_smooth` auf beiden Legs; zähle $\Gamma_E(p)$.
4. $\chi^2$ vs. Uniform $1/16$; Shuffle-Null (marginal-erhaltend).
5. **Zusatz:** mod-$3$-Defekt-Check (split/inert vs. $E\cup B$ / $A\cup C$).

---

## 7. Erwartung vs. Gauß (vorläufig)

| Größe | Gauß $\mathbb{Z}[i]$ | Eisenstein $\mathbb{Z}[\omega]$ |
|-------|----------------------|----------------------------------|
| Grobe Bipartition | exakt (mod $4$) | exakt (mod $3$) |
| $\Gamma$-Bias bei $X=10^5$ | $\chi^2$ erhöht, Shuffle-Null schwach | empirisch zu messen |
| Marginalen $\kappa(a')$, $\kappa(b')$ | nicht durch Split allein fixiert | analog |
| Strukturelle Schicht | $4n+1$ (tetrahedral) | $6n\pm 1$ (hexagonal) |

Die hexagonale Schicht ist **feiner** als die tetrahedrale; ob $\Gamma_E$-Biases **stärker**
oder **schwächer** als bei Gauß sind, ist eine **empirische** Frage (siehe JSON).

---

## 8. Zusammenfassung

Die **Eisenstein–EABC-Spaltungshypothese** testet, ob die Idealzerlegung in $\mathbb{Z}[\omega]$
nach glatt-EABC-Projektion eine messbare $\Gamma_E$-Verteilung erzeugt. Die mod-$3$-Bipartition
$E\cup B$ vs. $A\cup C$ ist arithmetisch trivial; der kanonische Nicht-Trivialitätstest ist
$\Gamma_E(p)=(\kappa(a'),\kappa(b'))$ auf $\{E,A,B,C\}^2$.

---

*Epistemische Einordnung: Eisenstein-Spaltung = Theorem; $\Gamma_E$-Bias = offene Conjecture.*
