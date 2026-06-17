# Gauß–EABC-Spaltungshypothese (glatt-EABC)

**Status:** Forschungshypothese + Experiment  
**Kanonsiche Erweiterung von:** `collatz_eabc_normabstieg_hypothese.md` (§8–§9)  
**Experiment:** `collatz_eabc_gauss_spaltung_test.py` → `collatz_eabc_gauss_spaltung.json`  
**Hurwitz-Orbit (4D):** `collatz_eabc_hurwitz_spaltung.md` → `collatz_eabc_hurwitz_orbit_test.py`  
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

**Lemma (Parität).** Da $p$ ungerade ist, sind $a$ und $b$ von **entgegengesetzter Parität**:
genau eine der beiden Legs ist gerade. Die $2$-adische Valuation $\nu_2\ge 1$ sitzt auf
**genau einer** Seite — dies ist die **glatte Paritätskomponente** und darf in $\Gamma$ nicht
verloren gehen.

**Label:** **Lemma** (elementar).

---

## 2. Definition (glatt-EABC-Zerlegung)

**Definition 2 (glatter Kern).** Für $n\in\mathbb{N}_{>0}$ schreibe eindeutig
$$n = 2^{\alpha}\,3^{\beta}\, n',\qquad \gcd(n',6)=1.$$
Schreibe $(\alpha,\beta,n') = \mathrm{strip\_smooth}(n)$ bzw. kompakt $(\nu_2(n),\nu_3(n),n')$.

**Definition 3 ($\kappa$-Klassifikation).** Für $n'\equiv 1,5,7,11\pmod{12}$:
$$\kappa(n') \in \{E,A,B,C\},\quad E\equiv 1,\; A\equiv 5,\; B\equiv 7,\; C\equiv 11 \pmod{12}$$
(via `eabc_from_lean.class_of`).

**Lemma.** Ist $\gcd(n',6)=1$, so ist $\kappa(n')$ stets definiert
($n'\bmod 12 \in \{1,5,7,11\}$).

**Label:** **Definition** (EABC-Programm); Lemma = elementare Kongruenzrechnung.

---

## 3. Definition (volle $\Gamma$-Signatur)

Für split $p=a^2+b^2$ mit glatter Zerlegung beider Legs:
$$a = 2^{\alpha_a}\,3^{\beta_a}\, a',\qquad b = 2^{\alpha_b}\,3^{\beta_b}\, b',\qquad
\gcd(a',6)=\gcd(b',6)=1,$$
setze die **volle Signatur**
$$\Gamma(p) := \bigl((\alpha_a,\beta_a,\kappa(a')),\;(\alpha_b,\beta_b,\kappa(b'))\bigr).$$

**Kompakte Schreibweise (6-Tupel):**
$$\Gamma(p) = \bigl(\nu_2(a),\,\nu_3(a),\,\kappa(a'),\,\nu_2(b),\,\nu_3(b),\,\kappa(b')\bigr).$$

**Prinzip $N=(N_{\mathrm{glatt}},N_{\mathrm{EABC}})$:** Die glatten Faktoren $2^\alpha 3^\beta$
werden **vor** der EABC-Orientierung abgetrennt; $\nu_2,\nu_3$ bleiben als **explizite**
glatte Koordinaten erhalten, $\kappa$ lebt nur auf den primitiven Kernen $a',b'$.

**Observablenräume:**

| Ebene | Objekt | Größe |
|-------|--------|-------|
| voll | $\Gamma(p)$ als 6-Tupel | endlich pro $X$, wächst mit $X$ (empirisch $\sim 10^3$ bei $X=10^6$) |
| glatt | $(\alpha_a,\beta_a,\alpha_b,\beta_b)$ | Paritätsmuster + $3$-Glätte |
| EABC-Kern | $(\kappa(a'),\kappa(b'))$ | $|\{E,A,B,C\}^2|=16$ |

Die **16er-EABC-Ebene** ist eine **Marginalprojektion**; der kanonische Test behält
$(\nu_2,\nu_3)$ pro Leg.

---

## 4. Forschungsfrage (Experiment)

Für $X\to\infty$, split-Primzahlen $p\le X$:
$$\mu_X(\gamma) := \frac{\#\{p\le X \text{ split} : \Gamma(p)=\gamma\}}{\#\{p\le X \text{ split}\}}.$$

**Fragen:**

1. Existiert $\lim_{X\to\infty}\mu_X$ auf dem beobachteten Signaturraum?
2. Ist die **bedingte** Verteilung von $(\kappa(a'),\kappa(b'))$ bei festem
   $(\alpha_a,\beta_a,\alpha_b,\beta_b)$ asymptotisch uniform ($1/16$)?
3. Sind die **Marginalen** $\kappa(a')$, $\kappa(b')$ nach glatt-strip uniform ($1/4$ je Klasse)?

**Falsifikation:** Wenn $E,A,B,C$ auf den primitiven Koordinaten nach glatt-strip
**marginal uniform** sind und keine stabile bedingte Kopplung über $(\alpha,\beta)$-Muster
besteht → **schwache** Evidenz für eine $\mathbb{Z}[i]\to$EABC-Brücke.
**Stabile Abweichungen** jenseits der Paritäts-/Glätte-Struktur → **reeller Anker**.

---

## 5. Conjecture (Spaltungs-Orientierung)

$$\boxed{
\text{Untersuche die EABC-Verteilung der primitiven Koordinaten in } p=a^2+b^2
}$$

$$\boxed{
\text{Die Gauß-Spaltung } p=(a+bi)(a-bi) \text{ trägt eine nichttriviale EABC-Orientierung }
\Gamma(p) \text{, die aus } p\equiv 1\pmod 4 \text{ und der Paritätsregel allein nicht folgt.}
}$$

**Lesart:** Wenn $\kappa(a'),\kappa(b')$ nach Stratifizierung nach $(\nu_2,\nu_3)$-Mustern
signifikant von Uniformität abweichen, ist $\mathbb{Z}[i]\to$EABC eine **reelle** Brücke
jenseits der trivialen bipartiten Zuordnung (§8 in `collatz_eabc_normabstieg_hypothese.md`).

**Label:** **Conjecture** — empirisch zu prüfen; Uniformität **schwächt** die Konjektur.

---

## 6. Experiment (Tao-Stil)

**Setup.** Für $X\in\{10^4,10^5,10^6\}$ (erweiterbar):

1. Alle $p\le X$, $p\equiv 1\pmod 4$, $p>3$.
2. Kanonisches $(a,b)$; `strip_smooth` auf beiden Legs → volle $\Gamma(p)$.
3. **Paritäts-Check:** genau eine Leg mit $\nu_2\ge 1$.
4. Zähle $\mu_X$ auf vollem 6-Tupel-Raum; berichte $|\Gamma|_{\mathrm{beobachtet}}$.
5. **Marginal:** $\kappa(a'),\kappa(b')$ vs.\ $1/4$; $\chi^2$.
6. **Bedingt:** $(\kappa(a'),\kappa(b'))$ bei festem $(\alpha_a,\beta_a,\alpha_b,\beta_b)$;
   $\chi^2$ vs.\ $1/16$ pro Muster (Schwelle $n\ge 8$).
7. **Shuffle-Nullen:** (i) marginal-erhaltend auf $\kappa$-Paaren;
   (ii) **bedingt** innerhalb jedes $(\alpha,\beta)$-Musters.
8. Top-$k$-Abweichungen auf voller $\Gamma$- und auf $\kappa$-Paar-Ebene.

**Ehrliche Auswertung:**

| Ergebnis | Epistemik |
|----------|-----------|
| $\kappa$-Marginalen $\approx 1/4$, bedingte $\chi^2$ klein | **Falsifikation schwach** |
| erhöhte $\chi^2$ nur durch Randverteilungen | **marginal getrieben**, vorsichtig |
| stabile Abweichung + signifikante bedingte $\chi^2$ vs. Shuffle-Null | **Anker** für $\mathbb{Z}[i]\to$EABC |

---

## 7. Verhältnis zu früheren Tests

| Test | Objekt | Problem |
|------|--------|---------|
| `collatz_eabc_gauss_defekt_test.py` | $\kappa(p)$ vs. split/inert | **arithmetisch trivial** (mod $4\times 12$) |
| `collatz_eabc_gauss_faktor_eabc_test.py` | rohe $(a,b)\bmod 12$ | höchstens **eine** EABC-Leg sichtbar |
| **dieser Test** | volle $\Gamma(p)$ mit $(\nu_2,\nu_3,\kappa)$ | Parität + glatt + beide primitive Legs |

---

## 8. Zusammenfassung

Die **Gauß–EABC-Spaltungshypothese** fragt, ob die Idealzerlegung in $\mathbb{Z}[i]$
nach glatt-EABC-Projektion eine messbare Orientierungsverteilung $\Gamma(p)$ erzeugt —
**nicht** nur das $\kappa$-Paar, sondern mit expliziter $(\nu_2,\nu_3)$-Paritätsstruktur.
Der kanonische Test ist `collatz_eabc_gauss_spaltung_test.py`.

---

*Epistemische Einordnung: Gauß-Spaltung = Theorem; Paritätslemma = Lemma; $\Gamma$-Bias = offene Conjecture.*
