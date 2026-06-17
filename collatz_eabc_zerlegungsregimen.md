# Primzahlen als Übergänge zwischen EABC-Zerlegungsregimen

**Status:** Forschungshypothese + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Experiment:** `collatz_eabc_Z_decomposition_test.py` → `collatz_eabc_Z_decomposition.json`  
**Tao-Labels:** Definition | Theorem | Conjecture | Heuristik | Experiment

> **Kanonsiche Plattenformulierung:** `collatz_eabc_plattenuebergang.md` — $(\Sigma_n)_{n\ge 1}$ als
> diskrete Normschichtfolge auf $\Lambda_{\mathbb{O}}$; Prim = produkt-irreduzible Übergangsschicht.
> Dieses Dokument enthält die **operativen** $Z(n)$-, $\Delta Z$- und $\mathcal{Z}_{\mathrm{regime}}$-Definitionen
> (Quaternion-Testbed).

**Querverweise:**
- `collatz_eabc_plattenuebergang.md` — **Haupt-Hypothese** (Plattenübergang, boxed Conjecture, scharfe Fragen)
- `collatz_eabc_quaternion_mass_hypothese.md` §1–§5 ($\Sigma_n$, $\mu_n$, $I(\mu_n)$) — **Quaternion-Referenz**
- `collatz_eabc_oktonion_singularitaet.md` §3 (Faktorisierung $n=ab$, $\Sigma_a\star\Sigma_b$) — **8D-Analogie**
- `collatz_eabc_shell_defekt_test.py` (Spektral-Defekt $D(n)$, Prim-Emergenz §12)
- `collatz_eabc_euklidische_hebung.md` (Norm-Multiplikativität $N(xy)=N(x)N(y)$)

---

## 1. Motivation: Primzahlen sind keine speziellen Punkte — aber Regimengrenzen

**Heuristik (Perspektivwechsel).** In der klassischen Lesart sind Primzahlen **irreduzible** Elemente:
keine nichttriviale Faktorisierung $n=ab$ mit $a,b>1$. In der EABC-Schalenlesart ist das fundamentale
Objekt nicht die Punktmenge $|\Sigma_n|$, sondern die **interne Zerlegbarkeit** des Schalenmaßes $\mu_n$.

**Theorem (Norm-Multiplikativität).** Auf Hurwitz-Quaternionen gilt $N(xy)=N(x)N(y)$.

**Heuristik (schematische Schalen-Faltung).** Für $n=ab$ entstehen Punkte in $\Sigma_n$ aus
**Produktstruktur** der Norm — nicht als kartesisches Produkt $\Sigma_a\times\Sigma_b$, aber
schematisch für EABC-Verteilungen:
$$\Sigma_{ab} \;\leadsto\; \mu_{ab} \approx \mu_a \otimes \mu_b \quad\text{(Vier-Bein-Tensorbild).}$$

Primzahlen wären dann **Übergänge zwischen Zerlegungsregimen**: Stellen, an denen die interne
EABC-Dekomposierbarkeit **nicht** aus kleineren Schalen rekonstruierbar ist.

**Label:** Norm-Multiplikativität = **Theorem**; $\mu_a\otimes\mu_b$-Schema = **Heuristik**.

---

## 2. Operational definitions (Hurwitz, $n\le 200$)

Alle Definitionen sind in `collatz_eabc_Z_decomposition_test.py` implementiert.

| Symbol | Definition | Label |
|--------|------------|-------|
| $\mathcal{Z}_{\mathrm{fact}}(n)$ | Anzahl **ungeordneter** Faktorisierungen $n=ab$ mit $a,b>1$ | **Definition** (klassische Baseline) |
| $\mathcal{Z}_{\mathrm{EABC}}(n)$ | Anzahl **verschiedener** Kanäle $(\sigma_a,\sigma_b)$ über alle Faktorisierungen; $\sigma$ = kompakte $\mu$-Signatur $(H,\chi,K_{\mathrm{tr}},|\Gamma|)$ | **Definition** + **Experiment** |
| $\mathcal{Z}_{\mathrm{regime}}(n)$ | $1$ wenn Prim **oder** kein Kanal $\mu_a\otimes\mu_b$ erklärt $I(\mu_n)$ (additiv + Marginal-Heuristik); sonst $0$ | **Definition** + **Heuristik** |
| $Z(n)$ | Primär $Z(n):=\mathcal{Z}_{\mathrm{EABC}}(n)$ | **Definition** |
| $\Delta Z(n)$ | $Z(n)-Z(n-1)$ | **Definition** |

**Kritische Falsifikation:** Wenn $Z(n)\approx\mathcal{Z}_{\mathrm{fact}}(n)$ und $\Delta Z$ nur $\omega(n)$
repackagiert, ist **nichts Neues** gegenüber klassischer Faktorisierung gewonnen.

---

## 3. Experiment: Springt $\Delta Z$ an Primzahlen?

**Experiment.** `collatz_eabc_Z_decomposition_test.py` berechnet für $n\in[2,N]$:
- $Z_{\mathrm{fact}}$, $Z_{\mathrm{EABC}}$, $Z_{\mathrm{regime}}$, $\Delta Z$
- Prim-vs.-zusammengesetzt: Mittel von $|\Delta Z|$ und $Z$
- Pearson-Korrelation mit $\omega(n)$
- Korrelation $Z_{\mathrm{fact}}$ vs.\ $Z_{\mathrm{EABC}}$

**Label:** **Experiment** — explorativ, kleines $n$.

---

## 4. Boxed Conjecture

> **Konjektur (EABC-Zerlegungsregimen).**
> Primzahlen sind **Übergänge zwischen EABC-Zerlegungsregimen**: an $n=p$ ändert sich die
> interne Dekomposierbarkeit von $\mu_n$ qualitativ — messbar durch Sprünge in $\Delta Z(n)$ oder
> $\mathcal{Z}_{\mathrm{regime}}(p)=1$ bei $\mathcal{Z}_{\mathrm{regime}}(ab)=0$ für zusammengesetztes $ab$.
>
> **Epistemische Bedingung:** Die Konjektur ist nur **nicht-trivial**, wenn $Z_{\mathrm{EABC}}$
> von $\mathcal{Z}_{\mathrm{fact}}$ und von $\omega(n)$ **systematisch abweicht**.

---

## 5. Epistemische Einordnung

| Aussage | Label |
|---------|-------|
| $N(xy)=N(x)N(y)$ auf $\mathbb{H}$ | **Theorem** |
| $\Sigma_{ab}\leadsto\mu_a\otimes\mu_b$ (schematisch) | **Heuristik** |
| $Z_{\mathrm{fact}}$, $Z_{\mathrm{EABC}}$, $Z_{\mathrm{regime}}$ | **Definition** |
| Prim = Regimengrenze in $\Delta Z$ | **Conjecture** |
| $Z_{\mathrm{EABC}}\approx Z_{\mathrm{fact}}$ $\Rightarrow$ Falsifikation | **Methodologie** |
| `collatz_eabc_Z_decomposition_test.py` | **Experiment** |

### Verbindung zu verwandten Hypothesen

| Dokument | Bezug |
|----------|-------|
| `collatz_eabc_plattenuebergang.md` | **Kanonsiche** Plattenhypothese ($\Sigma_n$-Folge, boxed Conjecture) |
| `collatz_eabc_quaternion_mass_hypothese.md` §14 | Prim-Emergenz via $D(n)=I(\mu_n)-I_{\mathrm{ref}}$ — **komplementär** (Spektral-Defekt vs. Zerlegungszahl) |
| `collatz_eabc_oktonion_singularitaet.md` §3 | $n=ab$: Faltung $\Sigma_a\star\Sigma_b$; Prim = irreduzibel in Maßstruktur — **8D-Programm** |
| `collatz_eabc_shell_defekt_test.py` | Rolling-$D(n)$: kein Prim-Überhang; $\Delta Z$-Test ist **unabhängiger** Kanal |

---

## 6. Ergebnis (automatisch aus JSON)

Siehe `collatz_eabc_Z_decomposition.json` → Felder `verdict`, `prime_vs_composite`, `omega_correlation`,
`z_fact_vs_z_eabc`.

*Kanonsiche Zerlegungsregimen-Hypothese: Primzahlen als Regimengrenzen der internen EABC-Dekomposierbarkeit von $\mu_n$ — testbar via $Z(n)$ und $\Delta Z(n)$; epistemisch hohl wenn $Z_{\mathrm{EABC}}$ nur $\omega(n)$ repackagiert.*
