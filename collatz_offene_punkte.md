# Collatz: offene Punkte (Stand Juni 2026)

Synthese aus `collatz_schlussartikel_arxiv.tex` (Epilog, Uniformität, Lean §9.2),
`collatz_z2_attraktor.lean`, `collatz_uniformity.lean`, `collatz_density_appendix.lean`
und `collatz_z2_attraktor.tex` §8 Roadmap.

**Ehrliche Einordnung:** Die Collatz-Vermutung ist **nicht** bewiesen. Dieses Dokument
listet die verbleibende Lücke zwischen gesicherter Strukturtheorie und punktweiser
Konvergenz.

---

## 1. Kernvermutung (Collatz selbst)

- **Klassische Form:** Jede Collatz-Bahn $n, \mathrm{Col}(n), \mathrm{Col}^2(n), \ldots$
  erreicht nach endlich vielen Schritten den Zyklus $(1,2,4)$.
- **Odd-to-odd-Form:** Für jedes ungerade $n\in\mathbb{N}$ existiert $K<\infty$ mit
  $U^{K}(n)=1$ (Abbildung $U$ aus `collatz_z2_attraktor.lean`).
- **Echte Ausnahmemenge:** $E_\infty=\{n\in\mathbb{N}_{\mathrm{odd}}:\forall K,\;U^K(n)\neq 1\}$;
  Collatz $\Leftrightarrow$ $E_\infty=\emptyset$ (`ExceptionSetInfinity` in `CollatzEabc.Open`).
- **Nicht äquivalent:** $E_{\mathrm{diag}}\cap\mathbb{N}=\emptyset$ mit
  $E_{\mathrm{diag}}=\overline{\bigcup_N E_{N,N}}$ (`ExceptionSet` / `ExceptionSetDiag` in Lean) —
  Gegenbeispiel $n=27$ vor Erreichen von $1$; siehe `collatz_equivalenz_e_infty.tex`.
- **Status:** Offen. Tao (2019) beweist nur logarithmische Dichte $1$ für
  „fast bounded“-Bahnen — **nicht** punktweise Konvergenz.
- **Zyklen $>1$:** Bekannte Ausschlüsse modulo $2^k$ für endliches $k$; „für alle $k$“
  ist nicht dasselbe wie „für alle $n$“.

---

## 1a. Abgrenzung: ABCE/CEAB-Drehung (bewiesen)

- **Status: bewiesen** — kein offener Punkt. Markierte Primvierlinge
  $Q(p)=(p,p+2,p+6,p+8)$ werden modulo $12$ eindeutig als $\mathrm{ABCE}$ bzw.
  $\mathrm{CEAB}$ klassifiziert ($p\bmod{12}\in\{5,11\}$).
- **CEAB** ist die zyklische Verschiebung von **ABCE** auf dem Flavor-Ring
  $E\to A\to B\to C\to E$ (Operator $T$ in `EABC.lean`, `chiralityOrder`).
- **Belege:** `EABC.lean` ($T^4=\mathrm{id}$), `Projektionszeuge.tex`
  (Projektionszeuge, Startkanten), `collatz_hurwitz_polytop_eabc.tex`
  (Satz mod-$12$-Vierling), `Miller_alt.tex` ($C_4$-Rotation der Primlücken).
- **Nicht vermischen:** Diese Drehung betrifft mod-$12$-Chiralität und
  Primvierlings-Kodierung — **nicht** die offene Collatz-Brücke
  „endliche EABC-Grammatik $\Rightarrow$ $E_\infty=\emptyset$“
  (vgl. `collatz_equivalenz_e_infty.tex`).

---

## 2. Uniformitätslücke (Stufe E: $\mathrm{dist}_2(T^k(n),E)\to 0$)

- **Zentrale offene Frage:** Negativer mittlerer Drift ($\Lambda\approx -0.830$) und
  mod-12-Mischung implizieren **nicht**, dass jede einzelne Trajektorie konvergiert.
- **Präzise Formulierung (Vermutung):** Für jedes $n\in\mathbb{N}$ existiert $K(n)<\infty$,
  sodass die odd-to-odd-Bahn den trivialen Attraktorbereich trifft; geometrisch:
  $\mathrm{dist}_2(U^k(n),E)\not\to 0$ entlang einer nicht-konvergierenden Bahn.
- **Lean-Objekt:** Stufe E in `collatz_z2_attraktor.lean` — Stufen A–D sind ohne `sorry`;
  globale Uniformität fehlt (`collatz_uniformity_conjecture`).
- **Zwischenstufe:** Uniformitäts-Vermutung ist strukturell **schwächer** als voller Collatz
  (aus Uniformität allein folgt Konvergenz nach $1$ nicht automatisch), aber **stärker** als
  Tao/Birkhoff auf Blockmodellen.
- **Widerlegte Heuristik:** Schranken der Form $\mathrm{dist}_2(n,1)\geq c\cdot 2^{-\log n}$
  scheitern an der LTE-Familie $n_0=2^{k+1}3^r-1$ (`dist_to_one_not_uniform_bound`).

---

## 3. Attraktor $E$ (noch zu beweisen: „richtig“, nicht nur finite Approximation)

- **Definition (Lean):** `ExceptionSet = closure (⋃_N ExceptionSetApprox N N)` —
  2-adische Hülle der Starts $n\leq N$, deren Bahn $U^k(n)\neq 1$ für alle $k\leq N$.
- **Abgrenzung:** Collatz $\Leftrightarrow$ $E_\infty=\emptyset$ (präzise); $E_{\mathrm{diag}}$ ist
  nur der 2-adische Beobachtungsschatten — nicht dieselbe Aussage.
- **Offene Fragen (TeX §8 F1–F5):**
  - Ist punktweise Uniformität $\mathrm{dist}_2(U^k(n),E_{\mathrm{diag}})\to 0$ äquivalent zu Collatz
    oder strikt schwächer?
  - Stimmt der Limes $E$ mit der kleinsten abgeschlossenen $U$-invarianten Hülle von
    $A_{\mathrm{triv}}$ überein?
  - Ist $E$ wirklich der „richtige“ Attraktor, oder nur ein endlich approximierter Limes
    schlechter Präfixe?
- **Kugel vs. Zyklus:** $\|x-1\|_2<1$ enthält **alle** ungeraden $n\in\mathbb{N}$
  (`odd_nat_in_trivial_ball`); daher verwendet die Approximation Erreichen von $U^k(n)=1$,
  nicht bloße Kugelnähe.
- **TeX-Konvention:** $\mathrm{dist}_2(n,\emptyset)=1$ vs. Mathlib `sInf ∅ = 0` —
  noch nicht eingebaut.

---

## 4. Numerisch / heuristisch ($\Lambda$, Mixing $W\leq 70$, RPF)

- **Negativer Drift:** $\mathbb{E}[\log_2 F(B_k)]=\Lambda<0$ ist Eigenschaft des
  **stationären** mod-12-Transfer-Operators — nicht der einzelnen Trajektorie.
- **$\Lambda\approx -0.830$:** Numerisch/heuristisch gestützt; kein Lean-Beweis der
  Stationarität oder Ergodizität der echten Collatz-Dynamik.
- **mod-12-Mischung:** $|\lambda_2|\approx 0.35$, $p=\min M_{st}>0$;
  $P(\text{kein E-Besuch in }n\text{ Blöcken})\leq(1-p)^{\lfloor n/12\rfloor}$
  (`exception_probability_tendsto_zero` in `collatz_uniformity.lean`).
- **Mixing-Zeit $W(n)\leq 70$:** Numerisch untersucht; nicht formalisiert.
- **RPF-Operator:** Irreduzibilität und Spektrallücke $|\lambda_2|<1$ — endliche
  Kombinatorik noch nicht als Lean-Theorem (`Finset`-Matrix).
- **Informationszerreibung:** Heuristik $I_n(k)\sim(0.35)^k$; rigorose Definition als
  bedingte Entropie $H(B_{k+1}\mid B_1,\ldots,B_k)$ offen.
- **Bernoulli-Normschale:** Lyapunov-Nutzung ist No-Go (bewiesen).
- **$\Phi_{\mathrm{pref}}$ (wohldefiniert):** Vorabbildung auf EABC-Wörtern
  $\Phi_{\mathrm{pref}}:\{E,A,B,C\}^{<\omega}\to\mathbb{C}\times\mathbb{R}$ mit
  $z(w)=z_0+2^{-|w|}u(w)$, $t(w)=|w|$; Lean-Kern in
  `collatz_eabc_core/CollatzEabc/PrefProjection.lean` (Phasen als $\mathbb{Z}\times\mathbb{Z}$).
  **Kein Collatz-Anspruch** — die dynamische Brücke $\Phi=\Phi_{\mathrm{pref}}\circ\kappa$
  bleibt offen (`collatz_kepler_gedankenexperiment.tex`, Abschnitt~Φ\_pref).
- **Geometrisches Gedankenexperiment:** Spekulative Skizze
  $\Phi:\mathbb{Z}_2\to\mathcal{M}=\mathbb{C}\times\mathbb{R}$,
  Zwölfer-Kreis $K_{12}$, Stirling--Bernoulli--Zeta-Anbindung und
  $\mathcal{E}_\infty\subset\mathcal{M}$ in
  `collatz_kepler_gedankenexperiment.tex` (kein Beweisanspruch).
- **Bernoulli-Uhr (Gedankenmodell):** Auf zwei chiralen Ellipsenachsen
  $\mathcal{E}_\pm$ (ABCE/CEAB) ordnet jede Zelle $m$ das Tripel
  $(B_{2m-2},B_{2m},B_{2m+2})$ mit $r_m=2^{-m}$, $\theta_m=\pi m/2$ und
  $B_{2m}=-2m\,\zeta(1-2m)$ zu; $U_\pm(m)$ ist eine chirale Gewichtsformel
  (TeX §Bernoulli-Uhr). Lean-Kern: `BernoulliClock.lean`
  (`BernoulliCell`, `bernoulliTriplet`). **Kein Collatz-Beweis**, keine
  Lyapunov-Nutzung (Bernoulli-Normschale bleibt No-Go).

---

## 5. Mathlib-Lücken (Kingman, Birkhoff, etc.)

Aus `collatz_uniformity.lean` (Kommentarblock) und Epilog §9.2:

| Thema | Status |
|-------|--------|
| `padicValNat`, `PadicInt 2`, Ultrametrik, `closure` | **In Mathlib, genutzt** |
| Geometrische Reihe / Tail-Summen | **Genutzt** (`tail_series_formula`) |
| `emultiplicity_*` (volles LTE $p=2$) | Mathlib vorhanden, **noch nicht eingebunden** |
| Kingman-Subadditivität | **Fehlt in Mathlib** |
| Punktweiser Birkhoff-Satz | **Fehlt** (nur fast-sicher) |
| Collatz-Ergodizität | **Fehlt** |
| Vollständiger Perron-Frobenius | **Fehlt** |
| von-Staudt-Clausen (v4.29) | **Fehlt / unvollständig** |
| Asymptotische `Nat.density` | **Fehlt** |
| mod-12-Irreduzibilität als endliche Matrix | **Priorisiert, offen** |

Wo Mathlib endet, endet der formalisierte Teil — dokumentiert, nicht verschwiegen.

---

## 6. Lean-technisch (Lake-CI, sorry=0 aber Stufe E fehlt)

| Datei | sorry | Inhalt |
|-------|-------|--------|
| `collatz_density_appendix.lean` | 0 | C-Ketten-Dichte, Tail-Reihe, `exception_probability_tendsto_zero` |
| `collatz_uniformity.lean` | 0 | LTE-Reset, Mischschranken, Dichte-Duplikat |
| `collatz_z2_attraktor.lean` | 0 (Stufen A–D) | Metrik, `ExceptionSetApprox`, Monotonie, $U$-Invarianz |
| `collatz_uniformity_e.lean` | gezielt | Stufe E: Vermutung + Beweisversuche mit kommentierten `sorry` |

- **Build:** `cd ptolemaeus-lean/Ptolemaeus && lake env lean ../../collatz_z2_attraktor.lean`
  (Mathlib 4.29 via Ptolemaeus-Lake).
- **CI:** Collatz-Lean-Dateien liegen im Repo-Root, nicht im Ptolemaeus-Modul —
  automatische Lake-CI prüft sie derzeit **nicht** bei jedem Push.
- **Paradox:** `sorry=0` in Kernfiles suggeriert Vollständigkeit, obwohl Stufe E —
  der eigentliche Collatz-relevante Schritt — offen ist.

---

## Proof Landscape (Epilog, Juni 2026)

```
[ Irreduzibler mod-12 RPF-Operator ]
                 ↓
     [ Stationäres Maß π ]
                 ↓
[ Negativer Drift Λ ≈ −0.830 ]  (numerisch/heuristisch)
                 ↓
[ Informationszerreibung / Mischung W(n)≤70 ]
                 ↓
          Brücke fehlt:
  Uniformität über E ⊂ ℤ₂ (nicht dist zu 1)
                 ↓
       [ Volle Konvergenz — OFFEN ]
```

---

## Priorisierte Lean-Roadmap (Anhang)

1. Irreduzibilität der mod-12-Matrix (`Finset`-Kombinatorik)
2. Spektrallücke $|\lambda_2|<1$ (exakte rationale Matrix)
3. Uniforme Zweiblock-Kompensation
4. Roth-Analogon in $\mathbb{Z}_2$
5. Punktweise Collatz-Konvergenz (äquivalent zur Vermutung)

---

## Referenzen

- `EABC.lean`, `Projektionszeuge.tex`, `collatz_hurwitz_polytop_eabc.tex`, `Miller_alt.tex` — ABCE/CEAB (bewiesen)
- `collatz_equivalenz_e_infty.tex` — $E_\infty$ vs. $E_{\mathrm{diag}}$, offene Brücke
- `collatz_schlussartikel_arxiv.tex` — §Uniformität, Epilog, Anhang Lean §9.2
- `collatz_z2_attraktor.tex` — §8 Roadmap, §Offene Fragen
- `collatz_z2_attraktor.lean` — Stufen A–E
- `collatz_uniformity_e.lean` — Beweisversuche Stufe E
- `collatz_kepler_gedankenexperiment.tex` — geometrisches Gedankenexperiment ($\Phi$, $K_{12}$, Stirling/Bernoulli; spekulativ)
