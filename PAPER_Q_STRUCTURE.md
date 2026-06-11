# Rationale Eigenraumstruktur von Π_Γ (Paper-Abschnitt)

Operator über **ℚ⁴** aus dem holonomen Rundweg **Γ = R·K** (Permutationsmatrizen der Klein-Vierer-Gruppe):

\[
\Pi_\Gamma = \frac{1}{4}\bigl(I + \Gamma + \Gamma^2 + \Gamma^3\bigr).
\]

Sage-Verifikation: `./run_rundweg.sh` (Block „ARITHMETISCHE STRUKTUR …“).

## Spektraltabelle über ℚ

| Eigenwert λ | Dimension | Basis (Spaltenkonvention) |
|-------------|-----------|---------------------------|
| 0 | 1 | (0, 1, 0, −1) |
| 1 | 3 | (1, 0, 0, 0), (0, 1, 0, 1), (0, 0, 1, 0) |

**Interpretation:** Der Eigenraum λ = 0 ist der antisymmetrische „Fluktuations“-Träger; λ = 1 ist der glatte Projektionsanteil (Kollaps auf **S_glatt**).

## Nächster Lift (Paper)

- **Pfad A:** Frobenius-Resonanz über **ℚ(ζ₅)** — natürlich zur 5-Falt-Symmetrie und zu `Global Lokal.py` / EABC-Klassen mod 12.
- **Pfad B:** Hurwitz-Ideale (5, 7) — siehe `Gabriel Hurwitz.tex`, `rh scholze.tex`, `divisionsalgebren.tex`.
  Verifikation: `./run_hurwitz_resonanz.sh`, Doku `PAPER_HURWITZ_RESONANZ.md`, LaTeX `paper_hurwitz_snippet.tex`.

Empfehlung im Haupttext: zuerst diese ℚ-Tabelle zitieren, dann Pfad A für den Frobenius-Lift oder Pfad B für Hurwitz-Resonanz (siehe Eltern-Antwort).

## Lift $\mathbb{Q}(\zeta_5)$ und Galois

Verifikation: `./run_lift_zeta5.sh` (`rundweg_lift_zeta5.sage`).

Lift $K = \mathbb{Q}(\zeta_5)$, $\zeta = e^{2\pi i/5}$. Die Rundweg-Matrizen $R,K,\Gamma$ und
$\Pi_\Gamma = \frac{1}{4}(I+\Gamma+\Gamma^2+\Gamma^3)$ werden koeffizientweise nach $K$ gehoben.
**Spektrum unverändert:** dieselben Eigenwerte $\lambda \in \{0,1\}$, dieselben Dimensionen $1$ und $3$;
die $\mathbb{Q}$-Basen liegen in den $K$-Eigenräumen (keine Spaltung beim ersten Lift).

Galois-Gruppe $\mathrm{Gal}(K/\mathbb{Q}) \cong (\mathbb{Z}/5\mathbb{Z})^\times$ (Ordnung $4$).
Ein erzeugendes $\sigma$ wirkt auf dem Körper nichttrivial: $\sigma(\zeta) = \zeta^2$.
(Hinweis: $p=5$ ist hier **kein** Frobenius-Test — $\zeta \mapsto \zeta^5 = 1$ ist kein Automorphismus von $K$.)

Auf dem **4D-Modul** $V = K^4$ (Koordinaten wie im Rundweg) ist die Sage-Wirkung koeffizientweise:
$\sigma(v_1,\ldots,v_4) = (\sigma(v_1),\ldots,\sigma(v_4))$. Dann gilt:

- $\sigma(\Pi_\Gamma) = \Pi_\Gamma$ (Matrix über $\mathbb{Q}$),
- $\sigma$ fixiert punktweise $\mathrm{im}(\Pi_\Gamma)$ und den $\lambda=0$-Eigenraum mit Basis $(0,1,0,-1)^\top$,
- alle verwendeten Eigenraum-Basen sind $\mathbb{Q}$-rational $\Rightarrow$ unter $\sigma$ **trivial auf Vektoren**.

Die nichttriviale Galois-Information steckt zunächst nur in $K$, nicht in der gewählten $\mathbb{Q}$-Darstellung von $V$.

## Arithmetischer Resonanzpunkt (Paper-Definition)

**Vorschlag (präzise):** Ein Punkt $x$ im erweiterten Modul $V_{\mathrm{arith}} \supset \mathbb{Q}^4$
ist ein *arithmetischer Resonanzpunkt*, wenn gleichzeitig

1. $x \in \mathrm{Fix}(\Pi_\Gamma)$ (glatter bzw. projektiver Anteil: $\Pi_\Gamma x = x$), und
2. $x$ liegt in einem Untermodul, auf dem $\mathrm{Gal}(K/\mathbb{Q})$ **nicht** punktweise trivial wirkt
   (echte $\zeta$-Mischung in mindestens einer Koordinate).

Äquivalent (operativ): Schnittmenge des glatten Eigenraums $\lambda=1$ mit dem Komplement der
$\mathbb{Q}$-rationalen Spalte — dort, wo Koordinaten als $\zeta$-Linearkombinationen identifiziert sind
und $\sigma$ nicht die Identität auf $x$ ist.

**Warum der aktuelle $\sigma$-Test noch „leer“ ist:** $V_{\mathrm{arith}} = K^4$ mit Standard-Lift
enthält nur $\mathbb{Q}$-Vektoren in den Eigenräumen. Für $v \in \mathbb{Q}^4 \hookrightarrow K^4$ gilt
$\sigma(v) = v$, weil $\sigma|_{\mathbb{Q}} = \mathrm{id}$. $\Pi_\Gamma$ und seine Eigenprojektionen
sind über $\mathbb{Q}$ definiert — Galois sieht die Rundweg-Symmetrie als **bereits arithmetisch abgeschlossen**
auf diesem Modul; Resonanz im Sinne von (2) tritt erst auf, wenn $V$ mit zyklotomischer Koordinatenwahl
angereichert wird.

**Nächster nichttrivialer Test (Pfad A, Schritt 2):**

- Identifikation der Achsen $E,A,B,C$ mit $\zeta^1,\zeta^2,\zeta^3,\zeta^4$ (oder $1,\zeta,\zeta^2,\zeta^3$)
  statt reiner $\mathbb{Q}$-Koordinaten; Modul z.\,B. $V_{\mathrm{arith}} \cong K$ oder $K^4$ mit
  Einbettung über Vandermonde / zyklotomische Spur.
- Frobenius bei ungeradem Prim $p \neq 5$: $\varphi_p(\zeta) = \zeta^p$ auf $K$; prüfen, ob
  $\varphi_p$ auf einer **gemischten** Basis von $\mathrm{Fix}(\Pi_\Gamma)$ nichttrivial wird und ob
  $\varphi_p(\Pi_\Gamma) = \Pi_\Gamma$ erhalten bleibt.
- Kriterium für „Resonanz sichtbar“: existiert $x \in \ker(\Pi_\Gamma - I)$ mit
  $\varphi_p(x) \neq x$ und $\varphi_p(\Pi_\Gamma) = \Pi_\Gamma$ (Galois-Äquivarianz des Operators,
  aber nicht punktweise Fixierung aller glatten Vektoren).

LaTeX zum Einfügen: `paper_resonanz_snippet.tex`.
