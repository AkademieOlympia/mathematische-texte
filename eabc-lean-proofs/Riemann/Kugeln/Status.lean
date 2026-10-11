import Mathlib.Analysis.Real.Sqrt
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
Eröffnung: Klärung der Kugeln.

Zweite Bibliothek des Pakets `Riemann`. Dieses Modul importiert keine
Datei aus `Riemann.*`. Gate 3A–3E bleibt eingefroren, und aus ihm folgt
hier nichts über BK-I2.7D oder BK-I2.7E.

Die Zahlen 240, 12 und 18 werden nur als verschiedene natürliche Zahlen
festgehalten. Die reelle Zielzahl `keplerDensity` ist
\(\pi/\sqrt{18}\). Sie ist kein Satz über ein Gitter
\(\Lambda_{\mathrm{BK},s}\).

BK-I2.7D ist drei Schritte, in dieser Reihenfolge. Schritt 2 und Schritt 3
sind erst definiert, wenn Schritt 1 eine Abbildung hingeschrieben hat.

1. Konstruktion von \(S\) und \(\Pi : S \to \mathbb{R}^3\).
2. Diskretheit des Bildes und eine Minimalschale aus zwölf Punkten.
3. Kongruenz der Gram-Matrix dieses Bildes mit der fcc-Kontaktschale.

Der Forschungsbedarf liegt bei Schritt 1. Eine orthogonale Abbildung kann
eine Ikosaederschale nicht in die fcc-Schale führen; das ist in
`Kugeln.GramSeparation` bewiesen. Dort stehen auch die Vielfachheiten der
fcc-Rohschale. Sie sind ein Abnahmezertifikat und keine Abbildung.
\(\Pi\) selbst ist nicht definiert. Der Kandidat \(\Pi_{\mathrm{Im}}\) auf der
binären Ikosaedergruppe ist in `Kugeln.IcosianImag` falsifiziert.

Kein weiterer BK-I2.7D-Claim ohne eine auswertbare Menge und eine
auswertbare Matrix. Die Spaltenzahl folgt aus der Dimension der Quelle.
fcc und hcp werden hier durch die Stapelfolgen \(ABCABC\ldots\) und
\(ABAB\ldots\) unterschieden, und dadurch, dass nur die erste ein
Bravais-Gitter ist. Beide haben lokal zwölf Kontakte und dieselbe
Packungsdichte. Eine Zwölfer-Schale ist damit noch nicht BK-I2.7E.
-/

noncomputable section

namespace Kugeln

/-- Kusszahl von \(E_8\) in Dimension 8 und Kontaktzahl der fcc-Schale
in Dimension 3 sind verschiedene Zahlen. -/
theorem e8Kissing_ne_fccContact : (240 : ℕ) ≠ 12 := by decide

/-- Die Kusszahl 240 ist nicht das Schalenverhältnis 18. -/
theorem e8Kissing_ne_icosianShellRatio : (240 : ℕ) ≠ 18 := by decide

/-- Die fcc-Kontaktzahl 12 ist nicht das Schalenverhältnis 18. -/
theorem fccContact_ne_icosianShellRatio : (12 : ℕ) ≠ 18 := by decide

/-- Zielzahl von BK-I2.7E. Kein Kovolumen und keine Minimaldistanz. -/
noncomputable def keplerDensity : ℝ := Real.pi / Real.sqrt 18

theorem keplerDensity_pos : 0 < keplerDensity := by
  unfold keplerDensity
  refine div_pos Real.pi_pos ?_
  refine Real.sqrt_pos.mpr ?_
  exact Nat.cast_pos.mpr (by decide : 0 < (18 : ℕ))

end Kugeln
