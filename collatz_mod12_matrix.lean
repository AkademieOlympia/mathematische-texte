/-
  Forschungsziel: Irreduzibilität der empirischen mod-12-Matrix aus
  `collatz_mixing_test.py` (B/C-Zeilen mit Nullen, aber stark zusammenhängend).

  Sorry erlaubt — Kernpaket `collatz_eabc_core` enthält den vollständigen Beweis
  für die vereinfachte Matrix in `CollatzEabc.Mod12Matrix`.
-/

import CollatzEabc.Mod12Matrix

open CollatzEabc

/-- Empirische odd-to-odd-Übergangsmatrix ist irreduzibel (Pfad auf `toQuiver`). -/
theorem eabc_mod12_empirical_matrix_irreducible :
    Matrix.IsIrreducible eabcMod12TransitionEmpirical := by
  sorry
