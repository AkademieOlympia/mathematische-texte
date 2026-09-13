/-
  EABC/AlphaEabcBounds.lean

  Formal enclosure of
    α⁻¹_EABC = 137 + √(D_r² - D_t²)
  with D_r = 8/(137 φ), D_t = 96/37675, φ = (1+√5)/2
  in the rational CODATA window (137.03599908, 137.03600000).

  Strategy: monotone rational certificates for √5, hence for φ and D_r,
  then compare radicand squares against ε_min² and ε_max² by `norm_num` over ℚ.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Data.Real.Sqrt
import Mathlib.Tactic

namespace EABC

noncomputable section

/-- Golden ratio. -/
def phi : ℝ := (1 + Real.sqrt 5) / 2

def rk_E8 : ℕ := 8
def mod12_period : ℕ := 12

/-- Trialitätsfaktor `c = 23/8 - 1/96 = 275/96`. -/
def c_tri : ℚ := 275 / 96

theorem c_tri_eq : c_tri = (23 : ℚ) / 8 - 1 / 96 := by
  unfold c_tri
  norm_num

theorem c_tri_den : (96 : ℕ) = rk_E8 * mod12_period := by
  unfold rk_E8 mod12_period
  rfl

def D_r : ℝ := (8 : ℝ) / (137 * phi)

/-- `96 / (275 * 137) = 96 / 37675`. -/
def D_t : ℝ := (96 : ℝ) / 37675

theorem D_t_from_c_tri : D_t = (1 : ℝ) / ((c_tri : ℝ) * 137) := by
  unfold D_t c_tri
  norm_num

def radikand : ℝ := D_r ^ 2 - D_t ^ 2

def alpha_inv_EABC : ℝ := 137 + Real.sqrt radikand

/-- CODATA lower edge `137.03599908`. -/
def codata_lo : ℝ := 137 + (3599908 : ℝ) / 100000000

/-- CODATA upper edge `137.036`. -/
def codata_hi : ℝ := 137 + (36 : ℝ) / 1000

def eps_lo : ℝ := (3599908 : ℝ) / 100000000
def eps_hi : ℝ := (36 : ℝ) / 1000

/-! ### Lemma 1: rational bounds for `√5` and `φ` -/

theorem sqrt5_bounds :
    (2236067977 : ℝ) / 1000000000 < Real.sqrt 5 ∧
      Real.sqrt 5 < (2236067978 : ℝ) / 1000000000 := by
  constructor
  · rw [Real.lt_sqrt (by positivity)]
    norm_num
  · rw [Real.sqrt_lt (by positivity) (by positivity)]
    norm_num

def phi_lo : ℝ := (1 + (2236067977 : ℝ) / 1000000000) / 2
def phi_hi : ℝ := (1 + (2236067978 : ℝ) / 1000000000) / 2

theorem phi_lo_pos : 0 < phi_lo := by
  unfold phi_lo
  positivity

theorem phi_bounds : phi_lo < phi ∧ phi < phi_hi := by
  have h := sqrt5_bounds
  unfold phi phi_lo phi_hi
  constructor <;> linarith [h.1, h.2]

/-! ### Lemma 2: bounds for the radial defect `D_r` -/

def Dr_lo : ℝ := 8 / (137 * phi_hi)
def Dr_hi : ℝ := 8 / (137 * phi_lo)

theorem Dr_lo_pos : 0 < Dr_lo := by
  unfold Dr_lo phi_hi
  positivity

theorem Dr_bounds : Dr_lo < D_r ∧ D_r < Dr_hi := by
  have hp := phi_bounds
  have h137 : (0 : ℝ) < 137 := by norm_num
  have h8 : (0 : ℝ) < 8 := by norm_num
  have hphi : 0 < phi := lt_trans phi_lo_pos hp.1
  unfold D_r Dr_lo Dr_hi
  constructor
  · exact div_lt_div_of_pos_left h8 (mul_pos h137 hphi)
      (mul_lt_mul_of_pos_left hp.2 h137)
  · exact div_lt_div_of_pos_left h8 (mul_pos h137 phi_lo_pos)
      (mul_lt_mul_of_pos_left hp.1 h137)

/-! ### Lemma 3: bounds for the radicand -/

def radikand_lo : ℝ := Dr_lo ^ 2 - ((96 : ℝ) / 37675) ^ 2
def radikand_hi : ℝ := Dr_hi ^ 2 - ((96 : ℝ) / 37675) ^ 2

theorem radikand_bounds : radikand_lo < radikand ∧ radikand < radikand_hi := by
  have hdr := Dr_bounds
  unfold radikand radikand_lo radikand_hi D_t
  constructor
  · nlinarith [hdr.1, Dr_lo_pos]
  · have hDr_pos : 0 < D_r := lt_trans Dr_lo_pos hdr.1
    nlinarith [hdr.2, hDr_pos]

/-! ### Main theorem: exact CODATA enclosure -/

theorem eps_lo_sq_lt_radikand_lo : eps_lo ^ 2 < radikand_lo := by
  unfold eps_lo radikand_lo Dr_lo phi_hi
  norm_num

theorem radikand_hi_lt_eps_hi_sq : radikand_hi < eps_hi ^ 2 := by
  unfold eps_hi radikand_hi Dr_hi phi_lo
  norm_num

theorem radikand_lo_pos : 0 < radikand_lo := by
  have : (0 : ℝ) < eps_lo ^ 2 := by
    unfold eps_lo
    positivity
  linarith [eps_lo_sq_lt_radikand_lo]

theorem alpha_inv_in_codata_window :
    codata_lo < alpha_inv_EABC ∧ alpha_inv_EABC < codata_hi := by
  have hr := radikand_bounds
  have h_rad_pos : 0 < radikand := lt_trans radikand_lo_pos hr.1
  have h_sqrt_lo : eps_lo < Real.sqrt radikand := by
    rw [Real.lt_sqrt (by unfold eps_lo; positivity)]
    linarith [hr.1, eps_lo_sq_lt_radikand_lo]
  have h_sqrt_hi : Real.sqrt radikand < eps_hi := by
    rw [Real.sqrt_lt (le_of_lt h_rad_pos) (by unfold eps_hi; positivity)]
    linarith [hr.2, radikand_hi_lt_eps_hi_sq]
  unfold eps_lo at h_sqrt_lo
  unfold eps_hi at h_sqrt_hi
  unfold alpha_inv_EABC codata_lo codata_hi
  constructor <;> linarith [h_sqrt_lo, h_sqrt_hi]

end

end EABC
