import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Basic.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.Adjugate
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormDet
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

open Matrix BigOperators

namespace E8SelfDuality

/-!
Konkrete E8-Basis in ℝ⁸, Gram-Matrix und die drei vom Compiler geprüften Blöcke
`e8_basis_matrix_det`, `e8_inner_integer`, `unimodular_coordinates_integral`.
-/

noncomputable def e8_basis_matrix : Matrix (Fin 8) (Fin 8) ℝ := !![
  (1 / 2 : ℝ), (-1 / 2 : ℝ), (-1 / 2 : ℝ), (-1 / 2 : ℝ), (-1 / 2 : ℝ), (-1 / 2 : ℝ), (-1 / 2 : ℝ), (1 / 2 : ℝ);
  1,  1,  0,  0,  0,  0,  0,  0;
  -1,  1,  0,  0,  0,  0,  0,  0;
  0, -1,  1,  0,  0,  0,  0,  0;
  0,  0, -1,  1,  0,  0,  0,  0;
  0,  0,  0, -1,  1,  0,  0,  0;
  0,  0,  0,  0, -1,  1,  0,  0;
  0,  0,  0,  0,  0, -1,  1,  0
]

/-- Zeile `i` ist der `i`-te Basisvektor. -/
noncomputable def e8_basis : Fin 8 → (Fin 8 → ℝ) :=
  fun i => e8_basis_matrix i

/-- Das E8-Gitter als ℤ-Spann der Basisvektoren. -/
def E8LatticeSet : Set (Fin 8 → ℝ) :=
  { x | ∃ n : Fin 8 → ℤ, x = ∑ i, (n i : ℝ) • e8_basis i }

/-- Das Dualgitter Λ* = { y | ∀ x ∈ Λ, ⟨x, y⟩ ∈ ℤ }. -/
def latticeDual (Λ : Set (Fin 8 → ℝ)) : Set (Fin 8 → ℝ) :=
  { y | ∀ x ∈ Λ, ∃ m : ℤ, ∑ i, x i * y i = (m : ℝ) }

def e8_gram_z : Matrix (Fin 8) (Fin 8) ℤ := !![
   2,  0, -1,  0,  0,  0,  0,  0;
   0,  2,  0, -1,  0,  0,  0,  0;
  -1,  0,  2, -1,  0,  0,  0,  0;
   0, -1, -1,  2, -1,  0,  0,  0;
   0,  0,  0, -1,  2, -1,  0,  0;
   0,  0,  0,  0, -1,  2, -1,  0;
   0,  0,  0,  0,  0, -1,  2, -1;
   0,  0,  0,  0,  0,  0, -1,  2
]

noncomputable def e8_gram_r : Matrix (Fin 8) (Fin 8) ℝ :=
  fun i j => ∑ k, e8_basis i k * e8_basis j k

theorem e8_gram_eq : e8_gram_r = e8_gram_z.map fun x : ℤ => (x : ℝ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [e8_gram_r, e8_gram_z, e8_basis, e8_basis_matrix, Fin.sum_univ_eight, Matrix.map_apply] <;>
    norm_num

theorem e8_gram_z_det_one : e8_gram_z.det = 1 := by
  unfold e8_gram_z
  eval_det

theorem e8_basis_gram_det_one : e8_gram_r.det = 1 := by
  rw [e8_gram_eq, ← Int.cast_det, e8_gram_z_det_one]
  norm_num

theorem e8_basis_matrix_det : e8_basis_matrix.det = -1 := by
  unfold e8_basis_matrix
  eval_det

lemma e8_basis_linearIndependent : LinearIndependent ℝ e8_basis := by
  unfold e8_basis
  exact Matrix.linearIndependent_rows_of_det_ne_zero (by rw [e8_basis_matrix_det]; norm_num)

lemma e8_inner_integer (nx ny : Fin 8 → ℤ) :
    ∃ m : ℤ,
      ∑ k, (∑ i, (nx i : ℝ) * e8_basis i k) * (∑ j, (ny j : ℝ) * e8_basis j k) = (m : ℝ) := by
  classical
  let m : ℤ := ∑ i, ∑ j, nx i * ny j * e8_gram_z i j
  refine ⟨m, ?_⟩
  have hG : ∀ i j, ∑ k, e8_basis i k * e8_basis j k = (e8_gram_z i j : ℝ) := by
    intro i j
    simpa [e8_gram_r, Matrix.map_apply] using congr_fun (congr_fun e8_gram_eq i) j
  calc
    ∑ k, (∑ i, (nx i : ℝ) * e8_basis i k) * (∑ j, (ny j : ℝ) * e8_basis j k)
        = ∑ k, ∑ i, ∑ j, ((nx i : ℝ) * e8_basis i k) * ((ny j : ℝ) * e8_basis j k) := by
          simp_rw [Finset.sum_mul, Finset.mul_sum]
      _ = ∑ i, ∑ j, ∑ k, (nx i : ℝ) * (ny j : ℝ) * (e8_basis i k * e8_basis j k) := by
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl fun i _ => ?_
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl fun j _ => ?_
          refine Finset.sum_congr rfl fun k _ => ?_
          ring
      _ = ∑ i, ∑ j, (nx i : ℝ) * (ny j : ℝ) * ∑ k, e8_basis i k * e8_basis j k := by
          simp_rw [Finset.mul_sum]
      _ = ∑ i, ∑ j, (nx i : ℝ) * (ny j : ℝ) * (e8_gram_z i j : ℝ) := by
          simp_rw [hG]
      _ = (m : ℝ) := by
          simp [m, Int.cast_sum, Int.cast_mul]

lemma e8_subset_dual : E8LatticeSet ⊆ latticeDual E8LatticeSet := by
  rintro y ⟨ny, rfl⟩ x ⟨nx, rfl⟩
  simpa [e8_basis, mul_comm, smul_eq_mul, Pi.smul_apply] using e8_inner_integer nx ny

/-- Wenn `G` ganzzahlig und unimodular ist und `G c` ganzzahlige Koordinaten hat, dann auch `c`. -/
lemma unimodular_coordinates_integral (c : Fin 8 → ℝ)
    (h : ∀ i, ∃ z : ℤ, ∑ j, (e8_gram_z i j : ℝ) * c j = (z : ℝ)) :
    ∃ n : Fin 8 → ℤ, ∀ i, c i = (n i : ℝ) := by
  classical
  choose z hz using h
  let G : Matrix (Fin 8) (Fin 8) ℤ := e8_gram_z
  let A : Matrix (Fin 8) (Fin 8) ℤ := G.adjugate
  let ι : ℤ →+* ℝ := Int.castRingHom ℝ
  have hGc : G.map ι *ᵥ c = fun i => (z i : ℝ) := by
    ext i
    rw [Matrix.mulVec_apply_eq_sum]
    simpa [Matrix.map_apply] using hz i
  have hAG : A.map ι * G.map ι = 1 := by
    rw [← Matrix.map_mul]
    have hdet : A * G = 1 := by
      rw [show A * G = G.adjugate * G by simp [A]]
      rw [Matrix.adjugate_mul, e8_gram_z_det_one]
      simp
    rw [hdet]
    exact Matrix.map_one ι (map_zero ι) (map_one ι)
  refine ⟨fun i => ∑ j, A i j * z j, ?_⟩
  intro i
  have hc : c = A.map ι *ᵥ G.map ι *ᵥ c := by
    rw [Matrix.mulVec_mulVec, hAG, Matrix.one_mulVec]
  have hi : ∑ j, ι (A i j) * (z j : ℝ) = c i := by
    have hi0 := congr_fun hc.symm i
    simp only [hGc, Matrix.mulVec_apply_eq_sum, Matrix.map_apply] at hi0
    exact hi0
  calc
    c i = ∑ j, ι (A i j * z j) := by
      rw [← hi]
      congr 1
      ext j
      exact (map_mul ι (A i j) (z j)).symm
    _ = ι (∑ j, A i j * z j) := (map_sum ι (fun j => A i j * z j) Finset.univ).symm

lemma e8_dual_subset_e8 : latticeDual E8LatticeSet ⊆ E8LatticeSet := by
  sorry

theorem e8_self_duality : latticeDual E8LatticeSet = E8LatticeSet := by
  exact Set.Subset.antisymm e8_dual_subset_e8 e8_subset_dual

end E8SelfDuality
