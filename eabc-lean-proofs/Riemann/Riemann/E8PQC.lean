import Mathlib

set_option linter.unusedSimpArgs false
set_option linter.unusedTactic false
set_option linter.unreachableTactic false

open Matrix MeasureTheory Module Set
open scoped BigOperators Real

noncomputable section

namespace E8PQC

def e8_basis : Fin 8 → (Fin 8 → ℝ) :=
  ![![ (1/2 : ℝ), (-1/2 : ℝ), (-1/2 : ℝ), (-1/2 : ℝ), (-1/2 : ℝ), (-1/2 : ℝ), (-1/2 : ℝ), (1/2 : ℝ)],
    ![ (1 : ℝ), 1, 0, 0, 0, 0, 0, 0],
    ![(-1 : ℝ), 1, 0, 0, 0, 0, 0, 0],
    ![ (0 : ℝ), -1, 1, 0, 0, 0, 0, 0],
    ![ (0 : ℝ), 0, -1, 1, 0, 0, 0, 0],
    ![ (0 : ℝ), 0, 0, -1, 1, 0, 0, 0],
    ![ (0 : ℝ), 0, 0, 0, -1, 1, 0, 0],
    ![ (0 : ℝ), 0, 0, 0, 0, -1, 1, 0]]

def e8_basis_matrix : Matrix (Fin 8) (Fin 8) ℝ :=
  fun i j => e8_basis i j

def E8LatticeSet : Set (Fin 8 → ℝ) :=
  { x | ∃ n : Fin 8 → ℤ, x = ∑ i, (n i : ℝ) • e8_basis i }

def latticeDual (Λ : Set (Fin 8 → ℝ)) : Set (Fin 8 → ℝ) :=
  { y | ∀ x ∈ Λ, ∃ m : ℤ, ∑ i, x i * y i = (m : ℝ) }

def e8_gram_z : Matrix (Fin 8) (Fin 8) ℤ :=
  !![2, 0, -1, 0, 0, 0, 0, 0;
     0, 2, 0, -1, 0, 0, 0, 0;
     -1, 0, 2, -1, 0, 0, 0, 0;
     0, -1, -1, 2, -1, 0, 0, 0;
     0, 0, 0, -1, 2, -1, 0, 0;
     0, 0, 0, 0, -1, 2, -1, 0;
     0, 0, 0, 0, 0, -1, 2, -1;
     0, 0, 0, 0, 0, 0, -1, 2]

def e8_basis_twice : Matrix (Fin 8) (Fin 8) ℤ :=
  !![1, -1, -1, -1, -1, -1, -1, 1;
     2, 2, 0, 0, 0, 0, 0, 0;
     -2, 2, 0, 0, 0, 0, 0, 0;
     0, -2, 2, 0, 0, 0, 0, 0;
     0, 0, -2, 2, 0, 0, 0, 0;
     0, 0, 0, -2, 2, 0, 0, 0;
     0, 0, 0, 0, -2, 2, 0, 0;
     0, 0, 0, 0, 0, -2, 2, 0]

def e8_gram_r : Matrix (Fin 8) (Fin 8) ℝ :=
  fun i j => ∑ k, e8_basis i k * e8_basis j k

def addSmulRow (A : Matrix (Fin 8) (Fin 8) ℚ) (i j : Fin 8) (c : ℚ) :
    Matrix (Fin 8) (Fin 8) ℚ :=
  A.updateRow i (A i + c • A j)

lemma det_addSmulRow (A B : Matrix (Fin 8) (Fin 8) ℚ) (i j : Fin 8) (c : ℚ)
    (hij : i ≠ j) (h : addSmulRow A i j c = B) : B.det = A.det := by
  rw [← h, addSmulRow]
  exact Matrix.det_updateRow_add_smul_self A hij c

def M0 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    -1, 0, 2, -1, 0, 0, 0, 0;
    0, -1, -1, 2, -1, 0, 0, 0;
    0, 0, 0, -1, 2, -1, 0, 0;
    0, 0, 0, 0, -1, 2, -1, 0;
    0, 0, 0, 0, 0, -1, 2, -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M1 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, -1, -1, 2, -1, 0, 0, 0;
    0, 0, 0, -1, 2, -1, 0, 0;
    0, 0, 0, 0, -1, 2, -1, 0;
    0, 0, 0, 0, 0, -1, 2, -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M2 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, 0, -1, (3/2), -1, 0, 0, 0;
    0, 0, 0, -1, 2, -1, 0, 0;
    0, 0, 0, 0, -1, 2, -1, 0;
    0, 0, 0, 0, 0, -1, 2, -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M3 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, 0, 0, (5/6), -1, 0, 0, 0;
    0, 0, 0, -1, 2, -1, 0, 0;
    0, 0, 0, 0, -1, 2, -1, 0;
    0, 0, 0, 0, 0, -1, 2, -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M4 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, 0, 0, (5/6), -1, 0, 0, 0;
    0, 0, 0, 0, (4/5), -1, 0, 0;
    0, 0, 0, 0, -1, 2, -1, 0;
    0, 0, 0, 0, 0, -1, 2, -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M5 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, 0, 0, (5/6), -1, 0, 0, 0;
    0, 0, 0, 0, (4/5), -1, 0, 0;
    0, 0, 0, 0, 0, (3/4), -1, 0;
    0, 0, 0, 0, 0, -1, 2, -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M6 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, 0, 0, (5/6), -1, 0, 0, 0;
    0, 0, 0, 0, (4/5), -1, 0, 0;
    0, 0, 0, 0, 0, (3/4), -1, 0;
    0, 0, 0, 0, 0, 0, (2/3), -1;
    0, 0, 0, 0, 0, 0, -1, 2]

def M7 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    2, 0, -1, 0, 0, 0, 0, 0;
    0, 2, 0, -1, 0, 0, 0, 0;
    0, 0, (3/2), -1, 0, 0, 0, 0;
    0, 0, 0, (5/6), -1, 0, 0, 0;
    0, 0, 0, 0, (4/5), -1, 0, 0;
    0, 0, 0, 0, 0, (3/4), -1, 0;
    0, 0, 0, 0, 0, 0, (2/3), -1;
    0, 0, 0, 0, 0, 0, 0, (1/2)]


lemma Mstep0 : addSmulRow M0 2 0 ((1/2) : ℚ) = M1 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M0, M1, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mstep1 : addSmulRow M1 3 1 ((1/2) : ℚ) = M2 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M1, M2, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mstep2 : addSmulRow M2 3 2 ((2/3) : ℚ) = M3 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M2, M3, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mstep3 : addSmulRow M3 4 3 ((6/5) : ℚ) = M4 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M3, M4, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mstep4 : addSmulRow M4 5 4 ((5/4) : ℚ) = M5 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M4, M5, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mstep5 : addSmulRow M5 6 5 ((4/3) : ℚ) = M6 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M5, M6, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mstep6 : addSmulRow M6 7 6 ((3/2) : ℚ) = M7 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, M6, M7, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Mupper : M7.IsUpperTriangular := by
  intro i j hij
  fin_cases i <;> fin_cases j <;> (try norm_num at hij) <;> simp [M7]

lemma Mfinal_det : M7.det = 1 := by
  rw [Matrix.det_of_isUpperTriangular Mupper]
  simp only [Fin.prod_univ_succ, Fin.prod_univ_zero, M7]
  norm_num

lemma Mstart_det : M0.det = 1 := by
  have h0 := det_addSmulRow M0 M1 2 0 ((1/2) : ℚ) (by decide) Mstep0
  have h1 := det_addSmulRow M1 M2 3 1 ((1/2) : ℚ) (by decide) Mstep1
  have h2 := det_addSmulRow M2 M3 3 2 ((2/3) : ℚ) (by decide) Mstep2
  have h3 := det_addSmulRow M3 M4 4 3 ((6/5) : ℚ) (by decide) Mstep3
  have h4 := det_addSmulRow M4 M5 5 4 ((5/4) : ℚ) (by decide) Mstep4
  have h5 := det_addSmulRow M5 M6 6 5 ((4/3) : ℚ) (by decide) Mstep5
  have h6 := det_addSmulRow M6 M7 7 6 ((3/2) : ℚ) (by decide) Mstep6
  rw [← h0, ← h1, ← h2, ← h3, ← h4, ← h5, ← h6, Mfinal_det]

lemma Mmap_eq : (Int.castRingHom ℚ).mapMatrix e8_gram_z = M0 := by
  simp only [RingHom.mapMatrix]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [e8_gram_z, M0, Matrix.map_apply]

theorem e8_gram_z_det_one : e8_gram_z.det = 1 := by
  generalize hA : e8_gram_z = A
  have hmap := (Int.castRingHom ℚ).map_det A
  have hM : (Int.castRingHom ℚ).mapMatrix A = M0 := hA ▸ Mmap_eq
  exact Int.cast_injective (hmap.trans ((congrArg Matrix.det hM).trans Mstart_det))

def T0 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    2, 2, 0, 0, 0, 0, 0, 0;
    -2, 2, 0, 0, 0, 0, 0, 0;
    0, -2, 2, 0, 0, 0, 0, 0;
    0, 0, -2, 2, 0, 0, 0, 0;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T1 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    -2, 2, 0, 0, 0, 0, 0, 0;
    0, -2, 2, 0, 0, 0, 0, 0;
    0, 0, -2, 2, 0, 0, 0, 0;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T2 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, -2, 2, 0, 0, 0, 0, 0;
    0, 0, -2, 2, 0, 0, 0, 0;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T3 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 3, 1, 1, 1, 1, -1;
    0, 0, -2, 2, 0, 0, 0, 0;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T4 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, -2, 2, 0, 0, 0, 0;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T5 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 4, 2, 2, 2, -2;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T6 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, -2, 2, 0, 0, 0;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T7 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, 0, 4, 2, 2, -2;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T8 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, 0, 0, -2, -2, 2;
    0, 0, 0, 0, -2, 2, 0, 0;
    0, 0, 0, 0, 0, -2, 2, 0]

def T9 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, 0, 0, -2, -2, 2;
    0, 0, 0, 0, 0, 4, 2, -2;
    0, 0, 0, 0, 0, -2, 2, 0]

def T10 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, 0, 0, -2, -2, 2;
    0, 0, 0, 0, 0, 0, -2, 2;
    0, 0, 0, 0, 0, -2, 2, 0]

def T11 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, 0, 0, -2, -2, 2;
    0, 0, 0, 0, 0, 0, -2, 2;
    0, 0, 0, 0, 0, 0, 4, -2]

def T12 : Matrix (Fin 8) (Fin 8) ℚ :=
  !![
    1, -1, -1, -1, -1, -1, -1, 1;
    0, 4, 2, 2, 2, 2, 2, -2;
    0, 0, -2, -2, -2, -2, -2, 2;
    0, 0, 0, -2, -2, -2, -2, 2;
    0, 0, 0, 0, -2, -2, -2, 2;
    0, 0, 0, 0, 0, -2, -2, 2;
    0, 0, 0, 0, 0, 0, -2, 2;
    0, 0, 0, 0, 0, 0, 0, 2]


lemma Tstep0 : addSmulRow T0 1 0 (-2 : ℚ) = T1 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T0, T1, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep1 : addSmulRow T1 2 0 (2 : ℚ) = T2 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T1, T2, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep2 : addSmulRow T2 3 1 ((1/2) : ℚ) = T3 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T2, T3, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep3 : addSmulRow T3 3 2 ((3/2) : ℚ) = T4 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T3, T4, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep4 : addSmulRow T4 4 2 (-1 : ℚ) = T5 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T4, T5, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep5 : addSmulRow T5 4 3 (2 : ℚ) = T6 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T5, T6, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep6 : addSmulRow T6 5 3 (-1 : ℚ) = T7 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T6, T7, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep7 : addSmulRow T7 5 4 (2 : ℚ) = T8 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T7, T8, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep8 : addSmulRow T8 6 4 (-1 : ℚ) = T9 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T8, T9, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep9 : addSmulRow T9 6 5 (2 : ℚ) = T10 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T9, T10, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep10 : addSmulRow T10 7 5 (-1 : ℚ) = T11 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T10, T11, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tstep11 : addSmulRow T11 7 6 (2 : ℚ) = T12 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    (simp [addSmulRow, Matrix.updateRow_apply, T11, T12, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul] <;> norm_num)


lemma Tupper : T12.IsUpperTriangular := by
  intro i j hij
  fin_cases i <;> fin_cases j <;> (try norm_num at hij) <;> simp [T12]

lemma Tfinal_det : T12.det = (-256) := by
  rw [Matrix.det_of_isUpperTriangular Tupper]
  simp only [Fin.prod_univ_succ, Fin.prod_univ_zero, T12]
  norm_num

lemma Tstart_det : T0.det = (-256) := by
  have h0 := det_addSmulRow T0 T1 1 0 (-2 : ℚ) (by decide) Tstep0
  have h1 := det_addSmulRow T1 T2 2 0 (2 : ℚ) (by decide) Tstep1
  have h2 := det_addSmulRow T2 T3 3 1 ((1/2) : ℚ) (by decide) Tstep2
  have h3 := det_addSmulRow T3 T4 3 2 ((3/2) : ℚ) (by decide) Tstep3
  have h4 := det_addSmulRow T4 T5 4 2 (-1 : ℚ) (by decide) Tstep4
  have h5 := det_addSmulRow T5 T6 4 3 (2 : ℚ) (by decide) Tstep5
  have h6 := det_addSmulRow T6 T7 5 3 (-1 : ℚ) (by decide) Tstep6
  have h7 := det_addSmulRow T7 T8 5 4 (2 : ℚ) (by decide) Tstep7
  have h8 := det_addSmulRow T8 T9 6 4 (-1 : ℚ) (by decide) Tstep8
  have h9 := det_addSmulRow T9 T10 6 5 (2 : ℚ) (by decide) Tstep9
  have h10 := det_addSmulRow T10 T11 7 5 (-1 : ℚ) (by decide) Tstep10
  have h11 := det_addSmulRow T11 T12 7 6 (2 : ℚ) (by decide) Tstep11
  rw [← h0, ← h1, ← h2, ← h3, ← h4, ← h5, ← h6, ← h7, ← h8, ← h9, ← h10, ← h11, Tfinal_det]

lemma Tmap_eq : (Int.castRingHom ℚ).mapMatrix e8_basis_twice = T0 := by
  simp only [RingHom.mapMatrix]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [e8_basis_twice, T0, Matrix.map_apply]

theorem e8_basis_twice_det : e8_basis_twice.det = (-256) := by
  generalize hA : e8_basis_twice = A
  have hmap := (Int.castRingHom ℚ).map_det A
  have hM : (Int.castRingHom ℚ).mapMatrix A = T0 := hA ▸ Tmap_eq
  exact Int.cast_injective (hmap.trans ((congrArg Matrix.det hM).trans Tstart_det))

theorem e8_gram_eq :
    e8_gram_r = e8_gram_z.map (Int.castRingHom ℝ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    (simp [e8_gram_r, e8_gram_z, e8_basis, Matrix.map_apply, Fin.sum_univ_succ,
      Fin.sum_univ_zero] <;> norm_num)

theorem e8_basis_gram_det_one : e8_gram_r.det = 1 := by
  rw [e8_gram_eq]
  generalize hA : e8_gram_z = A
  have hmap := (Int.castRingHom ℝ).map_det A
  have hdet : (A.det : ℝ) = 1 := by
    have : A.det = 1 := hA ▸ e8_gram_z_det_one
    simp [this]
  have hdef : (Int.castRingHom ℝ).mapMatrix A = A.map (Int.castRingHom ℝ) := rfl
  rw [← hdef]
  exact hmap.symm.trans hdet

lemma twice_map_eq :
    (Int.castRingHom ℝ).mapMatrix e8_basis_twice = (2 : ℝ) • e8_basis_matrix := by
  simp only [RingHom.mapMatrix]
  ext i j
  fin_cases i <;> fin_cases j <;>
    (simp [e8_basis_twice, e8_basis_matrix, e8_basis, Matrix.map_apply, Pi.smul_apply,
      smul_eq_mul] <;> norm_num)

theorem e8_basis_matrix_det : e8_basis_matrix.det = -1 := by
  have htwice : ((2 : ℝ) • e8_basis_matrix).det = -256 := by
    generalize hA : e8_basis_twice = A
    have hmap := (Int.castRingHom ℝ).map_det A
    have hM : (Int.castRingHom ℝ).mapMatrix A = (2 : ℝ) • e8_basis_matrix :=
      hA ▸ twice_map_eq
    have hdet : (A.det : ℝ) = -256 := by
      have : A.det = -256 := hA ▸ e8_basis_twice_det
      simp [this]
    exact (congrArg Matrix.det hM).symm.trans (hmap.symm.trans hdet)
  have hsmul := Matrix.det_smul e8_basis_matrix (2 : ℝ)
  rw [hsmul, Fintype.card_fin] at htwice
  have hpow : (2 : ℝ) ^ 8 = 256 := by norm_num
  rw [hpow] at htwice
  apply mul_left_cancel₀ (by norm_num : (256 : ℝ) ≠ 0)
  calc
    256 * e8_basis_matrix.det = -256 := htwice
    _ = 256 * (-1) := by norm_num

lemma e8_basis_linearIndependent : LinearIndependent ℝ e8_basis := by
  have hdet : e8_basis_matrix.det ≠ 0 := by
    rw [e8_basis_matrix_det]
    norm_num
  have h := Matrix.linearIndependent_rows_of_det_ne_zero (A := e8_basis_matrix) hdet
  have hfun : (fun i => e8_basis_matrix i) = e8_basis := by
    funext i
    rfl
  simpa [hfun] using h

noncomputable def e8Basis : Basis (Fin 8) ℝ (Fin 8 → ℝ) :=
  basisOfLinearIndependentOfCardEqFinrank (b := e8_basis) e8_basis_linearIndependent
    (by simp [Fintype.card_fin])

lemma e8_inner_integer (nx ny : Fin 8 → ℤ) :
    ∃ m : ℤ,
      ∑ k, (∑ i, (nx i : ℝ) * e8_basis i k) * (∑ j, (ny j : ℝ) * e8_basis j k) =
        (m : ℝ) := by
  refine ⟨∑ i, ∑ j, nx i * ny j * e8_gram_z i j, ?_⟩
  push_cast
  have hmul : ∀ i j k,
      (nx i : ℝ) * e8_basis i k * ((ny j : ℝ) * e8_basis j k) =
        ((nx i : ℝ) * (ny j : ℝ)) * (e8_basis i k * e8_basis j k) := by
    intro i j k
    ring
  simp_rw [Finset.sum_mul, Finset.mul_sum, hmul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl ?_
  intro i _
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl ?_
  intro j _
  rw [← Finset.mul_sum]
  have hg : ∑ k, e8_basis i k * e8_basis j k = (e8_gram_z i j : ℝ) := by
    simpa [e8_gram_r] using congrFun (congrFun e8_gram_eq i) j
  rw [hg]

lemma e8_subset_dual : E8LatticeSet ⊆ latticeDual E8LatticeSet := by
  rintro y ⟨ny, rfl⟩ x ⟨nx, rfl⟩
  exact e8_inner_integer nx ny

lemma unimodular_coordinates_integral
    (c : Fin 8 → ℝ)
    (h : ∀ i, ∃ z : ℤ, ∑ j, (e8_gram_z i j : ℝ) * c j = (z : ℝ)) :
    ∃ n : Fin 8 → ℤ, ∀ i, c i = (n i : ℝ) := by
  choose z hz using h
  let v : Fin 8 → ℤ := z
  generalize hA : e8_gram_z = A
  have h_adj : adjugate A * A = 1 := by
    have hmul := adjugate_mul A
    have hdet : A.det = 1 := hA ▸ e8_gram_z_det_one
    rw [hdet, one_smul] at hmul
    exact hmul
  have hzA : ∀ k, ∑ j, (A k j : ℝ) * c j = (v k : ℝ) := by
    intro k
    have hk := hz k
    rwa [hA] at hk
  refine ⟨adjugate A *ᵥ v, ?_⟩
  intro i
  have h_map_mul :
      (adjugate A).map (Int.castRingHom ℝ) * A.map (Int.castRingHom ℝ) = 1 := by
    rw [← Matrix.map_mul, h_adj]
    exact Matrix.map_one (Int.castRingHom ℝ) (map_zero _) (map_one _)
  have hGc : A.map (Int.castRingHom ℝ) *ᵥ c = fun k => (v k : ℝ) := by
    ext k
    rw [Matrix.mulVec_apply_eq_sum]
    simpa [Matrix.map_apply] using hzA k
  have h_cast :
      (((adjugate A) *ᵥ v) i : ℝ) =
        (((adjugate A).map (Int.castRingHom ℝ)) *ᵥ fun k => (v k : ℝ)) i :=
    RingHom.map_mulVec (Int.castRingHom ℝ) (adjugate A) v i
  rw [h_cast, ← hGc, Matrix.mulVec_mulVec, h_map_mul, Matrix.one_mulVec]

lemma e8_basis_mem (i : Fin 8) : e8_basis i ∈ E8LatticeSet := by
  refine ⟨fun j => if j = i then (1 : ℤ) else 0, ?_⟩
  ext k
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hj
    simp [hj]
  · intro hi
    exact (hi (Finset.mem_univ i)).elim

lemma e8_basis_dot_sum (i : Fin 8) (c : Fin 8 → ℝ) :
    ∑ k, e8_basis i k * (∑ j, c j * e8_basis j k) =
      ∑ j, (e8_gram_z i j : ℝ) * c j := by
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl ?_
  intro j _
  rw [show ∑ k, e8_basis i k * (c j * e8_basis j k) =
      c j * ∑ k, e8_basis i k * e8_basis j k by
    have hterm : ∀ k, e8_basis i k * (c j * e8_basis j k) =
        c j * (e8_basis i k * e8_basis j k) := by
      intro k
      ring
    simp_rw [hterm]
    rw [← Finset.mul_sum]]
  have hg : ∑ k, e8_basis i k * e8_basis j k = (e8_gram_z i j : ℝ) := by
    simpa [e8_gram_r] using congrFun (congrFun e8_gram_eq i) j
  rw [hg, mul_comm]

lemma e8Basis_coord (c : Fin 8 → ℝ) (i : Fin 8) :
    e8Basis.repr (∑ j, c j • e8Basis j) i = c i := by
  classical
  let v : Fin 8 →₀ ℝ := Finsupp.equivFunOnFinite.symm c
  have hv : Finsupp.linearCombination ℝ (⇑e8Basis) v = ∑ j, c j • e8Basis j := by
    rw [Finsupp.linearCombination_apply, Finsupp.sum_fintype]
    · simp [v, Finsupp.coe_equivFunOnFinite_symm]
    · intro j
      exact zero_smul ℝ _
  have hr : e8Basis.repr (∑ j, c j • e8Basis j) = v := by
    rw [← hv]
    exact e8Basis.repr_linearCombination v
  rw [hr]
  simp [v, Finsupp.coe_equivFunOnFinite_symm]

lemma e8_dual_subset_e8 : latticeDual E8LatticeSet ⊆ E8LatticeSet := by
  intro y hy
  let c : Fin 8 → ℝ := fun j => e8Basis.repr y j
  have hcoe : ∀ j, e8Basis j = e8_basis j := by
    intro j
    simp [e8Basis]
  have hy_expand : y = ∑ j, c j • e8_basis j := by
    have : y = ∑ j, c j • e8Basis j := (e8Basis.sum_repr y).symm
    simpa [hcoe] using this
  have h_int : ∀ i, ∃ z : ℤ, ∑ j, (e8_gram_z i j : ℝ) * c j = (z : ℝ) := by
    intro i
    obtain ⟨z, hz⟩ := hy (e8_basis i) (e8_basis_mem i)
    refine ⟨z, ?_⟩
    rw [← hz]
    have hy_comp : ∀ k, y k = ∑ j, c j * e8_basis j k := by
      intro k
      rw [hy_expand]
      simp [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    simp_rw [hy_comp]
    exact (e8_basis_dot_sum i c).symm
  obtain ⟨n, hn⟩ := unimodular_coordinates_integral c h_int
  refine ⟨n, ?_⟩
  rw [hy_expand]
  congr 1
  ext i
  rw [hn i]

theorem e8_self_duality : latticeDual E8LatticeSet = E8LatticeSet := by
  ext y
  constructor
  · exact fun hy => e8_dual_subset_e8 hy
  · exact fun hy => e8_subset_dual hy

def E8FundamentalDomain : Set (Fin 8 → ℝ) :=
  { x | ∃ c : Fin 8 → ℝ, (∀ i, 0 ≤ c i ∧ c i < 1) ∧ x = ∑ i, c i • e8_basis i }

def fourierKernel (k : Fin 8 → ℝ) (x : Fin 8 → ℝ) : ℂ :=
  Complex.exp (-2 * Real.pi * Complex.I * ((∑ i, k i * x i) : ℂ))

def fourierChar (k : Fin 8 → ℝ) (x : Fin 8 → ℝ) : ℂ :=
  Complex.exp (2 * Real.pi * Complex.I * ((∑ i, k i * x i) : ℂ))

lemma e8FundamentalDomain_eq_zspan :
    E8FundamentalDomain = ZSpan.fundamentalDomain e8Basis := by
  ext x
  constructor
  · rintro ⟨c, hc, rfl⟩
    rw [ZSpan.mem_fundamentalDomain]
    intro i
    have hsum : (∑ j, c j • e8_basis j) = ∑ j, c j • e8Basis j := by
      simp [e8Basis]
    rw [hsum, e8Basis_coord]
    exact hc i
  · intro hx
    rw [ZSpan.mem_fundamentalDomain] at hx
    refine ⟨fun i => e8Basis.repr x i, hx, ?_⟩
    have hsum :
        (∑ j, e8Basis.repr x j • e8_basis j) = ∑ j, e8Basis.repr x j • e8Basis j := by
      simp [e8Basis]
    rw [hsum]
    exact (e8Basis.sum_repr x).symm

theorem e8_fundamental_domain_volume : volume E8FundamentalDomain = 1 := by
  rw [e8FundamentalDomain_eq_zspan, ZSpan.volume_fundamentalDomain]
  have hM : Matrix.of (⇑e8Basis) = e8_basis_matrix := by
    ext i j
    simp [e8Basis, e8_basis_matrix, Matrix.of_apply]
  rw [hM, e8_basis_matrix_det]
  simp [ENNReal.ofReal_one]

lemma fourier_integral_nonzero (m : ℤ) (hm : m ≠ 0) :
    ∫ t in Icc (0 : ℝ) 1,
      Complex.exp (2 * Real.pi * Complex.I * (m : ℂ) * (t : ℂ)) = 0 := by
  let c : ℂ := 2 * Real.pi * Complex.I * (m : ℂ)
  have hc : c ≠ 0 := by
    intro h
    apply hm
    have h2 : (2 : ℂ) ≠ 0 := by norm_num
    have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
    have hI : Complex.I ≠ 0 := Complex.I_ne_zero
    rcases mul_eq_zero.mp h with h | hm0
    · rcases mul_eq_zero.mp h with h | hI0
      · rcases mul_eq_zero.mp h with h2' | hpi'
        · exact absurd h2' h2
        · exact absurd hpi' hpi
      · exact absurd hI0 hI
    · exact_mod_cast hm0
  rw [MeasureTheory.integral_Icc_eq_integral_Ioc]
  rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  have hfun :
      (fun t : ℝ => Complex.exp (2 * Real.pi * Complex.I * (m : ℂ) * (t : ℂ))) =
        fun t : ℝ => Complex.exp (c * t) := by
    ext t
    simp [c, mul_assoc]
  rw [hfun, integral_exp_mul_complex hc]
  have hexp1 : Complex.exp (c * (1 : ℝ)) = 1 := by
    have hc1 : c * (1 : ℝ) = (m : ℂ) * (2 * Real.pi * Complex.I) := by
      simp only [c, mul_one]
      ac_rfl
    rw [hc1]
    exact Complex.exp_int_mul_two_pi_mul_I m
  have hexp0 : Complex.exp (c * (0 : ℝ)) = 1 := by simp
  rw [hexp1, hexp0]
  simp

lemma fourier_orthogonality_1d (m : ℤ) :
    ∫ t in Icc (0 : ℝ) 1,
      Complex.exp (2 * Real.pi * Complex.I * (m : ℂ) * (t : ℂ)) =
        if m = 0 then 1 else 0 := by
  by_cases hm : m = 0
  · subst hm
    rw [MeasureTheory.integral_Icc_eq_integral_Ioc]
    rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    simp only [Int.cast_zero, mul_zero, Complex.exp_zero, intervalIntegral.integral_const]
    simp
  · simp [hm, fourier_integral_nonzero m hm]

lemma e8_phase_coord (k c : Fin 8 → ℝ) :
    ∑ r, k r * (∑ j, c j * e8_basis j r) =
      ∑ j, c j * (∑ r, k r * e8_basis j r) := by
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hterm : ∀ r, k r * (c j * e8_basis j r) = c j * (k r * e8_basis j r) := by
    intro r
    ring
  simp_rw [hterm]

lemma fourierChar_coord_factorization_all (q : Fin 8 → ℝ) (hq : q ∈ E8LatticeSet) :
    ∃ m : Fin 8 → ℤ,
      (∀ j, ∑ i, e8_basis j i * q i = (m j : ℝ)) ∧
        (∀ c : Fin 8 → ℝ, fourierChar q (∑ j, c j • e8_basis j) =
          ∏ j, Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ))) := by
  have hdual : q ∈ latticeDual E8LatticeSet := e8_subset_dual hq
  choose m hm using fun j : Fin 8 => hdual (e8_basis j) (e8_basis_mem j)
  refine ⟨m, hm, fun c : Fin 8 → ℝ => ?_⟩
  unfold fourierChar
  have hx : ∀ r, (∑ j, c j • e8_basis j) r = ∑ j, c j * e8_basis j r := by
    intro r
    simp [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  have hphase :
      ∑ r, q r * (∑ j, c j • e8_basis j) r = ∑ j, (m j : ℝ) * c j := by
    simp_rw [hx, e8_phase_coord]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [show ∑ r, q r * e8_basis j r = ∑ i, e8_basis j i * q i by
      refine Finset.sum_congr rfl fun r _ => mul_comm _ _]
    rw [hm j, mul_comm]
  have hdist :
      (∑ i, (q i : ℂ) * ((∑ j, c j • e8_basis j) i : ℂ)) =
        ((∑ r, q r * (∑ j, c j • e8_basis j) r : ℝ) : ℂ) := by
    calc
      ∑ i, (q i : ℂ) * ((∑ j, c j • e8_basis j) i : ℂ)
          = ∑ i, ((q i * (∑ j, c j • e8_basis j) i : ℝ) : ℂ) := by
            refine Finset.sum_congr rfl fun i _ => (RCLike.ofReal_mul (K := ℂ) _ _).symm
        _ = ((∑ r, q r * (∑ j, c j • e8_basis j) r : ℝ) : ℂ) :=
            (RCLike.ofReal_sum (K := ℂ) _ _).symm
  have harg :
      2 * Real.pi * Complex.I *
          ∑ i, (q i : ℂ) * ((∑ j, c j • e8_basis j) i : ℂ) =
        ∑ j, 2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ) := by
    rw [hdist, hphase]
    push_cast
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    ac_rfl
  rw [harg]
  exact Complex.exp_sum _ _

lemma fourierChar_coord_factorization
    (q c : Fin 8 → ℝ)
    (hq : q ∈ E8LatticeSet) :
    ∃ m : Fin 8 → ℤ,
      fourierChar q (∑ j, c j • e8_basis j) =
        ∏ j, Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ)) := by
  obtain ⟨m, -, hm⟩ := fourierChar_coord_factorization_all q hq
  exact ⟨m, hm c⟩

open scoped ComplexConjugate in
lemma star_fourierChar (k x : Fin 8 → ℝ) :
    star (fourierChar k x) = fourierChar (-k) x := by
  unfold fourierChar
  rw [RCLike.star_def, ← Complex.exp_conj]
  congr 1
  set s : ℂ := ∑ i, (k i : ℂ) * (x i : ℂ)
  have hs_neg : ∑ i, ((-k) i : ℂ) * (x i : ℂ) = -s := by
    simp_rw [Pi.neg_apply, Complex.ofReal_neg, neg_mul, Finset.sum_neg_distrib]
    rfl
  have hs_conj : conj s = s := by
    rw [map_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [map_mul, Complex.conj_ofReal, Complex.conj_ofReal]
  rw [hs_neg, map_mul, map_mul, map_mul, hs_conj, Complex.conj_I, Complex.conj_ofReal,
    Complex.conj_ofNat]
  ring

lemma fourierChar_add (k l x : Fin 8 → ℝ) :
    fourierChar k x * fourierChar l x = fourierChar (k + l) x := by
  unfold fourierChar
  rw [← Complex.exp_add]
  congr 1
  have hs :
      ∑ i, ((k + l) i : ℂ) * (x i : ℂ) =
        (∑ i, (k i : ℂ) * (x i : ℂ)) + ∑ i, (l i : ℂ) * (x i : ℂ) := by
    simp_rw [Pi.add_apply, Complex.ofReal_add, add_mul, Finset.sum_add_distrib]
  rw [hs, mul_add]

lemma fourierChar_mul_star (k l x : Fin 8 → ℝ) :
    fourierChar k x * star (fourierChar l x) = fourierChar (k - l) x := by
  rw [star_fourierChar]
  simpa [sub_eq_add_neg] using fourierChar_add k (-l) x

lemma E8LatticeSet.sub_mem {k l : Fin 8 → ℝ}
    (hk : k ∈ E8LatticeSet) (hl : l ∈ E8LatticeSet) : l - k ∈ E8LatticeSet := by
  rcases hk with ⟨nk, hk⟩
  rcases hl with ⟨nl, hl⟩
  refine ⟨nl - nk, ?_⟩
  rw [hk, hl]
  ext r
  simp [Finset.sum_apply, Pi.sub_apply, Pi.smul_apply, sub_smul, Int.cast_sub,
    Finset.sum_sub_distrib]

lemma e8_fourier_product_integral (q : Fin 8 → ℝ) (hq : q ∈ E8LatticeSet) :
    ∃ m : Fin 8 → ℤ,
      (∀ j, ∑ i, e8_basis j i * q i = (m j : ℝ)) ∧
        ∫ x in E8FundamentalDomain, fourierChar q x =
          ∏ j, ∫ t in Icc (0 : ℝ) 1,
            Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (t : ℂ)) := by
  obtain ⟨m, hm_inner, hm_prod⟩ := fourierChar_coord_factorization_all q hq
  let fL : (Fin 8 → ℝ) →ₗ[ℝ] Fin 8 → ℝ := Matrix.toLin' e8_basis_matrix.transpose
  have hdet : LinearMap.det fL = -1 := by
    rw [LinearMap.det_toLin', Matrix.det_transpose, e8_basis_matrix_det]
  have hf : LinearMap.det fL ≠ 0 := by
    rw [hdet]
    norm_num
  have hcoord (c : Fin 8 → ℝ) : fL c = ∑ j, c j • e8_basis j := by
    ext r
    rw [Matrix.toLin'_apply, Matrix.mulVec_apply_eq_sum]
    simp only [Matrix.transpose_apply, e8_basis_matrix, Finset.sum_apply, Pi.smul_apply,
      smul_eq_mul]
    refine Finset.sum_congr rfl fun j _ => mul_comm _ _
  have hpres : MeasurePreserving (⇑fL) := by
    refine ⟨(LinearMap.continuous_of_finiteDimensional fL).measurable, ?_⟩
    rw [Real.map_linearMap_volume_pi_eq_smul_volume_pi hf]
    have habs : |(LinearMap.det fL)⁻¹| = 1 := by
      rw [hdet]
      norm_num
    simp [habs, ENNReal.ofReal_one, one_smul]
  let fc := LinearMap.toContinuousLinearMap fL
  have hdetc : fc.det ≠ 0 := by
    change LinearMap.det (fc : (Fin 8 → ℝ) →ₗ[ℝ] Fin 8 → ℝ) ≠ 0
    rw [LinearMap.coe_toContinuousLinearMap]
    exact hf
  have hembed : MeasurableEmbedding (⇑fL) := by
    have hfun :
        (⇑fL) = (fc.toContinuousLinearEquivOfDetNeZero hdetc :
          (Fin 8 → ℝ) → Fin 8 → ℝ) := by
      funext c
      exact (fc.toContinuousLinearEquivOfDetNeZero_apply hdetc c).symm
    rw [hfun]
    exact (fc.toContinuousLinearEquivOfDetNeZero hdetc).toHomeomorph.toMeasurableEquiv.measurableEmbedding
  let cube : Set (Fin 8 → ℝ) := Set.univ.pi fun _ => Ico (0 : ℝ) 1
  have hdom : E8FundamentalDomain = ⇑fL '' cube := by
    ext x
    constructor
    · rintro ⟨c, hc, rfl⟩
      refine ⟨c, ?_, hcoord c⟩
      intro i _
      simpa [mem_Ico] using hc i
    · rintro ⟨c, hc, rfl⟩
      refine ⟨c, ?_, hcoord c⟩
      intro i
      simpa [mem_Ico] using hc i (mem_univ i)
  have hpi :
      ∫ c in cube, ∏ j, Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ)) =
        ∏ j, ∫ t in Ico (0 : ℝ) 1,
          Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (t : ℂ)) := by
    rw [volume_pi, Measure.restrict_pi_pi (μ := fun _ : Fin 8 => (volume : Measure ℝ))
      (s := fun _ => Ico (0 : ℝ) 1)]
    exact integral_fintype_prod_eq_prod
      (μ := fun _ : Fin 8 => (volume : Measure ℝ).restrict (Ico (0 : ℝ) 1))
      (fun j (t : ℝ) => Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (t : ℂ)))
  refine ⟨m, hm_inner, ?_⟩
  rw [hdom, hpres.setIntegral_image_emb hembed (fourierChar q) cube]
  have hcube :
      ∫ c in cube, fourierChar q (fL c) =
        ∫ c in cube, ∏ j, Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ)) := by
    refine setIntegral_congr_fun (MeasurableSet.univ_pi fun _ => measurableSet_Ico) ?_
    intro c _
    change fourierChar q (fL c) = _
    rw [hcoord c]
    exact hm_prod c
  rw [hcube, hpi]
  refine Finset.prod_congr rfl fun j _ => (MeasureTheory.integral_Icc_eq_integral_Ico).symm

lemma e8_fourierChar_integral (q : Fin 8 → ℝ) (hq : q ∈ E8LatticeSet) :
    ∫ x in E8FundamentalDomain, fourierChar q x =
      if q = 0 then (1 : ℂ) else 0 := by
  obtain ⟨m, hm_inner, hint⟩ := e8_fourier_product_integral q hq
  rw [hint]
  simp_rw [fourier_orthogonality_1d]
  have hq0 : (∀ j, m j = 0) ↔ q = 0 := by
    constructor
    · intro h
      apply Matrix.eq_zero_of_mulVec_eq_zero (M := e8_basis_matrix)
        (by rw [e8_basis_matrix_det]; norm_num)
      ext j
      rw [Matrix.mulVec_apply_eq_sum]
      simp only [e8_basis_matrix, Pi.zero_apply]
      rw [hm_inner j, h j, Int.cast_zero]
    · intro hq0 j
      apply Int.cast_injective (α := ℝ)
      rw [← hm_inner j, hq0]
      simp
  have hprod :
      (∏ j, if m j = 0 then (1 : ℂ) else 0) = if ∀ j, m j = 0 then 1 else 0 := by
    by_cases h : ∀ j, m j = 0
    · simp [h]
    · rw [ite_eq_right h]
      obtain ⟨j, hj⟩ := Classical.not_forall.mp h
      exact Finset.prod_eq_zero (Finset.mem_univ j) (ite_eq_right hj)
  rw [hprod]
  by_cases h : q = 0
  · simp [h, hq0.mpr h]
  · simp [h, hq0.not.mpr h]

theorem e8_fourier_orthogonality
    (k l : Fin 8 → ℝ)
    (hk : k ∈ E8LatticeSet)
    (hl : l ∈ E8LatticeSet) :
    ∫ x in E8FundamentalDomain, star (fourierChar k x) * fourierChar l x =
      if k = l then 1 else 0 := by
  have hq : l - k ∈ E8LatticeSet := E8LatticeSet.sub_mem hk hl
  have hint : ∀ x, star (fourierChar k x) * fourierChar l x = fourierChar (l - k) x := by
    intro x
    rw [mul_comm]
    exact fourierChar_mul_star l k x
  simp_rw [hint]
  simpa [sub_eq_zero, eq_comm] using e8_fourierChar_integral (l - k) hq

lemma e8Domain_finiteMeasure :
    IsFiniteMeasure (volume.restrict E8FundamentalDomain) :=
  isFiniteMeasure_restrict.mpr <| by
    simp [e8_fundamental_domain_volume]

lemma fourierChar_norm (k x : Fin 8 → ℝ) : ‖fourierChar k x‖ = 1 := by
  unfold fourierChar
  rw [Complex.norm_exp, Real.exp_eq_one_iff]
  set s : ℂ := ∑ i, (k i : ℂ) * (x i : ℂ)
  have hs : s.im = 0 := by
    unfold s
    rw [Complex.im_sum]
    refine Finset.sum_eq_zero fun i _ => ?_
    rw [Complex.mul_im]
    simp [Complex.ofReal_re, Complex.ofReal_im]
  have hπ : ((2 : ℂ) * Real.pi).im = 0 := by
    rw [Complex.mul_im]
    simp [Complex.ofReal_re, Complex.ofReal_im]
  rw [mul_assoc, Complex.mul_re, Complex.I_mul_re, Complex.I_mul_im, hs, hπ]
  simp

lemma fourierKernel_eq_star_fourierChar (k x : Fin 8 → ℝ) :
    fourierKernel k x = star (fourierChar k x) := by
  rw [star_fourierChar]
  unfold fourierKernel fourierChar
  congr 1
  have hsum :
      ∑ i, ((-k) i : ℂ) * (x i : ℂ) = -((∑ i, k i * x i : ℝ) : ℂ) := by
    simp only [Pi.neg_apply]
    have hterm : ∀ i, ((-(k i) : ℝ) : ℂ) * (x i : ℂ) = -((k i * x i : ℝ) : ℂ) := by
      intro i
      rw [Complex.ofReal_neg, neg_mul, Complex.ofReal_mul]
    rw [Finset.sum_congr rfl fun i _ => hterm i, Finset.sum_neg_distrib]
    congr 1
    exact (RCLike.ofReal_sum (K := ℂ) Finset.univ (fun i => k i * x i)).symm
  have hdist : ∑ i, (k i : ℂ) * (x i : ℂ) = ((∑ i, k i * x i : ℝ) : ℂ) := by
    calc
      ∑ i, (k i : ℂ) * (x i : ℂ) = ∑ i, ((k i * x i : ℝ) : ℂ) := by
        refine Finset.sum_congr rfl fun i _ => (Complex.ofReal_mul (k i) (x i)).symm
      _ = ((∑ i, k i * x i : ℝ) : ℂ) :=
        (RCLike.ofReal_sum (K := ℂ) Finset.univ (fun i => k i * x i)).symm
  rw [hsum, hdist]
  ring

lemma fourierChar_memLp (k : Fin 8 → ℝ) :
    MemLp (fourierChar k) 2 (volume.restrict E8FundamentalDomain) := by
  have := e8Domain_finiteMeasure
  refine MemLp.of_bound ?_ 1 ?_
  · refine Continuous.aestronglyMeasurable ?_
    unfold fourierChar
    fun_prop
  · filter_upwards with x
    exact (fourierChar_norm k x).le

def e8CharacterLp (k : Fin 8 → ℝ) : Lp ℂ 2 (volume.restrict E8FundamentalDomain) :=
  (fourierChar_memLp k).toLp (fourierChar k)

lemma e8CharacterLp_inner (k l : Fin 8 → ℝ) :
    inner ℂ (e8CharacterLp k) (e8CharacterLp l) =
      ∫ x in E8FundamentalDomain, star (fourierChar k x) * fourierChar l x := by
  rw [MeasureTheory.L2.inner_def]
  simp_rw [RCLike.inner_apply', ← RCLike.star_def]
  refine integral_congr_ae ?_
  have hk := (fourierChar_memLp k).coeFn_toLp
  have hl := (fourierChar_memLp l).coeFn_toLp
  filter_upwards [hk, hl] with x hxk hxl
  unfold e8CharacterLp
  rw [hxk, hxl]

theorem e8CharacterLp_orthonormal :
    Orthonormal ℂ (fun k : E8LatticeSet => e8CharacterLp k.val) := by
  rw [orthonormal_iff_ite]
  intro k l
  rw [e8CharacterLp_inner, e8_fourier_orthogonality k.val l.val k.property l.property]
  by_cases h : k = l
  · subst h
    simp
  · have hval : k.val ≠ l.val := Subtype.coe_ne_coe.mpr h
    rw [ite_eq_right hval, ite_eq_right h]

lemma poissonCoeff_eq_inner (f : Lp ℂ 2 (volume.restrict E8FundamentalDomain))
    (k : E8LatticeSet) :
    inner ℂ (e8CharacterLp k.val) f =
      ∫ x in E8FundamentalDomain, f x * fourierKernel k.val x := by
  rw [MeasureTheory.L2.inner_def]
  simp_rw [RCLike.inner_apply', ← RCLike.star_def]
  refine integral_congr_ae ?_
  filter_upwards [(fourierChar_memLp k.val).coeFn_toLp] with x hx
  unfold e8CharacterLp
  rw [hx, fourierKernel_eq_star_fourierChar, mul_comm]

lemma e8_lp_norm_sq (f : Lp ℂ 2 (volume.restrict E8FundamentalDomain)) :
    ‖f‖ ^ 2 = ∫ x in E8FundamentalDomain, Complex.normSq (f x) := by
  rw [InnerProductSpace.norm_sq_eq_re_inner (𝕜 := ℂ), MeasureTheory.L2.inner_def]
  have hint : Integrable (fun a => inner ℂ (f a) (f a))
      (volume.restrict E8FundamentalDomain) := by
    have hpow := (Lp.memLp f).integrable_norm_rpow (by simp) (by simp)
    have hpow' : Integrable (fun a => ‖(f a : ℂ)‖ ^ 2)
        (volume.restrict E8FundamentalDomain) := by
      simpa [ENNReal.toReal_ofNat] using hpow
    refine (Integrable.ofReal (𝕜 := ℂ) hpow').congr ?_
    filter_upwards with a
    rw [inner_self_eq_norm_sq_to_K]
    exact RCLike.ofReal_pow (K := ℂ) ‖f a‖ 2
  rw [← integral_re hint]
  have hre : (fun a => RCLike.re (inner ℂ (f a) (f a))) =
      fun a => Complex.normSq (f a) := by
    funext a
    rw [inner_self_eq_norm_sq_to_K, RCLike.re_ofReal_pow, Complex.normSq_eq_norm_sq]
  rw [hre]

theorem e8_fourier_bessel_inequality
    (f : Lp ℂ 2 (volume.restrict E8FundamentalDomain)) :
    ∑' k : E8LatticeSet, Complex.normSq
        (∫ x in E8FundamentalDomain, f x * fourierKernel k.val x) ≤
      ∫ x in E8FundamentalDomain, Complex.normSq (f x) := by
  have hterm : ∀ k : E8LatticeSet,
      Complex.normSq (∫ x in E8FundamentalDomain, f x * fourierKernel k.val x) =
        ‖inner ℂ (e8CharacterLp k.val) f‖ ^ 2 := by
    intro k
    rw [poissonCoeff_eq_inner, Complex.normSq_eq_norm_sq]
  simp_rw [hterm, ← e8_lp_norm_sq]
  exact e8CharacterLp_orthonormal.tsum_inner_products_le f

-- Gate 2C uses the probability Haar measure on the circle, the same normalisation as
-- `UnitAddTorus.mFourierLp`.
local instance (priority := 1001) : MeasureSpace UnitAddCircle := ⟨AddCircle.haarAddCircle⟩

open UnitAddTorus

def e8CoordMap : (Fin 8 → ℝ) →ₗ[ℝ] Fin 8 → ℝ :=
  Matrix.toLin' e8_basis_matrix.transpose

lemma e8CoordMap_det : LinearMap.det e8CoordMap = -1 := by
  unfold e8CoordMap
  rw [LinearMap.det_toLin', Matrix.det_transpose, e8_basis_matrix_det]

lemma e8CoordMap_det_ne_zero : LinearMap.det e8CoordMap ≠ 0 := by
  rw [e8CoordMap_det]
  norm_num

lemma e8CoordMap_apply (c : Fin 8 → ℝ) :
    e8CoordMap c = ∑ j, c j • e8_basis j := by
  unfold e8CoordMap
  ext r
  rw [Matrix.toLin'_apply, Matrix.mulVec_apply_eq_sum]
  simp only [Matrix.transpose_apply, e8_basis_matrix, Finset.sum_apply, Pi.smul_apply,
    smul_eq_mul]
  refine Finset.sum_congr rfl fun j _ => mul_comm _ _

lemma e8CoordMap_measurable : Measurable e8CoordMap :=
  (LinearMap.continuous_of_finiteDimensional e8CoordMap).measurable

noncomputable def e8CoordEquiv : (Fin 8 → ℝ) ≃L[ℝ] Fin 8 → ℝ :=
  let fc := e8CoordMap.toContinuousLinearMap
  have hdet : fc.det ≠ 0 := by
    change LinearMap.det (fc : (Fin 8 → ℝ) →ₗ[ℝ] Fin 8 → ℝ) ≠ 0
    rw [LinearMap.coe_toContinuousLinearMap]
    exact e8CoordMap_det_ne_zero
  fc.toContinuousLinearEquivOfDetNeZero hdet

lemma e8CoordEquiv_apply (c : Fin 8 → ℝ) : e8CoordEquiv c = e8CoordMap c :=
  e8CoordMap.toContinuousLinearMap.toContinuousLinearEquivOfDetNeZero_apply _ c

lemma e8CoordMap_measurePreserving : MeasurePreserving e8CoordMap := by
  refine ⟨e8CoordMap_measurable, ?_⟩
  rw [Real.map_linearMap_volume_pi_eq_smul_volume_pi e8CoordMap_det_ne_zero]
  have habs : |(LinearMap.det e8CoordMap)⁻¹| = 1 := by
    rw [e8CoordMap_det]
    norm_num
  simp [habs, ENNReal.ofReal_one, one_smul]

lemma e8CoordMap_embedding : MeasurableEmbedding e8CoordMap := by
  let fc := e8CoordMap.toContinuousLinearMap
  have hdet : fc.det ≠ 0 := by
    change LinearMap.det (fc : (Fin 8 → ℝ) →ₗ[ℝ] Fin 8 → ℝ) ≠ 0
    rw [LinearMap.coe_toContinuousLinearMap]
    exact e8CoordMap_det_ne_zero
  have hfun : (e8CoordMap : (Fin 8 → ℝ) → Fin 8 → ℝ) =
      (fc.toContinuousLinearEquivOfDetNeZero hdet : (Fin 8 → ℝ) → Fin 8 → ℝ) := by
    funext c
    exact (fc.toContinuousLinearEquivOfDetNeZero_apply hdet c).symm
  rw [hfun]
  exact (fc.toContinuousLinearEquivOfDetNeZero hdet).toHomeomorph.toMeasurableEquiv.measurableEmbedding

def e8IcoCube : Set (Fin 8 → ℝ) :=
  Set.univ.pi fun _ => Ico (0 : ℝ) 1

def e8IocCube : Set (Fin 8 → ℝ) :=
  Set.univ.pi fun _ => Ioc (0 : ℝ) 1

lemma E8FundamentalDomain_eq_image : E8FundamentalDomain = e8CoordMap '' e8IcoCube := by
  ext x
  constructor
  · rintro ⟨c, hc, rfl⟩
    refine ⟨c, ?_, ?_⟩
    · intro i _
      simpa [mem_Ico] using hc i
    · rw [e8CoordMap_apply]
  · rintro ⟨c, hc, rfl⟩
    refine ⟨c, ?_, e8CoordMap_apply c⟩
    intro i
    simpa [mem_Ico] using hc i (mem_univ i)

lemma E8FundamentalDomain_measurable : MeasurableSet E8FundamentalDomain := by
  rw [E8FundamentalDomain_eq_image]
  exact e8CoordMap_embedding.measurableSet_image'
    (MeasurableSet.univ_pi fun _ => measurableSet_Ico)

lemma e8IcoCube_ae_eq_iocCube : e8IcoCube =ᵐ[volume] e8IocCube := by
  rw [volume_pi, e8IcoCube, e8IocCube]
  exact (Measure.univ_pi_Ico_ae_eq_Icc (μ := fun _ : Fin 8 => (volume : Measure ℝ))
      (f := fun _ => (0 : ℝ)) (g := fun _ => (1 : ℝ))).trans
    (Measure.univ_pi_Ioc_ae_eq_Icc (μ := fun _ : Fin 8 => (volume : Measure ℝ))
      (f := fun _ => (0 : ℝ)) (g := fun _ => (1 : ℝ))).symm

lemma e8CoordMap_image_ae (s t : Set (Fin 8 → ℝ))
    (hs : MeasurableSet s) (ht : MeasurableSet t) (hst : s =ᵐ[volume] t) :
    e8CoordMap '' s =ᵐ[volume] e8CoordMap '' t := by
  rw [← measure_symmDiff_eq_zero_iff] at hst ⊢
  have hms : MeasurableSet (symmDiff s t) := hs.symmDiff ht
  have hnull : volume (e8CoordMap '' symmDiff s t) = 0 := by
    rw [← e8CoordMap_measurePreserving.map_eq, e8CoordMap_embedding.map_apply,
      Set.preimage_image_eq _ e8CoordMap_embedding.injective, hst]
  refine measure_mono_null ?_ hnull
  intro y hy
  simp only [mem_symmDiff, mem_image] at hy ⊢
  rcases hy with ⟨⟨c, hc, rfl⟩, hy⟩ | ⟨⟨c, hc, rfl⟩, hy⟩
  · exact ⟨c, Or.inl ⟨hc, fun hct => hy ⟨c, hct, rfl⟩⟩, rfl⟩
  · exact ⟨c, Or.inr ⟨hc, fun hcs => hy ⟨c, hcs, rfl⟩⟩, rfl⟩

lemma e8_iocImage_ae_domain :
    e8CoordMap '' e8IocCube =ᵐ[volume] E8FundamentalDomain := by
  rw [E8FundamentalDomain_eq_image]
  exact e8CoordMap_image_ae e8IocCube e8IcoCube
    (MeasurableSet.univ_pi fun _ => measurableSet_Ioc)
    (MeasurableSet.univ_pi fun _ => measurableSet_Ico)
    e8IcoCube_ae_eq_iocCube.symm

def e8IcoRep : Set (Fin 8 → ℝ) :=
  {x | ∀ i, x i ∈ Ico (0 : ℝ) ((0 : ℝ) + 1)}

def e8IocRep : Set (Fin 8 → ℝ) :=
  {x | ∀ i, x i ∈ Ioc (0 : ℝ) ((0 : ℝ) + 1)}

lemma e8IcoRep_eq_cube : e8IcoRep = e8IcoCube := by
  ext x
  simp [e8IcoRep, e8IcoCube, mem_univ_pi, zero_add, mem_Ico]

lemma e8IocRep_eq_cube : e8IocRep = e8IocCube := by
  ext x
  simp [e8IocRep, e8IocCube, mem_univ_pi, zero_add, mem_Ioc]

noncomputable def e8TorusToIcoEquiv :
    UnitAddTorus (Fin 8) ≃ᵐ {x : Fin 8 → ℝ // x ∈ e8IcoRep} :=
  (MeasurableEquiv.piCongrRight fun _ : Fin 8 => AddCircle.measurableEquivIco 1 (0 : ℝ)).trans
    MeasurableEquiv.subtypePiEquivPi.symm

def e8TorusToDomain (x : UnitAddTorus (Fin 8)) : Fin 8 → ℝ :=
  e8CoordMap (e8TorusToIcoEquiv x).1

def e8DomainToTorus (y : Fin 8 → ℝ) : UnitAddTorus (Fin 8) :=
  fun i => (e8CoordEquiv.symm y i : UnitAddCircle)

lemma e8DomainToTorus_measurable : Measurable e8DomainToTorus := by
  refine Measurable.of_eval fun i => ?_
  exact AddCircle.measurable_mk'.comp
    ((continuous_apply i).comp e8CoordEquiv.continuous_invFun).measurable

lemma e8TorusToDomain_measurable : Measurable e8TorusToDomain :=
  e8CoordMap_measurable.comp (measurable_subtype_coe.comp e8TorusToIcoEquiv.measurable)

lemma equivIco_eq_equivIoc_of_ne_zero {x : UnitAddCircle} (hx : x ≠ 0) :
    (AddCircle.equivIco 1 (0 : ℝ) x : ℝ) = (AddCircle.equivIoc 1 (0 : ℝ) x : ℝ) := by
  let r : ℝ := AddCircle.equivIco 1 0 x
  have hr : r ∈ Ico (0 : ℝ) 1 := by
    simpa [zero_add] using (AddCircle.equivIco 1 (0 : ℝ) x).2
  have hr0 : r ≠ 0 := by
    intro hr0
    apply hx
    have hx0 : (r : UnitAddCircle) = x := by
      simpa [r] using (AddCircle.coe_equivIco (p := (1 : ℝ)) (a := (0 : ℝ)) (y := x)).symm
    simpa [hr0] using hx0.symm
  have hrIoc : r ∈ Ioc (0 : ℝ) ((0 : ℝ) + 1) := by
    rw [zero_add]
    exact ⟨lt_of_le_of_ne hr.1 (Ne.symm hr0), le_of_lt hr.2⟩
  have hxcoe : (r : UnitAddCircle) = x := by
    simpa [r] using (AddCircle.coe_equivIco (p := (1 : ℝ)) (a := (0 : ℝ)) (y := x)).symm
  have heq : AddCircle.equivIoc 1 0 x = ⟨r, hrIoc⟩ := by
    rw [← hxcoe]
    exact AddCircle.equivIoc_coe_eq hrIoc
  exact (congrArg Subtype.val heq).symm

lemma e8TorusToIcoEquiv_val (x : UnitAddTorus (Fin 8)) :
    (e8TorusToIcoEquiv x).1 = fun i => (AddCircle.equivIco 1 (0 : ℝ) (x i)).1 := by
  rfl

lemma e8_torus_zero_null (i : Fin 8) :
    volume {x : UnitAddTorus (Fin 8) | x i = 0} = 0 := by
  let b : Fin 8 → ℝ := fun _ => 0
  let e := measurableEquivPiIoc b
  have hmp := measurePreserving_equivPiIoc b
  have hrep0 : AddCircle.equivIoc 1 (0 : ℝ) (0 : UnitAddCircle) =
      ⟨(1 : ℝ), by rw [zero_add]; simp [mem_Ioc]⟩ := by
    have hper : ((1 : ℝ) : UnitAddCircle) = 0 := AddCircle.coe_period 1
    rw [← hper]
    exact AddCircle.equivIoc_coe_eq (by rw [zero_add]; simp [mem_Ioc])
  let t : Set {x : Fin 8 → ℝ // ∀ j, x j ∈ Ioc (b j) (b j + 1)} :=
    {y | y.1 i = 1}
  have ht' : MeasurableSet t := by
    have hpre : t = ((fun z : Fin 8 → ℝ => z i) ∘ Subtype.val) ⁻¹' ({1} : Set ℝ) := by
      ext y
      simp [t]
    rw [hpre]
    exact measurableSet_preimage ((measurable_pi_apply i).comp measurable_subtype_coe)
      (measurableSet_singleton (1 : ℝ))
  have hface : volume {z : Fin 8 → ℝ | z i = 1} = 0 := by
    rw [volume_pi]
    exact Measure.pi_hyperplane (fun _ => (volume : Measure ℝ)) i (1 : ℝ)
  have hsub : Subtype.val '' t ⊆ {z : Fin 8 → ℝ | z i = 1} := by
    intro z hz
    rcases hz with ⟨y, hy, rfl⟩
    have hyi : y.1 i = 1 := by
      simpa only [t, Set.mem_setOf_eq] using hy
    exact hyi
  have hcomap : (Measure.comap Subtype.val volume) t = 0 := by
    refine le_antisymm ?_ (by simp)
    exact le_trans (Measure.comap_apply_le Subtype.val volume ht'.nullMeasurableSet)
      (le_trans (measure_mono hsub) (le_of_eq hface))
  have hpre : {x : UnitAddTorus (Fin 8) | x i = 0} ⊆ e ⁻¹' t := by
    intro x hx
    have hx0 : x i = 0 := by simpa [Set.mem_setOf_eq] using hx
    show (e x).1 i = 1
    rw [coe_measurableEquivPiIoc_apply]
    change (AddCircle.equivIoc 1 (b i) (x i)).1 = 1
    rw [show b i = 0 from rfl, hx0, hrep0]
  have hzero : volume (e ⁻¹' t) = 0 := by
    rw [hmp.measure_preimage_equiv t, hcomap]
  exact measure_mono_null hpre hzero

lemma e8_coordRep_ae :
    (fun x : UnitAddTorus (Fin 8) => (e8TorusToIcoEquiv x).1) =ᵐ[volume]
      fun x => (measurableEquivPiIoc (fun _ : Fin 8 => (0 : ℝ)) x).1 := by
  unfold Filter.EventuallyEq
  rw [ae_iff]
  refine measure_mono_null ?_ (measure_iUnion_null e8_torus_zero_null)
  intro x hx
  simp only [mem_setOf, mem_iUnion] at hx ⊢
  rw [e8TorusToIcoEquiv_val, coe_measurableEquivPiIoc_apply] at hx
  obtain ⟨i, hi⟩ := Function.ne_iff.mp hx
  refine ⟨i, ?_⟩
  by_contra hne
  exact hi (equivIco_eq_equivIoc_of_ne_zero hne)

lemma e8TorusToDomain_ae_ioc :
    e8TorusToDomain =ᵐ[volume] fun x : UnitAddTorus (Fin 8) =>
      e8CoordMap (measurableEquivPiIoc (fun _ : Fin 8 => (0 : ℝ)) x).1 := by
  filter_upwards [e8_coordRep_ae] with x hx
  simp [e8TorusToDomain, hx]

lemma e8_iocChart_measurePreserving :
    MeasurePreserving
      (fun x : UnitAddTorus (Fin 8) =>
        e8CoordMap (measurableEquivPiIoc (fun _ : Fin 8 => (0 : ℝ)) x).1)
      volume (volume.restrict (e8CoordMap '' e8IocRep)) := by
  have hEquiv := measurePreserving_equivPiIoc (fun _ : Fin 8 => (0 : ℝ))
  have hs : MeasurableSet e8IocRep := by
    rw [e8IocRep_eq_cube]
    exact MeasurableSet.univ_pi fun _ => measurableSet_Ioc
  have hVal : MeasurePreserving (Subtype.val : {x : Fin 8 → ℝ // x ∈ e8IocRep} → Fin 8 → ℝ)
      (Measure.comap Subtype.val (volume : Measure (Fin 8 → ℝ))) (volume.restrict e8IocRep) :=
    ⟨measurable_subtype_coe, map_comap_subtype_coe hs volume⟩
  have hLin := e8CoordMap_measurePreserving.restrict_image_emb e8CoordMap_embedding e8IocRep
  exact hLin.comp (hVal.comp hEquiv)

lemma e8TorusToDomain_measurePreserving :
    MeasurePreserving e8TorusToDomain volume (volume.restrict E8FundamentalDomain) := by
  have hCongr := e8_iocChart_measurePreserving.congr e8TorusToDomain_measurable
    e8TorusToDomain_ae_ioc.symm
  have hrest : volume.restrict (e8CoordMap '' e8IocRep) =
      volume.restrict E8FundamentalDomain := by
    apply Measure.restrict_congr_set
    rw [e8IocRep_eq_cube]
    exact e8_iocImage_ae_domain
  exact ⟨hCongr.measurable, hCongr.map_eq.trans hrest⟩

lemma e8_chart_left_inv (x : UnitAddTorus (Fin 8)) :
    e8DomainToTorus (e8TorusToDomain x) = x := by
  funext i
  have hcoord : e8CoordEquiv.symm (e8TorusToDomain x) = (e8TorusToIcoEquiv x).1 := by
    rw [e8TorusToDomain, ← e8CoordEquiv_apply]
    exact e8CoordEquiv.symm_apply_apply _
  rw [e8DomainToTorus, congr_fun hcoord i, e8TorusToIcoEquiv_val]
  simpa using (AddCircle.coe_equivIco (p := (1 : ℝ)) (a := (0 : ℝ)) (y := x i))

lemma e8DomainToTorus_measurePreserving :
    MeasurePreserving e8DomainToTorus (volume.restrict E8FundamentalDomain) volume := by
  have hcomp_fun : e8DomainToTorus ∘ e8TorusToDomain = id := funext e8_chart_left_inv
  have hcomp : MeasurePreserving (e8DomainToTorus ∘ e8TorusToDomain) volume volume := by
    rw [hcomp_fun]
    exact MeasurePreserving.id volume
  simpa [e8TorusToDomain_measurePreserving.map_eq] using
    hcomp.map_of_comp e8DomainToTorus_measurable e8TorusToDomain_measurable

lemma e8CoordMap_mem_domain {c : Fin 8 → ℝ} (hc : c ∈ e8IcoRep) :
    e8CoordMap c ∈ E8FundamentalDomain := by
  rw [e8IcoRep_eq_cube] at hc
  refine ⟨c, ?_, e8CoordMap_apply c⟩
  intro i
  simpa [mem_Ico] using hc i (mem_univ i)

lemma e8CoordEquiv_symm_mem_ico {y : Fin 8 → ℝ} (hy : y ∈ E8FundamentalDomain) :
    e8CoordEquiv.symm y ∈ e8IcoRep := by
  rcases hy with ⟨c, hc, hy⟩
  have hsymm : e8CoordEquiv.symm y = c := by
    have hy' : y = e8CoordMap c := by rw [hy, e8CoordMap_apply]
    rw [hy', ← e8CoordEquiv_apply]
    exact e8CoordEquiv.symm_apply_apply c
  rw [hsymm, e8IcoRep_eq_cube]
  intro i _
  simpa [mem_Ico] using hc i

noncomputable def e8CubeToDomainEquiv :
    {x : Fin 8 → ℝ // x ∈ e8IcoRep} ≃ᵐ ↥E8FundamentalDomain where
  toEquiv :=
    { toFun := fun c => ⟨e8CoordMap c.1, e8CoordMap_mem_domain c.2⟩
      invFun := fun y => ⟨e8CoordEquiv.symm y.1, e8CoordEquiv_symm_mem_ico y.2⟩
      left_inv := fun c => by
        apply Subtype.ext
        convert e8CoordEquiv.symm_apply_apply c.1
        rw [e8CoordEquiv_apply]
      right_inv := fun y => by
        apply Subtype.ext
        convert e8CoordEquiv.apply_symm_apply y.1
        rw [e8CoordEquiv_apply] }
  measurable_toFun := Measurable.subtype_mk (e8CoordMap_measurable.comp measurable_subtype_coe)
  measurable_invFun :=
    Measurable.subtype_mk (e8CoordEquiv.continuous_invFun.measurable.comp measurable_subtype_coe)

/-- G2C.1: Maßtreue messbare Äquivalenz zwischen dem Standardtorus und der Fundamentaldomäne. -/
noncomputable def e8FundamentalDomainMeasureEquiv :
    UnitAddTorus (Fin 8) ≃ᵐ ↥E8FundamentalDomain :=
  e8TorusToIcoEquiv.trans e8CubeToDomainEquiv

lemma e8_chart_right_inv {y : Fin 8 → ℝ} (hy : y ∈ E8FundamentalDomain) :
    e8TorusToDomain (e8DomainToTorus y) = y := by
  have hc : e8CoordEquiv.symm y ∈ e8IcoRep := e8CoordEquiv_symm_mem_ico hy
  have hBy : e8CoordMap (e8CoordEquiv.symm y) = y := by
    rw [← e8CoordEquiv_apply, e8CoordEquiv.apply_symm_apply]
  have htorus : e8DomainToTorus y = fun i => ((e8CoordEquiv.symm y) i : UnitAddCircle) := rfl
  have hrep : (e8TorusToIcoEquiv (e8DomainToTorus y)).1 = e8CoordEquiv.symm y := by
    rw [e8TorusToIcoEquiv_val, htorus]
    ext i
    have hmem : (e8CoordEquiv.symm y) i ∈ Ico (0 : ℝ) ((0 : ℝ) + 1) := by
      rw [e8IcoRep] at hc
      simpa using hc i
    simpa using congrArg Subtype.val (AddCircle.equivIco_coe_eq hmem)
  rw [e8TorusToDomain, hrep, hBy]

lemma e8_chart_right_inv_ae :
    e8TorusToDomain ∘ e8DomainToTorus =ᵐ[volume.restrict E8FundamentalDomain] id := by
  filter_upwards [ae_restrict_mem E8FundamentalDomain_measurable] with y hy
  exact e8_chart_right_inv hy

lemma e8_gram_z_mul_adjugate : e8_gram_z * adjugate e8_gram_z = 1 := by
  simpa [e8_gram_z_det_one] using mul_adjugate e8_gram_z

lemma e8_gram_z_adjugate_mul : adjugate e8_gram_z * e8_gram_z = 1 := by
  simpa [e8_gram_z_det_one] using adjugate_mul e8_gram_z

def z8ToE8Coord (m : Fin 8 → ℤ) : Fin 8 → ℤ :=
  adjugate e8_gram_z *ᵥ m

def z8ToE8Fun (m : Fin 8 → ℤ) : Fin 8 → ℝ :=
  ∑ i, (z8ToE8Coord m i : ℝ) • e8_basis i

def z8ToE8 (m : Fin 8 → ℤ) : E8LatticeSet :=
  ⟨z8ToE8Fun m, z8ToE8Coord m, rfl⟩

lemma z8ToE8_pairing (m : Fin 8 → ℤ) (j : Fin 8) :
    ∑ k, e8_basis j k * z8ToE8Fun m k = (m j : ℝ) := by
  have hcomp : ∀ k, z8ToE8Fun m k = ∑ i, (z8ToE8Coord m i : ℝ) * e8_basis i k := by
    intro k
    simp [z8ToE8Fun, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  simp_rw [hcomp, e8_basis_dot_sum]
  have hcoord : (fun i => (z8ToE8Coord m i : ℝ)) =
      (adjugate e8_gram_z).map (Int.castRingHom ℝ) *ᵥ fun i => (m i : ℝ) := by
    ext i
    rw [z8ToE8Coord]
    exact RingHom.map_mulVec (Int.castRingHom ℝ) (adjugate e8_gram_z) m i
  have hmul : e8_gram_z.map (Int.castRingHom ℝ) *
      (adjugate e8_gram_z).map (Int.castRingHom ℝ) = 1 := by
    rw [← Matrix.map_mul, e8_gram_z_mul_adjugate]
    exact Matrix.map_one (Int.castRingHom ℝ) (RingHom.map_zero _) (RingHom.map_one _)
  have hvec : (e8_gram_z.map (Int.castRingHom ℝ) *ᵥ fun i => (z8ToE8Coord m i : ℝ)) j =
      (m j : ℝ) := by
    rw [hcoord, Matrix.mulVec_mulVec, hmul, Matrix.one_mulVec]
  simpa [Matrix.mulVec_apply_eq_sum, Matrix.map_apply] using hvec

lemma z8ToE8Coord_left_inv (m : Fin 8 → ℤ) : e8_gram_z *ᵥ z8ToE8Coord m = m := by
  simp [z8ToE8Coord, Matrix.mulVec_mulVec, e8_gram_z_mul_adjugate, Matrix.one_mulVec]

lemma z8ToE8Coord_right_inv (k : Fin 8 → ℤ) : z8ToE8Coord (e8_gram_z *ᵥ k) = k := by
  simp [z8ToE8Coord, Matrix.mulVec_mulVec, e8_gram_z_adjugate_mul, Matrix.one_mulVec]

lemma z8ToE8Fun_injective : Function.Injective z8ToE8Fun := by
  intro m n h
  have h0 : ∑ i, ((z8ToE8Coord m i - z8ToE8Coord n i : ℤ) : ℝ) • e8_basis i = 0 := by
    ext r
    have hm : z8ToE8Fun m r = ∑ i, (z8ToE8Coord m i : ℝ) * e8_basis i r := by
      simp [z8ToE8Fun, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    have hn : z8ToE8Fun n r = ∑ i, (z8ToE8Coord n i : ℝ) * e8_basis i r := by
      simp [z8ToE8Fun, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    rw [Finset.sum_apply, Pi.zero_apply]
    simp_rw [Pi.smul_apply, smul_eq_mul, Int.cast_sub, sub_mul, Finset.sum_sub_distrib]
    rw [← hm, ← hn, congr_fun h r, sub_self]
  have hli := (Fintype.linearIndependent_iff.mp e8_basis_linearIndependent)
    (fun i => ((z8ToE8Coord m i - z8ToE8Coord n i : ℤ) : ℝ)) h0
  have hcoord : z8ToE8Coord m = z8ToE8Coord n := by
    ext i
    exact sub_eq_zero.mp (Int.cast_eq_zero.mp (hli i))
  rw [← z8ToE8Coord_left_inv m, ← z8ToE8Coord_left_inv n, hcoord]

lemma z8ToE8_bijective : Function.Bijective z8ToE8 := by
  constructor
  · intro m n h
    exact z8ToE8Fun_injective (congrArg Subtype.val h)
  · intro q
    rcases q.2 with ⟨k, hk⟩
    refine ⟨e8_gram_z *ᵥ k, ?_⟩
    apply Subtype.ext
    rw [hk]
    simp [z8ToE8, z8ToE8Fun, z8ToE8Coord_right_inv]

/-- Φ : E₈ ≃ ℤ⁸, so `e8FrequencyEquiv.symm = z8ToE8`. -/
noncomputable def e8FrequencyEquiv : E8LatticeSet ≃ (Fin 8 → ℤ) :=
  (Equiv.ofBijective z8ToE8 z8ToE8_bijective).symm

lemma e8FrequencyEquiv_symm_apply (m : Fin 8 → ℤ) :
    e8FrequencyEquiv.symm m = z8ToE8 m := rfl

lemma fourierChar_z8ToE8 (m : Fin 8 → ℤ) (c : Fin 8 → ℝ) :
    fourierChar (z8ToE8Fun m) (e8CoordMap c) =
      ∏ j, Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ)) := by
  obtain ⟨m', hm_pair, hm_prod⟩ :=
    fourierChar_coord_factorization_all (z8ToE8Fun m) (z8ToE8 m).2
  have hm' : m' = m := by
    ext j
    apply Int.cast_injective (α := ℝ)
    rw [← hm_pair j, z8ToE8_pairing]
  rw [e8CoordMap_apply, hm_prod, hm']

lemma mFourier_on_coords (m : Fin 8 → ℤ) (c : Fin 8 → ℝ) :
    mFourier m (fun i => (c i : UnitAddCircle)) =
      ∏ j, Complex.exp (2 * Real.pi * Complex.I * (m j : ℂ) * (c j : ℂ)) := by
  simp only [mFourier, ContinuousMap.coe_mk]
  refine Finset.prod_congr rfl fun j _ => ?_
  rw [fourier_coe_apply]
  simp [div_one]

lemma mFourier_eq_fourierChar_on_domain (m : Fin 8 → ℤ) {y : Fin 8 → ℝ}
    (hy : y ∈ E8FundamentalDomain) :
    mFourier m (e8DomainToTorus y) = fourierChar (z8ToE8Fun m) y := by
  have hc : e8CoordEquiv.symm y ∈ e8IcoRep := e8CoordEquiv_symm_mem_ico hy
  have hBy : e8CoordMap (e8CoordEquiv.symm y) = y := by
    rw [← e8CoordEquiv_apply, e8CoordEquiv.apply_symm_apply]
  have htorus : e8DomainToTorus y = fun i => ((e8CoordEquiv.symm y) i : UnitAddCircle) := rfl
  rw [htorus, mFourier_on_coords, ← fourierChar_z8ToE8, hBy]

noncomputable def e8L2Push :
    Lp ℂ 2 (volume : Measure (UnitAddTorus (Fin 8))) →ₗᵢ[ℂ]
      Lp ℂ 2 (volume.restrict E8FundamentalDomain) :=
  Lp.compMeasurePreservingₗᵢ ℂ e8DomainToTorus e8DomainToTorus_measurePreserving

noncomputable def e8L2Pull :
    Lp ℂ 2 (volume.restrict E8FundamentalDomain) →ₗᵢ[ℂ]
      Lp ℂ 2 (volume : Measure (UnitAddTorus (Fin 8))) :=
  Lp.compMeasurePreservingₗᵢ ℂ e8TorusToDomain e8TorusToDomain_measurePreserving

lemma e8L2Pull_push (f : Lp ℂ 2 (volume : Measure (UnitAddTorus (Fin 8)))) :
    e8L2Pull (e8L2Push f) = f := by
  unfold e8L2Pull e8L2Push
  change (Lp.compMeasurePreserving e8TorusToDomain e8TorusToDomain_measurePreserving)
      ((Lp.compMeasurePreserving e8DomainToTorus e8DomainToTorus_measurePreserving) f) = f
  rw [← Lp.compMeasurePreserving_comp_apply]
  apply Subtype.ext
  rw [Lp.compMeasurePreserving_val]
  have hfun : e8DomainToTorus ∘ e8TorusToDomain = id := funext e8_chart_left_inv
  have hae : e8DomainToTorus ∘ e8TorusToDomain =ᵐ[volume] id := by
    rw [hfun]
  have hmp := e8DomainToTorus_measurePreserving.comp e8TorusToDomain_measurePreserving
  rw [AEEqFun.compMeasurePreserving_congr f.1 hmp measurable_id hae]
  have hid : hmp.congr measurable_id hae = MeasurePreserving.id volume := Subsingleton.elim _ _
  rw [hid, AEEqFun.compMeasurePreserving_id]

lemma e8L2Push_pull (f : Lp ℂ 2 (volume.restrict E8FundamentalDomain)) :
    e8L2Push (e8L2Pull f) = f := by
  unfold e8L2Push e8L2Pull
  change (Lp.compMeasurePreserving e8DomainToTorus e8DomainToTorus_measurePreserving)
      ((Lp.compMeasurePreserving e8TorusToDomain e8TorusToDomain_measurePreserving) f) = f
  rw [← Lp.compMeasurePreserving_comp_apply]
  apply Subtype.ext
  rw [Lp.compMeasurePreserving_val]
  have hmp := e8TorusToDomain_measurePreserving.comp e8DomainToTorus_measurePreserving
  rw [AEEqFun.compMeasurePreserving_congr f.1 hmp measurable_id e8_chart_right_inv_ae]
  have hid : hmp.congr measurable_id e8_chart_right_inv_ae =
      MeasurePreserving.id (volume.restrict E8FundamentalDomain) := Subsingleton.elim _ _
  rw [hid, AEEqFun.compMeasurePreserving_id]

/-- G2C.2: Die von der Koordinatenabbildung induzierte L²-Isometrie. -/
noncomputable def e8L2CoordinateEquiv :
    Lp ℂ 2 (volume : Measure (UnitAddTorus (Fin 8))) ≃ₗᵢ[ℂ]
      Lp ℂ 2 (volume.restrict E8FundamentalDomain) :=
  LinearIsometryEquiv.ofLinearIsometry e8L2Push e8L2Pull.toLinearMap
    (by
      apply LinearMap.ext
      intro f
      exact e8L2Push_pull f)
    (by
      apply LinearMap.ext
      intro f
      exact e8L2Pull_push f)

/-- G2C.3: `U_B` schickt die Standardmonome auf die E₈-Charaktere. -/
lemma e8_mFourier_transport (m : Fin 8 → ℤ) :
    e8L2CoordinateEquiv (mFourierLp 2 m) = e8CharacterLp (z8ToE8 m).val := by
  apply Lp.ext
  have hpush : e8L2CoordinateEquiv (mFourierLp 2 m) =ᵐ[volume.restrict E8FundamentalDomain]
      fun y => (mFourierLp 2 m : UnitAddTorus (Fin 8) → ℂ) (e8DomainToTorus y) := by
    rw [e8L2CoordinateEquiv, LinearIsometryEquiv.coe_ofLinearIsometry, e8L2Push]
    change Lp.compMeasurePreserving e8DomainToTorus e8DomainToTorus_measurePreserving
        (mFourierLp 2 m) =ᵐ[volume.restrict E8FundamentalDomain]
      fun y => (mFourierLp 2 m : UnitAddTorus (Fin 8) → ℂ) (e8DomainToTorus y)
    exact Lp.coeFn_compMeasurePreserving (mFourierLp 2 m) e8DomainToTorus_measurePreserving
  have hcont : (fun y => (mFourierLp 2 m : UnitAddTorus (Fin 8) → ℂ) (e8DomainToTorus y))
      =ᵐ[volume.restrict E8FundamentalDomain] fun y => mFourier m (e8DomainToTorus y) :=
    e8DomainToTorus_measurePreserving.quasiMeasurePreserving.ae_eq_comp (coeFn_mFourierLp 2 m)
  have hchar : (fun y => mFourier m (e8DomainToTorus y)) =ᵐ[volume.restrict E8FundamentalDomain]
      fun y => fourierChar (z8ToE8Fun m) y := by
    filter_upwards [ae_restrict_mem E8FundamentalDomain_measurable] with y hy
    exact mFourier_eq_fourierChar_on_domain m hy
  have hLp : e8CharacterLp (z8ToE8 m).val =ᵐ[volume.restrict E8FundamentalDomain]
      fun y => fourierChar (z8ToE8Fun m) y := by
    unfold e8CharacterLp
    simpa [z8ToE8] using (fourierChar_memLp (z8ToE8Fun m)).coeFn_toLp
  exact hpush.trans (hcont.trans (hchar.trans hLp.symm))

end E8PQC

namespace HilbertBasis

variable {ι κ 𝕜 E : Type*} [RCLike 𝕜] [NormedAddCommGroup E] [InnerProductSpace 𝕜 E]
  [CompleteSpace E]

noncomputable def reindex (b : HilbertBasis ι 𝕜 E) (e : ι ≃ κ) : HilbertBasis κ 𝕜 E :=
  HilbertBasis.mk (b.orthonormal.comp e.symm e.symm.injective) <| by
    have hr : Set.range (fun k : κ => b (e.symm k)) = Set.range (⇑b) := by
      ext x
      constructor
      · rintro ⟨k, rfl⟩
        exact ⟨e.symm k, rfl⟩
      · rintro ⟨i, rfl⟩
        exact ⟨e i, by simp⟩
    simpa [hr, Function.comp] using b.dense_span.ge

lemma reindex_apply (b : HilbertBasis ι 𝕜 E) (e : ι ≃ κ) (k : κ) :
    b.reindex e k = b (e.symm k) := by
  simpa [reindex, Function.comp] using congr_fun
    (HilbertBasis.coe_mk (b.orthonormal.comp e.symm e.symm.injective) _) k

noncomputable def map {F : Type*} [NormedAddCommGroup F] [InnerProductSpace 𝕜 F]
    (b : HilbertBasis ι 𝕜 E) (L : E ≃ₗᵢ[𝕜] F) : HilbertBasis ι 𝕜 F :=
  HilbertBasis.ofRepr (L.symm.trans b.repr)

lemma map_apply {F : Type*} [NormedAddCommGroup F] [InnerProductSpace 𝕜 F]
    (b : HilbertBasis ι 𝕜 E) (L : E ≃ₗᵢ[𝕜] F) (i : ι) : b.map L i = L (b i) := by
  classical
  change (L.symm.trans b.repr).symm (lp.single 2 i (1 : 𝕜)) = L (b i)
  rw [LinearIsometryEquiv.symm_trans, LinearIsometryEquiv.trans_apply, b.repr_symm_single,
    LinearIsometryEquiv.symm_symm]

end HilbertBasis

namespace E8PQC

open UnitAddTorus

local instance (priority := 1001) : MeasureSpace UnitAddCircle := ⟨AddCircle.haarAddCircle⟩

/-- Die transportierte E₈-Hilbertbasis von `L²(P_{E8})`. -/
noncomputable def e8FourierBasis :
    HilbertBasis E8LatticeSet ℂ (Lp ℂ 2 (volume.restrict E8FundamentalDomain)) :=
  (mFourierBasis (d := Fin 8)).reindex e8FrequencyEquiv.symm |>.map e8L2CoordinateEquiv

lemma e8FourierBasis_apply (k : E8LatticeSet) :
    e8FourierBasis k = e8CharacterLp k.val := by
  rw [e8FourierBasis, HilbertBasis.map_apply, HilbertBasis.reindex_apply, coe_mFourierBasis,
    e8_mFourier_transport, ← e8FrequencyEquiv_symm_apply, Equiv.symm_symm,
    Equiv.symm_apply_apply]

/-- Gate 2, Parseval auf `L²(P_{E8})`. -/
theorem poisson_plancherel_e8_isometry
    (f : Lp ℂ 2 (volume.restrict E8FundamentalDomain)) :
    ∫ x in E8FundamentalDomain, Complex.normSq (f x) =
      ∑' k : E8LatticeSet, Complex.normSq
        (∫ x in E8FundamentalDomain, f x * fourierKernel k.val x) := by
  have hsum := e8FourierBasis.tsum_inner_mul_inner f f
  simp_rw [e8FourierBasis_apply] at hsum
  have hcoeff : ∀ k : E8LatticeSet,
      inner ℂ f (e8CharacterLp k.val) * inner ℂ (e8CharacterLp k.val) f =
        ((‖inner ℂ (e8CharacterLp k.val) f‖ ^ 2 : ℝ) : ℂ) := by
    intro k
    rw [← inner_conj_symm f (e8CharacterLp k.val), starRingEnd_apply, RCLike.star_def,
      RCLike.conj_mul, ← RCLike.ofReal_pow]
    rfl
  simp_rw [hcoeff] at hsum
  rw [← Complex.ofReal_tsum, inner_self_eq_norm_sq_to_K,
    ← RCLike.ofReal_pow (K := ℂ)] at hsum
  have hparseval : ∑' k : E8LatticeSet, ‖inner ℂ (e8CharacterLp k.val) f‖ ^ 2 = ‖f‖ ^ 2 :=
    Complex.ofReal_injective hsum
  have hterm : ∀ k : E8LatticeSet,
      Complex.normSq (∫ x in E8FundamentalDomain, f x * fourierKernel k.val x) =
        ‖inner ℂ (e8CharacterLp k.val) f‖ ^ 2 := by
    intro k
    rw [← poissonCoeff_eq_inner, Complex.normSq_eq_norm_sq]
  rw [← e8_lp_norm_sq, ← hparseval]
  simp_rw [hterm]

theorem poisson_plancherel_isometry
    (f : (Fin 8 → ℝ) → ℂ)
    (hf_L2 : IntegrableOn (fun x => Complex.normSq (f x)) E8FundamentalDomain) :
    ∫ x in E8FundamentalDomain, Complex.normSq (f x) =
      ∑' (k : E8LatticeSet), Complex.normSq
        (∫ x in E8FundamentalDomain, f x * fourierKernel k.val x) := by
  sorry

end E8PQC
