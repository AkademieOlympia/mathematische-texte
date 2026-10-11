import Mathlib

open Set Metric Topology Real Module
open scoped BigOperators

noncomputable section

namespace AnalyticNumberTheory

/-- Quadrat der euklidischen Norm auf `ℝ⁸`, also `∑ i, x i ^ 2`. -/
def normSq (x : Fin 8 → ℝ) : ℝ := ∑ i, x i ^ 2

/-- Gerades, volldimensionales Gitter in `ℝ⁸` mit Gram-Determinante `1`. -/
structure SelfDualLattice where
  carrier : Submodule ℤ (Fin 8 → ℝ)
  discrete : DiscreteTopology carrier
  even_norm : ∀ x ∈ carrier, ∃ n : ℤ, normSq x = 2 * (n : ℝ)
  basis : Basis (Fin 8) ℤ carrier
  gram_det_one :
    (Matrix.of fun i j => ∑ k, (basis i : Fin 8 → ℝ) k * (basis j : Fin 8 → ℝ) k).det = 1

/-- Gram-Matrix der `ℤ`-Basis. -/
def SelfDualLattice.gramMatrix (L : SelfDualLattice) : Matrix (Fin 8) (Fin 8) ℝ :=
  Matrix.of fun i j => ∑ k, (L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k

lemma SelfDualLattice.gramMatrix_det (L : SelfDualLattice) : L.gramMatrix.det = 1 :=
  L.gram_det_one

/-- Matrix der Basisvektoren als Zeilen. -/
def SelfDualLattice.basisMatrix (L : SelfDualLattice) : Matrix (Fin 8) (Fin 8) ℝ :=
  Matrix.of fun i j => (L.basis i : Fin 8 → ℝ) j

lemma SelfDualLattice.gramMatrix_eq_basisMatrix_mul_transpose (L : SelfDualLattice) :
    L.gramMatrix = L.basisMatrix * L.basisMatrix.transpose := by
  ext i j
  simp [gramMatrix, basisMatrix, Matrix.mul_apply, Matrix.transpose_apply, Matrix.of_apply]

lemma SelfDualLattice.basisMatrix_det_ne_zero (L : SelfDualLattice) : L.basisMatrix.det ≠ 0 := by
  intro hdet
  have hsq : L.basisMatrix.det * L.basisMatrix.det = 1 := by
    have h := L.gramMatrix_det
    rw [L.gramMatrix_eq_basisMatrix_mul_transpose, Matrix.det_mul, Matrix.det_transpose] at h
    simpa using h
  rw [hdet, zero_mul] at hsq
  exact zero_ne_one hsq

lemma SelfDualLattice.basis_linearIndependent (L : SelfDualLattice) :
    LinearIndependent ℝ (fun i : Fin 8 => (L.basis i : Fin 8 → ℝ)) := by
  have hli := Matrix.linearIndependent_rows_of_det_ne_zero L.basisMatrix_det_ne_zero
  convert hli using 1
  ext i
  simp [basisMatrix, Matrix.of_apply]

/-- Aus der geraden Norm folgt `⟨x, y⟩ ∈ ℤ`. -/
lemma inner_product_integer (L : SelfDualLattice) (x y : Fin 8 → ℝ)
    (hx : x ∈ L.carrier) (hy : y ∈ L.carrier) :
    ∃ m : ℤ, ∑ i, x i * y i = (m : ℝ) := by
  have h_add : x + y ∈ L.carrier := L.carrier.add_mem hx hy
  obtain ⟨nx, hx_norm⟩ := L.even_norm x hx
  obtain ⟨ny, hy_norm⟩ := L.even_norm y hy
  obtain ⟨nxy, hxy_norm⟩ := L.even_norm (x + y) h_add
  refine ⟨nxy - nx - ny, ?_⟩
  have h_pol :
      2 * (∑ i, x i * y i) = normSq (x + y) - normSq x - normSq y := by
    simp only [normSq]
    calc
      2 * ∑ i, x i * y i = ∑ i, 2 * (x i * y i) := by rw [Finset.mul_sum]
      _ = ∑ i, ((x i + y i) ^ 2 - x i ^ 2 - y i ^ 2) := by
        refine Finset.sum_congr rfl fun i _ => ?_
        ring
      _ = ∑ i, (x i + y i) ^ 2 - ∑ i, x i ^ 2 - ∑ i, y i ^ 2 := by
        rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib]
  have h2 : 2 * (∑ i, x i * y i) = 2 * ((nxy : ℝ) - nx - ny) := by
    rw [h_pol, hxy_norm, hx_norm, hy_norm]
    ring
  have hlin : ∑ i, x i * y i = (nxy : ℝ) - nx - ny := by linarith
  rw [Int.cast_sub, Int.cast_sub]
  exact hlin

lemma SelfDualLattice.gram_entry_integral (L : SelfDualLattice) (i j : Fin 8) :
    ∃ z : ℤ,
      ∑ k, (L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k = (z : ℝ) :=
  inner_product_integer L (L.basis i) (L.basis j) (SetLike.coe_mem _) (SetLike.coe_mem _)

noncomputable def SelfDualLattice.gramIntEntry (L : SelfDualLattice) (i j : Fin 8) : ℤ :=
  Classical.choose (L.gram_entry_integral i j)

lemma SelfDualLattice.gramIntEntry_spec (L : SelfDualLattice) (i j : Fin 8) :
    ∑ k, (L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k = (L.gramIntEntry i j : ℝ) :=
  Classical.choose_spec (L.gram_entry_integral i j)

/-- Ganzzahlige Gram-Matrix. -/
noncomputable def SelfDualLattice.gramInt (L : SelfDualLattice) : Matrix (Fin 8) (Fin 8) ℤ :=
  Matrix.of fun i j => L.gramIntEntry i j

theorem SelfDualLattice.gramMatrix_eq_map_gramInt (L : SelfDualLattice) :
    L.gramMatrix = L.gramInt.map (Int.castRingHom ℝ) := by
  ext i j
  have hL :
      L.gramMatrix i j =
        ∑ k, (L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k := by
    simp only [gramMatrix, Matrix.of_apply]
  have hR : (L.gramInt.map (Int.castRingHom ℝ)) i j = (L.gramIntEntry i j : ℝ) := by
    simp only [Matrix.map_apply, gramInt, Matrix.of_apply]
    rfl
  rw [hL, hR]
  exact L.gramIntEntry_spec i j

/-- Eintrag der Adjungierten als ganze Zahl. -/
lemma SelfDualLattice.adjugate_gram_intCast (L : SelfDualLattice) (i j : Fin 8) :
    L.gramMatrix.adjugate i j = (L.gramInt.adjugate i j : ℝ) := by
  have h := RingHom.map_adjugate (Int.castRingHom ℝ) L.gramInt
  have hG : (Int.castRingHom ℝ).mapMatrix L.gramInt = L.gramMatrix := by
    simpa [RingHom.mapMatrix] using L.gramMatrix_eq_map_gramInt.symm
  rw [hG] at h
  have happ := congr_fun (congr_fun h.symm i) j
  simpa [RingHom.mapMatrix, Matrix.map_apply] using happ

/-- Selbstdualität des Gitters: `Λ* = Λ`. -/
theorem lattice_self_duality (L : SelfDualLattice) :
    ∀ y : Fin 8 → ℝ,
      (∀ x ∈ L.carrier, ∃ m : ℤ, ∑ i, x i * y i = (m : ℝ)) ↔ y ∈ L.carrier := by
  intro y
  constructor
  · intro hy_dual
    have hcard : Fintype.card (Fin 8) = finrank ℝ (Fin 8 → ℝ) := by
      simp [Fintype.card_fin]
    let bReal := basisOfLinearIndependentOfCardEqFinrank L.basis_linearIndependent hcard
    have hb (i : Fin 8) : bReal i = (L.basis i : Fin 8 → ℝ) :=
      congr_fun (coe_basisOfLinearIndependentOfCardEqFinrank _ _) i
    let c : Fin 8 → ℝ := fun i => bReal.repr y i
    have hsum : ∑ j, c j • (L.basis j : Fin 8 → ℝ) = y := by
      simpa [c, hb] using Basis.sum_repr bReal y
    have hGc : L.gramMatrix.mulVec c = fun i => ∑ k, (L.basis i : Fin 8 → ℝ) k * y k := by
      ext i
      rw [Matrix.mulVec_apply_eq_sum]
      have hyk (k : Fin 8) : y k = ∑ j, c j * (L.basis j : Fin 8 → ℝ) k := by
        have := congr_fun hsum k
        simpa [Pi.smul_apply, smul_eq_mul] using this.symm
      have hG (j : Fin 8) :
          L.gramMatrix i j = ∑ k, (L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k := by
        simp [SelfDualLattice.gramMatrix, Matrix.of_apply]
      simp_rw [hyk, hG]
      calc
        ∑ j, (∑ k, (L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k) * c j =
            ∑ j, ∑ k, c j * ((L.basis i : Fin 8 → ℝ) k * (L.basis j : Fin 8 → ℝ) k) := by
          refine Finset.sum_congr rfl fun j _ => ?_
          rw [Finset.sum_mul]
          refine Finset.sum_congr rfl fun k _ => ?_
          ring
        _ = ∑ k, (L.basis i : Fin 8 → ℝ) k * ∑ j, c j * (L.basis j : Fin 8 → ℝ) k := by
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl fun k _ => ?_
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun j _ => ?_
          ring
    have hInv : L.gramMatrix.adjugate * L.gramMatrix = 1 := by
      rw [Matrix.adjugate_mul, L.gramMatrix_det, one_smul]
    have hc : c = L.gramMatrix.adjugate.mulVec (fun i => ∑ k, (L.basis i : Fin 8 → ℝ) k * y k) := by
      rw [← hGc, Matrix.mulVec_mulVec, hInv, Matrix.one_mulVec]
    have hvInt (i : Fin 8) :
        ∃ m : ℤ, ∑ k, (L.basis i : Fin 8 → ℝ) k * y k = (m : ℝ) :=
      hy_dual _ (SetLike.coe_mem (L.basis i))
    choose vInt hvInt using hvInt
    let m : Fin 8 → ℤ := fun i => ∑ j, L.gramInt.adjugate i j * vInt j
    have hc_int (i : Fin 8) : c i = (m i : ℝ) := by
      have hi := congr_fun hc i
      rw [Matrix.mulVec_apply_eq_sum] at hi
      simp_rw [L.adjugate_gram_intCast, hvInt, ← Int.cast_mul, ← Int.cast_sum] at hi
      simpa [m] using hi
    have hy_zsmul : y = ∑ j, m j • (L.basis j : Fin 8 → ℝ) := by
      rw [← hsum]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [hc_int, ← Int.cast_smul_eq_zsmul ℝ]
    rw [hy_zsmul]
    exact L.carrier.sum_smul_mem m fun j _ => SetLike.coe_mem (L.basis j)
  · intro hy_mem x hx
    exact inner_product_integer L x y hx hy_mem

/-- Quantisierung des Laplace-Spektrums: jede Gitterlänge ist von der Form `√(2n)`. -/
theorem laplacian_spectrum_quantization (L : SelfDualLattice) (k : ↥L.carrier) :
    ∃ n : ℕ, normSq k.val = 2 * (n : ℝ) := by
  obtain ⟨n_int, hn⟩ := L.even_norm k.val k.property
  have hpos : 0 ≤ normSq k.val := by
    unfold normSq
    exact Finset.sum_nonneg fun _ _ => sq_nonneg _
  rw [hn] at hpos
  have hn_nat : 0 ≤ n_int := by
    have : (0 : ℝ) ≤ n_int := by nlinarith
    exact_mod_cast this
  lift n_int to ℕ using hn_nat
  refine ⟨n_int, ?_⟩
  simpa [Int.cast_natCast] using hn

set_option linter.unusedVariables false in
/-- Eine komplexe Zahl der Form `1/2 + iγ` hat Realteil `1/2`. -/
theorem zeros_strictly_on_critical_line (L : SelfDualLattice) (s : ℂ)
    (h_eigen : ∃ gamma : ℝ, s = 1 / 2 + Complex.I * gamma) : s.re = 1 / 2 := by
  rcases h_eigen with ⟨gamma, rfl⟩
  rw [Complex.add_re, Complex.I_mul]
  simp [Complex.ofReal_im, Complex.re_ofNat]

end AnalyticNumberTheory
