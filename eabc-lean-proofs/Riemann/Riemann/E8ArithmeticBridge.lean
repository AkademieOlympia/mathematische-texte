import Mathlib
import Riemann.E8PQC
import Riemann.E8ModularForms

open scoped BigOperators ArithmeticFunction ArithmeticFunction.sigma ArithmeticFunction.zeta
  Complex LSeries.notation

noncomputable section

namespace E8ArithmeticBridge

open E8PQC ArithmeticFunction LSeries

/-!
Gate 3 beginnt hier. Gate 1 und Gate 2 bleiben in `Riemann/E8PQC.lean`.
Der Torus-Plancherel auf `P_{E8}` ist nicht die Theta-Transformation:
die braucht die Fourier-Transformation der Gauß-Funktion auf `ℝ⁸`.

Gate 3A, in dieser Reihenfolge:

* `gaussian_fourier_8` aus Mathlibs Fourier-Transformation auf euklidischen Räumen.
  Dieselbe Formel für `Re a > 0` steht in `Riemann/E8ModularForms.lean`.
* `gaussian_poisson_z8`: Produktformel für `exp(-π t ∑ n_i²)` auf `ℤ⁸`.
  Das ist die Standardnorm, nicht die Quadratnorm der Gitterpunkte von `E₈`.
* `poisson_summation_zd`: Poisson-Summe auf `ℤ^d` durch Induktion aus der 1D-Formel.
  Basis, `Fin.snoc`-Aufspaltung, die 1D-Formel bei `0`, `rapidDecay_summable`,
  `fourier_iterate_snoc`, Stetigkeit und schneller Abfall der partiellen Transformierten,
  `coordSlice_fourier_decay_2` und die gemeinsame Summierbarkeit sind bewiesen.
  `poisson_summation_z8` ist `d = 8`. Der Übergang auf `E₈` ist `e8LatticeCoord`
  zusammen mit `z8ToE8`.
* `poisson_summation_e8`: Summe über `E₈` gleich Summe der Fourier-Transformierten
  über das duale Gitter. Kovolumen 1, weil `|det B| = 1`.
* `e8_poisson_gaussian`: daraus mit `gaussian_fourier_8` und `e8_self_duality`.
  Keine Basismatrix `U` mit `Uᵀ G U = G⁻¹`.
* `e8_theta_S`: daraus `Θ(it) = t^{-4} Θ(i/t)` auf der imaginären Achse.
* `e8Theta_S`: `Θ(-1/τ) = τ^4 Θ(τ)` auf ganz `ℍ`, mit `thetaGaussian_fourier_8`.
* `e8_theta_T`: `Θ(τ+1) = Θ(τ)`, weil `‖v‖²` gerade ist.
* `e8Theta_differentiableOn` in `Riemann/E8ThetaHolomorphy.lean`:
  Holomorphie auf `{z | 0 < Im z}` über `hasDerivAt_tsum_of_isPreconnected`
  und die Gauß-Majorante der Ableitung.
* `e8Theta_limit_at_cusp` in `Riemann/E8ThetaCusp.lean`:
  `Θ(iy) → 1` für `y → +∞`, aus `‖v‖² ≥ 2` für `v ≠ 0`.
* `e8Theta_slash_action_eq` und `e8Theta_modularForm` in `Riemann/E8ThetaSlash.lean`:
  Slash-Invarianz vom Gewicht 4 auf ganz `SL(2, ℤ)` über `S` und `T`,
  und die Bündelung als `ModularForm 𝒮ℒ 4`.
* `e8_theta_eq_E4` in `Riemann/E8ThetaIdentification.lean`:
  `Θ_{E8} = E₄`, aus `dim M₄ = 1` und den beiden Spitzenwerten `1`.
* `e8Theta_qExpansion_coeff` und `e8_representation_numbers` in
  `Riemann/E8ThetaQExpansion.lean`: der Koeffizient von `q^n` ist `r_{E8}(2n)`,
  und für `n ≥ 1` folgt `r_{E8}(2n) = 240 σ₃(n)`.
* `e8_mellin_product_identity`: die Dirichlet-Faltung, bereits ohne das Gitter bewiesen
* `e8RepresentationDirichlet_eq_sigma`, `sigma_three_dirichlet_summable`,
  `e8RepresentationDirichlet_summable` und `e8RepresentationDirichlet_eq_zeta`
  in `Riemann/E8Dirichlet.lean`: Koeffizienten für jedes `s`, absolute
  Konvergenz für `Re s > 4`, und `∑ r_{E8}(2n) n^{-s} = 240 ζ(s) ζ(s-3)`.

Gate 3B isoliert danach einen einzelnen `ζ`-Faktor. G3-ZF scheidet `ζ(s-3)` nicht aus.
-/

/-- Quadrat der euklidischen Norm. -/
def e8NormSq (x : Fin 8 → ℝ) : ℝ := ∑ i, x i ^ 2

/-- Schalenzahl `r_{E8}(m) = # {v ∈ E₈ | ‖v‖² = m}`. -/
def e8_representationNumber (m : ℕ) : ℕ :=
  Nat.card {v : E8LatticeSet // e8NormSq (v : Fin 8 → ℝ) = (m : ℝ)}

/-- `Θ_{E8}(q) = ∑_{v ∈ E₈} q^{‖v‖²/2}`. -/
def e8_theta (q : ℂ) : ℂ :=
  ∑' v : E8LatticeSet, q ^ (((e8NormSq (v : Fin 8 → ℝ)) / 2 : ℝ) : ℂ)

/-- Dieselbe euklidische Norm auf `EuclideanSpace ℝ (Fin 8)`. -/
def e8ToEuclidean (x : Fin 8 → ℝ) : EuclideanSpace ℝ (Fin 8) :=
  WithLp.toLp 2 x

lemma e8NormSq_eq_norm_sq (x : Fin 8 → ℝ) :
    e8NormSq x = ‖e8ToEuclidean x‖ ^ 2 := by
  unfold e8NormSq e8ToEuclidean
  rw [PiLp.norm_sq_eq_of_L2]
  simp [Real.norm_eq_abs, sq_abs]

open UpperHalfPlane
open scoped Real

/-- Theta-Reihe auf der oberen Halbebene: `Θ(τ) = ∑_v exp(π i τ ‖v‖²)`. -/
def e8Theta (τ : ℍ) : ℂ :=
  ∑' v : E8LatticeSet,
    cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ))

/-- Gauß-Funktion `exp(-π t ‖x‖²)` auf `ℝ⁸`. -/
def e8Gaussian (t : ℝ) (x : EuclideanSpace ℝ (Fin 8)) : ℂ :=
  cexp (-((π * t : ℝ) : ℂ) * ‖x‖ ^ 2)

open scoped FourierTransform RealInnerProductSpace Asymptotics

/-- Fourier-Transformation der achtdimensionalen Gauß-Funktion.
`b = π t` in `fourier_gaussian_innerProductSpace` und `finrank = 8` ergeben den Faktor `t^{-4}`. -/
theorem gaussian_fourier_8 {t : ℝ} (ht : 0 < t) (ξ : EuclideanSpace ℝ (Fin 8)) :
    𝓕 (e8Gaussian t) ξ = (t : ℂ) ^ (-4 : ℂ) * cexp (-π * ‖ξ‖ ^ 2 / t) := by
  have hb : 0 < ((π * t : ℝ) : ℂ).re := by
    rw [Complex.ofReal_re]
    exact mul_pos Real.pi_pos ht
  unfold e8Gaussian
  rw [fourier_gaussian_innerProductSpace hb, finrank_euclideanSpace_fin]
  have hexp : ((8 : ℕ) / 2 : ℂ) = 4 := by norm_num
  have hdiv : (π : ℂ) / ((π * t : ℝ) : ℂ) = (t : ℂ)⁻¹ := by
    rw [Complex.ofReal_mul]
    field_simp
  have hargπ : (t : ℂ).arg ≠ π := by
    rw [Complex.arg_ofReal_of_nonneg ht.le]
    exact Real.pi_ne_zero.symm
  have harg : -π ^ 2 * ‖ξ‖ ^ 2 / ((π * t : ℝ) : ℂ) = -π * ‖ξ‖ ^ 2 / t := by
    rw [Complex.ofReal_mul]
    field_simp
  rw [hexp, hdiv, harg]
  congr 1
  rw [Complex.inv_cpow (t : ℂ) 4 hargπ, ← Complex.cpow_neg]

/-- `B^{-T} ℤ⁸` trifft genau `E₈`. Das ist die Bijektion `z8ToE8` aus Gate 1. -/
lemma e8_dual_reindex (x : Fin 8 → ℝ) :
    x ∈ E8LatticeSet ↔ ∃ m : Fin 8 → ℤ, x = z8ToE8Fun m := by
  constructor
  · intro hx
    refine ⟨e8FrequencyEquiv ⟨x, hx⟩, ?_⟩
    have h := e8FrequencyEquiv.symm_apply_apply ⟨x, hx⟩
    rw [e8FrequencyEquiv_symm_apply] at h
    exact (congrArg Subtype.val h).symm
  · rintro ⟨m, rfl⟩
    exact (z8ToE8 m).property

open MeasureTheory

/-- Die Basis `B` als lineare Abbildung des euklidischen Raums. -/
def e8CoordLinear : EuclideanSpace ℝ (Fin 8) →ₗ[ℝ] EuclideanSpace ℝ (Fin 8) :=
  (WithLp.linearEquiv 2 ℝ (Fin 8 → ℝ)).symm.toLinearMap ∘ₗ
    e8CoordMap ∘ₗ (WithLp.linearEquiv 2 ℝ (Fin 8 → ℝ)).toLinearMap

/-- Dieselbe Abbildung als Funktion. -/
def e8CoordMapOnEuclidean (c : EuclideanSpace ℝ (Fin 8)) : EuclideanSpace ℝ (Fin 8) :=
  e8CoordLinear c

lemma e8CoordLinear_det : LinearMap.det e8CoordLinear = -1 := by
  simpa [e8CoordLinear, e8CoordMap_det] using
    LinearMap.det_conj e8CoordMap (WithLp.linearEquiv 2 ℝ (Fin 8 → ℝ)).symm

lemma e8CoordLinear_measurePreserving :
    Measure.map (e8CoordLinear : EuclideanSpace ℝ (Fin 8) → EuclideanSpace ℝ (Fin 8)) volume =
      volume := by
  rw [Measure.map_linearMap_addHaar_eq_smul_addHaar volume (f := e8CoordLinear)
      (by rw [e8CoordLinear_det]; norm_num), e8CoordLinear_det]
  norm_num

/-- `g ∘ B` ist integrierbar genau dann, wenn `g` es ist: `|det B| = 1`. -/
lemma e8CoordLinear_integrable_comp (g : EuclideanSpace ℝ (Fin 8) → ℂ) :
    Integrable (g ∘ e8CoordMapOnEuclidean) ↔ Integrable g := by
  let f : EuclideanSpace ℝ (Fin 8) →L[ℝ] EuclideanSpace ℝ (Fin 8) :=
    e8CoordLinear.toContinuousLinearMap
  have hf : f.det ≠ 0 := by
    rw [LinearMap.det_toContinuousLinearMap, e8CoordLinear_det]
    norm_num
  let e := f.toContinuousLinearEquivOfDetNeZero hf
  have hfun : (e : EuclideanSpace ℝ (Fin 8) → EuclideanSpace ℝ (Fin 8)) = e8CoordMapOnEuclidean := by
    ext c
    simp [e, f, e8CoordMapOnEuclidean,
      ContinuousLinearMap.toContinuousLinearEquivOfDetNeZero_apply]
  rw [← hfun]
  have hembed : MeasurableEmbedding (e : EuclideanSpace ℝ (Fin 8) → EuclideanSpace ℝ (Fin 8)) :=
    e.toHomeomorph.toMeasurableEquiv.measurableEmbedding
  have hmap : Measure.map (e : EuclideanSpace ℝ (Fin 8) → EuclideanSpace ℝ (Fin 8)) volume =
      volume := by
    rw [hfun]
    exact e8CoordLinear_measurePreserving
  refine (hembed.integrable_map_iff (μ := volume) (g := g)).symm.trans ?_
  rw [hmap]

/-- `⟪c, m⟫ = ⟪B c, B^{-T} m⟫`. Die rechte Frequenz ist `z8ToE8Fun m`. -/
lemma e8_fourier_phase (m : Fin 8 → ℤ) (c : EuclideanSpace ℝ (Fin 8)) :
    ⟪c, e8ToEuclidean fun i => (m i : ℝ)⟫ =
      ⟪e8CoordMapOnEuclidean c, e8ToEuclidean (z8ToE8Fun m)⟫ := by
  rw [EuclideanSpace.inner_eq_star_dotProduct, EuclideanSpace.inner_eq_star_dotProduct]
  simp only [e8ToEuclidean, WithLp.ofLp_toLp, star_trivial, dotProduct]
  unfold e8CoordMapOnEuclidean e8CoordLinear
  rw [LinearMap.comp_apply, LinearMap.comp_apply, LinearEquiv.coe_toLinearMap,
    LinearEquiv.coe_toLinearMap, WithLp.coe_linearEquiv, WithLp.coe_symm_linearEquiv,
    WithLp.ofLp_toLp, e8CoordMap_apply]
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hpull : ∑ x, z8ToE8Fun m x * (c.ofLp j * e8_basis j x) =
      c.ofLp j * ∑ x, e8_basis j x * z8ToE8Fun m x := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun x _ => ?_
    ring
  rw [hpull, z8ToE8_pairing, mul_comm]

/-- Die Reihe `∑_{n∈ℤ} exp(-a n²)` konvergiert für `a > 0`. -/
lemma summable_gaussian_coord {a : ℝ} (ha : 0 < a) :
    Summable fun n : ℤ => Real.exp (-a * (n : ℝ) ^ 2) := by
  have hc : -a < 0 := neg_lt_zero.mpr ha
  have hsq (i : ℕ) : (i : ℝ) ≤ (i : ℝ) ^ 2 := by
    cases i with
    | zero => norm_num
    | succ i =>
      have h1 : (1 : ℝ) ≤ (i.succ : ℝ) := by exact_mod_cast Nat.succ_le_succ (Nat.zero_le i)
      nlinarith
  have hnat : Summable fun n : ℕ => Real.exp (-a * (n : ℝ) ^ 2) := by
    simpa [mul_assoc, neg_mul] using Real.summable_exp_nat_mul_of_ge hc hsq
  apply Summable.of_nat_of_neg
  · exact hnat
  · simpa [Int.cast_neg, neg_sq] using hnat

/-- Summe eines koordinatenweisen Produkts über `Fin n → ℤ`. -/
lemma tsum_coord_prod {n : ℕ} (f : Fin n → ℤ → ℝ)
    (hf : ∀ i, Summable (f i)) (h0 : ∀ i, 0 ≤ f i) :
    Summable (fun p : Fin n → ℤ => ∏ i, f i (p i)) ∧
      ∑' p : Fin n → ℤ, ∏ i, f i (p i) = ∏ i, ∑' k : ℤ, f i k := by
  induction n with
  | zero =>
    refine ⟨Summable.of_finite, ?_⟩
    rw [tsum_fintype, Finset.univ_unique, Finset.sum_singleton, Fin.prod_univ_zero,
      Fin.prod_univ_zero]
  | succ n ih =>
    obtain ⟨hs, heq⟩ := ih (fun i => f i.succ) (fun i => hf i.succ) (fun i => h0 i.succ)
    let e : ℤ × (Fin n → ℤ) ≃ (Fin (n + 1) → ℤ) := Fin.consEquiv fun _ => ℤ
    have hfun (a : ℤ) (q : Fin n → ℤ) :
        (∏ i, f i (e (a, q) i)) = f 0 a * ∏ i, f i.succ (q i) := by
      have hcons : e (a, q) = Fin.cons a q := by
        unfold e
        rfl
      rw [hcons, Fin.prod_univ_succ]
      simp only [Fin.cons_zero, Fin.cons_succ]
    have hg0 : 0 ≤ fun q : Fin n → ℤ => ∏ i, f i.succ (q i) :=
      fun q => Finset.prod_nonneg fun i _ => h0 i.succ (q i)
    have hmul : Summable fun p : ℤ × (Fin n → ℤ) =>
        f 0 p.1 * ∏ i, f i.succ (p.2 i) := by
      refine (summable_prod_of_nonneg fun p => mul_nonneg (h0 0 p.1) (hg0 p.2)).2 ⟨?_, ?_⟩
      · intro a
        exact (hs.mul_left (f 0 a))
      · have hpull (a : ℤ) :
            ∑' q : Fin n → ℤ, f 0 a * ∏ i : Fin n, f i.succ (q i) =
              f 0 a * ∑' q : Fin n → ℤ, ∏ i : Fin n, f i.succ (q i) :=
          hs.tsum_mul_left (f 0 a)
        simp_rw [hpull]
        exact (hf 0).mul_right _
    have hcomp : Summable fun p : ℤ × (Fin n → ℤ) => ∏ i, f i (e p i) :=
      hmul.congr fun p => (hfun p.1 p.2).symm
    refine ⟨e.summable_iff.mp hcomp, ?_⟩
    rw [← e.tsum_eq]
    simp_rw [hfun]
    rw [← (hf 0).tsum_mul_tsum hs hmul, heq, Fin.prod_univ_succ]

/-- Eindimensionale Gauß-Poisson-Summation: Mathlib `Real.tsum_exp_neg_mul_int_sq`. -/
lemma gaussian_poisson_1d {t : ℝ} (ht : 0 < t) :
    ∑' n : ℤ, Real.exp (-Real.pi * t * (n : ℝ) ^ 2) =
      t ^ (-1 / 2 : ℝ) * ∑' m : ℤ, Real.exp (-(Real.pi / t) * (m : ℝ) ^ 2) := by
  rw [Real.tsum_exp_neg_mul_int_sq ht]
  have hcoeff : (1 : ℝ) / t ^ (1 / 2 : ℝ) = t ^ (-1 / 2 : ℝ) := by
    rw [one_div, ← Real.rpow_neg ht.le]
    congr 1
    ring
  rw [hcoeff]
  congr 1
  refine tsum_congr fun m => ?_
  congr 1
  ring

/-- Achtfache Produktformel auf `ℤ⁸` für die Standardnorm `∑ n_i²`.
Das ist nicht die Theta-Reihe von `E₈`: dort steht `‖B n‖² = nᵀ G n`. -/
theorem gaussian_poisson_z8 {t : ℝ} (ht : 0 < t) :
    ∑' n : Fin 8 → ℤ, Real.exp (-Real.pi * t * ∑ i, (n i : ℝ) ^ 2) =
      t ^ (-4 : ℝ) *
        ∑' m : Fin 8 → ℤ, Real.exp (-(Real.pi / t) * ∑ i, (m i : ℝ) ^ 2) := by
  have hsplit (c : ℝ) (n : Fin 8 → ℤ) :
      Real.exp (-c * ∑ i, (n i : ℝ) ^ 2) = ∏ i, Real.exp (-c * (n i : ℝ) ^ 2) := by
    rw [Finset.mul_sum, Real.exp_sum]
  have hleft (n : Fin 8 → ℤ) :
      Real.exp (-Real.pi * t * ∑ i, (n i : ℝ) ^ 2) =
        ∏ i, Real.exp (-Real.pi * t * (n i : ℝ) ^ 2) := by
    simpa [mul_assoc] using hsplit (Real.pi * t) n
  have hright (m : Fin 8 → ℤ) :
      Real.exp (-(Real.pi / t) * ∑ i, (m i : ℝ) ^ 2) =
        ∏ i, Real.exp (-(Real.pi / t) * (m i : ℝ) ^ 2) := by
    simpa using hsplit (Real.pi / t) m
  simp_rw [hleft]
  let φ : ℤ → ℝ := fun k => Real.exp (-Real.pi * t * (k : ℝ) ^ 2)
  let ψ : ℤ → ℝ := fun k => Real.exp (-(Real.pi / t) * (k : ℝ) ^ 2)
  have hφ : Summable φ := by
    simpa [φ, mul_assoc, neg_mul] using summable_gaussian_coord (mul_pos Real.pi_pos ht)
  have hψ : Summable ψ := by
    simpa [ψ] using summable_gaussian_coord (div_pos Real.pi_pos ht)
  have hφ0 : 0 ≤ φ := fun _ => Real.exp_nonneg _
  have hψ0 : 0 ≤ ψ := fun _ => Real.exp_nonneg _
  obtain ⟨-, hprod⟩ := tsum_coord_prod (fun _ : Fin 8 => φ) (fun _ => hφ) (fun _ => hφ0)
  rw [hprod, Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  have h1 : ∑' k : ℤ, φ k = t ^ (-1 / 2 : ℝ) * ∑' k : ℤ, ψ k := gaussian_poisson_1d ht
  rw [h1]
  obtain ⟨-, hψprod⟩ := tsum_coord_prod (fun _ : Fin 8 => ψ) (fun _ => hψ) (fun _ => hψ0)
  have hψeq : (∑' k : ℤ, ψ k) ^ 8 =
      ∑' m : Fin 8 → ℤ, Real.exp (-(Real.pi / t) * ∑ i, (m i : ℝ) ^ 2) := by
    have h8 : (∑' k : ℤ, ψ k) ^ 8 = ∏ i : Fin 8, ∑' k : ℤ, ψ k := by
      have hc : (Finset.univ : Finset (Fin 8)).card = 8 := by simp
      rw [← hc]
      exact (Finset.prod_const (∑' k : ℤ, ψ k)).symm
    rw [h8, ← hψprod]
    refine tsum_congr fun m => ?_
    exact (hright m).symm
  rw [mul_pow, hψeq]
  congr 1
  rw [← Real.rpow_natCast, ← Real.rpow_mul ht.le]
  norm_num

/-- Ein Gitterpunkt von `ℤ^d` im euklidischen Raum. -/
def zdPoint {d : ℕ} (n : Fin d → ℤ) : EuclideanSpace ℝ (Fin d) :=
  WithLp.toLp 2 fun i => (n i : ℝ)

/-- Schneller Abfall: zu jeder Potenz `k` ist `F` am Unendlichen `O(‖x‖^{-k})`. -/
def RapidDecay {d : ℕ} (F : EuclideanSpace ℝ (Fin d) → ℂ) : Prop :=
  ∀ k : ℕ, F =O[Filter.cocompact (EuclideanSpace ℝ (Fin d))] fun x => ‖x‖ ^ (-(k : ℝ))

/-- Letzte Koordinate abtrennen: `(x', t) ↦ (x', t)` in `ℝ^{d+1}`. -/
def snocPoint {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) :
    EuclideanSpace ℝ (Fin (d + 1)) :=
  WithLp.toLp 2 (Fin.snoc (WithLp.ofLp x') t)

/-- Schnitt von `F` bei festem `x'` längs der letzten Koordinate. -/
def coordSlice {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) : ℂ :=
  F (snocPoint x' t)

lemma snocPoint_norm_sq {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) :
    ‖snocPoint x' t‖ ^ 2 = ‖x'‖ ^ 2 + t ^ 2 := by
  rw [PiLp.norm_sq_eq_of_L2, PiLp.norm_sq_eq_of_L2]
  unfold snocPoint
  simp only [Real.norm_eq_abs, sq_abs]
  rw [Fin.sum_univ_castSucc]
  simp only [Fin.snoc_castSucc, Fin.snoc_last]

lemma snocPoint_norm_ge {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) :
    |t| ≤ ‖snocPoint x' t‖ := by
  have hsq : t ^ 2 ≤ ‖snocPoint x' t‖ ^ 2 := by
    rw [snocPoint_norm_sq]
    nlinarith [sq_nonneg ‖x'‖]
  exact abs_le_of_sq_le_sq hsq (norm_nonneg _)

lemma continuous_snocPoint {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) :
    Continuous (snocPoint x') := by
  unfold snocPoint
  have hlin : Continuous (⇑(WithLp.linearEquiv 2 ℝ (Fin (d + 1) → ℝ)).symm) :=
    (WithLp.linearEquiv 2 ℝ (Fin (d + 1) → ℝ)).symm.toLinearMap.continuous_of_finiteDimensional
  let x0 : Fin d → ℝ := WithLp.ofLp x'
  have hpi : Continuous fun t : ℝ => (Fin.snoc x0 t : Fin (d + 1) → ℝ) := by
    refine continuous_pi fun i => ?_
    induction i using Fin.lastCases with
    | last =>
      simp only [Fin.snoc_last]
      exact continuous_id
    | cast _ => simpa [Fin.snoc_castSucc, x0] using continuous_const
  convert hlin.comp hpi using 2 with t
  simp [x0, WithLp.coe_symm_linearEquiv]

lemma coordSlice_continuous {d : ℕ} {F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ}
    (hF : Continuous F) (x' : EuclideanSpace ℝ (Fin d)) :
    Continuous (coordSlice F x') :=
  hF.comp (continuous_snocPoint x')

lemma tendsto_snocPoint_cocompact {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) :
    Filter.Tendsto (snocPoint x') (Filter.cocompact ℝ)
      (Filter.cocompact (EuclideanSpace ℝ (Fin (d + 1)))) := by
  refine tendsto_cocompact_of_tendsto_dist_comp_atTop (0 : EuclideanSpace ℝ (Fin (d + 1))) ?_
  simp_rw [dist_zero_right]
  exact Filter.tendsto_atTop_mono (snocPoint_norm_ge x')
    (tendsto_norm_cocompact_atTop (E := ℝ))

/-- Der räumliche Schnitt fällt wie `|t|^{-2}`. Das genügt für die 1D-Poisson-Formel. -/
lemma coordSlice_decay_2 {d : ℕ} {F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ}
    (hdecay : RapidDecay F) (x' : EuclideanSpace ℝ (Fin d)) :
    coordSlice F x' =O[Filter.cocompact ℝ] (|·| ^ (-2 : ℝ)) := by
  have hcomp := (hdecay 2).comp_tendsto (tendsto_snocPoint_cocompact x')
  refine hcomp.trans ?_
  rw [Asymptotics.isBigO_iff]
  refine ⟨1, ?_⟩
  filter_upwards [(tendsto_norm_cocompact_atTop (E := ℝ)).eventually
      (Filter.eventually_ge_atTop (1 : ℝ))] with t ht
  simp only [Function.comp_apply, Real.norm_eq_abs,
    abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _), one_mul,
    abs_of_nonneg (Real.rpow_nonneg (abs_nonneg _) _)]
  have ht0 : (0 : ℝ) < |t| := lt_of_lt_of_le zero_lt_one ht
  exact Real.rpow_le_rpow_of_nonpos ht0 (snocPoint_norm_ge x' t) (by norm_num)

lemma zdPoint_snoc {d : ℕ} (n' : Fin d → ℤ) (k : ℤ) :
    zdPoint (Fin.snoc n' k) = snocPoint (zdPoint n') (k : ℝ) := by
  apply (WithLp.linearEquiv 2 ℝ (Fin (d + 1) → ℝ)).injective
  ext i
  simp only [zdPoint, snocPoint, WithLp.coe_linearEquiv, WithLp.ofLp_toLp]
  induction i using Fin.lastCases with
  | last => simp [Fin.snoc_last]
  | cast _ => simp [Fin.snoc_castSucc]

lemma zd_snoc_summable {d : ℕ} {g : (Fin (d + 1) → ℤ) → ℂ} (hg : Summable g) :
    Summable fun p : (Fin d → ℤ) × ℤ => g (Fin.snoc p.1 p.2) := by
  let e : ℤ × (Fin d → ℤ) ≃ (Fin (d + 1) → ℤ) := Fin.snocEquiv fun _ => ℤ
  have hs : Summable fun p : ℤ × (Fin d → ℤ) => g (e p) := e.summable_iff.mpr hg
  have hsnoc (k : ℤ) (n' : Fin d → ℤ) : e (k, n') = Fin.snoc n' k := by
    ext i
    simp [e, Fin.snocEquiv_apply]
  have hswap : Summable fun p : ℤ × (Fin d → ℤ) => g (Fin.snoc p.2 p.1) := by
    refine hs.congr fun p => ?_
    simp [hsnoc]
  exact hswap.prod_symm.congr fun p => by simp [Prod.swap]

/-- Aufspaltung der Summe längs der letzten Koordinate.
`Fin.snocEquiv` ist `ℤ × (Fin d → ℤ) ≃ Fin (d + 1) → ℤ`. -/
lemma zd_tsum_snoc {d : ℕ} {g : (Fin (d + 1) → ℤ) → ℂ} (hg : Summable g) :
    ∑' n : Fin (d + 1) → ℤ, g n =
      ∑' (n' : Fin d → ℤ) (k : ℤ), g (Fin.snoc n' k) := by
  let e : ℤ × (Fin d → ℤ) ≃ (Fin (d + 1) → ℤ) := Fin.snocEquiv fun _ => ℤ
  have hs : Summable fun p : ℤ × (Fin d → ℤ) => g (e p) := e.summable_iff.mpr hg
  have hsnoc (k : ℤ) (n' : Fin d → ℤ) : e (k, n') = Fin.snoc n' k := by
    ext i
    simp [e, Fin.snocEquiv_apply]
  have hprod : Summable fun p : (Fin d → ℤ) × ℤ => g (Fin.snoc p.1 p.2) := by
    have hswap : Summable fun p : ℤ × (Fin d → ℤ) => g (Fin.snoc p.2 p.1) := by
      simpa [hsnoc] using hs
    simpa [Prod.swap] using hswap.prod_symm
  calc
    ∑' n, g n = ∑' p, g (e p) := (e.tsum_eq g).symm
    _ = ∑' (k : ℤ) (n' : Fin d → ℤ), g (e (k, n')) := hs.tsum_prod
    _ = ∑' (k : ℤ) (n' : Fin d → ℤ), g (Fin.snoc n' k) := by simp_rw [hsnoc]
    _ = ∑' (n' : Fin d → ℤ) (k : ℤ), g (Fin.snoc n' k) := hprod.tsum_comm

lemma zd_tsum_snoc_outer {d : ℕ} {g : (Fin (d + 1) → ℤ) → ℂ} (hg : Summable g) :
    ∑' n : Fin (d + 1) → ℤ, g n =
      ∑' (k : ℤ) (n' : Fin d → ℤ), g (Fin.snoc n' k) := by
  rw [zd_tsum_snoc hg]
  exact (zd_snoc_summable hg).tsum_comm.symm

/-- Poisson-Formel auf `ℤ` bei `0`. Der Kreisfaktor `fourier n 0` ist `1`. -/
lemma poisson_1d_at_zero {f : ℝ → ℂ} (hf : Continuous f)
    (hdecay : f =O[Filter.cocompact ℝ] (|·| ^ (-2 : ℝ)))
    (hFdecay : (𝓕 f) =O[Filter.cocompact ℝ] (|·| ^ (-2 : ℝ))) :
    ∑' n : ℤ, f (n : ℝ) = ∑' m : ℤ, (𝓕 f) (m : ℝ) := by
  have h := Real.tsum_eq_tsum_fourier_of_rpow_decay hf one_lt_two hdecay hFdecay 0
  simp only [zero_add] at h
  refine h.trans (tsum_congr fun n => ?_)
  rw [show (fourier n) ((0 : ℝ) : UnitAddCircle) = 1 by
    simp, mul_one]

/-- Eindimensionale Vergleichsreihe `(1 + |n|)^{-2}`. -/
lemma summable_one_add_abs_inv_sq :
    Summable fun n : ℤ => (1 + |(n : ℝ)|) ^ (-(2 : ℝ)) := by
  refine (Real.summable_abs_int_rpow one_lt_two).of_norm_bounded_eventually ?_
  filter_upwards [Filter.eventually_cofinite_ne (0 : ℤ)] with n hn
  have hpos : (0 : ℝ) < |(n : ℝ)| := abs_pos.mpr (by simpa [Int.cast_eq_zero] using hn)
  rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (by positivity) _)]
  exact Real.rpow_le_rpow_of_nonpos hpos (le_add_of_nonneg_left zero_le_one) (by norm_num)

/-- Produkt der eindimensionalen Reihen über `Fin d → ℤ`. -/
lemma summable_coord_decay (d : ℕ) :
    Summable fun n : Fin d → ℤ => ∏ i, (1 + |(n i : ℝ)|) ^ (-(2 : ℝ)) :=
  (tsum_coord_prod (fun _ : Fin d => fun k : ℤ => (1 + |(k : ℝ)|) ^ (-(2 : ℝ)))
    (fun _ => summable_one_add_abs_inv_sq)
    (fun _ _ => Real.rpow_nonneg (by positivity) _)).1

/-- In jeder Kugel liegen nur endlich viele Punkte von `ℤ^d`. -/
lemma finite_zdPoint_norm_le {d : ℕ} (R : ℝ) :
    {n : Fin d → ℤ | ‖zdPoint n‖ ≤ R}.Finite := by
  classical
  by_cases hR : 0 ≤ R
  · refine (Set.Finite.pi fun _ : Fin d => Set.finite_Icc (⌈-R⌉ : ℤ) ⌊R⌋).subset ?_
    intro n hn
    rw [Set.mem_univ_pi]
    intro i
    have habs : |(n i : ℝ)| ≤ R := by
      have hle : |(n i : ℝ)| ≤ ‖zdPoint n‖ := by
        simpa [zdPoint, PiLp.toLp_apply, Real.norm_eq_abs] using
          PiLp.norm_apply_le (zdPoint n) i
      exact hle.trans hn
    rw [abs_le] at habs
    exact ⟨Int.ceil_le.mpr habs.1, Int.le_floor.mpr habs.2⟩
  · refine Set.finite_empty.subset ?_
    intro n hn
    exact absurd (hn.trans' (norm_nonneg _)) (not_le.mpr (lt_of_not_ge hR))

lemma tendsto_zdPoint_norm_atTop {d : ℕ} :
    Filter.Tendsto (fun n : Fin d → ℤ => ‖zdPoint n‖) Filter.cofinite Filter.atTop := by
  rw [Filter.tendsto_atTop]
  intro R
  rw [Filter.eventually_cofinite]
  refine (finite_zdPoint_norm_le R).subset ?_
  intro n hn
  exact le_of_lt (lt_of_not_ge hn)

lemma tendsto_zdPoint_cocompact {d : ℕ} :
    Filter.Tendsto zdPoint Filter.cofinite (Filter.cocompact (EuclideanSpace ℝ (Fin d))) :=
  tendsto_cocompact_of_tendsto_dist_comp_atTop (0 : EuclideanSpace ℝ (Fin d)) <| by
    simp_rw [dist_zero_right]
    exact tendsto_zdPoint_norm_atTop

/-- Für `‖x‖ ≥ 1` majorisiert das Koordinatenprodukt die euklidische Potenz `‖x‖^{-2d}`. -/
lemma zdPoint_norm_le_coordProd {d : ℕ} (n : Fin d → ℤ) (hn : 1 ≤ ‖zdPoint n‖) :
    ‖zdPoint n‖ ^ (-(2 * d : ℝ)) ≤
      (2 : ℝ) ^ (2 * (d : ℝ)) * ∏ i, (1 + |(n i : ℝ)|) ^ (-(2 : ℝ)) := by
  set N : ℝ := ‖zdPoint n‖
  have hN : (0 : ℝ) < N := lt_of_lt_of_le zero_lt_one hn
  have hcoord (i : Fin d) : 1 + |(n i : ℝ)| ≤ 2 * N := by
    have habs : |(n i : ℝ)| ≤ N := by
      simpa [zdPoint, PiLp.toLp_apply, Real.norm_eq_abs, N] using
        PiLp.norm_apply_le (zdPoint n) i
    linarith
  have hpos (i : Fin d) : (0 : ℝ) ≤ 1 + |(n i : ℝ)| := by positivity
  have hprod_le : ∏ i, (1 + |(n i : ℝ)|) ≤ (2 * N) ^ d := by
    calc
      ∏ i, (1 + |(n i : ℝ)|) ≤ ∏ i : Fin d, 2 * N :=
        Finset.prod_le_prod₀ (fun i _ => hpos i) (fun i _ => hcoord i)
      _ = (2 * N) ^ d := by rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  have hprod_pos : (0 : ℝ) < ∏ i, (1 + |(n i : ℝ)|) :=
    Finset.prod_pos fun i _ => by positivity
  have hexp : -(2 * (d : ℝ)) = (d : ℝ) * (-2) := by ring
  calc
    N ^ (-(2 * (d : ℝ)))
        = (2 : ℝ) ^ (2 * (d : ℝ)) *
            ((2 : ℝ) ^ (-(2 * (d : ℝ))) * N ^ (-(2 * (d : ℝ)))) := by
          have h2 : (2 : ℝ) ^ (2 * (d : ℝ)) * (2 : ℝ) ^ (-(2 * (d : ℝ))) = 1 := by
            rw [← Real.rpow_add (by norm_num : (0 : ℝ) < 2)]
            simp
          calc
            N ^ (-(2 * (d : ℝ))) = 1 * N ^ (-(2 * (d : ℝ))) := by rw [one_mul]
            _ = ((2 : ℝ) ^ (2 * (d : ℝ)) * (2 : ℝ) ^ (-(2 * (d : ℝ)))) *
                  N ^ (-(2 * (d : ℝ))) := by rw [← h2]
            _ = (2 : ℝ) ^ (2 * (d : ℝ)) *
                  ((2 : ℝ) ^ (-(2 * (d : ℝ))) * N ^ (-(2 * (d : ℝ)))) := by rw [mul_assoc]
      _ = (2 : ℝ) ^ (2 * (d : ℝ)) * ((2 * N) ^ (-(2 * (d : ℝ)))) := by
          rw [Real.mul_rpow (by norm_num) (le_of_lt hN)]
      _ = (2 : ℝ) ^ (2 * (d : ℝ)) * (((2 * N) ^ d) ^ (-(2 : ℝ))) := by
          rw [hexp, Real.rpow_mul (by positivity : (0 : ℝ) ≤ 2 * N), Real.rpow_natCast]
      _ ≤ (2 : ℝ) ^ (2 * (d : ℝ)) *
            ((∏ i, (1 + |(n i : ℝ)|)) ^ (-(2 : ℝ))) := by
          refine mul_le_mul_of_nonneg_left ?_ (by positivity)
          exact Real.rpow_le_rpow_of_nonpos hprod_pos hprod_le (by norm_num)
      _ = (2 : ℝ) ^ (2 * (d : ℝ)) * ∏ i, (1 + |(n i : ℝ)|) ^ (-(2 : ℝ)) := by
          rw [← Real.finsetProd_rpow Finset.univ (fun i => 1 + |(n i : ℝ)|) (fun i _ => hpos i)]

/-- Schneller Abfall jeder Potenz macht die Gittersumme absolut konvergent.
Die Majorante ist `∏ (1 + |n_i|)^{-2}`; die euklidische Potenz `‖n‖^{-2d}` liegt darunter. -/
lemma rapidDecay_summable {d : ℕ} {F : EuclideanSpace ℝ (Fin d) → ℂ} (hF : RapidDecay F) :
    Summable fun n : Fin d → ℤ => F (zdPoint n) := by
  cases d with
  | zero => exact Summable.of_finite
  | succ d =>
    let k : ℕ := 2 * (d + 1)
    have hk : (k : ℝ) = 2 * ((d + 1 : ℕ) : ℝ) := by simp [k, Nat.cast_mul]
    have hO := (hF k).comp_tendsto tendsto_zdPoint_cocompact
    rw [Asymptotics.isBigO_iff] at hO
    obtain ⟨C, hC⟩ := hO
    have hg : Summable fun n : Fin (d + 1) → ℤ =>
        max C 0 * (2 : ℝ) ^ (k : ℝ) * ∏ i, (1 + |(n i : ℝ)|) ^ (-(2 : ℝ)) :=
      (summable_coord_decay (d + 1)).mul_left _
    refine hg.of_norm_bounded_eventually ?_
    filter_upwards [hC,
      tendsto_zdPoint_norm_atTop.eventually (Filter.eventually_ge_atTop (1 : ℝ))] with n hn hnorm
    have hrpow : 0 ≤ ‖zdPoint n‖ ^ (-(k : ℝ)) := Real.rpow_nonneg (norm_nonneg _) _
    calc
      ‖F (zdPoint n)‖ ≤ C * ‖‖zdPoint n‖ ^ (-(k : ℝ))‖ := hn
      _ = C * ‖zdPoint n‖ ^ (-(k : ℝ)) := by
          rw [Real.norm_eq_abs, abs_of_nonneg hrpow]
      _ ≤ max C 0 * ‖zdPoint n‖ ^ (-(k : ℝ)) :=
          mul_le_mul_of_nonneg_right (le_max_left _ _) hrpow
      _ ≤ max C 0 * ((2 : ℝ) ^ (k : ℝ) * ∏ i, (1 + |(n i : ℝ)|) ^ (-(2 : ℝ))) := by
          refine mul_le_mul_of_nonneg_left ?_ (le_max_right _ _)
          rw [hk]
          simpa using zdPoint_norm_le_coordProd n hnorm
      _ = max C 0 * (2 : ℝ) ^ (k : ℝ) * ∏ i, (1 + |(n i : ℝ)|) ^ (-(2 : ℝ)) := by
          ring

/-- `Fin.snoc` auf den Koordinatenfunktionen. -/
def snocPiLinear (d : ℕ) : (Fin d → ℝ) × ℝ →ₗ[ℝ] (Fin (d + 1) → ℝ) where
  toFun p := Fin.snoc p.1 p.2
  map_add' p q := by
    ext i
    induction i using Fin.lastCases with
    | last => simp [Fin.snoc_last]
    | cast i => simp [Fin.snoc_castSucc]
  map_smul' c p := by
    ext i
    induction i using Fin.lastCases with
    | last => simp [Fin.snoc_last]
    | cast i => simp [Fin.snoc_castSucc]

def snocPiInv (d : ℕ) : (Fin (d + 1) → ℝ) →ₗ[ℝ] (Fin d → ℝ) × ℝ where
  toFun x := (fun i => x i.castSucc, x (Fin.last d))
  map_add' x y := by ext <;> simp
  map_smul' c x := by ext <;> simp

def snocPiLinearEquiv (d : ℕ) : ((Fin d → ℝ) × ℝ) ≃ₗ[ℝ] (Fin (d + 1) → ℝ) :=
  LinearEquiv.ofLinearMap (snocPiLinear d) (snocPiInv d)
    (LinearMap.ext fun x => by
      ext i
      induction i using Fin.lastCases with
      | last => simp [snocPiLinear, snocPiInv, Fin.snoc_last]
      | cast j => simp [snocPiLinear, snocPiInv, Fin.snoc_castSucc])
    (LinearMap.ext fun p => by
      ext i
      · simp [snocPiLinear, snocPiInv, Fin.snoc_castSucc]
      · simp [snocPiLinear, snocPiInv, Fin.snoc_last])

/-- `ℝ^d × ℝ` mit der `L²`-Norm ist linear isometrisch zu `ℝ^{d+1}`. -/
noncomputable def snocLinearEquiv (d : ℕ) :
    WithLp 2 (EuclideanSpace ℝ (Fin d) × ℝ) ≃ₗ[ℝ] EuclideanSpace ℝ (Fin (d + 1)) :=
  (WithLp.linearEquiv 2 ℝ (EuclideanSpace ℝ (Fin d) × ℝ)).trans <|
    (LinearEquiv.prodCongr (WithLp.linearEquiv 2 ℝ (Fin d → ℝ)) (LinearEquiv.refl ℝ ℝ)).trans <|
      (snocPiLinearEquiv d).trans (WithLp.linearEquiv 2 ℝ (Fin (d + 1) → ℝ)).symm

lemma snocLinearEquiv_apply {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) :
    snocLinearEquiv d (WithLp.toLp 2 (x', t)) = snocPoint x' t := by
  simp [snocLinearEquiv, snocPiLinearEquiv, snocPiLinear, snocPoint,
    LinearEquiv.prodCongr_apply]

noncomputable def snocLinearIsometry (d : ℕ) :
    WithLp 2 (EuclideanSpace ℝ (Fin d) × ℝ) ≃ₗᵢ[ℝ] EuclideanSpace ℝ (Fin (d + 1)) where
  toLinearEquiv := snocLinearEquiv d
  norm_map' x := by
    rw [← sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)]
    have happly : snocLinearEquiv d x = snocPoint (WithLp.fst x) (WithLp.snd x) := by
      have hx : x = WithLp.toLp 2 (WithLp.fst x, WithLp.snd x) := by
        simp only [WithLp.fst, WithLp.snd]
      rw [hx]
      exact snocLinearEquiv_apply _ _
    rw [happly, snocPoint_norm_sq, WithLp.prod_norm_sq_eq_of_L2, Real.norm_eq_abs, sq_abs]

/-- Maßerhaltung: `snocPoint` schiebt `volume.prod volume` auf `volume`. -/
lemma measurePreserving_snocPoint (d : ℕ) :
    MeasurePreserving (fun p : EuclideanSpace ℝ (Fin d) × ℝ => snocPoint p.1 p.2)
      (volume.prod volume) := by
  have hcomp := (snocLinearIsometry d).measurePreserving.comp
    (WithLp.volume_preserving_toLp (EuclideanSpace ℝ (Fin d)) ℝ)
  have heq : (fun p : EuclideanSpace ℝ (Fin d) × ℝ => snocPoint p.1 p.2) =
      (snocLinearIsometry d) ∘ WithLp.toLp 2 := by
    funext p
    exact (snocLinearEquiv_apply p.1 p.2).symm
  rwa [heq]

/-- Stetigkeit und schneller Abfall jeder Potenz ergeben `L¹`. -/
lemma rapidDecay_integrable {d : ℕ} {F : EuclideanSpace ℝ (Fin d) → ℂ}
    (hF : Continuous F) (hdecay : RapidDecay F) : Integrable F := by
  let r : ℝ := (d : ℝ) + 1
  have hr : (Module.finrank ℝ (EuclideanSpace ℝ (Fin d)) : ℝ) < r := by
    simp [r]
  have hrcast : r = ((d + 1 : ℕ) : ℝ) := by simp [r]
  have hint : Integrable (fun x : EuclideanSpace ℝ (Fin d) => (1 + ‖x‖) ^ (-r)) :=
    integrable_one_add_norm hr
  obtain ⟨C, hC⟩ := (Asymptotics.isBigO_iff.mp (hdecay (d + 1)))
  obtain ⟨R, hRout⟩ :=
    (Metric.mem_cocompact_iff_closedBall_compl_subset (0 : EuclideanSpace ℝ (Fin d))).mp hC
  let R' : ℝ := max R 1
  have hsub : (Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R')ᶜ ⊆
      {x | ‖F x‖ ≤ C * ‖‖x‖ ^ (-r)‖} := by
    intro x hx
    have hxR : x ∈ (Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R)ᶜ := by
      intro hxball
      exact hx (Metric.closedBall_subset_closedBall (le_max_left R 1) hxball)
    simpa [hrcast] using hRout hxR
  have hball : IntegrableOn F (Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R') :=
    hF.locallyIntegrable.integrableOn_isCompact (isCompact_closedBall _ _)
  set C' : ℝ := max C 0 * (2 : ℝ) ^ r
  have hrest : IntegrableOn F (Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R')ᶜ := by
    have hbase : IntegrableOn (fun x : EuclideanSpace ℝ (Fin d) =>
        C' * (1 + ‖x‖) ^ (-r)) (Metric.closedBall 0 R')ᶜ :=
      (hint.const_mul C').integrableOn
    rw [IntegrableOn] at hbase ⊢
    refine Integrable.mono' hbase (hF.aestronglyMeasurable.restrict) ?_
    filter_upwards [ae_restrict_mem (measurableSet_closedBall.compl)] with x hx
    have hx1 : (1 : ℝ) < ‖x‖ := by
      rw [Set.mem_compl_iff, Metric.mem_closedBall, dist_zero_right, not_le] at hx
      exact lt_of_le_of_lt (le_max_right R 1) hx
    have hpow : ‖x‖ ^ (-r) ≤ (2 : ℝ) ^ r * (1 + ‖x‖) ^ (-r) := by
      have hle : 1 + ‖x‖ ≤ 2 * ‖x‖ := by linarith
      have hcmp := Real.rpow_le_rpow_of_nonpos (by linarith : (0 : ℝ) < 1 + ‖x‖) hle
        (by linarith : -r ≤ 0)
      have hsplit : (2 * ‖x‖) ^ (-r) = (2 : ℝ) ^ (-r) * ‖x‖ ^ (-r) :=
        Real.mul_rpow (by norm_num) (norm_nonneg _)
      have h2 : (2 : ℝ) ^ r * (2 : ℝ) ^ (-r) = 1 := by
        rw [← Real.rpow_add (by norm_num : (0 : ℝ) < 2)]
        simp
      calc
        ‖x‖ ^ (-r) = 1 * ‖x‖ ^ (-r) := by rw [one_mul]
        _ = ((2 : ℝ) ^ r * (2 : ℝ) ^ (-r)) * ‖x‖ ^ (-r) := by rw [h2]
        _ = (2 : ℝ) ^ r * ((2 : ℝ) ^ (-r) * ‖x‖ ^ (-r)) := by rw [mul_assoc]
        _ = (2 : ℝ) ^ r * (2 * ‖x‖) ^ (-r) := by rw [← hsplit]
        _ ≤ (2 : ℝ) ^ r * (1 + ‖x‖) ^ (-r) :=
          mul_le_mul_of_nonneg_left hcmp (by positivity)
    have hnorm : ‖‖x‖ ^ (-r)‖ = ‖x‖ ^ (-r) := by
      rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _)]
    have hC' : C * ‖x‖ ^ (-r) ≤ C' * (1 + ‖x‖) ^ (-r) := by
      calc
        C * ‖x‖ ^ (-r) ≤ max C 0 * ‖x‖ ^ (-r) :=
          mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.rpow_nonneg (norm_nonneg _) _)
        _ ≤ max C 0 * ((2 : ℝ) ^ r * (1 + ‖x‖) ^ (-r)) :=
          mul_le_mul_of_nonneg_left hpow (le_max_right _ _)
        _ = C' * (1 + ‖x‖) ^ (-r) := by ring
    calc
      ‖F x‖ ≤ C * ‖‖x‖ ^ (-r)‖ := hsub hx
      _ = C * ‖x‖ ^ (-r) := by rw [hnorm]
      _ ≤ C' * (1 + ‖x‖) ^ (-r) := hC'
  rw [← integrableOn_univ]
  have hcover : Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R' ∪
      (Metric.closedBall 0 R')ᶜ = Set.univ := Set.union_compl_self _
  rw [← hcover]
  exact hball.union hrest

/-- Globale Schranke: schneller Abfall und Stetigkeit ergeben `‖F x‖ ≤ C (1+‖x‖)^{-N}`. -/
lemma rapidDecay_one_add_le {d : ℕ} {F : EuclideanSpace ℝ (Fin d) → ℂ}
    (hF : Continuous F) (hdecay : RapidDecay F) (N : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x, ‖F x‖ ≤ C * (1 + ‖x‖) ^ (-(N : ℝ)) := by
  obtain ⟨C0, hC0⟩ := Asymptotics.isBigO_iff.mp (hdecay N)
  obtain ⟨R, hRout⟩ :=
    (Metric.mem_cocompact_iff_closedBall_compl_subset (0 : EuclideanSpace ℝ (Fin d))).mp hC0
  let R' : ℝ := max R 1
  have hbd : BddAbove ((fun x => ‖F x‖) '' Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R') :=
    (isCompact_closedBall 0 R').bddAbove_image hF.norm.continuousOn
  obtain ⟨M, hM⟩ := hbd
  let Cball : ℝ := max M 0 * (1 + R') ^ (N : ℝ)
  let Cout : ℝ := max C0 0 * (2 : ℝ) ^ (N : ℝ)
  have hCnonneg : 0 ≤ max Cball Cout :=
    le_trans (mul_nonneg (le_max_right M 0) (Real.rpow_nonneg (by positivity) _)) (le_max_left _ _)
  refine ⟨max Cball Cout, hCnonneg, ?_⟩
  intro x
  by_cases hx : x ∈ Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R'
  · have hFx : ‖F x‖ ≤ max M 0 :=
      ((mem_upperBounds.mp hM) _ (Set.mem_image_of_mem _ hx)).trans (le_max_left _ _)
    have hfactor : (1 : ℝ) ≤ (1 + R') ^ (N : ℝ) * (1 + ‖x‖) ^ (-(N : ℝ)) := by
      have hnx : ‖x‖ ≤ R' := by simpa [Metric.mem_closedBall, dist_zero_right] using hx
      have hbase : 1 + ‖x‖ ≤ 1 + R' := by linarith
      have hpow := Real.rpow_le_rpow_of_nonpos (by positivity : (0 : ℝ) < 1 + ‖x‖) hbase
        (neg_nonpos.mpr (Nat.cast_nonneg N))
      calc
        (1 : ℝ) = (1 + R') ^ (N : ℝ) * (1 + R') ^ (-(N : ℝ)) := by
          rw [← Real.rpow_add (by positivity : (0 : ℝ) < 1 + R')]
          simp
        _ ≤ (1 + R') ^ (N : ℝ) * (1 + ‖x‖) ^ (-(N : ℝ)) :=
          mul_le_mul_of_nonneg_left hpow (by positivity)
    calc
      ‖F x‖ ≤ max M 0 * 1 := by simpa using hFx
      _ ≤ max M 0 * ((1 + R') ^ (N : ℝ) * (1 + ‖x‖) ^ (-(N : ℝ))) :=
        mul_le_mul_of_nonneg_left hfactor (le_max_right _ _)
      _ = Cball * (1 + ‖x‖) ^ (-(N : ℝ)) := by ring
      _ ≤ max Cball Cout * (1 + ‖x‖) ^ (-(N : ℝ)) :=
        mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.rpow_nonneg (by positivity) _)
  · have hxR : x ∈ (Metric.closedBall (0 : EuclideanSpace ℝ (Fin d)) R)ᶜ := by
      intro hxball
      exact hx (Metric.closedBall_subset_closedBall (le_max_left R 1) hxball)
    have hx1 : (1 : ℝ) < ‖x‖ := by
      have hnot : ¬ ‖x‖ ≤ R' := by simpa [Metric.mem_closedBall, dist_zero_right] using hx
      exact lt_of_le_of_lt (le_max_right R 1) (lt_of_not_ge hnot)
    have hpow : ‖x‖ ^ (-(N : ℝ)) ≤ (2 : ℝ) ^ (N : ℝ) * (1 + ‖x‖) ^ (-(N : ℝ)) := by
      have hle : 1 + ‖x‖ ≤ 2 * ‖x‖ := by linarith
      have hcmp := Real.rpow_le_rpow_of_nonpos (by linarith : (0 : ℝ) < 1 + ‖x‖) hle
        (neg_nonpos.mpr (Nat.cast_nonneg N))
      have hsplit : (2 * ‖x‖) ^ (-(N : ℝ)) = (2 : ℝ) ^ (-(N : ℝ)) * ‖x‖ ^ (-(N : ℝ)) :=
        Real.mul_rpow (by norm_num) (norm_nonneg _)
      have h2 : (2 : ℝ) ^ (N : ℝ) * (2 : ℝ) ^ (-(N : ℝ)) = 1 := by
        rw [← Real.rpow_add (by norm_num : (0 : ℝ) < 2)]
        simp
      calc
        ‖x‖ ^ (-(N : ℝ)) = 1 * ‖x‖ ^ (-(N : ℝ)) := by rw [one_mul]
        _ = ((2 : ℝ) ^ (N : ℝ) * (2 : ℝ) ^ (-(N : ℝ))) * ‖x‖ ^ (-(N : ℝ)) := by rw [h2]
        _ = (2 : ℝ) ^ (N : ℝ) * ((2 : ℝ) ^ (-(N : ℝ)) * ‖x‖ ^ (-(N : ℝ))) := by rw [mul_assoc]
        _ = (2 : ℝ) ^ (N : ℝ) * (2 * ‖x‖) ^ (-(N : ℝ)) := by rw [← hsplit]
        _ ≤ (2 : ℝ) ^ (N : ℝ) * (1 + ‖x‖) ^ (-(N : ℝ)) :=
          mul_le_mul_of_nonneg_left hcmp (by positivity)
    have hnorm : ‖‖x‖ ^ (-(N : ℝ))‖ = ‖x‖ ^ (-(N : ℝ)) := by
      rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _)]
    calc
      ‖F x‖ ≤ C0 * ‖‖x‖ ^ (-(N : ℝ))‖ := hRout hxR
      _ = C0 * ‖x‖ ^ (-(N : ℝ)) := by rw [hnorm]
      _ ≤ max C0 0 * ‖x‖ ^ (-(N : ℝ)) :=
        mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.rpow_nonneg (norm_nonneg _) _)
      _ ≤ max C0 0 * ((2 : ℝ) ^ (N : ℝ) * (1 + ‖x‖) ^ (-(N : ℝ))) :=
        mul_le_mul_of_nonneg_left hpow (le_max_right _ _)
      _ = Cout * (1 + ‖x‖) ^ (-(N : ℝ)) := by ring
      _ ≤ max Cball Cout * (1 + ‖x‖) ^ (-(N : ℝ)) :=
        mul_le_mul_of_nonneg_right (le_max_right _ _) (Real.rpow_nonneg (by positivity) _)

lemma inner_snocPoint {d : ℕ} (y x' : EuclideanSpace ℝ (Fin d)) (s t : ℝ) :
    ⟪snocPoint y s, snocPoint x' t⟫ = ⟪y, x'⟫ + ⟪s, t⟫ := by
  rw [PiLp.inner_apply, PiLp.inner_apply]
  unfold snocPoint
  rw [Fin.sum_univ_castSucc]
  simp [Fin.snoc_castSucc, Fin.snoc_last, RCLike.inner_apply]

/-- Fubini: Fourier in der letzten Koordinate, danach in den ersten `d`, ist `𝓕` auf `ℝ^{d+1}`.
Die Phase hat Betrag 1, und `snocPoint` erhält das Lebesgue-Maß. -/
lemma fourier_iterate_snoc {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F)
    (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) :
    𝓕 (fun y => 𝓕 (coordSlice F y) t) x' = 𝓕 F (snocPoint x' t) := by
  have hFint : Integrable F := rapidDecay_integrable hF hdecay
  let ξ : EuclideanSpace ℝ (Fin (d + 1)) := snocPoint x' t
  simp_rw [Real.fourier_eq, coordSlice]
  have hslice (y : EuclideanSpace ℝ (Fin d)) :
      𝐞 (-⟪y, x'⟫) • (∫ s : ℝ, 𝐞 (-⟪s, t⟫) • F (snocPoint y s)) =
        ∫ s : ℝ, 𝐞 (-(⟪y, x'⟫ + ⟪s, t⟫)) • F (snocPoint y s) := by
    suffices (𝐞 (-⟪y, x'⟫) : ℂ) •
        (∫ s : ℝ, (𝐞 (-⟪s, t⟫) : ℂ) • F (snocPoint y s)) =
        ∫ s : ℝ, (𝐞 (-(⟪y, x'⟫ + ⟪s, t⟫)) : ℂ) • F (snocPoint y s) by
      exact this
    rw [← integral_smul]
    refine integral_congr_ae (Filter.Eventually.of_forall fun s => ?_)
    change (𝐞 (-⟪y, x'⟫) : ℂ) • ((𝐞 (-⟪s, t⟫) : ℂ) • F (snocPoint y s)) =
      (𝐞 (-(⟪y, x'⟫ + ⟪s, t⟫)) : ℂ) • F (snocPoint y s)
    rw [← mul_smul]
    refine congr_arg (fun c : ℂ => c • F (snocPoint y s)) ?_
    rw [← Circle.coe_mul, ← AddChar.map_add_eq_mul, ← neg_add]
  simp_rw [hslice, ← inner_snocPoint]
  have hg : Integrable
      (fun p : EuclideanSpace ℝ (Fin d) × ℝ =>
        𝐞 (-⟪snocPoint p.1 p.2, ξ⟫) • F (snocPoint p.1 p.2)) (volume.prod volume) :=
    (measurePreserving_snocPoint d).integrable_comp_of_integrable
      ((Real.fourierIntegral_convergent_iff ξ).2 hFint)
  rw [← integral_prod _ hg]
  let e : (EuclideanSpace ℝ (Fin d) × ℝ) ≃ᵐ EuclideanSpace ℝ (Fin (d + 1)) :=
    (MeasurableEquiv.toLp 2 (EuclideanSpace ℝ (Fin d) × ℝ)).trans
      (snocLinearIsometry d).toMeasurableEquiv
  have hembed : MeasurableEmbedding (fun p : EuclideanSpace ℝ (Fin d) × ℝ => snocPoint p.1 p.2) := by
    have hcoe : (e : (EuclideanSpace ℝ (Fin d) × ℝ) → EuclideanSpace ℝ (Fin (d + 1))) =
        fun p => snocPoint p.1 p.2 := by
      funext p
      simp only [e, MeasurableEquiv.coe_trans, LinearIsometryEquiv.coe_toMeasurableEquiv]
      exact (snocLinearEquiv_apply p.1 p.2).symm
    rw [← hcoe]
    exact e.measurableEmbedding
  exact (measurePreserving_snocPoint d).integral_comp hembed
    (fun v => 𝐞 (-⟪v, ξ⟫) • F v)

lemma snocPoint_norm_ge_fst {d : ℕ} (x' : EuclideanSpace ℝ (Fin d)) (t : ℝ) :
    ‖x'‖ ≤ ‖snocPoint x' t‖ := by
  have hsq : ‖x'‖ ^ 2 ≤ ‖snocPoint x' t‖ ^ 2 := by
    rw [snocPoint_norm_sq]
    exact le_add_of_nonneg_right (sq_nonneg t)
  have habs := sq_le_sq.mp hsq
  rwa [abs_of_nonneg (norm_nonneg _), abs_of_nonneg (norm_nonneg _)] at habs

lemma tendsto_snocPoint_fst_cocompact {d : ℕ} (t : ℝ) :
    Filter.Tendsto (fun x' : EuclideanSpace ℝ (Fin d) => snocPoint x' t)
      (Filter.cocompact _) (Filter.cocompact _) := by
  refine tendsto_cocompact_of_tendsto_dist_comp_atTop (0 : EuclideanSpace ℝ (Fin (d + 1))) ?_
  simp_rw [dist_zero_right]
  exact Filter.tendsto_atTop_mono (fun x' => snocPoint_norm_ge_fst x' t)
    (tendsto_norm_cocompact_atTop (E := EuclideanSpace ℝ (Fin d)))

/-- Transversaler Abfall: `ξ' ↦ 𝓕 F (snoc ξ' t)` erbt `RapidDecay` von `𝓕 F`. -/
lemma partialFourier_transverse_decay {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (hFdecay : RapidDecay (𝓕 F)) (t : ℝ) :
    RapidDecay (𝓕 (fun x' : EuclideanSpace ℝ (Fin d) => 𝓕 (coordSlice F x') t)) := by
  intro k
  have hid : (fun x' => 𝓕 (fun y => 𝓕 (coordSlice F y) t) x') =
      fun x' => (𝓕 F) (snocPoint x' t) := by
    ext x'
    exact fourier_iterate_snoc F hF hdecay x' t
  change (fun x' => 𝓕 (fun y => 𝓕 (coordSlice F y) t) x') =O[Filter.cocompact
    (EuclideanSpace ℝ (Fin d))] fun x => ‖x‖ ^ (-(k : ℝ))
  rw [hid]
  refine ((hFdecay k).comp_tendsto (tendsto_snocPoint_fst_cocompact t)).trans ?_
  rw [Asymptotics.isBigO_iff]
  refine ⟨1, ?_⟩
  filter_upwards [(tendsto_norm_cocompact_atTop
      (E := EuclideanSpace ℝ (Fin d))).eventually (Filter.eventually_ge_atTop (1 : ℝ))] with x' hx'
  simp only [Function.comp_apply, Real.norm_eq_abs, one_mul,
    abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _)]
  have hx0 : (0 : ℝ) < ‖x'‖ := lt_of_lt_of_le zero_lt_one hx'
  exact Real.rpow_le_rpow_of_nonpos hx0 (snocPoint_norm_ge_fst x' t)
    (neg_nonpos.mpr (Nat.cast_nonneg k))

lemma continuous_snocPoint_fst {d : ℕ} (t : ℝ) :
    Continuous fun y : EuclideanSpace ℝ (Fin d) => snocPoint y t := by
  unfold snocPoint
  have hlin : Continuous (⇑(WithLp.linearEquiv 2 ℝ (Fin (d + 1) → ℝ)).symm) :=
    (WithLp.linearEquiv 2 ℝ (Fin (d + 1) → ℝ)).symm.toLinearMap.continuous_of_finiteDimensional
  have hpi : Continuous fun y : EuclideanSpace ℝ (Fin d) =>
      (Fin.snoc (WithLp.ofLp y) t : Fin (d + 1) → ℝ) := by
    refine continuous_pi fun i => ?_
    induction i using Fin.lastCases with
    | last =>
      simp only [Fin.snoc_last]
      exact continuous_const
    | cast j =>
      simpa [Fin.snoc_castSucc] using (PiLp.continuous_apply 2 (fun _ : Fin d => ℝ) j)
  convert hlin.comp hpi using 2 with y
  simp [WithLp.coe_symm_linearEquiv]

/-- `(1 + ‖snoc y s‖)^{-(a+b)} ≤ (1+‖y‖)^{-a} (1+‖s‖)^{-b}`. -/
lemma one_add_snoc_rpow_le {d : ℕ} (y : EuclideanSpace ℝ (Fin d)) (s : ℝ) (a b : ℕ) :
    (1 + ‖snocPoint y s‖) ^ (-((a + b : ℕ) : ℝ)) ≤
      (1 + ‖y‖) ^ (-(a : ℝ)) * (1 + ‖s‖) ^ (-(b : ℝ)) := by
  have hr : (0 : ℝ) < 1 + ‖snocPoint y s‖ := by positivity
  have hcast : (-((a + b : ℕ) : ℝ)) = (-(a : ℝ)) + (-(b : ℝ)) := by
    simp [Nat.cast_add]
    ring
  rw [hcast, Real.rpow_add hr]
  refine mul_le_mul ?_ ?_ (Real.rpow_nonneg (by positivity) _) (Real.rpow_nonneg (by positivity) _)
  · exact Real.rpow_le_rpow_of_nonpos (by positivity)
      (by linarith [snocPoint_norm_ge_fst y s]) (neg_nonpos.mpr (Nat.cast_nonneg a))
  · have hs : ‖s‖ ≤ ‖snocPoint y s‖ := by
      simpa [Real.norm_eq_abs] using snocPoint_norm_ge y s
    exact Real.rpow_le_rpow_of_nonpos (by positivity) (by linarith)
      (neg_nonpos.mpr (Nat.cast_nonneg b))

lemma integrable_one_add_abs_inv_sq :
    Integrable fun s : ℝ => (1 + ‖s‖) ^ (-(2 : ℝ)) :=
  integrable_one_add_norm (by simp [Module.finrank_self])

lemma continuous_fourier_of_integrable {d : ℕ} {F : EuclideanSpace ℝ (Fin d) → ℂ}
    (hF : Integrable F) : Continuous (𝓕 F) :=
  VectorFourier.fourierIntegral_continuous Real.continuous_fourierChar
    (innerSL ℝ).continuous₂ hF

/-- `y ↦ 𝓕(Schnitt bei y)(t)` ist stetig. Die Majorante `(1+‖s‖)^{-2}` hängt nicht von `y` ab. -/
lemma partialFourier_continuous {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (t : ℝ) :
    Continuous fun y : EuclideanSpace ℝ (Fin d) => 𝓕 (coordSlice F y) t := by
  obtain ⟨C, hC0, hC⟩ := rapidDecay_one_add_le hF hdecay 2
  refine continuous_of_dominated
      (F := fun y s => 𝐞 (-⟪s, t⟫) • F (snocPoint y s))
      (bound := fun s => C * (1 + ‖s‖) ^ (-(2 : ℝ))) ?_ ?_
      (integrable_one_add_abs_inv_sq.const_mul C) ?_
  · intro y
    have hphase : Continuous fun s : ℝ => ⟪s, t⟫ := by
      simp_rw [RCLike.inner_apply]
      exact (continuous_const.mul continuous_id : Continuous fun s : ℝ => t * s)
    exact (Continuous.smul (Real.continuous_fourierChar.comp (continuous_neg.comp hphase))
      (hF.comp (continuous_snocPoint y))).aestronglyMeasurable
  · intro y
    filter_upwards with s
    rw [Circle.norm_smul]
    calc
      ‖F (snocPoint y s)‖ ≤ C * (1 + ‖snocPoint y s‖) ^ (-(2 : ℝ)) := hC _
      _ ≤ C * (1 + ‖s‖) ^ (-(2 : ℝ)) :=
        mul_le_mul_of_nonneg_left (by simpa using one_add_snoc_rpow_le y s 0 2) hC0
  · filter_upwards with s
    exact ((hF.comp (continuous_snocPoint_fst s)).const_smul (𝐞 (-⟪s, t⟫)))

/-- Gleichmäßige Schranke in der Frequenz: `‖G_t(y)‖ ≤ C_k (1+‖y‖)^{-k}`. -/
lemma partialFourier_norm_le {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (k : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (t : ℝ) (y : EuclideanSpace ℝ (Fin d)),
      ‖𝓕 (coordSlice F y) t‖ ≤ C * (1 + ‖y‖) ^ (-(k : ℝ)) := by
  obtain ⟨C0, hC0, hFbound⟩ := rapidDecay_one_add_le hF hdecay (k + 2)
  refine ⟨C0 * ∫ s : ℝ, (1 + ‖s‖) ^ (-(2 : ℝ)),
      mul_nonneg hC0 (integral_nonneg fun _ => Real.rpow_nonneg (by positivity) _), ?_⟩
  intro t y
  rw [Real.fourier_eq]
  refine (norm_integral_le_integral_norm _).trans ?_
  have hpoint : ∀ s : ℝ, ‖𝐞 (-⟪s, t⟫) • F (snocPoint y s)‖ ≤
      C0 * (1 + ‖y‖) ^ (-(k : ℝ)) * (1 + ‖s‖) ^ (-(2 : ℝ)) := by
    intro s
    rw [Circle.norm_smul]
    calc
      ‖F (snocPoint y s)‖
          ≤ C0 * (1 + ‖snocPoint y s‖) ^ (-((k + 2 : ℕ) : ℝ)) := hFbound _
      _ ≤ C0 * ((1 + ‖y‖) ^ (-(k : ℝ)) * (1 + ‖s‖) ^ (-(2 : ℝ))) := by
        exact mul_le_mul_of_nonneg_left (one_add_snoc_rpow_le y s k 2) hC0
      _ = C0 * (1 + ‖y‖) ^ (-(k : ℝ)) * (1 + ‖s‖) ^ (-(2 : ℝ)) := by ring
  let K : ℝ := C0 * (1 + ‖y‖) ^ (-(k : ℝ))
  have hint : Integrable fun s : ℝ => K * (1 + ‖s‖) ^ (-(2 : ℝ)) :=
    integrable_one_add_abs_inv_sq.const_mul K
  have hcomplex : Integrable fun s : ℝ => 𝐞 (-⟪s, t⟫) • F (snocPoint y s) := by
    refine hint.mono' ?_ (Filter.Eventually.of_forall fun s => by simpa [K] using hpoint s)
    have hphase : Continuous fun s : ℝ => ⟪s, t⟫ := by
      simp_rw [RCLike.inner_apply]
      exact (continuous_const.mul continuous_id : Continuous fun s : ℝ => t * s)
    exact (Continuous.smul (Real.continuous_fourierChar.comp (continuous_neg.comp hphase))
      (hF.comp (continuous_snocPoint y))).aestronglyMeasurable
  have hpoint' : (fun s : ℝ => ‖𝐞 (-⟪s, t⟫) • F (snocPoint y s)‖) ≤
      fun s => K * (1 + ‖s‖) ^ (-(2 : ℝ)) := by
    intro s
    simpa [K] using hpoint s
  calc
    ∫ s, ‖𝐞 (-⟪s, t⟫) • F (snocPoint y s)‖
        ≤ ∫ s, K * (1 + ‖s‖) ^ (-(2 : ℝ)) :=
      integral_mono hcomplex.norm hint hpoint'
    _ = K * ∫ s, (1 + ‖s‖) ^ (-(2 : ℝ)) := integral_const_mul _ _
    _ = (C0 * ∫ s : ℝ, (1 + ‖s‖) ^ (-(2 : ℝ))) * (1 + ‖y‖) ^ (-(k : ℝ)) := by
      simp only [K]
      ring

lemma partialFourier_rapidDecay {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (t : ℝ) :
    RapidDecay fun y : EuclideanSpace ℝ (Fin d) => 𝓕 (coordSlice F y) t := by
  intro k
  obtain ⟨C, hC0, hC⟩ := partialFourier_norm_le F hF hdecay k
  rw [Asymptotics.isBigO_iff]
  refine ⟨C, ?_⟩
  filter_upwards [(tendsto_norm_cocompact_atTop
      (E := EuclideanSpace ℝ (Fin d))).eventually (Filter.eventually_ge_atTop (1 : ℝ))] with y hy
  have hcmp : (1 + ‖y‖) ^ (-(k : ℝ)) ≤ ‖y‖ ^ (-(k : ℝ)) :=
    Real.rpow_le_rpow_of_nonpos (by linarith : (0 : ℝ) < ‖y‖) (by linarith)
      (neg_nonpos.mpr (Nat.cast_nonneg k))
  calc
    ‖𝓕 (coordSlice F y) t‖ ≤ C * (1 + ‖y‖) ^ (-(k : ℝ)) := hC t y
    _ ≤ C * ‖y‖ ^ (-(k : ℝ)) :=
      mul_le_mul_of_nonneg_left hcmp hC0
    _ = C * ‖‖y‖ ^ (-(k : ℝ))‖ := by
      rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _)]

/-- Die drei Eingaben der Induktion für `G_t`. -/
lemma partialFourier_poisson_input {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (hFdecay : RapidDecay (𝓕 F)) (t : ℤ) :
    Continuous (fun x' : EuclideanSpace ℝ (Fin d) => 𝓕 (coordSlice F x') (t : ℝ)) ∧
      RapidDecay (fun x' => 𝓕 (coordSlice F x') (t : ℝ)) ∧
      RapidDecay (𝓕 (fun x' => 𝓕 (coordSlice F x') (t : ℝ))) :=
  ⟨partialFourier_continuous F hF hdecay (t : ℝ),
    partialFourier_rapidDecay F hF hdecay (t : ℝ),
    partialFourier_transverse_decay F hF hdecay hFdecay (t : ℝ)⟩

/-- Gleichmäßige Schranke im Ort: `‖G_t(y)‖ ≤ C_k (1+|t|)^{-k}`. -/
lemma partialFourier_norm_le_t {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (hFdecay : RapidDecay (𝓕 F)) (k : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (y : EuclideanSpace ℝ (Fin d)) (t : ℝ),
      ‖𝓕 (coordSlice F y) t‖ ≤ C * (1 + ‖t‖) ^ (-(k : ℝ)) := by
  have hFint : Integrable F := rapidDecay_integrable hF hdecay
  have hFcont : Continuous (𝓕 F) := continuous_fourier_of_integrable hFint
  obtain ⟨C0, hC0, hFbound⟩ := rapidDecay_one_add_le hFcont hFdecay ((d + 1) + k)
  have hint : Integrable fun ξ : EuclideanSpace ℝ (Fin d) =>
      (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) :=
    integrable_one_add_norm (by simp)
  refine ⟨C0 * ∫ ξ : EuclideanSpace ℝ (Fin d), (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)),
      mul_nonneg hC0 (integral_nonneg fun _ => Real.rpow_nonneg (by positivity) _), ?_⟩
  intro y t
  let H : EuclideanSpace ℝ (Fin d) → ℂ := fun z => 𝓕 (coordSlice F z) t
  have hHcont : Continuous H := partialFourier_continuous F hF hdecay t
  have hHdecay : RapidDecay H := partialFourier_rapidDecay F hF hdecay t
  have hHint : Integrable H := rapidDecay_integrable hHcont hHdecay
  have hFH : 𝓕 H = fun ξ => (𝓕 F) (snocPoint ξ t) := by
    ext ξ
    simpa [H] using fourier_iterate_snoc F hF hdecay ξ t
  have hFHint : Integrable (𝓕 H) := by
    rw [hFH]
    have hslice_cont : Continuous fun ξ => (𝓕 F) (snocPoint ξ t) :=
      hFcont.comp (continuous_snocPoint_fst t)
    obtain ⟨Cslice, hCslice, hslice⟩ := rapidDecay_one_add_le hFcont hFdecay (d + 1)
    have hdom : Integrable fun ξ : EuclideanSpace ℝ (Fin d) =>
        Cslice * (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) := hint.const_mul _
    refine hdom.mono' hslice_cont.aestronglyMeasurable (Filter.Eventually.of_forall fun ξ => ?_)
    calc
      ‖(𝓕 F) (snocPoint ξ t)‖ ≤ Cslice * (1 + ‖snocPoint ξ t‖) ^ (-((d + 1 : ℕ) : ℝ)) :=
        hslice _
      _ ≤ Cslice * (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) := by
        exact mul_le_mul_of_nonneg_left (by simpa using one_add_snoc_rpow_le ξ t (d + 1) 0) hCslice
  have hinv : 𝓕⁻ (𝓕 H) = H := hHcont.fourierInv_fourier_eq hHint hFHint
  change ‖H y‖ ≤ _
  rw [← congr_fun hinv y, Real.fourierInv_eq, hFH]
  refine (norm_integral_le_integral_norm _).trans ?_
  have hpoint : ∀ ξ : EuclideanSpace ℝ (Fin d),
      ‖𝐞 (⟪ξ, y⟫) • (𝓕 F) (snocPoint ξ t)‖ ≤
        C0 * (1 + ‖t‖) ^ (-(k : ℝ)) * (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) := by
    intro ξ
    rw [Circle.norm_smul]
    calc
      ‖(𝓕 F) (snocPoint ξ t)‖
          ≤ C0 * (1 + ‖snocPoint ξ t‖) ^ (-(((d + 1) + k : ℕ) : ℝ)) := hFbound _
      _ ≤ C0 * ((1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) * (1 + ‖t‖) ^ (-(k : ℝ))) := by
        have hs := one_add_snoc_rpow_le ξ t (d + 1) k
        simpa [mul_comm] using mul_le_mul_of_nonneg_left hs hC0
      _ = C0 * (1 + ‖t‖) ^ (-(k : ℝ)) * (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) := by ring
  let K : ℝ := C0 * (1 + ‖t‖) ^ (-(k : ℝ))
  have hintK : Integrable fun ξ : EuclideanSpace ℝ (Fin d) =>
      K * (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) := hint.const_mul K
  have hcomplex : Integrable fun ξ : EuclideanSpace ℝ (Fin d) =>
      𝐞 (⟪ξ, y⟫) • (𝓕 F) (snocPoint ξ t) := by
    refine hintK.mono' ?_ (Filter.Eventually.of_forall fun ξ => by simpa [K] using hpoint ξ)
    have hinner : Continuous fun ξ : EuclideanSpace ℝ (Fin d) => ⟪ξ, y⟫ :=
      continuous_inner.comp (continuous_id.prodMk continuous_const)
    exact (Continuous.smul (Real.continuous_fourierChar.comp hinner)
      (hFcont.comp (continuous_snocPoint_fst t))).aestronglyMeasurable
  calc
    ∫ ξ, ‖𝐞 (⟪ξ, y⟫) • (𝓕 F) (snocPoint ξ t)‖
        ≤ ∫ ξ, K * (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) :=
      integral_mono hcomplex.norm hintK (by
        intro ξ
        simpa [K] using hpoint ξ)
    _ = K * ∫ ξ, (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ)) := integral_const_mul _ _
    _ = (C0 * ∫ ξ : EuclideanSpace ℝ (Fin d), (1 + ‖ξ‖) ^ (-((d + 1 : ℕ) : ℝ))) *
          (1 + ‖t‖) ^ (-(k : ℝ)) := by
      simp only [K]
      ring

/-- Der 1D-Schnitt in der Frequenz fällt wie `|t|^{-2}`, gleichmäßig im Ort. -/
lemma coordSlice_fourier_decay_2 {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (hFdecay : RapidDecay (𝓕 F))
    (x' : EuclideanSpace ℝ (Fin d)) :
    (𝓕 (coordSlice F x')) =O[Filter.cocompact ℝ] (|·| ^ (-2 : ℝ)) := by
  obtain ⟨C, hC0, hC⟩ := partialFourier_norm_le_t F hF hdecay hFdecay 2
  rw [Asymptotics.isBigO_iff]
  refine ⟨C, ?_⟩
  filter_upwards [(tendsto_norm_cocompact_atTop (E := ℝ)).eventually
      (Filter.eventually_ge_atTop (1 : ℝ))] with t ht
  have hcmp : (1 + ‖t‖) ^ (-(2 : ℝ)) ≤ |t| ^ (-(2 : ℝ)) := by
    have hpow : (1 + ‖t‖) ^ (-(2 : ℝ)) ≤ ‖t‖ ^ (-(2 : ℝ)) :=
      Real.rpow_le_rpow_of_nonpos (by linarith : (0 : ℝ) < ‖t‖)
        (by linarith : ‖t‖ ≤ 1 + ‖t‖) (by norm_num : (-(2 : ℝ)) ≤ 0)
    simpa [Real.norm_eq_abs] using hpow
  calc
    ‖𝓕 (coordSlice F x') t‖ ≤ C * (1 + ‖t‖) ^ (-(2 : ℝ)) := hC x' t
    _ ≤ C * |t| ^ (-(2 : ℝ)) := mul_le_mul_of_nonneg_left hcmp hC0
    _ = C * ‖|t| ^ (-(2 : ℝ))‖ := by
      rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (abs_nonneg _) _)]

lemma partial_fourier_joint_summable {d : ℕ} (F : EuclideanSpace ℝ (Fin (d + 1)) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (hFdecay : RapidDecay (𝓕 F)) :
    Summable fun p : (Fin d → ℤ) × ℤ =>
      𝓕 (coordSlice F (zdPoint p.1)) (p.2 : ℝ) := by
  obtain ⟨Cy, hCy, hy⟩ := partialFourier_norm_le F hF hdecay (4 * d)
  obtain ⟨Ct, hCt, ht⟩ := partialFourier_norm_le_t F hF hdecay hFdecay 4
  let K : ℝ := Real.sqrt (Cy * Ct)
  have hK : 0 ≤ K := Real.sqrt_nonneg _
  have hySum : Summable fun n : Fin d → ℤ => (1 + ‖zdPoint n‖) ^ (-(2 * (d : ℝ))) := by
    cases d with
    | zero => exact Summable.of_finite
    | succ d =>
      refine ((summable_coord_decay (d + 1)).mul_left
          ((2 : ℝ) ^ (2 * ((d + 1 : ℕ) : ℝ)))).of_norm_bounded_eventually ?_
      filter_upwards [tendsto_zdPoint_norm_atTop.eventually
          (Filter.eventually_ge_atTop (1 : ℝ))] with n hn
      rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (by positivity) _)]
      have hstep : (1 + ‖zdPoint n‖) ^ (-(2 * ((d + 1 : ℕ) : ℝ))) ≤
          ‖zdPoint n‖ ^ (-(2 * ((d + 1 : ℕ) : ℝ))) :=
        Real.rpow_le_rpow_of_nonpos (by linarith) (le_add_of_nonneg_left zero_le_one)
          (by simp; linarith)
      exact hstep.trans (by simpa using zdPoint_norm_le_coordProd n hn)
  have hprod := hySum.mul_of_nonneg summable_one_add_abs_inv_sq
    (fun _ => Real.rpow_nonneg (by positivity) _) (fun _ => Real.rpow_nonneg (by positivity) _)
  refine (hprod.mul_left K).of_norm_bounded fun p => ?_
  let Ay : ℝ := (1 + ‖zdPoint p.1‖) ^ (-((4 * d : ℕ) : ℝ))
  let Bt : ℝ := (1 + ‖(p.2 : ℝ)‖) ^ (-(4 : ℝ))
  have hu := hy (p.2 : ℝ) (zdPoint p.1)
  have hv := ht (zdPoint p.1) (p.2 : ℝ)
  have hnn : 0 ≤ ‖𝓕 (coordSlice F (zdPoint p.1)) (p.2 : ℝ)‖ := norm_nonneg _
  have hAy : 0 ≤ Cy * Ay := mul_nonneg hCy (Real.rpow_nonneg (by positivity) _)
  have hBt : 0 ≤ Ct * Bt := mul_nonneg hCt (Real.rpow_nonneg (by positivity) _)
  have hmul : ‖𝓕 (coordSlice F (zdPoint p.1)) (p.2 : ℝ)‖ *
      ‖𝓕 (coordSlice F (zdPoint p.1)) (p.2 : ℝ)‖ ≤ (Cy * Ay) * (Ct * Bt) :=
    mul_le_mul hu hv hnn hAy
  have hroot : ‖𝓕 (coordSlice F (zdPoint p.1)) (p.2 : ℝ)‖ ≤
      Real.sqrt ((Cy * Ay) * (Ct * Bt)) :=
    (Real.le_sqrt hnn (mul_nonneg hAy hBt)).2 (by simpa [sq] using hmul)
  have hsplit : Real.sqrt ((Cy * Ay) * (Ct * Bt)) =
      K * (1 + ‖zdPoint p.1‖) ^ (-(2 * (d : ℝ))) * (1 + |(p.2 : ℝ)|) ^ (-(2 : ℝ)) := by
    rw [Real.sqrt_mul hAy (Ct * Bt), Real.sqrt_mul hCy Ay, Real.sqrt_mul hCt Bt]
    have hre : (Real.sqrt Cy * Real.sqrt Ay) * (Real.sqrt Ct * Real.sqrt Bt) =
        (Real.sqrt Cy * Real.sqrt Ct) * (Real.sqrt Ay * Real.sqrt Bt) := by ring
    rw [hre, ← Real.sqrt_mul hCy Ct]
    have hKeq : Real.sqrt (Cy * Ct) = K := rfl
    have hAsqrt : Real.sqrt Ay = (1 + ‖zdPoint p.1‖) ^ (-(2 * (d : ℝ))) := by
      simp only [Ay, Real.sqrt_eq_rpow]
      rw [← Real.rpow_mul (by positivity : (0 : ℝ) ≤ 1 + ‖zdPoint p.1‖)]
      congr 1
      simp [Nat.cast_mul]
      ring
    have hBsqrt : Real.sqrt Bt = (1 + |(p.2 : ℝ)|) ^ (-(2 : ℝ)) := by
      simp only [Bt, Real.sqrt_eq_rpow, Real.norm_eq_abs]
      rw [← Real.rpow_mul (by positivity : (0 : ℝ) ≤ 1 + |(p.2 : ℝ)|)]
      norm_num
    rw [hKeq, hAsqrt, hBsqrt]
    ring
  calc
    ‖𝓕 (coordSlice F (zdPoint p.1)) (p.2 : ℝ)‖ ≤ Real.sqrt ((Cy * Ay) * (Ct * Bt)) := hroot
    _ = K * ((1 + ‖zdPoint p.1‖) ^ (-(2 * (d : ℝ))) *
          (1 + |(p.2 : ℝ)|) ^ (-(2 : ℝ))) := by
      rw [hsplit]
      ring

/-- Poisson-Summation auf `ℤ^d`.
Die Induktion benutzt `Real.tsum_eq_tsum_fourier_of_rpow_decay` in der letzten Koordinate.
Beide Kanäle sind geschlossen: der transversale Abfall von `G_t` aus dem schnellen Abfall
von `F`, der longitudinale Abfall aus der Fourier-Inversion von `𝓕 F`. -/
theorem poisson_summation_zd :
    ∀ d : ℕ, ∀ F : EuclideanSpace ℝ (Fin d) → ℂ,
      Continuous F → RapidDecay F → RapidDecay (𝓕 F) →
        ∑' n : Fin d → ℤ, F (zdPoint n) =
          ∑' m : Fin d → ℤ, (𝓕 F) (zdPoint m) := by
  intro d
  induction d with
  | zero =>
    intro F hF _ _
    have hpt : zdPoint (default : Fin 0 → ℤ) = 0 := Subsingleton.elim _ _
    have hvol : (volume : Measure (EuclideanSpace ℝ (Fin 0))) = Measure.dirac 0 := by
      have hpara : parallelepiped (EuclideanSpace.basisFun (Fin 0) ℝ) =
          {(0 : EuclideanSpace ℝ (Fin 0))} := by
        ext x
        rw [mem_parallelepiped_iff, Set.mem_singleton_iff]
        constructor
        · rintro ⟨_, _, ht⟩
          simpa [Finset.univ_eq_empty, Finset.sum_empty] using ht
        · intro hx
          refine ⟨fun i : Fin 0 => (i.elim0 : ℝ), ?_, ?_⟩
          · rw [Set.mem_Icc]
            constructor <;> intro i <;> exact i.elim0
          · simp [hx, Finset.univ_eq_empty, Finset.sum_empty]
      have hμ0 : volume ({(0 : EuclideanSpace ℝ (Fin 0))} : Set _) = 1 := by
        simpa [hpara] using
          (EuclideanSpace.basisFun (Fin 0) ℝ).volume_parallelepiped
      apply Measure.ext
      intro s hs
      by_cases h0 : (0 : EuclideanSpace ℝ (Fin 0)) ∈ s
      · have hsU : s = Set.univ := by
          ext x
          simp [h0, Subsingleton.elim x (0 : EuclideanSpace ℝ (Fin 0))]
        have huniv : (Set.univ : Set (EuclideanSpace ℝ (Fin 0))) = {0} := by
          ext x
          simp [Subsingleton.elim x (0 : EuclideanSpace ℝ (Fin 0))]
        rw [hsU, huniv, hμ0, Measure.dirac_apply_of_mem (by simp)]
      · have hsE : s = ∅ := by
          ext x
          simp [Subsingleton.elim x (0 : EuclideanSpace ℝ (Fin 0)), h0]
        simp [hsE]
    have hFourier : 𝓕 F (0 : EuclideanSpace ℝ (Fin 0)) = F 0 := by
      rw [Real.fourier_eq, hvol, integral_dirac]
      simp [one_smul]
    have hzero (n : Fin 0 → ℤ) : zdPoint n = 0 := Subsingleton.elim _ _
    simp_rw [hzero, hFourier]
  | succ d ih =>
    intro F hcont hdecay hFdecay
    have hsplit := zd_tsum_snoc (rapidDecay_summable hdecay)
    have h1d (n' : Fin d → ℤ) :
        ∑' k : ℤ, F (zdPoint (Fin.snoc n' k)) =
          ∑' t : ℤ, 𝓕 (coordSlice F (zdPoint n')) (t : ℝ) := by
      simpa [zdPoint_snoc, coordSlice] using
        poisson_1d_at_zero (coordSlice_continuous hcont (zdPoint n'))
          (coordSlice_decay_2 hdecay (zdPoint n'))
          (coordSlice_fourier_decay_2 F hcont hdecay hFdecay (zdPoint n'))
    rw [hsplit]
    simp_rw [h1d]
    have hjoint := partial_fourier_joint_summable F hcont hdecay hFdecay
    rw [← hjoint.tsum_comm]
    have hsumFF := rapidDecay_summable hFdecay
    rw [zd_tsum_snoc_outer hsumFF]
    refine tsum_congr fun t => ?_
    rcases partialFourier_poisson_input F hcont hdecay hFdecay t with ⟨hGc, hGd, hGF⟩
    rw [ih (fun x' => 𝓕 (coordSlice F x') (t : ℝ)) hGc hGd hGF]
    refine tsum_congr fun m' => ?_
    rw [fourier_iterate_snoc F hcont hdecay (zdPoint m') (t : ℝ), ← zdPoint_snoc]

/-- Spezialfall `d = 8` der induktiven Formel. Der Schritt nach `E₈` bleibt `poisson_summation_e8`. -/
theorem poisson_summation_z8 (F : EuclideanSpace ℝ (Fin 8) → ℂ)
    (hF : Continuous F) (hdecay : RapidDecay F) (hFdecay : RapidDecay (𝓕 F)) :
    ∑' n : Fin 8 → ℤ, F (zdPoint n) =
      ∑' m : Fin 8 → ℤ, (𝓕 F) (zdPoint m) :=
  poisson_summation_zd 8 F hF hdecay hFdecay

/-- Poisson-Summation auf dem Standardgitter `ℤ⁸` mit Abfall `|x|^{-2}`.
In Dimension 8 reicht diese eine Potenz nicht. Die induktive Fassung mit schnellem Abfall
ist `poisson_summation_zd`. Für `E₈` wird keine der beiden Aussagen gebraucht,
solange `poisson_summation_e8` offen ist. -/
theorem integer_lattice_poisson_8 (f : EuclideanSpace ℝ (Fin 8) → ℂ) (hf : Continuous f)
    (hdecay : f =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))] fun x => ‖x‖ ^ (-2 : ℝ))
    (hFdecay : (𝓕 f) =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))] fun x => ‖x‖ ^ (-2 : ℝ)) :
    ∑' n : Fin 8 → ℤ, f (e8ToEuclidean fun i => (n i : ℝ)) =
      ∑' m : Fin 8 → ℤ, (𝓕 f) (e8ToEuclidean fun i => (m i : ℝ)) := by
  sorry

/-- Fourier-Transformation unter `x ↦ B x`.
`LinearMap.det e8CoordMap = -1`, also `|det B|⁻¹ = 1`, und
`⟪B c, z8ToE8Fun m⟫ = ⟪c, m⟫` nach `z8ToE8_pairing`.
Mathlib hat `fourier_comp_linearIsometry` nur für Isometrien, nicht für dieses `B`. -/
theorem fourier_comp_e8CoordMap (g : EuclideanSpace ℝ (Fin 8) → ℂ) (m : Fin 8 → ℤ) :
    𝓕 (g ∘ e8CoordMapOnEuclidean) (e8ToEuclidean fun i => (m i : ℝ)) =
      𝓕 g (e8ToEuclidean (z8ToE8Fun m)) := by
  set freq : EuclideanSpace ℝ (Fin 8) := e8ToEuclidean fun i => (m i : ℝ)
  set dualPt : EuclideanSpace ℝ (Fin 8) := e8ToEuclidean (z8ToE8Fun m)
  let ψ : EuclideanSpace ℝ (Fin 8) → ℂ :=
    fun y => Complex.exp (↑(-2 * π * ⟪y, dualPt⟫) * Complex.I) • g y
  by_cases hg : Integrable g
  · rw [Real.fourier_eq', Real.fourier_eq']
    have hψ : Integrable ψ := by
      have hchar := (Real.fourierIntegral_convergent_iff dualPt).mpr hg
      refine hchar.congr (ae_of_all _ fun y => ?_)
      simp [ψ, Circle.smul_def, Real.fourierChar_apply, mul_neg, neg_mul]
    have hmeas : AEStronglyMeasurable ψ (Measure.map (⇑e8CoordLinear) volume) := by
      rw [e8CoordLinear_measurePreserving]
      exact hψ.aestronglyMeasurable
    have htransfer : ∫ v, ψ v = ∫ c, ψ (e8CoordLinear c) := by
      conv_lhs => rw [← e8CoordLinear_measurePreserving]
      exact integral_map (μ := volume) e8CoordLinear.continuous_of_finiteDimensional.aemeasurable
        hmeas
    rw [htransfer]
    refine integral_congr_ae (ae_of_all _ fun c => ?_)
    dsimp only [freq, ψ, e8CoordMapOnEuclidean]
    rw [e8_fourier_phase m c]
    dsimp only [dualPt]
    unfold e8CoordMapOnEuclidean
    rfl
  · have hzero_g : 𝓕 g dualPt = 0 := by
      rw [Real.fourier_eq, integral_undef]
      exact (Real.fourierIntegral_convergent_iff dualPt).not.mpr hg
    have hzero_B : 𝓕 (g ∘ e8CoordMapOnEuclidean) freq = 0 := by
      rw [Real.fourier_eq, integral_undef]
      exact (Real.fourierIntegral_convergent_iff freq).not.mpr
        ((e8CoordLinear_integrable_comp g).not.mpr hg)
    rw [hzero_g, hzero_B]

/-- `‖x‖^s exp(-c ‖x‖²)` verschwindet im Unendlichen. -/
lemma tendsto_norm_rpow_gaussian {c : ℝ} (hc : 0 < c) (s : ℝ) :
    Filter.Tendsto (fun x : EuclideanSpace ℝ (Fin 8) => ‖x‖ ^ s * Real.exp (-c * ‖x‖ ^ 2))
      (Filter.cocompact _) (nhds 0) := by
  convert Filter.Tendsto.comp
      ((rpow_mul_exp_neg_mul_sq_isLittleO_exp_neg hc s).tendsto_zero_of_tendsto (y := 0)
        (Real.tendsto_exp_atBot.comp
          (Filter.tendsto_id.const_mul_atTop_of_neg (show -(1 / 2 : ℝ) < 0 by norm_num))))
      (tendsto_norm_cocompact_atTop (E := EuclideanSpace ℝ (Fin 8))) using 1
  ext x
  rfl

/-- Die Gauß-Funktion und jede negative Potenz der Norm. -/
lemma gaussian_rapid_decay {c : ℝ} (hc : 0 < c) (k : ℕ) :
    (fun x : EuclideanSpace ℝ (Fin 8) =>
        (Real.exp (-c * ‖x‖ ^ 2) : ℂ)) =O[Filter.cocompact _]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) := by
  rw [Asymptotics.isBigO_iff]
  refine ⟨1, ?_⟩
  have hbig := (tendsto_norm_rpow_gaussian hc k).eventually
    (Metric.closedBall_mem_nhds (0 : ℝ) zero_lt_one)
  have hnorm := (tendsto_norm_cocompact_atTop (E := EuclideanSpace ℝ (Fin 8))).eventually
    (Filter.eventually_ge_atTop (1 : ℝ))
  filter_upwards [hbig, hnorm] with x hx h1
  have hpos : 0 < ‖x‖ := lt_of_lt_of_le zero_lt_one h1
  have hmul : ‖x‖ ^ (k : ℝ) * Real.exp (-c * ‖x‖ ^ 2) ≤ 1 := by
    have habs : |‖x‖ ^ (k : ℝ) * Real.exp (-c * ‖x‖ ^ 2)| ≤ 1 := by
      simpa [Real.dist_eq, abs_of_nonneg (mul_nonneg (Real.rpow_nonneg (norm_nonneg _) _)
        (Real.exp_nonneg _))] using hx
    exact (abs_le.mp habs).2
  have hexp : Real.exp (-c * ‖x‖ ^ 2) ≤ ‖x‖ ^ (-(k : ℝ)) := by
    rw [← mul_le_mul_iff_of_pos_left (Real.rpow_pos_of_pos hpos k), ← Real.rpow_add hpos]
    simpa [Real.rpow_neg (norm_nonneg _), Real.rpow_zero] using hmul
  rw [Complex.norm_real, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs, one_mul,
    abs_of_nonneg (Real.exp_nonneg _), abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _)]
  exact hexp

lemma e8Gaussian_continuous (t : ℝ) : Continuous (e8Gaussian t) := by
  have hreal : Continuous fun x : EuclideanSpace ℝ (Fin 8) => -(Real.pi * t) * ‖x‖ ^ 2 :=
    continuous_const.mul (continuous_norm.pow 2)
  have hcexp : Continuous fun x : EuclideanSpace ℝ (Fin 8) =>
      cexp (↑(-(Real.pi * t) * ‖x‖ ^ 2)) :=
    Complex.continuous_exp.comp (Complex.continuous_ofReal.comp hreal)
  convert hcexp using 2 with x
  unfold e8Gaussian
  push_cast
  ring

lemma e8Gaussian_norm (t : ℝ) (x : EuclideanSpace ℝ (Fin 8)) :
    ‖e8Gaussian t x‖ = Real.exp (-(Real.pi * t) * ‖x‖ ^ 2) := by
  unfold e8Gaussian
  rw [show -((Real.pi * t : ℝ) : ℂ) * ‖x‖ ^ 2 = (↑(-(Real.pi * t) * ‖x‖ ^ 2) : ℂ) by
    push_cast; ring]
  rw [Complex.norm_exp, Complex.ofReal_re]

lemma e8Gaussian_decay {t : ℝ} (ht : 0 < t) (k : ℕ) :
    e8Gaussian t =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) := by
  have hreal := gaussian_rapid_decay (mul_pos Real.pi_pos ht) k
  refine hreal.congr' (Filter.Eventually.of_forall fun x => ?_) Filter.EventuallyEq.rfl
  unfold e8Gaussian
  rw [show -((Real.pi * t : ℝ) : ℂ) * ‖x‖ ^ 2 = (↑(-(Real.pi * t) * ‖x‖ ^ 2) : ℂ) by
    push_cast; ring]
  rw [← Complex.ofReal_exp]

lemma e8Gaussian_fourier_norm {t : ℝ} (ht : 0 < t) (ξ : EuclideanSpace ℝ (Fin 8)) :
    ‖𝓕 (e8Gaussian t) ξ‖ =
      ‖(t : ℂ) ^ (-4 : ℂ)‖ * Real.exp (-(Real.pi / t) * ‖ξ‖ ^ 2) := by
  rw [gaussian_fourier_8 ht, norm_mul, Complex.norm_exp]
  congr 1
  rw [show -π * ‖ξ‖ ^ 2 / t = (↑(-(Real.pi / t) * ‖ξ‖ ^ 2) : ℂ) by
    push_cast; ring]
  rw [Complex.ofReal_re]

lemma e8Gaussian_fourier_decay {t : ℝ} (ht : 0 < t) (k : ℕ) :
    𝓕 (e8Gaussian t) =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) := by
  have hreal := gaussian_rapid_decay (div_pos Real.pi_pos ht) k
  rw [Asymptotics.isBigO_iff] at hreal ⊢
  obtain ⟨C, hC⟩ := hreal
  refine ⟨‖(t : ℂ) ^ (-4 : ℂ)‖ * C, ?_⟩
  filter_upwards [hC] with x hx
  rw [e8Gaussian_fourier_norm ht]
  have hnorm : ‖(Real.exp (-(Real.pi / t) * ‖x‖ ^ 2) : ℂ)‖ ≤ C * ‖Complex.ofReal (‖x‖ ^ (-(k : ℝ)))‖ := by
    simpa using hx
  have habs : ‖(Real.exp (-(Real.pi / t) * ‖x‖ ^ 2) : ℂ)‖ =
      Real.exp (-(Real.pi / t) * ‖x‖ ^ 2) := by
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Real.exp_nonneg _)]
  calc
    ‖(t : ℂ) ^ (-4 : ℂ)‖ * Real.exp (-(Real.pi / t) * ‖x‖ ^ 2)
      = ‖(t : ℂ) ^ (-4 : ℂ)‖ * ‖(Real.exp (-(Real.pi / t) * ‖x‖ ^ 2) : ℂ)‖ := by
        rw [habs]
    _ ≤ ‖(t : ℂ) ^ (-4 : ℂ)‖ * (C * ‖Complex.ofReal (‖x‖ ^ (-(k : ℝ)))‖) := by
        gcongr
    _ = ‖(t : ℂ) ^ (-4 : ℂ)‖ * C * ‖Complex.ofReal (‖x‖ ^ (-(k : ℝ)))‖ := by ring

/-- `B` als stetiger linearer Automorphismus. `|det B| = 1`. -/
noncomputable def e8CoordCLE :
    EuclideanSpace ℝ (Fin 8) ≃L[ℝ] EuclideanSpace ℝ (Fin 8) :=
  e8CoordLinear.toContinuousLinearMap.toContinuousLinearEquivOfDetNeZero
    (by rw [LinearMap.det_toContinuousLinearMap, e8CoordLinear_det]; norm_num)

lemma e8CoordCLE_apply (c : EuclideanSpace ℝ (Fin 8)) :
    e8CoordCLE c = e8CoordMapOnEuclidean c := by
  simp [e8CoordCLE, e8CoordMapOnEuclidean,
    ContinuousLinearMap.toContinuousLinearEquivOfDetNeZero_apply]

lemma e8CoordCLE_symm_comp :
    e8CoordCLE.symm.toContinuousLinearMap.comp e8CoordCLE.toContinuousLinearMap =
      ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin 8)) := by
  apply ContinuousLinearMap.ext
  intro x
  exact e8CoordCLE.symm_apply_apply x

lemma e8CoordCLE_comp_symm :
    e8CoordCLE.toContinuousLinearMap.comp e8CoordCLE.symm.toContinuousLinearMap =
      ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin 8)) := by
  apply ContinuousLinearMap.ext
  intro x
  exact e8CoordCLE.apply_symm_apply x

/-- `n ↦ ∑ nᵢ eᵢ` trifft jeden Punkt von `E₈` genau einmal. -/
def e8LatticeCoord (n : Fin 8 → ℤ) : E8LatticeSet :=
  ⟨∑ i, (n i : ℝ) • e8_basis i, n, rfl⟩

lemma e8LatticeCoord_injective : Function.Injective e8LatticeCoord := by
  intro m n h
  have hval : ∑ i, (m i : ℝ) • e8_basis i = ∑ i, (n i : ℝ) • e8_basis i :=
    congrArg Subtype.val h
  have h0 : ∑ i, ((m i - n i : ℤ) : ℝ) • e8_basis i = 0 := by
    ext r
    simp only [Finset.sum_apply, Pi.smul_apply, Pi.zero_apply, smul_eq_mul, Int.cast_sub, sub_mul,
      Finset.sum_sub_distrib]
    have hm : (∑ i, (m i : ℝ) • e8_basis i) r = ∑ i, (m i : ℝ) * e8_basis i r := by
      simp [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    have hn : (∑ i, (n i : ℝ) • e8_basis i) r = ∑ i, (n i : ℝ) * e8_basis i r := by
      simp [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    rw [← hm, ← hn, congr_fun hval r, sub_self]
  have hli := (Fintype.linearIndependent_iff.mp e8_basis_linearIndependent)
    (fun i => ((m i - n i : ℤ) : ℝ)) h0
  ext i
  exact sub_eq_zero.mp (Int.cast_eq_zero.mp (hli i))

lemma e8LatticeCoord_surjective : Function.Surjective e8LatticeCoord := by
  intro v
  rcases v.property with ⟨n, hn⟩
  exact ⟨n, Subtype.ext hn.symm⟩

lemma e8CoordMapOnEuclidean_zdPoint (n : Fin 8 → ℤ) :
    e8CoordMapOnEuclidean (zdPoint n) =
      e8ToEuclidean (∑ i, (n i : ℝ) • e8_basis i) := by
  unfold e8CoordMapOnEuclidean e8CoordLinear e8ToEuclidean zdPoint
  simp only [LinearMap.comp_apply, LinearEquiv.coe_toLinearMap, WithLp.coe_linearEquiv,
    WithLp.coe_symm_linearEquiv, e8CoordMap_apply]

/-- Der Abfall `O(‖x‖ ^ (-k))` überträgt sich längs eines linearen Linksinversen. -/
lemma rapidDecay_comp_leftInverse
    (fwd bwd : ContinuousLinearMap (RingHom.id ℝ)
      (EuclideanSpace ℝ (Fin 8)) (EuclideanSpace ℝ (Fin 8)))
    (hLR : bwd.comp fwd = ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin 8)))
    (hfwd : Filter.Tendsto fwd (Filter.cocompact _) (Filter.cocompact _))
    {F : EuclideanSpace ℝ (Fin 8) → ℂ} (hF : RapidDecay F) :
    RapidDecay (F ∘ fwd) := by
  intro k
  refine ((hF k).comp_tendsto hfwd).trans ?_
  rw [Asymptotics.isBigO_iff]
  let K := ‖bwd‖
  have hK : 0 < K := by
    refine lt_of_le_of_ne (norm_nonneg _) ?_
    intro h0
    have hbwd : bwd = 0 := norm_eq_zero.mp h0.symm
    have h1 : (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin 8))) = 0 := by
      rw [← hLR, hbwd, ContinuousLinearMap.zero_comp]
    let v : EuclideanSpace ℝ (Fin 8) := WithLp.toLp 2 fun _ => (1 : ℝ)
    have hv : v ≠ 0 := by
      intro hv0
      have hcoord := congrArg (fun x : EuclideanSpace ℝ (Fin 8) => (WithLp.ofLp x) 0) hv0
      simp [v] at hcoord
    have : v = 0 := by
      have happly := congr_fun (congrArg DFunLike.coe h1) v
      simpa using happly
    exact hv this
  refine ⟨K ^ (k : ℝ), ?_⟩
  filter_upwards [(tendsto_norm_cocompact_atTop
      (E := EuclideanSpace ℝ (Fin 8))).eventually (Filter.eventually_ge_atTop (1 : ℝ))] with x hx
  have hx0 : (0 : ℝ) < ‖x‖ := lt_of_lt_of_le zero_lt_one hx
  have hlower : ‖x‖ ≤ K * ‖fwd x‖ := by
    have hid : bwd (fwd x) = x := by
      have hfun := congr_fun (congrArg DFunLike.coe hLR) x
      simpa using hfun
    calc
      ‖x‖ = ‖bwd (fwd x)‖ := by rw [hid]
      _ ≤ ‖bwd‖ * ‖fwd x‖ := bwd.le_opNorm (fwd x)
      _ = K * ‖fwd x‖ := rfl
  have hbase : ‖x‖ / K ≤ ‖fwd x‖ := (div_le_iff₀ hK).2 (by simpa [mul_comm] using hlower)
  have hpow : ‖fwd x‖ ^ (-(k : ℝ)) ≤ (‖x‖ / K) ^ (-(k : ℝ)) :=
    Real.rpow_le_rpow_of_nonpos (div_pos hx0 hK) hbase (neg_nonpos.mpr (Nat.cast_nonneg k))
  have hdiv : (‖x‖ / K) ^ (-(k : ℝ)) = K ^ (k : ℝ) * ‖x‖ ^ (-(k : ℝ)) := by
    rw [Real.div_rpow hx0.le hK.le, div_eq_mul_inv, Real.rpow_neg hx0.le, Real.rpow_neg hK.le,
      inv_inv, mul_comm]
  simp only [Function.comp_apply, Real.norm_eq_abs,
    abs_of_nonneg (Real.rpow_nonneg (norm_nonneg _) _)]
  exact hpow.trans (le_of_eq hdiv)

lemma rapidDecay_of_ofReal {F : EuclideanSpace ℝ (Fin 8) → ℂ}
    (hF : ∀ k : ℕ, F =O[Filter.cocompact _] fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ)))) :
    RapidDecay F := by
  intro k
  refine (hF k).trans ?_
  rw [Asymptotics.isBigO_iff]
  refine ⟨1, Filter.Eventually.of_forall fun x => ?_⟩
  simp [Complex.norm_real, one_mul]

/-- `ℱ(g ∘ B)(ξ) = ℱ g((B⁻¹)† ξ)`. Die Rechnung ist der Gitterfall ohne Ganzzahligkeit. -/
lemma fourier_comp_e8CoordLinear (g : EuclideanSpace ℝ (Fin 8) → ℂ)
    (ξ : EuclideanSpace ℝ (Fin 8)) :
    𝓕 (g ∘ e8CoordMapOnEuclidean) ξ =
      𝓕 g (ContinuousLinearMap.adjoint e8CoordCLE.symm.toContinuousLinearMap ξ) := by
  set A := ContinuousLinearMap.adjoint e8CoordCLE.symm.toContinuousLinearMap
  by_cases hg : Integrable g
  · rw [Real.fourier_eq', Real.fourier_eq']
    set ψ : EuclideanSpace ℝ (Fin 8) → ℂ :=
      fun y => Complex.exp (↑(-2 * π * ⟪y, A ξ⟫) * Complex.I) • g y
    have hψ : Integrable ψ := by
      have hchar := (Real.fourierIntegral_convergent_iff (A ξ)).mpr hg
      refine hchar.congr (ae_of_all _ fun y => ?_)
      simp [ψ, Circle.smul_def, Real.fourierChar_apply, mul_neg, neg_mul]
    have hmeas : AEStronglyMeasurable ψ (Measure.map (⇑e8CoordLinear) volume) := by
      rw [e8CoordLinear_measurePreserving]
      exact hψ.aestronglyMeasurable
    have htransfer : ∫ v, ψ v = ∫ c, ψ (e8CoordLinear c) := by
      conv_lhs => rw [← e8CoordLinear_measurePreserving]
      exact integral_map (μ := volume) e8CoordLinear.continuous_of_finiteDimensional.aemeasurable
        hmeas
    rw [htransfer]
    refine integral_congr_ae (ae_of_all _ fun c => ?_)
    have hinner : ⟪e8CoordLinear c, A ξ⟫ = ⟪c, ξ⟫ := by
      have h := ContinuousLinearMap.adjoint_inner_right e8CoordCLE.symm.toContinuousLinearMap
        (e8CoordCLE.toContinuousLinearMap c) ξ
      have hsym := congr_fun (congrArg DFunLike.coe e8CoordCLE_symm_comp) c
      simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.id_apply] at hsym
      rw [hsym] at h
      have happ : e8CoordCLE.toContinuousLinearMap c = e8CoordLinear c := by
        simpa [e8CoordMapOnEuclidean] using e8CoordCLE_apply c
      rw [happ] at h
      exact h
    simp only [ψ, Function.comp_apply, e8CoordMapOnEuclidean, hinner]
  · have hzero_g : 𝓕 g (A ξ) = 0 := by
      rw [Real.fourier_eq, integral_undef]
      exact (Real.fourierIntegral_convergent_iff (A ξ)).not.mpr hg
    have hzero_B : 𝓕 (g ∘ e8CoordMapOnEuclidean) ξ = 0 := by
      rw [Real.fourier_eq, integral_undef]
      exact (Real.fourierIntegral_convergent_iff ξ).not.mpr
        ((e8CoordLinear_integrable_comp g).not.mpr hg)
    rw [hzero_g, hzero_B]

lemma rapidDecay_fourier_comp (f : EuclideanSpace ℝ (Fin 8) → ℂ)
    (hFdecay : ∀ k : ℕ, (𝓕 f) =O[Filter.cocompact _]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ)))) :
    RapidDecay (𝓕 (f ∘ e8CoordMapOnEuclidean)) := by
  let A := ContinuousLinearMap.adjoint e8CoordCLE.symm.toContinuousLinearMap
  have hfun : 𝓕 (f ∘ e8CoordMapOnEuclidean) = (𝓕 f) ∘ A := by
    ext ξ
    exact fourier_comp_e8CoordLinear f ξ
  rw [hfun]
  have hcomp : (ContinuousLinearMap.adjoint e8CoordCLE.toContinuousLinearMap).comp A =
      ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin 8)) := by
    rw [← ContinuousLinearMap.adjoint_comp, e8CoordCLE_symm_comp, ContinuousLinearMap.adjoint_id]
  have hinj : LinearMap.ker A.toLinearMap = ⊥ := by
    rw [LinearMap.ker_eq_bot]
    intro ξ η hξη
    apply_fun ContinuousLinearMap.adjoint e8CoordCLE.toContinuousLinearMap at hξη
    have hid (x : EuclideanSpace ℝ (Fin 8)) :
        (ContinuousLinearMap.adjoint e8CoordCLE.toContinuousLinearMap) (A x) = x := by
      simpa using congr_fun (congrArg DFunLike.coe hcomp) x
    exact (hid ξ).symm.trans (hξη.trans (hid η))
  have htend := (LinearMap.isClosedEmbedding_of_injective hinj).tendsto_cocompact
  exact rapidDecay_comp_leftInverse A
    (ContinuousLinearMap.adjoint e8CoordCLE.toContinuousLinearMap) hcomp htend
    (rapidDecay_of_ofReal hFdecay)

/-- `E₈^* = E₈` als Gleichheit der Summationsindizes. -/
lemma tsum_latticeDual_eq (g : EuclideanSpace ℝ (Fin 8) → ℂ) :
    ∑' w : latticeDual E8LatticeSet, g (e8ToEuclidean (w : Fin 8 → ℝ)) =
      ∑' v : E8LatticeSet, g (e8ToEuclidean (v : Fin 8 → ℝ)) := by
  let e : latticeDual E8LatticeSet ≃ E8LatticeSet := Set.equivOfEq e8_self_duality
  have hcoe (w : latticeDual E8LatticeSet) : (e w : Fin 8 → ℝ) = w := by
    simp [e, Set.equivOfEq_apply]
  simp_rw [← hcoe]
  exact e.tsum_eq fun v => g (e8ToEuclidean (v : Fin 8 → ℝ))

/-- Poisson-Summation für dieses `E₈`.
Rechts steht das duale Gitter. Der Faktor `1 / vol` fehlt, weil `|det B| = 1`.
Die Summe über `E₈` ist die Summe über `B ℤ⁸`, `poisson_summation_z8` vertauscht sie mit der
Fourier-Summe, und `z8ToE8` zählt `E₈^*`. -/
theorem poisson_summation_e8 (f : EuclideanSpace ℝ (Fin 8) → ℂ) (hf : Continuous f)
    (hdecay : ∀ k : ℕ, f =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))))
    (hFdecay : ∀ k : ℕ, (𝓕 f) =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ)))) :
    ∑' v : E8LatticeSet, f (e8ToEuclidean (v : Fin 8 → ℝ)) =
      ∑' w : latticeDual E8LatticeSet, (𝓕 f) (e8ToEuclidean (w : Fin 8 → ℝ)) := by
  let F : EuclideanSpace ℝ (Fin 8) → ℂ := f ∘ e8CoordMapOnEuclidean
  have hcont : Continuous F :=
    hf.comp e8CoordLinear.continuous_of_finiteDimensional
  have hdecayF : RapidDecay (f ∘ e8CoordCLE.toContinuousLinearMap) :=
    rapidDecay_comp_leftInverse e8CoordCLE.toContinuousLinearMap
      e8CoordCLE.symm.toContinuousLinearMap e8CoordCLE_symm_comp
      e8CoordCLE.toHomeomorph.isClosedEmbedding.tendsto_cocompact
      (rapidDecay_of_ofReal hdecay)
  have hFeq : f ∘ e8CoordCLE.toContinuousLinearMap = F := by
    ext x
    simp [F, e8CoordCLE_apply]
  have hdecayF' : RapidDecay F := by rwa [hFeq] at hdecayF
  have hFdecayF : RapidDecay (𝓕 F) := rapidDecay_fourier_comp f hFdecay
  have hbij : Function.Bijective e8LatticeCoord :=
    ⟨e8LatticeCoord_injective, e8LatticeCoord_surjective⟩
  let eLattice : (Fin 8 → ℤ) ≃ E8LatticeSet := Equiv.ofBijective e8LatticeCoord hbij
  have hspatial := eLattice.tsum_eq (fun v => f (e8ToEuclidean (v : Fin 8 → ℝ)))
  have hpoint : ∀ n : Fin 8 → ℤ, f (e8ToEuclidean (eLattice n : Fin 8 → ℝ)) = F (zdPoint n) := by
    intro n
    rw [Equiv.ofBijective_apply]
    simp [F, e8LatticeCoord, e8CoordMapOnEuclidean_zdPoint]
  rw [← hspatial, tsum_congr hpoint, poisson_summation_z8 F hcont hdecayF' hFdecayF]
  have hfreq : ∀ m : Fin 8 → ℤ,
      (𝓕 F) (zdPoint m) = (𝓕 f) (e8ToEuclidean (z8ToE8Fun m)) :=
    fun m => fourier_comp_e8CoordMap f m
  rw [tsum_congr hfreq]
  have hdual := (Equiv.ofBijective z8ToE8 z8ToE8_bijective).tsum_eq
    (fun v => (𝓕 f) (e8ToEuclidean (v : Fin 8 → ℝ)))
  simp only [Equiv.ofBijective_apply] at hdual
  have hfun : ∀ m : Fin 8 → ℤ,
      (𝓕 f) (e8ToEuclidean (z8ToE8Fun m)) =
        (𝓕 f) (e8ToEuclidean (z8ToE8 m : Fin 8 → ℝ)) := by
    intro m
    simp [z8ToE8]
  rw [tsum_congr hfun, hdual, (tsum_latticeDual_eq (𝓕 f)).symm]

/-- Poisson-Summation der radialen Gauß-Funktion über `E₈`.
`poisson_summation_e8` liefert die Summe über das duale Gitter,
`e8_self_duality` setzt sie mit `E₈` gleich, `gaussian_fourier_8` wertet die Transformierte aus. -/
theorem e8_poisson_gaussian {t : ℝ} (ht : 0 < t) :
    ∑' v : E8LatticeSet, e8Gaussian t (e8ToEuclidean (v : Fin 8 → ℝ)) =
      ∑' v : E8LatticeSet, 𝓕 (e8Gaussian t) (e8ToEuclidean (v : Fin 8 → ℝ)) := by
  have hpoisson := poisson_summation_e8 (e8Gaussian t) (e8Gaussian_continuous t)
    (fun k => e8Gaussian_decay ht k) (fun k => e8Gaussian_fourier_decay ht k)
  exact hpoisson.trans (tsum_latticeDual_eq (𝓕 (e8Gaussian t)))

private lemma upperHalfPlane_of_pos_mul_I {t : ℝ} (ht : 0 < t) :
    0 < ((t : ℂ) * Complex.I).im := by
  rw [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
  simpa using ht

/-- Funktionalgleichung auf der imaginären Achse: `Θ(it) = t^{-4} Θ(i/t)`.
Das ist `e8_poisson_gaussian` plus `gaussian_fourier_8`. -/
theorem e8_theta_S {t : ℝ} (ht : 0 < t) :
    e8Theta ⟨(t : ℂ) * Complex.I, upperHalfPlane_of_pos_mul_I ht⟩ =
      (t : ℂ) ^ (-4 : ℂ) *
        e8Theta ⟨((t⁻¹ : ℝ) : ℂ) * Complex.I,
          upperHalfPlane_of_pos_mul_I (inv_pos.mpr ht)⟩ := by
  have hterm₁ (v : E8LatticeSet) :
      cexp (π * Complex.I * ((t : ℂ) * Complex.I) * (e8NormSq (v : Fin 8 → ℝ) : ℂ)) =
        e8Gaussian t (e8ToEuclidean (v : Fin 8 → ℝ)) := by
    rw [e8Gaussian, e8NormSq_eq_norm_sq]
    congr 1
    ring_nf
    rw [Complex.I_sq, Complex.ofReal_mul, Complex.ofReal_pow]
    ring
  have hterm₂ (v : E8LatticeSet) :
      cexp (π * Complex.I * (((t⁻¹ : ℝ) : ℂ) * Complex.I) *
          (e8NormSq (v : Fin 8 → ℝ) : ℂ)) =
        cexp (-π * ‖e8ToEuclidean (v : Fin 8 → ℝ)‖ ^ 2 / t) := by
    rw [e8NormSq_eq_norm_sq]
    congr 1
    ring_nf
    rw [Complex.I_sq, Complex.ofReal_inv, Complex.ofReal_pow]
    field_simp [ht.ne']
  simp only [e8Theta, hterm₁, hterm₂]
  rw [e8_poisson_gaussian ht]
  simp_rw [gaussian_fourier_8 ht]
  rw [tsum_mul_left]

open E8ModularForms

lemma thetaPhase_re (τ : ℍ) (r : ℝ) :
    (π * Complex.I * (τ : ℂ) * (r : ℂ)).re = -(Real.pi * τ.im) * r := by
  have hIre : (Complex.I * (τ : ℂ)).re = -τ.im := by
    simp [Complex.mul_re, Complex.I_re, Complex.I_im, coe_re, coe_im]
  have hIim : (Complex.I * (τ : ℂ)).im = τ.re := by
    simp [Complex.mul_im, Complex.I_re, Complex.I_im, coe_re, coe_im]
  have hcast : ((Real.pi * r : ℝ) : ℂ) = π * (r : ℂ) := by
    push_cast
    rfl
  rw [show π * Complex.I * (τ : ℂ) * (r : ℂ) =
      ((Real.pi * r : ℝ) : ℂ) * (Complex.I * (τ : ℂ)) by
    rw [hcast]; ac_rfl]
  simp [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, hIre, hIim]
  ring

lemma negInvPhase_re (τ : ℍ) (r : ℝ) :
    (-π * Complex.I * (r : ℂ) / (τ : ℂ)).re =
      -(Real.pi * (τ.im / ‖(τ : ℂ)‖ ^ 2)) * r := by
  have hzre : (-π * Complex.I * (r : ℂ)).re = 0 := by
    simp [Complex.mul_re, Complex.neg_re, Complex.ofReal_re, Complex.ofReal_im,
      Complex.I_re, Complex.I_im]
  have hzim : (-π * Complex.I * (r : ℂ)).im = -(Real.pi * r) := by
    simp [Complex.mul_im, Complex.mul_re, Complex.neg_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im]
  rw [Complex.div_re, hzre, hzim, Complex.normSq_eq_norm_sq, coe_im]
  ring

lemma thetaGaussian_continuous (τ : ℍ) :
    Continuous (fun x : EuclideanSpace ℝ (Fin 8) =>
      cexp (π * Complex.I * (τ : ℂ) * ‖x‖ ^ 2)) := by
  refine Complex.continuous_exp.comp (continuous_const.mul ?_)
  exact ((Complex.continuous_ofReal.comp continuous_norm).pow 2)

lemma thetaGaussian_decay (τ : ℍ) (k : ℕ) :
    (fun x : EuclideanSpace ℝ (Fin 8) => cexp (π * Complex.I * (τ : ℂ) * ‖x‖ ^ 2))
      =O[Filter.cocompact _] fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) := by
  have hg := e8Gaussian_decay τ.im_pos k
  have hf :
      (fun x : EuclideanSpace ℝ (Fin 8) => cexp (π * Complex.I * (τ : ℂ) * ‖x‖ ^ 2))
        =O[Filter.cocompact _] e8Gaussian τ.im := by
    rw [Asymptotics.isBigO_iff]
    refine ⟨1, Filter.Eventually.of_forall fun x => ?_⟩
    rw [one_mul, e8Gaussian_norm, Complex.norm_exp,
      show (‖x‖ : ℂ) ^ 2 = ((‖x‖ ^ 2 : ℝ) : ℂ) from (Complex.ofReal_pow ‖x‖ 2).symm,
      thetaPhase_re]
  exact hf.trans hg

lemma thetaGaussian_fourier_decay (τ : ℍ) (k : ℕ) :
    (𝓕 fun x : EuclideanSpace ℝ (Fin 8) => cexp (π * Complex.I * (τ : ℂ) * ‖x‖ ^ 2))
      =O[Filter.cocompact _] fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) := by
  have hrate : 0 < τ.im / ‖(τ : ℂ)‖ ^ 2 :=
    div_pos τ.im_pos (sq_pos_of_pos (norm_pos_iff.mpr τ.ne_zero))
  have hg := e8Gaussian_decay hrate k
  have hnorm (ξ : EuclideanSpace ℝ (Fin 8)) :
      ‖(𝓕 fun x : EuclideanSpace ℝ (Fin 8) =>
          cexp (π * Complex.I * (τ : ℂ) * ‖x‖ ^ 2)) ξ‖ =
        ‖(τ : ℂ) ^ (-4 : ℤ)‖ * ‖e8Gaussian (τ.im / ‖(τ : ℂ)‖ ^ 2) ξ‖ := by
    rw [thetaGaussian_fourier_8, norm_mul, e8Gaussian_norm, Complex.norm_exp,
      show (‖ξ‖ : ℂ) ^ 2 = ((‖ξ‖ ^ 2 : ℝ) : ℂ) from (Complex.ofReal_pow ‖ξ‖ 2).symm,
      negInvPhase_re]
  rw [Asymptotics.isBigO_iff] at hg ⊢
  obtain ⟨C, hC⟩ := hg
  refine ⟨‖(τ : ℂ) ^ (-4 : ℤ)‖ * C, ?_⟩
  filter_upwards [hC] with ξ hξ
  rw [hnorm]
  calc
    ‖(τ : ℂ) ^ (-4 : ℤ)‖ * ‖e8Gaussian (τ.im / ‖(τ : ℂ)‖ ^ 2) ξ‖
      ≤ ‖(τ : ℂ) ^ (-4 : ℤ)‖ * (C * ‖Complex.ofReal (‖ξ‖ ^ (-(k : ℝ)))‖) := by
        gcongr
    _ = ‖(τ : ℂ) ^ (-4 : ℤ)‖ * C * ‖Complex.ofReal (‖ξ‖ ^ (-(k : ℝ)))‖ := by ring

/-- S-Transformation auf der ganzen oberen Halbebene: `Θ(-1/τ) = τ^4 Θ(τ)`. -/
theorem e8Theta_S (τ : ℍ) :
    e8Theta ⟨(-↑τ)⁻¹, τ.im_inv_neg_coe_pos⟩ = (τ : ℂ) ^ (4 : ℤ) * e8Theta τ := by
  let f : EuclideanSpace ℝ (Fin 8) → ℂ :=
    fun x => cexp (π * Complex.I * (τ : ℂ) * ‖x‖ ^ 2)
  have hterm (v : E8LatticeSet) :
      cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ)) =
        f (e8ToEuclidean (v : Fin 8 → ℝ)) := by
    simp only [f]
    rw [e8NormSq_eq_norm_sq, Complex.ofReal_pow]
  have hdecay : ∀ k : ℕ, f =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) :=
    fun k => thetaGaussian_decay τ k
  have hFdecay : ∀ k : ℕ, (𝓕 f) =O[Filter.cocompact (EuclideanSpace ℝ (Fin 8))]
      fun x => Complex.ofReal (‖x‖ ^ (-(k : ℝ))) :=
    fun k => thetaGaussian_fourier_decay τ k
  have hpoisson := poisson_summation_e8 f (thetaGaussian_continuous τ) hdecay hFdecay
  rw [← tsum_congr hterm] at hpoisson
  rw [tsum_latticeDual_eq (𝓕 f)] at hpoisson
  simp_rw [show ∀ ξ, (𝓕 f) ξ =
      (τ : ℂ) ^ (-4 : ℤ) * cexp (-π * Complex.I * ‖ξ‖ ^ 2 / (τ : ℂ)) from
      fun ξ => thetaGaussian_fourier_8 τ ξ] at hpoisson
  rw [tsum_mul_left] at hpoisson
  have hterm₂ (v : E8LatticeSet) :
      cexp (-π * Complex.I * ‖e8ToEuclidean (v : Fin 8 → ℝ)‖ ^ 2 / (τ : ℂ)) =
        cexp (π * Complex.I * ((-↑τ)⁻¹) * (e8NormSq (v : Fin 8 → ℝ) : ℂ)) := by
    rw [e8NormSq_eq_norm_sq, Complex.ofReal_pow]
    congr 1
    field_simp [τ.ne_zero]
  simp_rw [hterm₂] at hpoisson
  have hleft : ∑' b : E8LatticeSet,
      cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (b : Fin 8 → ℝ) : ℂ)) = e8Theta τ := by
    simp [e8Theta]
  have hright : ∑' x : E8LatticeSet,
      cexp (π * Complex.I * ((-↑τ)⁻¹) * (e8NormSq (x : Fin 8 → ℝ) : ℂ)) =
        e8Theta ⟨(-↑τ)⁻¹, τ.im_inv_neg_coe_pos⟩ := by
    simp [e8Theta]
  rw [hleft, hright] at hpoisson
  have hpow : (τ : ℂ) ^ (4 : ℤ) * (τ : ℂ) ^ (-4 : ℤ) = 1 := by
    rw [← zpow_add₀ τ.ne_zero, add_neg_cancel, zpow_zero]
  have hmul := congrArg (fun z : ℂ => (τ : ℂ) ^ (4 : ℤ) * z) hpoisson
  rw [← mul_assoc, hpow, one_mul] at hmul
  exact hmul.symm

lemma e8_normSq_eq_gram (n : Fin 8 → ℤ) :
    e8NormSq (∑ i, (n i : ℝ) • e8_basis i) =
      ((∑ i, ∑ j, n i * n j * e8_gram_z i j : ℤ) : ℝ) := by
  unfold e8NormSq
  have hcoord (k : Fin 8) :
      (∑ i, (n i : ℝ) • e8_basis i) k = ∑ i, (n i : ℝ) * e8_basis i k := by
    simp [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  simp_rw [hcoord]
  have hsq : ∑ k, (∑ i, (n i : ℝ) * e8_basis i k) ^ 2 =
      ∑ k, (∑ i, (n i : ℝ) * e8_basis i k) * (∑ j, (n j : ℝ) * e8_basis j k) := by
    simp_rw [sq]
  rw [hsq]
  have hmul : ∀ i j k,
      (n i : ℝ) * e8_basis i k * ((n j : ℝ) * e8_basis j k) =
        ((n i : ℝ) * (n j : ℝ)) * (e8_basis i k * e8_basis j k) := by
    intro i j k
    ring
  have hrearr :
      ∑ k, (∑ i, (n i : ℝ) * e8_basis i k) * (∑ j, (n j : ℝ) * e8_basis j k) =
        ∑ i, ∑ j, ((n i : ℝ) * (n j : ℝ)) * ∑ k, e8_basis i k * e8_basis j k := by
    simp_rw [Finset.sum_mul, Finset.mul_sum, hmul]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.sum_comm]
  rw [hrearr]
  simp_rw [show ∀ i j, ∑ k, e8_basis i k * e8_basis j k = (e8_gram_z i j : ℝ) from fun i j => by
    simpa [e8_gram_r] using congrFun (congrFun e8_gram_eq i) j]
  push_cast
  rfl

/-- `nᵀ G n` ist gerade: die Diagonale von `G` ist `2`, und die Nebendiagonale tritt paarweise auf. -/
lemma e8_gram_quadratic_even (n : Fin 8 → ℤ) :
    ∃ k : ℤ, ∑ i, ∑ j, n i * n j * e8_gram_z i j = 2 * k := by
  classical
  have hdiag (i : Fin 8) : e8_gram_z i i = 2 := by
    fin_cases i <;> rfl
  have hsym (i j : Fin 8) : e8_gram_z i j = e8_gram_z j i := by
    fin_cases i <;> fin_cases j <;> rfl
  let a : Fin 8 → Fin 8 → ℤ := fun i j => n i * n j * e8_gram_z i j
  have ha (i j : Fin 8) : a i j = a j i := by
    simp only [a, hsym i j, mul_comm (n i) (n j)]
  have hdiagSum : ∑ i, a i i = 2 * ∑ i, n i * n i := by
    have hterm (i : Fin 8) : a i i = 2 * (n i * n i) := by
      simp only [a, hdiag i]
      ring
    simp_rw [hterm, ← Finset.mul_sum]
  have hsplit : ∑ i, ∑ j, a i j =
      ∑ i, (a i i + ∑ j ∈ Finset.univ.erase i, a i j) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    exact (Finset.add_sum_erase Finset.univ (a i) (Finset.mem_univ i)).symm
  have hpart (i : Fin 8) :
      ∑ j ∈ Finset.univ.erase i, a i j =
        ∑ j ∈ Finset.univ.filter (fun j => i < j), a i j +
          ∑ j ∈ Finset.univ.filter (fun j => j < i), a i j := by
    have hdisj : Disjoint (Finset.univ.filter (fun j : Fin 8 => i < j))
        (Finset.univ.filter (fun j : Fin 8 => j < i)) := by
      rw [Finset.disjoint_filter]
      intro j _ h1 h2
      exact lt_asymm h1 h2
    have hunion : (Finset.univ.filter (fun j : Fin 8 => i < j)) ∪
        (Finset.univ.filter (fun j : Fin 8 => j < i)) = Finset.univ.erase i := by
      ext j
      simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_erase, Finset.mem_univ, true_and]
      constructor
      · rintro (hij | hji)
        · exact ⟨(ne_of_lt hij).symm, trivial⟩
        · exact ⟨ne_of_lt hji, trivial⟩
      · intro hne
        rcases lt_trichotomy i j with h | h | h
        · exact Or.inl h
        · exact absurd h hne.1.symm
        · exact Or.inr h
    rw [← hunion, Finset.sum_union hdisj]
  have hrepr (pred : Fin 8 → Fin 8 → Prop) [DecidableRel pred] :
      ∑ i, ∑ j ∈ Finset.univ.filter (pred i), a i j =
        ∑ p ∈ ((Finset.univ : Finset (Fin 8)) ×ˢ Finset.univ).filter (fun p => pred p.1 p.2),
          a p.1 p.2 := by
    rw [Finset.sum_filter, Finset.sum_product]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_filter]
  have hlow : ∑ i, ∑ j ∈ Finset.univ.filter (fun j => j < i), a i j =
      ∑ i, ∑ j ∈ Finset.univ.filter (fun j => i < j), a i j := by
    rw [hrepr (fun i j => j < i), hrepr (fun i j => i < j)]
    refine Finset.sum_bij' (fun p _ => (p.2, p.1)) (fun p _ => (p.2, p.1)) ?_ ?_ ?_ ?_ ?_
    · intro p hp
      simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_univ, true_and] at hp ⊢
      exact hp
    · intro p hp
      simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_univ, true_and] at hp ⊢
      exact hp
    · intro p hp
      rfl
    · intro p hp
      rfl
    · intro p hp
      exact ha p.1 p.2
  have hoff : ∑ i, ∑ j ∈ Finset.univ.erase i, a i j =
      2 * ∑ i, ∑ j ∈ Finset.univ.filter (fun j => i < j), a i j := by
    simp_rw [hpart]
    rw [Finset.sum_add_distrib, hlow, ← two_mul]
  have hQ : ∑ i, ∑ j, a i j =
      2 * (∑ i, n i * n i + ∑ i, ∑ j ∈ Finset.univ.filter (fun j => i < j), a i j) := by
    rw [hsplit, Finset.sum_add_distrib, hdiagSum, hoff]
    ring
  refine ⟨∑ i, n i * n i +
      ∑ i, ∑ j ∈ Finset.univ.filter (fun j => i < j), n i * n j * e8_gram_z i j, ?_⟩
  simpa [a] using hQ

lemma e8_even_norm (v : E8LatticeSet) :
    ∃ k : ℤ, e8NormSq (v : Fin 8 → ℝ) = ((2 * k : ℤ) : ℝ) := by
  obtain ⟨n, hn⟩ := v.property
  obtain ⟨k, hk⟩ := e8_gram_quadratic_even n
  refine ⟨k, ?_⟩
  rw [hn, e8_normSq_eq_gram, hk]

lemma cexp_pi_I_even (k : ℤ) : cexp (π * Complex.I * ((2 * k : ℤ) : ℂ)) = 1 := by
  have h2 : cexp (2 * π * Complex.I) = 1 := by
    rw [show (2 : ℂ) * π * Complex.I = π * Complex.I + π * Complex.I by ring,
      Complex.exp_add, Complex.exp_pi_mul_I]
    norm_num
  rw [show π * Complex.I * ((2 * k : ℤ) : ℂ) = (k : ℤ) * (2 * π * Complex.I) by
    push_cast; ring, Complex.exp_int_mul, h2, one_zpow]

/-- T-Transformation: `Θ(τ+1) = Θ(τ)`, weil `‖v‖²` für `v ∈ E₈` gerade ist. -/
theorem e8_theta_T (τ : ℍ) :
    e8Theta ⟨(τ : ℂ) + 1, by simpa [Complex.add_im, Complex.one_im] using τ.im_pos⟩ =
      e8Theta τ := by
  refine tsum_congr fun v => ?_
  obtain ⟨k, hk⟩ := e8_even_norm v
  rw [hk]
  have hexp : π * Complex.I * ((τ : ℂ) + 1) * (((2 * k : ℤ) : ℝ) : ℂ) =
      π * Complex.I * (τ : ℂ) * (((2 * k : ℤ) : ℝ) : ℂ) +
        π * Complex.I * (((2 * k : ℤ) : ℝ) : ℂ) := by
    ring
  rw [hexp, Complex.exp_add]
  have h1 : cexp (π * Complex.I * (((2 * k : ℤ) : ℝ) : ℂ)) = 1 := by
    rw [show (((2 * k : ℤ) : ℝ) : ℂ) = ((2 * k : ℤ) : ℂ) by norm_cast]
    exact cexp_pi_I_even k
  rw [h1, mul_one]

open ModularForm EisensteinSeries UpperHalfPlane

/-- Koeffizienten der normierten Eisenstein-Reihe: für `n ≥ 1` ist
`(q-Entwicklung von E₄)_n = 240 σ₃(n)`. Das ist die modulare Seite, nicht die Gitterseite. -/
theorem e4_qExpansion_coeff_eq_sigma (n : ℕ+) :
    (qExpansion 1 (E₄ : UpperHalfPlane → ℂ)).coeff (n : ℕ) =
      240 * (σ 3 (n : ℕ) : ℂ) := by
  rw [E_qExpansion_coeff (k := 4) (by norm_num) ⟨2, rfl⟩ n]
  have hn : (n : ℕ) ≠ 0 := n.ne_zero
  simp only [hn, ite_false]
  have hB : bernoulli 4 = -1 / 30 := by decide +kernel
  rw [hB]
  simp only [show (4 : ℕ) - 1 = 3 by norm_num]
  have hscale : -(((2 : ℂ) * (4 : ℕ)) / ((-1 / 30 : ℚ) : ℂ)) = 240 := by norm_num
  rw [hscale]

/-- Dirichlet-Reihe der Teilersummen: `∑_{n≥1} σ₃(n) n^{-s} = ζ(s) ζ(s-3)` für `Re s > 4`.
Das ist die Faltung `σ₃ = ζ * n^3` und benutzt das `E₈`-Gitter nicht. -/
theorem e8_mellin_product_identity (s : ℂ) (hs : 4 < s.re) :
    ∑' n : ℕ+, ((σ 3 (n : ℕ) : ℕ) : ℂ) / (n : ℂ) ^ s =
      riemannZeta s * riemannZeta (s - 3) := by
  have hs1 : 1 < s.re := by linarith
  have hs3 : 1 < (s - 3).re := by
    rw [Complex.sub_re]
    norm_num
    linarith
  let σc : ArithmeticFunction ℂ := σ 3
  let ζc : ArithmeticFunction ℂ := ζ
  let pk : ArithmeticFunction ℂ := pow 3
  have hconv : σc = ζc * pk := by
    ext n
    rw [show σc = ((σ 3 : ArithmeticFunction ℕ) : ArithmeticFunction ℂ) from rfl,
      show ζc = ((ζ : ArithmeticFunction ℕ) : ArithmeticFunction ℂ) from rfl,
      show pk = ((pow 3 : ArithmeticFunction ℕ) : ArithmeticFunction ℂ) from rfl,
      ← zeta_mul_pow_eq_sigma]
    simp [natCoe_apply, mul_apply]
  have hterm_pow : ∀ n, term (↗pk) s n = term (1 : ℕ → ℂ) (s - 3) n := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp
    · rw [term_of_ne_zero hn, term_of_ne_zero hn]
      have hnC : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn
      simp only [pk, pow_apply, hn, and_false, ite_false, Pi.one_apply, natCoe_apply,
        Nat.cast_pow]
      have hpow3 : (n : ℂ) ^ (3 : ℕ) = (n : ℂ) ^ (3 : ℂ) :=
        (Complex.cpow_natCast (n : ℂ) 3).symm
      rw [hpow3, (Complex.cpow_sub (x := (n : ℂ)) (y := (3 : ℂ)) (z := s) hnC).symm,
        show (3 : ℂ) - s = -(s - 3) by ring, Complex.cpow_neg, inv_eq_one_div]
  have hsum_pow : LSeriesSummable (↗pk) s := by
    rw [LSeriesSummable, funext hterm_pow]
    exact LSeriesSummable_one_iff.mpr hs3
  have hsum_zeta : LSeriesSummable (↗ζc) s :=
    LSeriesSummable_zeta_iff.mpr hs1
  have hL : L ↗σc s = riemannZeta s * riemannZeta (s - 3) := by
    rw [hconv, LSeries_mul' hsum_zeta hsum_pow]
    have hζ : L ↗ζc s = L ↗ζ s :=
      LSeries_congr (fun {n} hn => by simp [ζc, natCoe_apply]) s
    rw [hζ, LSeries_zeta_eq_riemannZeta hs1]
    have hpow : L ↗pk s = LSeries (1 : ℕ → ℂ) (s - 3) :=
      tsum_congr hterm_pow
    rw [hpow, LSeries_one_eq_riemannZeta hs3]
  have hseries : L ↗σc s = ∑' n : ℕ, ((σ 3 n : ℕ) : ℂ) / (n : ℂ) ^ s := by
    rw [LSeries_def₀ (by simp [σc])]
    refine tsum_congr fun n => ?_
    simp [σc, natCoe_apply]
  have hg0 : (fun n : ℕ => ((σ 3 n : ℕ) : ℂ) / (n : ℂ) ^ s) 0 = 0 := by
    simp
  have hcoe : ∑' n : ℕ+, ((σ 3 (n : ℕ) : ℕ) : ℂ) / (n : ℂ) ^ s =
      ∑' n : ℕ+, (fun k : ℕ => ((σ 3 k : ℕ) : ℂ) / (k : ℂ) ^ s) n := by
    refine tsum_congr fun n => ?_
    simp
  rw [← hL, hseries, hcoe]
  exact tsum_pnat_eq_tsum_of_eq_zero
    (f := fun n : ℕ => ((σ 3 n : ℕ) : ℂ) / (n : ℂ) ^ s) hg0

/-- Für `0 < Re s < 1` liegt `s - 3` im Streifen `-3 < Re < -2`.
Dort hat `ζ` weder den Pol `1`, noch eine triviale Nullstelle `-2, -4, …`,
noch eine Nullstelle aus `Re ≥ 1`. Das ist kein Spezialfall der Riemannschen Vermutung. -/
theorem riemannZeta_sub_three_ne_zero {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    riemannZeta (s - 3) ≠ 0 := by
  set z := s - 3
  have hz_re : z.re = s.re - 3 := by simp [z, Complex.sub_re]
  have hz0 : -3 < z.re := by linarith
  have hz1 : z.re < -2 := by linarith
  set w := 1 - z
  have hw_re : w.re = 1 - z.re := by simp [w, Complex.sub_re, Complex.one_re]
  have hw3 : 3 < w.re := by linarith
  have hw4 : w.re < 4 := by linarith
  have hw_ne_neg : ∀ n : ℕ, w ≠ -n := by
    intro n hwn
    have hre := congrArg Complex.re hwn
    simp only [Complex.neg_re, Complex.natCast_re] at hre
    linarith [hw3, (Nat.cast_nonneg n : 0 ≤ (n : ℝ))]
  have hw_ne_one : w ≠ 1 := by
    intro hw1
    have : w.re = 1 := by rw [hw1]; simp
    linarith
  have hfe := riemannZeta_one_sub hw_ne_neg hw_ne_one
  have hzw : 1 - w = z := by simp [w]
  rw [hzw] at hfe
  have hw_zeta : riemannZeta w ≠ 0 := riemannZeta_ne_zero_of_one_lt_re (by linarith)
  have hGamma : Complex.Gamma w ≠ 0 := Complex.Gamma_ne_zero_of_re_pos (by linarith)
  have hπ : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hbase : (2 * (Real.pi : ℂ)) ≠ 0 := mul_ne_zero (by norm_num) hπ
  have hpow : (2 * (Real.pi : ℂ)) ^ (-w) ≠ 0 := by
    intro h0
    exact hbase ((Complex.cpow_eq_zero_iff _ _).mp h0).1
  have hcos : Complex.cos ((Real.pi : ℂ) * w / 2) ≠ 0 := by
    intro hzero
    obtain ⟨k, hk⟩ := Complex.cos_eq_zero_iff.mp hzero
    have h2 : (2 : ℂ) ≠ 0 := by norm_num
    apply_fun (· * 2) at hk
    rw [div_mul_cancel₀ _ h2, div_mul_cancel₀ _ h2] at hk
    have hw_int : w = (2 * k + 1 : ℤ) := by
      apply mul_left_cancel₀ hπ
      simpa [mul_comm] using hk
    have hre_int : w.re = (2 * (k : ℤ) + 1 : ℝ) := by
      rw [hw_int]
      simp
    have hZlt : (3 : ℤ) < 2 * k + 1 := by
      have : (3 : ℝ) < (2 * (k : ℤ) + 1 : ℝ) := by linarith
      exact_mod_cast this
    have hZgt : 2 * k + 1 < 4 := by
      have : (2 * (k : ℤ) + 1 : ℝ) < 4 := by linarith
      exact_mod_cast this
    omega
  rw [hfe]
  exact mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (by norm_num) hpow) hGamma) hcos) hw_zeta

/-!
### Gate 3B

Ein Kandidat liefert `A` und `C` mit `A(s) = C(s) ζ(s)` auf `0 < Re s < 1`.
Die Struktur prüft G3-ZF und die Faktorisierung. G3-E8 und G3-NC bleiben Auflagen
an die Konstruktion: `ζ(s-3)` besteht den Nullstellentest und ist trotzdem unzulässig,
weil es weder aus `E₈` stammt noch ohne Division durch die klassische Zeta-Funktion
eingeführt wird.
-/

/-- Isolationskandidat für Gate 3B. Keine Instanz: die Isolation ist offen. -/
structure Gate3Candidate where
  A : ℂ → ℂ
  C : ℂ → ℂ
  zero_free_in_strip : ∀ s : ℂ, 0 < s.re ∧ s.re < 1 → C s ≠ 0
  factorization : ∀ s : ℂ, 0 < s.re ∧ s.re < 1 → A s = C s * riemannZeta s

end E8ArithmeticBridge
