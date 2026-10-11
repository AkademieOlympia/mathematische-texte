import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.Meromorphic.Basic
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Riemann.E8Dirichlet

open Complex Filter Function
open scoped Topology

noncomputable section

namespace E8ArithmeticBridge

/-!
Gate 3C.

`completedRiemannZeta` ist Mathlibs \(\Lambda\). Die Symmetrie
`completedRiemannZeta_one_sub` ergibt \(\Lambda_{E_8}(4-s)=\Lambda_{E_8}(s)\).

Das Dirichlet-Objekt \(240\,\zeta(s)\zeta(s-3)\) ist auf ganz \(\mathbb{C}\)
meromorph. Bei \(s=1\) kompensiert die triviale Nullstelle \(\zeta(-2)=0\)
den Pol von \(\zeta\), bei \(s=4\) bleibt der einfache Pol von \(\zeta(s-3)\).
-/

/-- \(\Lambda_{E_8}(s) = 240\,\Lambda(s)\Lambda(s-3)\). -/
def e8CompletedDirichlet (s : ℂ) : ℂ :=
  240 * completedRiemannZeta s * completedRiemannZeta (s - 3)

/-- Funktionalgleichung: Spiegelung an \(\operatorname{Re}s = 2\). -/
theorem e8CompletedDirichlet_functional_equation (s : ℂ) :
    e8CompletedDirichlet (4 - s) = e8CompletedDirichlet s := by
  unfold e8CompletedDirichlet
  rw [show (4 : ℂ) - s = 1 - (s - 3) by ring, completedRiemannZeta_one_sub]
  rw [show (1 - (s - 3) - 3 : ℂ) = 1 - s by ring, completedRiemannZeta_one_sub]
  ring

/-- \(D_{E_8}(s) = 240\,\zeta(s)\zeta(s-3)\). -/
def e8MeromorphicDirichlet (s : ℂ) : ℂ :=
  240 * riemannZeta s * riemannZeta (s - 3)

/-- Auf \(\operatorname{Re}s > 4\) ist das die Dirichlet-Reihe aus Gate 3B. -/
theorem e8MeromorphicDirichlet_eq_series {s : ℂ} (hs : 4 < s.re) :
    e8MeromorphicDirichlet s =
      ∑' n : ℕ+, (e8_representationNumber (2 * (n : ℕ)) : ℂ) / (n : ℂ) ^ s := by
  rw [e8RepresentationDirichlet_eq_zeta hs]
  rfl

lemma analyticAt_riemannZeta_sq_factor :
    AnalyticAt ℂ (fun z => (z - 1) ^ 2 * riemannZeta z) 1 := by
  classical
  let φ : ℂ → ℂ := fun z => (z - 1) * riemannZeta z
  let Φ : ℂ → ℂ := update φ 1 1
  have hcont : ContinuousAt Φ 1 :=
    continuousAt_update_same.mpr riemannZeta_residue_one
  have hdiff : ∀ᶠ z in 𝓝[≠] 1, DifferentiableAt ℂ Φ z := by
    filter_upwards [self_mem_nhdsWithin] with z hz
    have hz1 : z ≠ 1 := hz
    refine (DifferentiableAt.congr_of_eventuallyEq
      ((differentiableAt_id.sub_const (1 : ℂ)).mul (differentiableAt_riemannZeta hz1)) ?_)
    filter_upwards [eventually_ne_nhds hz1] with w hw
    simp [Φ, φ, update_of_ne hw]
  have hA : AnalyticAt ℂ Φ 1 :=
    analyticAt_of_differentiable_on_punctured_nhds_of_continuousAt hdiff hcont
  have heq : (fun z => (z - 1) ^ 2 * riemannZeta z) = fun z => (z - 1) * Φ z := by
    ext z
    by_cases hz : z = 1
    · simp [hz, Φ, φ]
    · simp [Φ, φ, update_of_ne hz, pow_two, mul_assoc]
  rw [heq]
  exact (analyticAt_id.sub analyticAt_const).mul hA

lemma meromorphicAt_riemannZeta (s : ℂ) : MeromorphicAt riemannZeta s := by
  by_cases hs : s = 1
  · subst hs
    exact ⟨2, analyticAt_riemannZeta_sq_factor⟩
  · refine (analyticAt_iff_eventually_differentiableAt.mpr ?_).meromorphicAt
    filter_upwards [eventually_ne_nhds hs] with z hz
    exact differentiableAt_riemannZeta hz

lemma meromorphicAt_riemannZeta_sub_three (s : ℂ) :
    MeromorphicAt (fun z => riemannZeta (z - 3)) s := by
  have hcomp : MeromorphicAt (riemannZeta ∘ fun z => z - 3) s :=
    (meromorphicAt_comp_sub_const_iff_meromorphicAt (c := (3 : ℂ))
      (f := riemannZeta) (x := s)).mpr (meromorphicAt_riemannZeta (s - 3))
  exact hcomp.congr <| Filter.EventuallyEq.of_eq (by ext z; rfl)

theorem meromorphicAt_e8MeromorphicDirichlet (s : ℂ) :
    MeromorphicAt e8MeromorphicDirichlet s := by
  have h := (((MeromorphicAt.const (240 : ℂ) s).mul (meromorphicAt_riemannZeta s)).mul
    (meromorphicAt_riemannZeta_sub_three s))
  refine h.congr ?_
  filter_upwards with z
  simp [e8MeromorphicDirichlet, Pi.mul_apply]

theorem meromorphicOn_e8MeromorphicDirichlet :
    MeromorphicOn e8MeromorphicDirichlet Set.univ :=
  fun s _ => meromorphicAt_e8MeromorphicDirichlet s

lemma tendsto_riemannZeta_slope_at_neg_two :
    Tendsto (fun s => riemannZeta (s - 3) / (s - 1)) (𝓝[≠] (1 : ℂ))
      (𝓝 (deriv riemannZeta (-2))) := by
  have hzero : riemannZeta (-2) = 0 := by
    simpa using riemannZeta_neg_two_mul_nat_add_one 0
  have hslope :=
    (differentiableAt_riemannZeta (by norm_num : (-2 : ℂ) ≠ 1)).hasDerivAt.tendsto_slope
  have hmap : Tendsto (fun s : ℂ => s - 3) (𝓝[≠] (1 : ℂ)) (𝓝[≠] (-2)) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within (fun s => s - 3) ?_ ?_
    · have hsub : Tendsto (fun s : ℂ => s - 3) (𝓝 (1 : ℂ)) (𝓝 ((1 : ℂ) - 3)) :=
        tendsto_id.sub_const _
      rw [show ((1 : ℂ) - 3) = -2 by norm_num] at hsub
      exact hsub.mono_left (nhdsWithin_le_nhds (a := (1 : ℂ)))
    · filter_upwards [self_mem_nhdsWithin] with s hs
      rw [Set.mem_compl_iff, Set.mem_singleton_iff] at hs ⊢
      intro h
      apply hs
      linear_combination h
  refine (hslope.comp hmap).congr' ?_
  filter_upwards [self_mem_nhdsWithin] with s hs
  simp only [slope, Function.comp_apply, hzero, vsub_eq_sub, sub_zero, smul_eq_mul]
  have hsub : (s - 3) - (-2) = s - 1 := by ring
  rw [hsub, div_eq_inv_mul]

/-- Bei \(s=1\) hat \(240\,\zeta(s)\zeta(s-3)\) den endlichen Grenzwert \(240\,\zeta'(-2)\). -/
theorem e8MeromorphicDirichlet_tendsto_one :
    Tendsto e8MeromorphicDirichlet (𝓝[≠] (1 : ℂ))
      (𝓝 (240 * deriv riemannZeta (-2))) := by
  have hres : Tendsto (fun s => (s - 1) * riemannZeta s) (𝓝[≠] (1 : ℂ)) (𝓝 1) :=
    riemannZeta_residue_one
  have hprod : Tendsto (fun s => riemannZeta s * riemannZeta (s - 3)) (𝓝[≠] (1 : ℂ))
      (𝓝 (deriv riemannZeta (-2))) := by
    have hmul := hres.mul tendsto_riemannZeta_slope_at_neg_two
    simp only [one_mul] at hmul
    refine hmul.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with s hs
    have hs0 : s - 1 ≠ 0 := sub_ne_zero.mpr hs
    field_simp [hs0]
  refine (tendsto_const_nhds.mul hprod).congr' ?_
  filter_upwards with s
  simp [e8MeromorphicDirichlet, mul_assoc]

/-- Residuum bei \(s=4\): \((s-4) D_{E_8}(s) \to 240\,\zeta(4)\). -/
theorem e8MeromorphicDirichlet_residue_four :
    Tendsto (fun s => (s - 4) * e8MeromorphicDirichlet s) (𝓝[≠] (4 : ℂ))
      (𝓝 (240 * riemannZeta 4)) := by
  have hcont : Tendsto riemannZeta (𝓝 (4 : ℂ)) (𝓝 (riemannZeta 4)) :=
    (differentiableAt_riemannZeta (by norm_num)).continuousAt.tendsto
  have hmap : Tendsto (fun s : ℂ => s - 3) (𝓝[≠] (4 : ℂ)) (𝓝[≠] (1 : ℂ)) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within (fun s => s - 3) ?_ ?_
    · have hsub : Tendsto (fun s : ℂ => s - 3) (𝓝 (4 : ℂ)) (𝓝 ((4 : ℂ) - 3)) :=
        tendsto_id.sub_const _
      rw [show ((4 : ℂ) - 3) = 1 by norm_num] at hsub
      exact hsub.mono_left (nhdsWithin_le_nhds (a := (4 : ℂ)))
    · filter_upwards [self_mem_nhdsWithin] with s hs
      rw [Set.mem_compl_iff, Set.mem_singleton_iff] at hs ⊢
      intro h
      apply hs
      linear_combination h
  have hres : Tendsto (fun s => (s - 4) * riemannZeta (s - 3)) (𝓝[≠] (4 : ℂ)) (𝓝 1) := by
    refine (riemannZeta_residue_one.comp hmap).congr' ?_
    filter_upwards with s
    simp only [Function.comp_apply]
    have : (s - 3) - 1 = s - 4 := by ring
    rw [this]
  have hmul : Tendsto (fun s => riemannZeta s * ((s - 4) * riemannZeta (s - 3)))
      (𝓝[≠] (4 : ℂ)) (𝓝 (riemannZeta 4)) := by
    simpa [mul_one] using (hcont.mono_left nhdsWithin_le_nhds).mul hres
  refine (tendsto_const_nhds.mul hmul).congr' ?_
  filter_upwards with s
  simp only [e8MeromorphicDirichlet]
  ring

theorem e8MeromorphicDirichlet_residue_four_ne_zero :
    240 * riemannZeta 4 ≠ 0 :=
  mul_ne_zero (by norm_num) (riemannZeta_ne_zero_of_one_lt_re (by norm_num))

end E8ArithmeticBridge
