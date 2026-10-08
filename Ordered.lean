import Mathlib.Geometry.Convex.Cone.Pointed
import Mathlib.Algebra.Order.BigOperators.Group.List
import Operational.Raw

universe u v

variable {R : Type u} [LinearOrderedField R]
variable {V : Type v} [AddCommGroup V] [Module R V]

/--
An ordered operational system extends the raw operational system
by a raw positive state cone and positivity of measurable effects.

No topological closure, generating assumption, local tomography,
Hilbert-space structure, or quantum positivity is built in here.
-/
structure OrderedOperationalSystem
    (R : Type u) [LinearOrderedField R]
    (V : Type v) [AddCommGroup V] [Module R V]
    extends RawOperationalSystem R V where

  /-- Raw positive state set. -/
  VPlusRaw : Set V

  /-- The raw positive set contains zero. -/
  zero_mem_VPlusRaw :
    (0 : V) ∈ VPlusRaw

  /-- The raw positive set is closed under addition. -/
  add_mem_VPlusRaw :
    ∀ {x y : V},
      x ∈ VPlusRaw →
      y ∈ VPlusRaw →
      x + y ∈ VPlusRaw

  /-- The raw positive set is closed under multiplication by nonnegative scalars. -/
  smul_mem_VPlusRaw :
    ∀ {c : R} {x : V},
      0 ≤ c →
      x ∈ VPlusRaw →
      c • x ∈ VPlusRaw

  /--
  Every effect that occurs in an admissible raw measurement
  is nonnegative on raw positive states.
  -/
  effect_nonneg_raw :
    ∀ e : V →ₗ[R] R,
      (∃ M ∈ MRaw, e ∈ M) →
      ∀ x ∈ VPlusRaw,
        0 ≤ e x


namespace OrderedOperationalSystem

variable (sys : OrderedOperationalSystem R V)

/--
The positive cone on the operational quotient:
`V₊ := q_A '' V₊^raw`.

This is deliberately only an image set. No closure or completion is inserted.
-/
def VPlus : Set sys.toRawOperationalSystem.VQuot :=
  sys.toRawOperationalSystem.qA '' sys.VPlusRaw


/-- Zero belongs to the quotient positive cone. -/
lemma zero_mem_VPlus :
    (0 : sys.toRawOperationalSystem.VQuot) ∈ sys.VPlus := by
  refine ⟨0, sys.zero_mem_VPlusRaw, ?_⟩
  simp [RawOperationalSystem.qA]


/-- The quotient positive cone is closed under addition. -/
lemma add_mem_VPlus
    {xBar yBar : sys.toRawOperationalSystem.VQuot}
    (hx : xBar ∈ sys.VPlus)
    (hy : yBar ∈ sys.VPlus) :
    xBar + yBar ∈ sys.VPlus := by
  rcases hx with ⟨x, hx, rfl⟩
  rcases hy with ⟨y, hy, rfl⟩
  refine ⟨x + y, sys.add_mem_VPlusRaw hx hy, ?_⟩
  simpa using sys.toRawOperationalSystem.qA.map_add x y


/-- The quotient positive cone is closed under nonnegative scalar multiplication. -/
lemma smul_mem_VPlus
    {c : R}
    {xBar : sys.toRawOperationalSystem.VQuot}
    (hc : 0 ≤ c)
    (hx : xBar ∈ sys.VPlus) :
    c • xBar ∈ sys.VPlus := by
  rcases hx with ⟨x, hx, rfl⟩
  refine ⟨c • x, sys.smul_mem_VPlusRaw hc hx, ?_⟩
  simpa using sys.toRawOperationalSystem.qA.map_smul c x


/--
Mathlib bundle for the algebraic quotient cone.

Important: in Mathlib `PointedCone` means a convex cone containing `0`.
It does NOT encode salience (`C ∩ (-C) = {0}`), generatingness, or closedness.
-/
def vPlusPointedCone :
    PointedCone R sys.toRawOperationalSystem.VQuot where
  carrier := sys.VPlus

  zero_mem' :=
    sys.zero_mem_VPlus

  add_mem' := by
    intro x hx y hy
    exact sys.add_mem_VPlus hx hy

  smul_mem' := by
    intro c x hx
    simpa using sys.smul_mem_VPlus c.property hx


/-- The associated Mathlib `ConvexCone`. -/
def vPlusConvexCone :
    ConvexCone R sys.toRawOperationalSystem.VQuot :=
  sys.vPlusPointedCone


/--
H0.2b is derived, not assumed.

If a raw positive state has zero unit effect, then it is operationally invisible.
This follows from measurement normalization plus positivity of all measurement effects.
-/
theorem nullspace_h02b_derived
    {x : V}
    (hx : x ∈ sys.VPlusRaw)
    (hu : sys.uRaw x = 0) :
    x ∈ sys.toRawOperationalSystem.nullspace := by

  change ∀ e,
    e ∈ sys.toRawOperationalSystem.ERaw →
    e x = 0

  intro e he

  change ∃ M ∈ sys.MRaw, e ∈ M at he
  rcases he with ⟨M, hM, heM⟩

  let values : List R := M.map (fun f => f x)

  have hnonneg :
      ∀ a ∈ values, 0 ≤ a := by
    intro a ha
    rcases List.mem_map.mp ha with ⟨f, hf, rfl⟩
    exact sys.effect_nonneg_raw f ⟨M, hM, hf⟩ x hx

  have hsum : values.sum = 0 := by
    have happly :
        M.sum x = sys.uRaw x := by
      exact congrArg
        (fun f : V →ₗ[R] R => f x)
        (sys.measurement_sum M hM)
    simpa [values, hu] using happly

  exact
    List.all_zero_of_le_zero_le_of_sum_eq_zero
      hnonneg hsum
      (List.mem_map.mpr ⟨e, heM, rfl⟩)


/-- Every descended measurable effect is nonnegative on the quotient positive cone. -/
lemma effectQuot_nonneg
    {e : V →ₗ[R] R}
    (he : e ∈ sys.toRawOperationalSystem.ERaw)
    {xBar : sys.toRawOperationalSystem.VQuot}
    (hx : xBar ∈ sys.VPlus) :
    0 ≤ sys.toRawOperationalSystem.effectQuot e he xBar := by

  rcases hx with ⟨x, hx_raw, rfl⟩

  rw [RawOperationalSystem.effectQuot_mk]

  change ∃ M ∈ sys.MRaw, e ∈ M at he
  exact sys.effect_nonneg_raw e he x hx_raw


/-- The descended unit effect is nonnegative on the quotient positive cone. -/
lemma uQuot_nonneg
    {xBar : sys.toRawOperationalSystem.VQuot}
    (hx : xBar ∈ sys.VPlus) :
    0 ≤ sys.toRawOperationalSystem.uQuot xBar := by

  rcases hx with ⟨x, hx_raw, rfl⟩

  rw [RawOperationalSystem.uQuot_mk]

  exact sys.effect_nonneg_raw
    sys.uRaw
    ⟨[sys.uRaw], sys.unit_measurement, by simp⟩
    x hx_raw


/--
If the descended unit effect vanishes on a positive quotient state,
then that quotient state is zero.
-/
lemma uQuot_eq_zero_imp_eq_zero
    {xBar : sys.toRawOperationalSystem.VQuot}
    (hx : xBar ∈ sys.VPlus)
    (hu : sys.toRawOperationalSystem.uQuot xBar = 0) :
    xBar = 0 := by

  rcases hx with ⟨x, hx_raw, rfl⟩

  rw [RawOperationalSystem.uQuot_mk] at hu

  have h_ns :
      x ∈ sys.toRawOperationalSystem.nullspace :=
    sys.nullspace_h02b_derived hx_raw hu

  change
    (Submodule.Quotient.mk x :
      sys.toRawOperationalSystem.VQuot) = 0

  exact
    (Submodule.Quotient.mk_eq_zero
      sys.toRawOperationalSystem.nullspace).2 h_ns


/--
H0.1-salient (derived):
`V₊ ∩ (-V₊) = {0}`.

This is not an independent cone axiom: it follows from
measurement normalization and raw effect positivity through H0.2b.
-/
theorem VPlus_salient
    {xBar : sys.toRawOperationalSystem.VQuot}
    (hx : xBar ∈ sys.VPlus)
    (hnx : -xBar ∈ sys.VPlus) :
    xBar = 0 := by

  have hx0 :
      0 ≤ sys.toRawOperationalSystem.uQuot xBar :=
    sys.uQuot_nonneg hx

  have hnx0 :
      0 ≤ sys.toRawOperationalSystem.uQuot (-xBar) :=
    sys.uQuot_nonneg hnx

  have hneg :
      sys.toRawOperationalSystem.uQuot (-xBar) =
        -sys.toRawOperationalSystem.uQuot xBar := by
    simpa using sys.toRawOperationalSystem.uQuot.map_neg xBar

  rw [hneg] at hnx0

  have hx_le :
      sys.toRawOperationalSystem.uQuot xBar ≤ 0 :=
    neg_nonneg.mp hnx0

  have hu :
      sys.toRawOperationalSystem.uQuot xBar = 0 :=
    le_antisymm hx_le hx0

  exact sys.uQuot_eq_zero_imp_eq_zero hx hu


/-- Strict positivity of the descended unit effect on `V₊ \ {0}`. -/
lemma uQuot_strictly_positive
    {xBar : sys.toRawOperationalSystem.VQuot}
    (hx : xBar ∈ sys.VPlus)
    (hneq : xBar ≠ 0) :
    0 < sys.toRawOperationalSystem.uQuot xBar := by

  rcases lt_or_eq_of_le (sys.uQuot_nonneg hx) with hlt | heq

  · exact hlt

  · exfalso
    apply hneq
    exact sys.uQuot_eq_zero_imp_eq_zero hx heq.symm


end OrderedOperationalSystem
