import Mathlib

/-!
# Streaming Factorisation as a Monoid Homomorphism

This file formalises the streaming evaluation map as a `MonoidHom`:

    F : (FreeMonoid α, *, 1) →* (M, *, 1)

proves the universal property of the free monoid `List α = FreeMonoid α`,
and develops the EABC-specific corollaries.

All definitions from `streaming_factorization.lean` are reproduced here
so that this file is fully self-contained with a single `import Mathlib`.

## Main Results

| Name                              | Statement                                                              |
|-----------------------------------|------------------------------------------------------------------------|
| `streamEvalHom`                   | `streamEval φ : FreeMonoid α →* M` packaged as `MonoidHom`            |
| `streamEvalHom_is_hom`            | `streamEvalHom φ (s * t) = streamEvalHom φ s * streamEvalHom φ t`    |
| `streamEvalHom_eq_lift`           | `streamEvalHom φ = FreeMonoid.lift φ`                                 |
| `universal_property`              | Unique `MonoidHom` extending any `φ : α → M` via generators           |
| `streamEvalAddHom`                | `streamEvalAdd φ : FreeAddMonoid α →+ M` packaged as `AddMonoidHom`   |
| `eabcHom`                         | Explicit EABC `AddMonoidHom` for `PrimeParity` streams                |
| `D_is_corollary`                  | `D (eabcHom s) = (eabcHom s).1 - (eabcHom s).2`                      |
| `Q_is_corollary`                  | `Q (eabcHom s) = (eabcHom s).1 + (eabcHom s).2`                      |
| `W_E_not_monoid_hom`              | `W_E ∘ eabcHom` has no `AddMonoidHom` extension to `ℚ`               |
-/

-- ===========================================================================
-- Prerequisites  (self-contained reproduction from streaming_factorization.lean)
-- ===========================================================================

section StreamingPrereqs

variable {α : Type*}

-- ---------------------------------------------------------------------------
-- Multiplicative streaming evaluation
-- ---------------------------------------------------------------------------

section Mult

variable {M : Type*} [Monoid M]

/-- `streamEval φ s` is the monoid product of the `φ`-images of all blocks
    in the stream `s`: `F(s) = ∏ i, φ(s[i])`. -/
def streamEval (φ : α → M) (s : List α) : M :=
  (s.map φ).prod

@[simp]
theorem streamEval_nil (φ : α → M) : streamEval φ [] = 1 := by
  simp [streamEval]

@[simp]
theorem streamEval_singleton (φ : α → M) (b : α) : streamEval φ [b] = φ b := by
  simp [streamEval]

/-- Streaming Factorisation Theorem (multiplicative form):
    `F(s₁ ‖ s₂) = F(s₁) · F(s₂)`. -/
theorem StreamingFactorization (φ : α → M) (s₁ s₂ : List α) :
    streamEval φ (s₁ ++ s₂) = streamEval φ s₁ * streamEval φ s₂ := by
  simp [streamEval, List.map_append, List.prod_append]

end Mult

-- ---------------------------------------------------------------------------
-- Additive streaming evaluation
-- ---------------------------------------------------------------------------

section Add

variable {M : Type*} [AddMonoid M]

/-- `streamEvalAdd φ s` is the additive monoid sum of the `φ`-images of all
    blocks in the stream `s`: `F(s) = ∑ i, φ(s[i])`. -/
def streamEvalAdd (φ : α → M) (s : List α) : M :=
  (s.map φ).sum

@[simp]
theorem streamEvalAdd_nil (φ : α → M) : streamEvalAdd φ [] = 0 := by
  simp [streamEvalAdd]

@[simp]
theorem streamEvalAdd_singleton (φ : α → M) (b : α) : streamEvalAdd φ [b] = φ b := by
  simp [streamEvalAdd]

/-- Streaming Factorisation Theorem (additive form):
    `F(s₁ ‖ s₂) = F(s₁) + F(s₂)`. -/
theorem AdditiveStreamingFactorization (φ : α → M) (s₁ s₂ : List α) :
    streamEvalAdd φ (s₁ ++ s₂) = streamEvalAdd φ s₁ + streamEvalAdd φ s₂ := by
  simp [streamEvalAdd, List.map_append, List.sum_append]

end Add

-- ---------------------------------------------------------------------------
-- EABC state space
-- ---------------------------------------------------------------------------

/-- EABC state space `(N₊, N₋) ∈ ℤ × ℤ` with componentwise addition. -/
abbrev EABCState := ℤ × ℤ

/-- Signed imbalance: `D(N₊, N₋) := N₊ − N₋`. -/
def D (s : EABCState) : ℤ := s.1 - s.2

/-- Total density: `Q(N₊, N₋) := N₊ + N₋`. -/
def Q (s : EABCState) : ℤ := s.1 + s.2

/-- Efficiency observable: `W_E := D / Q` over ℚ.
    Declared `noncomputable` because rational division requires choice. -/
noncomputable def W_E (s : EABCState) : ℚ := (D s : ℚ) / (Q s : ℚ)

end StreamingPrereqs

-- ===========================================================================
-- §A  F as a Monoid Homomorphism
-- ===========================================================================

section StreamingHom

variable {α : Type*} {M : Type*} [Monoid M]

/-!
### §A – `streamEval` as a `MonoidHom`

`FreeMonoid α` is definitionally equal to `List α` (concatenation as
multiplication, `[]` as `1`).  Packaging `streamEval φ` as a `MonoidHom`
makes the algebraic structure explicit for Mathlib's homomorphism API.
-/

/-- `streamEvalHom φ` packages `streamEval φ` as a `MonoidHom` from the free
    monoid `FreeMonoid α` into any monoid `M`:
    - `map_one'`: empty stream → neutral element `1`.
    - `map_mul'`: stream concatenation → monoid product. -/
def streamEvalHom (φ : α → M) : FreeMonoid α →* M where
  toFun    := streamEval φ
  map_one' := streamEval_nil φ
  map_mul' := StreamingFactorization φ

@[simp]
theorem streamEvalHom_apply (φ : α → M) (s : FreeMonoid α) :
    streamEvalHom φ s = streamEval φ s := rfl

/-- The `MonoidHom` homomorphism property: evaluation respects stream
    concatenation (stated as an explicit theorem for readability). -/
theorem streamEvalHom_is_hom (φ : α → M) (s t : FreeMonoid α) :
    streamEvalHom φ (s * t) = streamEvalHom φ s * streamEvalHom φ t :=
  (streamEvalHom φ).map_mul s t

/-- `streamEvalHom φ` is definitionally the canonical extension
    `FreeMonoid.lift φ` provided by the universal property. -/
theorem streamEvalHom_eq_lift (φ : α → M) :
    streamEvalHom φ = FreeMonoid.lift φ :=
  FreeMonoid.hom_eq fun a => by
    simp [streamEvalHom_apply,
          show (FreeMonoid.of a : FreeMonoid α) = [a] from rfl,
          streamEval_singleton, FreeMonoid.lift_eval_of]

end StreamingHom

-- ===========================================================================
-- §B  Universal Property of the Free Monoid
-- ===========================================================================

section UniversalProperty

variable {α : Type*} {M : Type*} [Monoid M]

/-!
### §B – Universal Property of `FreeMonoid α`

`FreeMonoid α` is (up to isomorphism) the *unique* monoid generated freely
by the alphabet `α`.  Concretely: every function `φ : α → M` extends to a
unique `MonoidHom` `FreeMonoid α →* M`.

Mathlib provides `FreeMonoid.lift` for the construction and
`FreeMonoid.hom_eq` for uniqueness.
-/

/-- **Universal property of `FreeMonoid α`.**
    For every function `φ : α → M` there exists a *unique* `MonoidHom`
    `f : FreeMonoid α →* M` satisfying `f (FreeMonoid.of a) = φ a`. -/
theorem universal_property (φ : α → M) :
    ∃! (f : FreeMonoid α →* M), ∀ a : α, f (FreeMonoid.of a) = φ a := by
  refine ⟨FreeMonoid.lift φ, FreeMonoid.lift_eval_of φ, fun f hf => ?_⟩
  exact FreeMonoid.hom_eq fun a =>
    (hf a).trans (FreeMonoid.lift_eval_of φ a).symm

/-- `streamEvalHom φ` satisfies the universal property:
    on single-element streams it recovers the generator map `φ`. -/
theorem streamEvalHom_satisfies_universal (φ : α → M) (a : α) :
    streamEvalHom φ (FreeMonoid.of a) = φ a := by
  rw [streamEvalHom_apply,
      show (FreeMonoid.of a : FreeMonoid α) = [a] from rfl,
      streamEval_singleton]

end UniversalProperty

-- ===========================================================================
-- §C  EABC State Space as an AddMonoidHom
-- ===========================================================================

section EABCHom

/-!
### §C – EABC Zustandsraum als `AddMonoidHom`

`FreeAddMonoid PrimeParity` (= `List PrimeParity` with `++` as `+` and `[]`
as `0`) maps into `EABCState = ℤ × ℤ` via componentwise accumulation.
-/

/-- `PrimeParity` encodes the sign contribution of a prime to the EABC state:
    `pos` adds `(+1, 0)` and `neg` adds `(0, +1)`. -/
inductive PrimeParity : Type
  | pos : PrimeParity
  | neg : PrimeParity

/-- Map a parity bit to the corresponding EABC state contribution. -/
def parityToState : PrimeParity → EABCState
  | .pos => (1, 0)
  | .neg => (0, 1)

/-- `streamEvalAddHom φ` packages `streamEvalAdd φ` as an `AddMonoidHom` from
    `FreeAddMonoid α` (= `List α` with concatenation as `+`) into any
    `AddMonoid M`:
    - `map_zero'`: empty stream → zero element `0`.
    - `map_add'`: stream concatenation → additive sum. -/
def streamEvalAddHom {α : Type*} {M : Type*} [AddMonoid M]
    (φ : α → M) : FreeAddMonoid α →+ M where
  toFun    := streamEvalAdd φ
  map_zero' := streamEvalAdd_nil φ
  map_add' := AdditiveStreamingFactorization φ

/-- The EABC homomorphism: accumulates a `PrimeParity` stream into the EABC
    state `(N₊, N₋) ∈ ℤ × ℤ` via additiver Streaming-Faktorisierung. -/
def eabcHom : FreeAddMonoid PrimeParity →+ EABCState :=
  streamEvalAddHom parityToState

@[simp]
theorem eabcHom_apply (s : FreeAddMonoid PrimeParity) :
    eabcHom s = streamEvalAdd parityToState s := rfl

/-- `D` is determined by the state components. -/
theorem D_is_corollary (s : FreeAddMonoid PrimeParity) :
    D (eabcHom s) = (eabcHom s).1 - (eabcHom s).2 := rfl

/-- `Q` is determined by the state components. -/
theorem Q_is_corollary (s : FreeAddMonoid PrimeParity) :
    Q (eabcHom s) = (eabcHom s).1 + (eabcHom s).2 := rfl

end EABCHom

-- ===========================================================================
-- §D  W_E does not live in the Monoid
-- ===========================================================================

section WENotHom

/-!
### §D – `W_E` ist kein Monoid-Homomorphismus

`W_E = D/Q` is a *ratio* observable and therefore non-linear.
We show that no `AddMonoidHom f : FreeAddMonoid PrimeParity →+ ℚ` can
satisfy `f s = W_E (eabcHom s)` for all `s`.

**Counterexample**:

| stream        | `eabcHom` | `W_E` |
|---------------|-----------|-------|
| `[pos]`       | `(1, 0)`  |  1    |
| `[pos, pos]`  | `(2, 0)`  |  1    |

Additivity would force `f [pos, pos] = f [pos] + f [pos] = 1 + 1 = 2`,
but `W_E (eabcHom [pos, pos]) = 1`.  Contradiction: `2 = 1`.
-/

-- Helper: W_E on concrete states
private lemma W_E_one_zero : W_E (1, 0) = 1 := by
  unfold W_E D Q
  push_cast
  norm_num

private lemma W_E_two_zero : W_E (2, 0) = 1 := by
  unfold W_E D Q
  push_cast
  norm_num

-- Helper: eabcHom on concrete singleton and doubleton streams
private lemma eabcHom_pos_singleton :
    eabcHom ([.pos] : FreeAddMonoid PrimeParity) = (1, 0) := by
  simp [eabcHom_apply, streamEvalAdd_singleton, parityToState]

private lemma eabcHom_pos_doubleton :
    eabcHom ([.pos, .pos] : FreeAddMonoid PrimeParity) = (2, 0) := by
  show streamEvalAdd parityToState ([.pos] ++ [.pos]) = (2, 0)
  rw [AdditiveStreamingFactorization, streamEvalAdd_singleton, streamEvalAdd_singleton,
      show parityToState PrimeParity.pos = (1, 0) from rfl]
  ext <;> norm_num

/-- **`W_E` cannot be extended to an `AddMonoidHom`.**
    Proof by contradiction via the counterexample `s = [pos]`, `t = [pos]`. -/
theorem W_E_not_monoid_hom :
    ¬∃ (f : FreeAddMonoid PrimeParity →+ ℚ),
        ∀ s : FreeAddMonoid PrimeParity, f s = W_E (eabcHom s) := by
  intro ⟨f, hf⟩
  -- Step 1: the hypothesis pins f on the singleton stream.
  have h_single : f ([.pos] : FreeAddMonoid PrimeParity) = 1 := by
    rw [hf, eabcHom_pos_singleton, W_E_one_zero]
  -- Step 2: additivity of f forces f [pos, pos] = 2.
  have h_add : f ([.pos, .pos] : FreeAddMonoid PrimeParity) = 2 := by
    have : ([.pos, .pos] : FreeAddMonoid PrimeParity) =
           ([.pos] : FreeAddMonoid PrimeParity) + ([.pos] : FreeAddMonoid PrimeParity) := rfl
    rw [this, map_add, h_single]
    norm_num
  -- Step 3: the hypothesis also pins f [pos, pos] = 1.
  have h_pair : f ([.pos, .pos] : FreeAddMonoid PrimeParity) = 1 := by
    rw [hf, eabcHom_pos_doubleton, W_E_two_zero]
  -- Contradiction: 2 = 1.
  linarith

end WENotHom
