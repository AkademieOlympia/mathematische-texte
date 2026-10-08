import Mathlib

/-!
# Ψ₉: EABC–Walsh–Galois arithmetic

A robust formal core for the nine-component arithmetic discussed in the
EABC / Walsh / Galois construction.

This file deliberately separates:

* purely linear/algebraic statements, which are proved here;
* the arithmetic encoding of a natural number via its prime factorization,
  which should be added as a later layer;
* the finite modulo-12 Goldbach switching rule, proved here on the residue
  classes {1,5,7,11}.

No Goldbach existence statement is assumed.
-/

namespace Psi9

noncomputable section

/-- The golden ratio φ. -/
def phi : ℝ := (1 + Real.sqrt 5) / 2

/-- A four-component EABC vector. -/
abbrev EABCVec := Fin 4 → ℝ

/-- Nine-component container: one scale coordinate and two four-vectors. -/
structure Vector9D where
  tau : ℝ
  U   : EABCVec
  S   : EABCVec

deriving Repr

instance : Zero Vector9D where
  zero := ⟨0, 0, 0⟩

instance : Add Vector9D where
  add v w := ⟨v.tau + w.tau, v.U + w.U, v.S + w.S⟩

/-!
## Walsh transform

We define the four Walsh coordinates explicitly.  This is equivalent to
multiplication by

  1/2 * [[ 1,  1,  1,  1],
         [ 1, -1,  1, -1],
         [ 1,  1, -1, -1],
         [ 1, -1, -1,  1]].

The explicit definition makes the basic proofs stable under Mathlib changes.
-/

/-- Walsh-Hadamard transform W₄ x. -/
def walshVector (x : EABCVec) : EABCVec :=
  ![
    (x 0 + x 1 + x 2 + x 3) / 2,
    (x 0 - x 1 + x 2 - x 3) / 2,
    (x 0 + x 1 - x 2 - x 3) / 2,
    (x 0 - x 1 - x 2 + x 3) / 2
  ]

/-- Ur-vector U = φ s. -/
def urVector (s : EABCVec) : EABCVec :=
  phi • s

/-- Shadow-vector S = (1 - φ) s. -/
def shadowVector (s : EABCVec) : EABCVec :=
  (1 - phi) • s

/-- Linear Ψ₉-construction from an external scale τ and an EABC vector x. -/
def psi9Of (tau : ℝ) (x : EABCVec) : Vector9D :=
  let s := walshVector x
  ⟨tau, urVector s, shadowVector s⟩

--------------------------------------------------------------------------------
-- SATZ 1: Linear homomorphy and Galois decoupling
--------------------------------------------------------------------------------

/-- Walsh transform is additive. -/
theorem walsh_add (x y : EABCVec) :
    walshVector (x + y) = walshVector x + walshVector y := by
  ext i
  fin_cases i <;>
    simp [walshVector] <;>
    ring

/-- The linear Ψ₉ construction is additive. -/
theorem psi9Of_add (tau₁ tau₂ : ℝ) (x y : EABCVec) :
    psi9Of (tau₁ + tau₂) (x + y) = psi9Of tau₁ x + psi9Of tau₂ y := by
  cases h₁ : psi9Of tau₁ x with
  | mk t₁ U₁ S₁ =>
    cases h₂ : psi9Of tau₂ y with
    | mk t₂ U₂ S₂ =>
      simp [psi9Of, walsh_add, urVector, shadowVector] at h₁ h₂ ⊢
      ext i <;> simp [walshVector, urVector, shadowVector] <;> ring

/-- Exact Galois decoupling: U - S = √5 · s. -/
theorem galois_decoupling (s : EABCVec) :
    urVector s - shadowVector s = (Real.sqrt 5) • s := by
  ext i
  simp [urVector, shadowVector, phi]
  ring

/-- The first Ur/Shadow construction is linearly dependent:
    φ S = (1-φ) U.  Thus the 4+4 block has real rank at most four. -/
theorem galois_linear_dependence (s : EABCVec) :
    phi • shadowVector s = (1 - phi) • urVector s := by
  ext i
  simp [urVector, shadowVector]
  ring

--------------------------------------------------------------------------------
-- SATZ 2: Neutrality of the special factors 2 and 3 at EABC level
--------------------------------------------------------------------------------

/-- If the EABC vector vanishes, both Ur and Shadow sectors vanish. -/
theorem zero_eabc_is_neutral :
    let s := walshVector (0 : EABCVec)
    urVector s = 0 ∧ shadowVector s = 0 := by
  simp [walshVector, urVector, shadowVector]

/-- A version with an arbitrary vector known to be zero. -/
theorem neutral_special_factors
    (x23 : EABCVec) (hzero : x23 = 0) :
    let s := walshVector x23
    urVector s = 0 ∧ shadowVector s = 0 := by
  subst x23
  simpa using zero_eabc_is_neutral

--------------------------------------------------------------------------------
-- SATZ 3: Isotropy collapse of an equal EABC occupation
--------------------------------------------------------------------------------

/-- Constant EABC occupation is a pure Walsh singlet. -/
theorem walsh_constant_vector (k : ℝ) :
    walshVector (fun _ : Fin 4 => k) = ![2 * k, 0, 0, 0] := by
  ext i
  fin_cases i <;>
    simp [walshVector] <;>
    ring

/-- Equal EABC occupation kills all three anisotropy channels. -/
theorem quadruplet_isotropy_collapse (k : ℝ) :
    let x : EABCVec := fun _ => k
    let s := walshVector x
    s 1 = 0 ∧ s 2 = 0 ∧ s 3 = 0 := by
  dsimp
  have hs := walsh_constant_vector k
  rw [hs]
  simp

--------------------------------------------------------------------------------
-- SATZ 4: Finite modulo-12 Goldbach/Walsh switch
--------------------------------------------------------------------------------

/-- The s_x Walsh character on the four units modulo 12.

    E=1 and B=7 have character +1.
    A=5 and C=11 have character -1.
-/
def sxResidue (r : ℕ) : ℝ :=
  if r % 12 = 1 ∨ r % 12 = 7 then 1
  else if r % 12 = 5 ∨ r % 12 = 11 then -1
  else 0

/-- Normalized s_x value of a pair of residue classes. -/
def sxPairResidues (r s : ℕ) : ℝ :=
  (sxResidue r + sxResidue s) / 2

/-- Deterministic s_x switch as a function of the even sum modulo 12. -/
def goldbachSxSwitch (Nmod12 : ℕ) : ℝ :=
  if Nmod12 % 12 = 2 ∨ Nmod12 % 12 = 8 then 1
  else if Nmod12 % 12 = 4 ∨ Nmod12 % 12 = 10 then -1
  else 0

/-- Pure finite residue theorem.

For r,s ∈ {1,5,7,11}, the normalized Walsh s_x value of the pair is
completely determined by (r+s) mod 12:

  2,8  ↦ +1
  4,10 ↦ -1
  0,6  ↦ 0.

This theorem is independent of the Goldbach conjecture; it only concerns
the possible residue classes of primes > 3.
-/
theorem goldbach_walsh_residue_switching
    (r s : ℕ)
    (hr : r = 1 ∨ r = 5 ∨ r = 7 ∨ r = 11)
    (hs : s = 1 ∨ s = 5 ∨ s = 7 ∨ s = 11) :
    sxPairResidues r s = goldbachSxSwitch ((r + s) % 12) := by
  rcases hr with rfl | rfl | rfl | rfl <;>
  rcases hs with rfl | rfl | rfl | rfl <;>
  norm_num [sxPairResidues, sxResidue, goldbachSxSwitch]

/-- The same finite rule, phrased for residues already reduced modulo 12. -/
theorem goldbach_walsh_switch_table :
    sxPairResidues 1 1 = 1 ∧
    sxPairResidues 1 5 = 0 ∧
    sxPairResidues 1 7 = 1 ∧
    sxPairResidues 1 11 = 0 ∧
    sxPairResidues 5 5 = -1 ∧
    sxPairResidues 5 7 = 0 ∧
    sxPairResidues 5 11 = -1 ∧
    sxPairResidues 7 7 = 1 ∧
    sxPairResidues 7 11 = 0 ∧
    sxPairResidues 11 11 = -1 := by
  norm_num [sxPairResidues, sxResidue]

end

end Psi9
