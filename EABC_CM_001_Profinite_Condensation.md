# EABC-CM-001: Profinite Condensation of the Four-Phase System

**Status:** Foundational specification · Layer [B→C] · No analytic or physical claims

## 0. Purpose

This module gives the finite EABC residue system a functorial and profinite formulation. It does not prove new statements about prime distributions, Kerr geometry, Collatz dynamics, or spectral asymptotics. Its role is to define the mathematical object on which such later structures may act.

The core object is the finite set
\[
Q_{\mathrm{EABC}}=\{E,A,B,C\},
\]
identified with the reduced residue classes modulo 12 by
\[
E\leftrightarrow 1,\qquad A\leftrightarrow 5,\qquad B\leftrightarrow 7,\qquad C\leftrightarrow 11.
\]

For each profinite test space \(S\), the EABC states over \(S\) are the continuous maps
\[
\mathcal E(S)=\operatorname{Cont}(S,Q_{\mathrm{EABC}}).
\]

Because \(Q_{\mathrm{EABC}}\) is finite and discrete, these maps are exactly the locally constant EABC labellings of \(S\).

---

## 1. Claim wall

This module establishes only:

1. the finite EABC state object;
2. its interpretation as a representable condensed set;
3. the functorial pullback of EABC states;
4. a natural four-phase labelling;
5. finite-quotient factorisation of each individual state;
6. descent along finite jointly surjective covers;
7. the non-uniformity boundary between pointwise finite factorisation and a universal finite resolution.

It does **not** establish:

- any statistical statement about primes;
- any spectral gap;
- any Mellin continuation;
- any physical interpretation;
- any \(E_8\)-regularisation;
- any equivalence with Kerr geometry;
- any uniform finite quotient through which all EABC states factor.

---

## 2. Two distinct four-element structures

The reduced unit group modulo 12 is
\[
U_{12}=(\mathbb Z/12\mathbb Z)^\times=\{1,5,7,11\}.
\]
Under multiplication modulo 12, every nontrivial element has order two:
\[
5^2\equiv7^2\equiv11^2\equiv1\pmod{12}.
\]
Hence
\[
U_{12}\cong C_2\times C_2,
\]
not \(C_4\).

The cyclic order
\[
C\to E\to A\to B\to C
\]
and the phase assignment
\[
C\mapsto1,\qquad E\mapsto i,\qquad A\mapsto-1,\qquad B\mapsto-i
\]
therefore define an additional **cyclic phase labelling**, not a group character of \(U_{12}\).

This distinction is structural:

- \(U_{12}\) carries the arithmetic multiplication law;
- the imposed \(C_4\)-order carries the EABC phase geometry.

They coexist on the same four-element set but must not be identified.

---

## 3. The finite EABC object

### Definition 3.1 — EABC state set

Let
\[
Q:=\{E,A,B,C\}
\]
with the discrete topology.

### Definition 3.2 — Residue realisation

Define the bijection
\[
\rho:Q\xrightarrow{\sim}U_{12}
\]
by
\[
\rho(E)=1,\quad \rho(A)=5,\quad \rho(B)=7,\quad \rho(C)=11.
\]

Transporting multiplication along \(\rho\) equips \(Q\) with a Klein-four-group structure.

### Definition 3.3 — Phase code

Define
\[
\phi:Q\to\mathbb Z/4\mathbb Z
\]
by
\[
\phi(C)=0,\qquad \phi(E)=1,\qquad \phi(A)=2,\qquad \phi(B)=3.
\]

### Definition 3.4 — Complex phase map

Let
\[
\vartheta:Q\to\mu_4=\{1,i,-1,-i\}
\]
be
\[
\vartheta(q)=i^{\phi(q)}.
\]

Thus
\[
\vartheta(C)=1,\quad
\vartheta(E)=i,\quad
\vartheta(A)=-1,\quad
\vartheta(B)=-i.
\]

The map \(\vartheta\) is a bijection of finite sets. It is not a homomorphism from the multiplicative group \(U_{12}\cong C_2\times C_2\) to \(\mu_4\cong C_4\).

---

## 4. Profinite test spaces

Let \(\mathbf{ProFin}\) denote the category of profinite spaces and continuous maps.

### Definition 4.1 — EABC states over a test space

For \(S\in\mathbf{ProFin}\), define
\[
\mathcal E(S):=\operatorname{Cont}(S,Q).
\]

For a continuous map \(f:T\to S\), define pullback
\[
f^*:\mathcal E(S)\to\mathcal E(T),
\qquad
f^*(x)=x\circ f.
\]

Hence
\[
\mathcal E:\mathbf{ProFin}^{\mathrm{op}}\to\mathbf{Set}
\]
is a contravariant functor.

### Proposition 4.2 — Functoriality

For continuous maps \(U\xrightarrow{g}T\xrightarrow{f}S\),
\[
(f\circ g)^*=g^*\circ f^*,
\]
and
\[
(\operatorname{id}_S)^*=\operatorname{id}_{\mathcal E(S)}.
\]

**Proof.** Both identities follow from associativity and unitality of composition. ∎

---

## 5. Condensed interpretation

The functor \(\mathcal E\) is the representable presheaf
\[
h_Q(S)=\operatorname{Cont}(S,Q).
\]
Representable presheaves on the profinite site satisfy the sheaf condition for the condensed topology. Therefore \(\mathcal E\) is a condensed set.

### Theorem 5.1 — Representable condensation

The EABC functor
\[
\mathcal E=h_Q
\]
is a condensed set represented by the finite discrete profinite space \(Q\).

**Interpretation.** The elementary four-state set is the value at a point,
\[
\mathcal E(*)\cong Q,
\]
but the complete object consists of all continuously parametrised EABC families \(\mathcal E(S)\).

---

## 6. Natural phase transformation

Define the phase functor
\[
\mathcal P(S):=\operatorname{Cont}(S,\mu_4).
\]

For every profinite \(S\), define
\[
\Theta_S:\mathcal E(S)\to\mathcal P(S),
\qquad
\Theta_S(x)=\vartheta\circ x.
\]

### Theorem 6.1 — Naturality

The family \((\Theta_S)_S\) defines a natural transformation
\[
\Theta:\mathcal E\Rightarrow\mathcal P.
\]
That is, for every continuous \(f:T\to S\),
\[
f^*\circ\Theta_S=\Theta_T\circ f^*.
\]

**Proof.** For \(x\in\mathcal E(S)\),
\[
(f^*\circ\Theta_S)(x)
=(\vartheta\circ x)\circ f
=\vartheta\circ(x\circ f)
=(\Theta_T\circ f^*)(x).
\]
∎

### Corollary 6.2 — Phase coding commutes with restriction

Every restriction, localisation, or finite pullback of an EABC state may be phase-coded before or after pullback with the same result.

---

## 7. Finite-quotient factorisation

### Theorem 7.1 — Individual finite-resolution theorem

Let \(S\) be profinite and let
\[
x:S\to Q
\]
be continuous. Then there exist

- a finite discrete set \(S_x\),
- a continuous surjection \(\pi_x:S\twoheadrightarrow S_x\),
- a map \(\bar x:S_x\to Q\),

such that
\[
x=\bar x\circ\pi_x.
\]

**Proof.** Since \(Q\) is finite discrete, each fibre \(x^{-1}(q)\) is clopen. The nonempty fibres form a finite clopen partition of \(S\). Let \(S_x\) be the finite set of blocks of this partition and let \(\pi_x\) send a point to its block. Then \(x\) is constant on each block and induces \(\bar x\). ∎

### Corollary 7.2 — Finite visibility

Every profinite EABC state is already visible on some finite quotient.

The quantifier order is
\[
\boxed{\forall x\in\mathcal E(S)\;\exists S_x\text{ finite such that }x\text{ factors through }S_x.}
\]

This theorem does not produce one quotient valid for all states.

---

## 8. The non-uniformity wall

### Proposition 8.1 — No universal finite resolution on an infinite profinite source

Let \(S\) be an infinite profinite space with arbitrarily fine nontrivial clopen partitions, and let \(|Q|\ge2\). Then there is generally no fixed finite quotient
\[
\pi:S\twoheadrightarrow S_0
\]
through which every continuous map \(S\to Q\) factors.

**Proof idea.** A fixed finite quotient sees only unions of its finitely many fibres. Since \(S\) has strictly finer clopen partitions, choose a clopen subset \(U\subset S\) that is not a union of \(\pi\)-fibres. Choose distinct \(q_0,q_1\in Q\) and define
\[
x_U(s)=
\begin{cases}
q_1,&s\in U,\\
q_0,&s\notin U.
\end{cases}
\]
Then \(x_U\) is continuous but cannot factor through \(\pi\). ∎

### Explicit form for \(\widehat{\mathbb Z}\)

For every fixed modulus \(m\), choose a prime \(p\nmid m\). The clopen condition
\[
n\equiv0\pmod p
\]
defines a continuous two-valued EABC state on \(\widehat{\mathbb Z}\) that does not factor through reduction modulo \(m\).

Hence
\[
\boxed{
\forall x\;\exists m(x)
}
\]
does not imply
\[
\boxed{
\exists m\;\forall x.
}
\]

This is the condensed analogue of the existing EABC freeze boundary
\[
\forall F\,\exists N_0(F)
\quad\not\Rightarrow\quad
\exists N_0\,\forall F.
\]

---

## 9. Descent

Let
\[
p:S'\twoheadrightarrow S
\]
be a continuous surjection of profinite spaces.

The two projections from the fibre product are
\[
\operatorname{pr}_1,\operatorname{pr}_2:S'\times_S S'\rightrightarrows S'.
\]

### Theorem 9.1 — Effective descent for EABC states

The diagram
\[
\mathcal E(S)
\longrightarrow
\mathcal E(S')
\rightrightarrows
\mathcal E(S'\times_S S')
\]
is an equalizer.

Equivalently, an EABC state \(x'\in\mathcal E(S')\) descends uniquely to \(S\) iff
\[
x'\circ\operatorname{pr}_1
=
x'\circ\operatorname{pr}_2.
\]

**Proof.** If \(x'=x\circ p\), compatibility is immediate. Conversely, compatibility says that \(x'\) is constant on every fibre of \(p\). Therefore it defines a unique set map \(x:S\to Q\). Since \(p\) is a quotient map and \(x\circ p=x'\) is continuous, \(x\) is continuous. ∎

### Interpretation

Local EABC data may be promoted to a global state only after overlap compatibility has been verified. In project governance terms:

\[
\boxed{
\text{A local claim crosses the claim wall only after descent data are checked.}
}
\]

---

## 10. Arithmetic source object

There is a canonical continuous projection
\[
\pi_{12}:\widehat{\mathbb Z}^{\,\times}\to U_{12}.
\]
Composing with \(\rho^{-1}\) gives the universal mod-12 EABC labelling
\[
\varepsilon_{12}:\widehat{\mathbb Z}^{\,\times}\to Q.
\]

Its phase lift is
\[
\theta_{12}
=
\vartheta\circ\varepsilon_{12}:
\widehat{\mathbb Z}^{\,\times}\to\mu_4.
\]

This object records only the mod-12 EABC class. Any finer arithmetic structure must be supplied by a larger quotient, a tower of quotients, or an additional sheaf/module over the same test space.

---

## 11. Fourier layer

Let
\[
V_Q=\mathbb C^Q\cong\mathbb C^4.
\]
Choose the cyclic phase order
\[
(C,E,A,B).
\]
The associated unnormalised Fourier matrix is
\[
F_4=
\begin{pmatrix}
1&1&1&1\\
1&i&-1&-i\\
1&-1&1&-1\\
1&-i&-1&i
\end{pmatrix}.
\]

The normalised transform
\[
U_4=\frac12F_4
\]
is unitary.

This Fourier structure belongs to the imposed cyclic phase labelling, not to the multiplicative Klein-four-group law on \(U_{12}\).

For every profinite \(S\), it acts pointwise on continuous families
\[
\operatorname{Cont}(S,V_Q)
\to
\operatorname{Cont}(S,V_Q).
\]
Thus the finite Fourier transform extends naturally to all profinite EABC families.

---

## 12. Minimal theorem stack

The module is considered mathematically complete once the following stack is proved:

1. `eabcState_equiv_unitsMod12`
2. `phaseCode_bijective`
3. `phaseMap_not_groupHom`
4. `eabcPullback_id`
5. `eabcPullback_comp`
6. `phase_natural`
7. `continuous_eabc_factors_finite`
8. `eabc_descent_equalizer`
9. `no_uniform_finite_resolution_hatZ`
10. `fourier_unitary`

---

## 13. Lean-oriented finite core

The finite core can be implemented independently of condensed-mathematics APIs:

```lean
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Fintype.Card
import Mathlib.Data.Complex.Exponential
import Mathlib.LinearAlgebra.Matrix.Fourier

namespace EABC

inductive State
  | E | A | B | C
  deriving DecidableEq, Fintype, Repr

open State

/-- Residue representative modulo 12. -/
def residue : State → ZMod 12
  | E => 1
  | A => 5
  | B => 7
  | C => 11

/-- Imposed cyclic phase code; this is not the multiplicative law of `(ZMod 12)ˣ`. -/
def phaseCode : State → ZMod 4
  | C => 0
  | E => 1
  | A => 2
  | B => 3

/-- Concrete fourth-root phase. -/
def phase : State → ℂ
  | C => 1
  | E => Complex.I
  | A => -1
  | B => -Complex.I

@[simp] theorem residue_E : residue E = 1 := rfl
@[simp] theorem residue_A : residue A = 5 := rfl
@[simp] theorem residue_B : residue B = 7 := rfl
@[simp] theorem residue_C : residue C = 11 := rfl

@[simp] theorem phase_C : phase C = 1 := rfl
@[simp] theorem phase_E : phase E = Complex.I := rfl
@[simp] theorem phase_A : phase A = -1 := rfl
@[simp] theorem phase_B : phase B = -Complex.I := rfl

/-- Pointwise pullback of an EABC field. -/
def pullback {S T : Type*} (f : T → S) (x : S → State) : T → State :=
  x ∘ f

@[simp] theorem pullback_id {S : Type*} (x : S → State) :
    pullback id x = x := by
  funext s
  rfl

@[simp] theorem pullback_comp {R S T : Type*}
    (f : S → R) (g : T → S) (x : R → State) :
    pullback g (pullback f x) = pullback (f ∘ g) x := by
  rfl

/-- Phase coding commutes with pullback. -/
theorem phase_natural {S T : Type*} (f : T → S) (x : S → State) :
    phase ∘ pullback f x = pullback f (phase ∘ x) := by
  rfl

end EABC
```

The topological/condensed layer should then be added in a separate file so that the finite kernel remains stable if the category API changes.

---

## 14. Lean-oriented topological scaffold

A later file should introduce:

```lean
namespace EABC.Condensed

-- Intended object:
-- E(S) = ContinuousMap S EABC.State

-- pullback along a continuous map
-- natural phase transformation
-- finite clopen quotient generated by the fibres of a state
-- descent along a continuous surjection

end EABC.Condensed
```

Recommended separation:

```text
KeplerHurwitz/EABC/FiniteState.lean
KeplerHurwitz/EABC/PhaseFourier.lean
KeplerHurwitz/EABC/ProfiniteState.lean
KeplerHurwitz/EABC/CondensedDescent.lean
```

The exact imports and category-level declarations should be fixed against the repository's current Mathlib revision before claiming build-green status.

---

## 15. Governance status

Suggested status string:

```text
EABC-CM-001 · PROFINITE FOUR-PHASE FOUNDATION
DEFINITIONS LOCKED · FINITE CORE READY FOR FORMALISATION
CONDENSED/DESCENT LAYER SPECIFIED · NOT YET BUILD-VERIFIED
NO ANALYTIC, PRIME-DISTRIBUTION, OR PHYSICAL CLAIMS
```

---

## 16. Next module

The natural successor is:

### EABC-CM-002 — Condensed Complement Exact Sequence

Target sequence:
\[
0\to\mathcal V_{2,3}
\to\mathcal A
\to\mathcal E
\to0.
\]

Before this can be stated, the objects \(\mathcal A\), \(\mathcal V_{2,3}\), and the quotient morphism must be defined as condensed abelian groups or condensed modules. The word “split” should be reserved for the existence of a section
\[
s:\mathcal E\to\mathcal A,
\qquad
\pi\circ s=\operatorname{id}_{\mathcal E}.
\]

Until such a section is constructed, the current complement remains a set-theoretic or residue-theoretic decomposition, not a categorical direct-sum theorem.

---

## 17. Final theorem-level summary

The foundational result of EABC-CM-001 is:

\[
\boxed{
\mathcal E(S)=\operatorname{Cont}(S,Q_{\mathrm{EABC}})
}
\]
forms a representable condensed set, the imposed four-phase code defines a natural transformation
\[
\Theta:\mathcal E\Rightarrow\underline{\mu_4},
\]
every individual EABC state factors through a finite quotient, and compatible local states descend uniquely.

The associated non-uniformity wall is
\[
\boxed{
\forall x\;\exists S_x
\quad\not\Rightarrow\quad
\exists S_0\;\forall x.
}
\]

This is the exact profinite counterpart of the established EABC quantifier boundary.
