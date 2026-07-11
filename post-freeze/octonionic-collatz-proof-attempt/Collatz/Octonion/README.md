# Collatz / Octonion — Lean scaffold (post-freeze)

**Status:** placeholder for octonionic lift formalization. No proof files yet.

## Intended interface (Module O3)

```text
π_ℕ(F_O(Q(n))) = S(n)
```

where `S` is the accelerated odd Collatz map from `CollatzEabc` (repo: `collatz_eabc_core/`).

## Do not merge into audit-freeze

This directory is isolated under `post-freeze/octonionic-collatz-proof-attempt/`.
Sorry-free baseline remains in `collatz_eabc_core/`; new lemmas here may use `sorry` until Variante C is formalized.

## Next formalization targets

1. `OddCoreLogCocycle` — O1 as Mathlib `Real.log` identity over `ℕ`
2. `ValuationOneLadder` — O2 closed form for `n0 = 2^(L+1)-1`
3. `NoGoBoundedPotential` — O2 corollary as structured `Prop`
4. `OctonionLift` — scaffold only; no closure claim
