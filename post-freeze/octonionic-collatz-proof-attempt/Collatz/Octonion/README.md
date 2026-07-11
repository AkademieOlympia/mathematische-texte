# Collatz / Octonion — Lean scaffold (post-freeze)

**Status:** Module O1 multiplicative cocycle **proved** (sorry-free); log cocycle as `Prop` + Python audit.

## Build

```bash
cd post-freeze/octonionic-collatz-proof-attempt
lake update   # first time only
lake build
grep -r "sorry" Collatz/   # should be empty
```

Depends on `../../collatz_eabc_core` (`CollatzZ2.collatzU`, `iterateU`).

## Module O1 (`OddCoreLogCocycle.lean`)

| Symbol | Status |
|---|---|
| `collatzU_mul_two_pow` | **Theorem** |
| `iterateU_succ_eq` | **Theorem** |
| `iterateU_cocycle_mulStatement` | **Prop** (multiplicative O1; Python audit) |
| `oddCoreLogCocycleStatement` | **Prop** (log form; Python audit) |

## Module O3 interface

```text
π_ℕ(F_O(Q(n))) = S(n)   — scaffold in VariantC.lean
```

## Module O4/O5

`VariantC.lean`: open `Prop` targets (Variante C, compensation bridge).

No Collatz proof closure.
