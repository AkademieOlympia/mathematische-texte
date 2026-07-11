# Post-Freeze Branch: Octonionic Collatz Proof Attempt

```
Status: Layer 5 research program (post-freeze).
Not part of audit-freeze-2026.
No Collatz proof claim.
Outside Synthesis_2026 / Bridges2026 audit dossier.
```

## Scope

| Module | Content | Status |
|---|---|---|
| **O1** | Odd-core log cocycle | Verified invariant (`verify_o1_o2.py`) |
| **O2** | Long $\nu_2=1$ runs at $n_0=2^{L+1}-1$ ($L$ steps, $j=0..L-1$) | No-go lemma (`verify_o1_o2.py`) |
| **O3** | Octonionic lift schema | Scaffold only |
| **O4–O5** | Compensation bridge, Variante C | Open boundary |

## Relation to `collatz_eabc_core`

Sorry-free Lean package at repo root: `collatz_eabc_core/` (density, uniformity, ℤ₂ attractor).
This branch adds **post-freeze** octonionic lift scaffolding under `Collatz/Octonion/` — no proof closure.

## Reproduce (O1/O2)

```bash
python3 post-freeze/octonionic-collatz-proof-attempt/verify_o1_o2.py
python3 post-freeze/octonionic-collatz-proof-attempt/verify_o1_o2.py --L-max 64
```

## Dossier

`notes/octonionic_collatz_proof_attempt_v1.md`
