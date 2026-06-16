/-
  CollatzEabc.Open — Stufe E: dokumentierte offene Uniformitätsvermutung.

  Kein `sorry`, kein Beweis — nur Prop-Definitionen und logische Äquivalenzen.
  Siehe `collatz_offene_punkte.md` und `collatz_equivalenz_e_infty.tex` im Repo-Root.
-/

import CollatzEabc.Z2Attraktor

namespace CollatzZ2

/-!
### E_∞ vs. E_diag

- `ExceptionSetInfinity` (E_∞): echte Collatz-Ausnahmen in ℕ,
  `E_∞ = { n | Odd n ∧ ∀ K, iterateU n K ≠ 1 }`. Collatz ⇔ E_∞ = ∅.
- `ExceptionSetDiag` / `ExceptionSet` (E_diag): 2-adischer Beobachtungsschatten
  `closure (⋃_N E_{N,N})` in ℤ₂ — **nicht** Collatz-äquivalent (z. B. n = 27).
-/

/-- E_∞: ungerade Starts, die nie 1 erreichen (`collatz_equivalenz_e_infty.tex`). -/
def ExceptionSetInfinity : Set ℕ := { n | Odd n ∧ ∀ K, iterateU n K ≠ 1 }

theorem collatz_iff_exceptionSetInfinity_empty :
    (∀ n, Odd n → ∃ K, iterateU n K = 1) ↔ ExceptionSetInfinity = ∅ := by
  constructor
  · intro h
    by_contra hne
    rcases (Set.nonempty_iff_ne_empty.mpr hne) with ⟨n, hn⟩
    simp only [ExceptionSetInfinity, Set.mem_setOf_eq] at hn
    rcases h n hn.1 with ⟨K, heq⟩
    exact hn.2 K heq
  · intro hempty n ho
    by_contra hnot
    push Not at hnot
    have hn : n ∈ ExceptionSetInfinity := ⟨ho, hnot⟩
    rw [hempty] at hn
    exact hn

/-- Uniformitätsvermutung (OFFEN): E_∞ = ∅, äquivalent zur Collatz-Vermutung. -/
def collatzUniformityConjecture : Prop := ExceptionSetInfinity = ∅

/-- Diagonaler Schatten E_diag ⊂ ℤ₂ (Alias für `ExceptionSet`). -/
abbrev ExceptionSetDiag : Set Z2 := ExceptionSet

/-- Äquivalente Formulierung über endliche Attraktor-Approximation. -/
def collatzUniformityConjectureApprox : Prop :=
  ∀ n : ℕ, Odd n → ∃ K : ℕ, inTrivialAttraktorApprox n K

theorem collatzUniformityConjecture_iff_iterate :
    collatzUniformityConjecture ↔ (∀ n, Odd n → ∃ K, iterateU n K = 1) := by
  unfold collatzUniformityConjecture
  exact collatz_iff_exceptionSetInfinity_empty.symm

/-- Collatz für ungerade Starts (odd-to-odd-Formulierung). -/
theorem collatzUniformityConjecture_iff :
    collatzUniformityConjecture ↔ collatzUniformityConjectureApprox := by
  rw [collatzUniformityConjecture_iff_iterate]
  constructor
  · intro h n ho
    rcases h n ho with ⟨K, heq⟩
    exact ⟨K, ⟨K, Nat.le_refl _, heq⟩⟩
  · intro h n ho
    rcases h n ho with ⟨K, hK⟩
    rcases hK with ⟨k, _, heq⟩
    exact ⟨k, heq⟩

theorem collatzUniformityConjectureApprox_iff_iterate :
    collatzUniformityConjectureApprox ↔ collatzUniformityConjecture :=
  collatzUniformityConjecture_iff.symm

end CollatzZ2
