/-
  CollatzEabc.Primvierling.Mod12 — mod-12-Klassifikation markierter Primvierlinge (Phase 2).

  Satz A–C aus `Projektionszeuge_Primvierling.tex` (Ebene A): echte Primzahl-Vierlinge,
  nicht bloße Ringarithmetik auf beliebigen Nat-Startwerten.

  **Phase 1:** Dieses Modul blieb absichtlich außerhalb von `Core.lean`, solange der
  Beweis von Satz A noch `sorry` trug — ein falscher Satz mit nur `Nat.Prime p` wäre
  in den Build geraten. Jetzt sorry-frei und importierbar.

  Siebstatistiken (`sieve_statistics` in Python) verifizieren die Implementierung nur;
  sie entdecken die Signaturen ABCE/CEAB nicht — das liefert die Bijektion
  `chiralityWord` ↔ `{5, 11}` (Satz B, Reverse).
-/

import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Tactic
import CollatzEabc.Primvierling.Basic
import CollatzEabc.Primvierling.Chirality

namespace CollatzEabc
namespace Primvierling

private lemma prime_gt_three_not_div_three {p : ℕ} (hp : Nat.Prime p) (h : 3 < p) :
    p % 3 ≠ 0 := by
  intro h0
  have h3dvd : 3 ∣ p := Nat.dvd_of_mod_eq_zero h0
  rcases hp.eq_one_or_self_of_dvd 3 h3dvd with h1 | hp3
  · omega
  · omega

private lemma prime_gt_three_odd {p : ℕ} (hp : Nat.Prime p) (h : 3 < p) : p % 2 = 1 := by
  rcases Nat.Prime.eq_two_or_odd hp with rfl | hodd
  · omega
  · exact hodd

private lemma mod12_offsets_of_five (p : ℕ) (h : p % 12 = 5) :
    (p + 2) % 12 = 7 ∧ (p + 6) % 12 = 11 ∧ (p + 8) % 12 = 1 := by
  omega

private lemma mod12_offsets_of_eleven (p : ℕ) (h : p % 12 = 11) :
    (p + 2) % 12 = 1 ∧ (p + 6) % 12 = 5 ∧ (p + 8) % 12 = 7 := by
  omega

theorem mod12_eq_five_of_chiralityWord_ABCE {p : ℕ}
    (h : chiralityWord p = some Chirality.ABCE) : p % 12 = 5 := by
  by_cases hp5 : p % 12 = 5
  · exact hp5
  · exfalso
    by_cases hp11 : p % 12 = 11
    · simp [chiralityWord, hp5, hp11] at h
    · simp [chiralityWord, hp5, hp11] at h

theorem mod12_eq_eleven_of_chiralityWord_CEAB {p : ℕ}
    (h : chiralityWord p = some Chirality.CEAB) : p % 12 = 11 := by
  by_cases hp11 : p % 12 = 11
  · exact hp11
  · exfalso
    by_cases hp5 : p % 12 = 5
    · simp [chiralityWord, hp5, hp11] at h
    · simp [chiralityWord, hp5, hp11] at h

/-- Satz A: zulässige Startrestklassen eines markierten Primvierlings.

  Beweisidee (Tao/Maynard-honest): `p` und `p+2` sind prim und `> 3`, also
  `p ≢ 0,1 (mod 3)`; mit `p` ungerade folgt `p ≡ 5 (mod 6)`, hence `p ≡ 5 ∨ 11 (mod 12)`. -/
theorem prime_quadruplet_start_mod12 (p : ℕ) (hq : IsPrimeQuadruplet p) (h : 3 < p) :
    p % 12 = 5 ∨ p % 12 = 11 := by
  rcases hq with ⟨hp, hp2, _, _⟩
  have hp3 := prime_gt_three_not_div_three hp h
  have hp3ne1 : p % 3 ≠ 1 := by
    intro h1
    have h02 : (p + 2) % 3 = 0 := by omega
    have h3dvd : 3 ∣ p + 2 := Nat.dvd_of_mod_eq_zero h02
    rcases hp2.eq_one_or_self_of_dvd 3 h3dvd with h1' | hp2eq
    · omega
    · omega
  have hp3two : p % 3 = 2 := by omega
  have hpodd := prime_gt_three_odd hp h
  omega

/-- Alias für den früheren Stub-Namen (Hypothesen jetzt korrekt: Primvierling). -/
theorem prime_gt3_mod12 (p : ℕ) (hq : IsPrimeQuadruplet p) (h : 3 < p) :
    p % 12 = 5 ∨ p % 12 = 11 :=
  prime_quadruplet_start_mod12 p hq h

private theorem quadrupletResidues_mod5 (p : ℕ) (h : p % 12 = 5) :
    quadrupletResidues p = [5, 7, 11, 1] := by
  rcases mod12_offsets_of_five p h with ⟨h2, h6, h8⟩
  rw [quadrupletResidues_eq, h, h2, h6, h8]

private theorem quadrupletResidues_mod11 (p : ℕ) (h : p % 12 = 11) :
    quadrupletResidues p = [11, 1, 5, 7] := by
  rcases mod12_offsets_of_eleven p h with ⟨h2, h6, h8⟩
  rw [quadrupletResidues_eq, h, h2, h6, h8]

private theorem quadrupletFlavors_mod5 (p : ℕ) (h : p % 12 = 5) :
    quadrupletFlavors p =
      [some EClass.A, some EClass.B, some EClass.C, some EClass.E] := by
  rcases mod12_offsets_of_five p h with ⟨h2, h6, h8⟩
  simp [quadrupletFlavors, Q, classOf, h, h2, h6, h8]

private theorem quadrupletFlavors_mod11 (p : ℕ) (h : p % 12 = 11) :
    quadrupletFlavors p =
      [some EClass.C, some EClass.E, some EClass.A, some EClass.B] := by
  rcases mod12_offsets_of_eleven p h with ⟨h2, h6, h8⟩
  simp [quadrupletFlavors, Q, classOf, h, h2, h6, h8]

/-- Satz B: Startrestklasse bestimmt das Chiralitätswort. -/
theorem chiralityWord_of_prime_quadruplet (p : ℕ) (hq : IsPrimeQuadruplet p) (h : 3 < p) :
    chiralityWord p = some Chirality.ABCE ∨ chiralityWord p = some Chirality.CEAB := by
  rcases prime_quadruplet_start_mod12 p hq h with h5 | h11
  · left; simp [chiralityWord, h5]
  · right; simp [chiralityWord, h11]

theorem chiralityWord_eq_ABCE_of_mod5 (p : ℕ) (h : p % 12 = 5) :
    chiralityWord p = some Chirality.ABCE := by
  simp [chiralityWord, h]

theorem chiralityWord_eq_CEAB_of_mod11 (p : ℕ) (h : p % 12 = 11) :
    chiralityWord p = some Chirality.CEAB := by
  simp [chiralityWord, h]

/-- Reverse von Satz B: ABCE bestimmt die Startrestklasse. -/
theorem chiralityWord_eq_ABCE_iff_mod5 (p : ℕ) :
    chiralityWord p = some Chirality.ABCE ↔ p % 12 = 5 :=
  ⟨mod12_eq_five_of_chiralityWord_ABCE, chiralityWord_eq_ABCE_of_mod5 p⟩

/-- Reverse von Satz B: CEAB bestimmt die Startrestklasse. -/
theorem chiralityWord_eq_CEAB_iff_mod11 (p : ℕ) :
    chiralityWord p = some Chirality.CEAB ↔ p % 12 = 11 :=
  ⟨mod12_eq_eleven_of_chiralityWord_CEAB, chiralityWord_eq_CEAB_of_mod11 p⟩

/-- Bijektion `{5, 11} ↔ {ABCE, CEAB}` via `chiralityWord` (für beliebiges `p`). -/
theorem chiralityWord_bijection (p : ℕ) :
    (chiralityWord p = some Chirality.ABCE ↔ p % 12 = 5) ∧
    (chiralityWord p = some Chirality.CEAB ↔ p % 12 = 11) :=
  ⟨chiralityWord_eq_ABCE_iff_mod5 p, chiralityWord_eq_CEAB_iff_mod11 p⟩

/-- Für markierte Primvierlinge (`p > 3`): Bijektion Startrestklasse ↔ Chiralitätswort. -/
theorem chiralityWord_prime_quadruplet_bijection (p : ℕ) (hq : IsPrimeQuadruplet p) (h : 3 < p) :
    (chiralityWord p = some Chirality.ABCE ↔ p % 12 = 5) ∧
    (chiralityWord p = some Chirality.CEAB ↔ p % 12 = 11) ∧
    (chiralityWord p = some Chirality.ABCE ∨ chiralityWord p = some Chirality.CEAB) := by
  refine ⟨chiralityWord_eq_ABCE_iff_mod5 p, chiralityWord_eq_CEAB_iff_mod11 p, ?_⟩
  exact chiralityWord_of_prime_quadruplet p hq h

/-- Satz C: ABCE und CEAB sind die einzigen chiralen Vierlings-Signaturen (markierter Start). -/
theorem only_chiral_quadruplet_signatures (p : ℕ) (hq : IsPrimeQuadruplet p) (h : 3 < p) :
    quadrupletFlavors p =
        [some EClass.A, some EClass.B, some EClass.C, some EClass.E] ∨
      quadrupletFlavors p =
        [some EClass.C, some EClass.E, some EClass.A, some EClass.B] := by
  rcases prime_quadruplet_start_mod12 p hq h with h5 | h11
  · exact Or.inl (quadrupletFlavors_mod5 p h5)
  · exact Or.inr (quadrupletFlavors_mod11 p h11)

/-- Witness-Brücke (minimal): Chiralitätswort ↔ mod-12-Restfolge auf Q(p). -/
theorem chiralityWord_residue_bridge (p : ℕ) (ch : Chirality) (h : chiralityWord p = some ch) :
    quadrupletResidues p = chiralityResidues ch := by
  cases ch
  · simpa [chiralityResidues_ABCE] using
      quadrupletResidues_mod5 p (mod12_eq_five_of_chiralityWord_ABCE h)
  · simpa [chiralityResidues_CEAB] using
      quadrupletResidues_mod11 p (mod12_eq_eleven_of_chiralityWord_CEAB h)

theorem chiralityWord_flavor_bridge (p : ℕ) (ch : Chirality) (h : chiralityWord p = some ch) :
    quadrupletFlavors p = (chiralityOrder ch).map (fun c => some c) := by
  cases ch
  · simpa [chiralityOrder] using
      quadrupletFlavors_mod5 p (mod12_eq_five_of_chiralityWord_ABCE h)
  · simpa [chiralityOrder] using
      quadrupletFlavors_mod11 p (mod12_eq_eleven_of_chiralityWord_CEAB h)

end Primvierling
end CollatzEabc
