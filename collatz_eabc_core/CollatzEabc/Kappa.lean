/-
  CollatzEabc.Kappa — Schnittstelle für treue Kodierung κ (Stufe 1 Generalangriff).

  **Kein Collatz-Beweis.** Nur Prop-Struktur, naive Präfix-Kodierung entlang
  odd-to-odd-Bahnen und dokumentierte offene Vermutung `kappaConjecture`.

  Referenz: `collatz_generalangriff_2026.md` (Stufe 1), `collatz_kepler_gedankenexperiment.tex`.
-/

import CollatzEabc.PrefProjection
import CollatzEabc.Z2Attraktor

namespace CollatzEabc

open CollatzZ2

/-- EABC-Wort als Buchstabenliste (vgl. `PrefProjection`). -/
abbrev EabcWord := List EabcLetter

/-- Mod-12-EABC-Klassifikation als `Option` (nur Reste 1,5,7,11). -/
def classOfLetter (n : ℕ) : Option EabcLetter :=
  match n % 12 with
  | 1 => some ⟨0, by decide⟩
  | 5 => some ⟨1, by decide⟩
  | 7 => some ⟨2, by decide⟩
  | 11 => some ⟨3, by decide⟩
  | _ => none

/-- Buchstabe aus `classOfLetter`, falls definiert. -/
def letterAt (n : ℕ) (h : (classOfLetter n).isSome) : EabcLetter :=
  (classOfLetter n).get h

/-- `collatzU` erhält Ungeradheit (`iterateU_odd` bei Schritt 1). -/
theorem collatzU_preserves_odd (n : ℕ) (h : Odd n) : Odd (collatzU n h) := by
  simpa [iterateU_succ n 0 h] using iterateU_odd n 1 h

/-- Erste `K` EABC-Schritte der odd-to-odd-Bahn ab ungeradem `n`.

Jeder Eintrag ist `classOfLetter (iterateU n i)` für `i < K`.
Für manche `n` (z. B. `n ≡ 3,9 (mod 12)`) ist der Eintrag `none` — das dokumentiert
Informationsverlust der naiven Paritäts→EABC-Abbildung. -/
def kappaPrefix (n : ℕ) (K : ℕ) (_h : Odd n) : List (Option EabcLetter) :=
  List.ofFn fun i : Fin K => classOfLetter (iterateU n i.val)

/-- Alle `K` Einträge von `kappaPrefix` sind definiert. -/
def kappaPrefixDefined (n : ℕ) (K : ℕ) (h : Odd n) : Prop :=
  ∀ i : Fin K, (classOfLetter (iterateU n i.val)).isSome

/-- Vollständiges Präfix als `EabcWord` (nur unter `kappaPrefixDefined`). -/
noncomputable def kappaPrefixWord (n : ℕ) (K : ℕ) (h : Odd n)
    (hdef : kappaPrefixDefined n K h) : EabcWord :=
  List.ofFn fun i : Fin K =>
    letterAt (iterateU n i.val) (by simpa [kappaPrefixDefined] using hdef i)

/-- Einzelindex-Form der Shift+Append-Dynamik (sorry-frei). -/
theorem kappaPrefix_get_shift (n : ℕ) (K : ℕ) (h : Odd n) (i : Fin K) :
    (kappaPrefix (collatzU n h) K (collatzU_preserves_odd n h))[i]! =
      if hi : i.val + 1 < K then
        (kappaPrefix n K h)[(⟨i.val + 1, hi⟩ : Fin K)]!
      else
        classOfLetter (iterateU n K) := by
  simp only [kappaPrefix, List.getElem_ofFn]
  by_cases hi : i.val + 1 < K
  · simp [hi, iterateU_succ n i.val h]
  · simp [hi]
    have hiEq : i.val = K - 1 := by omega
    have hiFin : i = ⟨K - 1, by omega⟩ := Fin.ext hiEq
    rw [hiFin]
    rw [← iterateU_succ n (K - 1) h]
    have hK : 1 ≤ K := by have := i.isLt; omega
    simp [Nat.sub_add_cancel hK]

/-!
### Treue Kodierung κ — Schnittstelle (offen)

Eine **treue** Kodierung liefert für jedes ungerade `n` ein EABC-Wort der Länge `K`,
dessen Buchstaben mit `classOfLetter ∘ iterateU` übereinstimmen, unter `U` die
Index-Shift-Dynamik erfüllen und auf Starts injektiv ist.
-/

/-- Daten einer hypothetischen treuen κ-Kodierung für feste Präfixlänge `K`. -/
structure FaithfulKappa (K : ℕ) where
  /-- Kodierung ungerader Starts als volles EABC-Wort der Länge `K`. -/
  encode : ∀ n : ℕ, Odd n → EabcWord
  /-- Länge stimmt. -/
  encode_len : ∀ n h, (encode n h).length = K
  /-- Jeder Buchstabe ist die mod-12-Klasse des entsprechenden Bahn-Glieds. -/
  encode_is_class :
    ∀ n h (i : Fin K), classOfLetter (iterateU n i) = some ((encode n h)[i]!)
  /-- Dynamiktreue (Index-Shift): `κ(U(n))_i = κ(n)_{i+1}` für `i < K-1`. -/
  dynamicsShift :
    ∀ n h (i : Fin K) (hi : i.val + 1 < K),
      (encode (collatzU n h) (collatzU_preserves_odd n h))[i]! =
        (encode n h)[(⟨i.val + 1, hi⟩ : Fin K)]!
  /-- Dynamiktreue (Append): letzter Buchstabe von `κ(U(n))` ist Klasse von `iterateU n K`. -/
  dynamicsAppend :
    ∀ n h (_hK : 0 < K),
      classOfLetter (iterateU n K) =
        some ((encode (collatzU n h) (collatzU_preserves_odd n h))[Fin.last (by omega)]!)
  /-- Injektivität auf ungeraden Starts (für festes `K`). -/
  injectiveOnStarts : ∀ n m (hn : Odd n) (hm : Odd m), n ≠ m → encode n hn ≠ encode m hm

/-- Existenz einer treuen κ-Kodierung der Länge `K` (offene Vermutung für `K ≥ 1`). -/
def faithfulKappaExists (K : ℕ) : Prop := Nonempty (FaithfulKappa K)

/-- **κ-Vermutung (Stufe 1):** Für jedes `K ≥ 1` existiert eine treue Kodierung.

Dies ist **schwächer** als Collatz und **stärker** als die naive `kappaPrefix`-Abbildung
(die `none`-Lücken und Kollisionen hat — siehe `collatz_kappa_test.py`). -/
def kappaConjecture : Prop := ∀ K, 0 < K → faithfulKappaExists K

/-- Dokumentierte Shift-Dynamik der naiven Präfix-Kodierung (Indexform). -/
abbrev naiveKappaDynamics (n : ℕ) (K : ℕ) (h : Odd n) (i : Fin K) : Prop :=
  (kappaPrefix (collatzU n h) K (collatzU_preserves_odd n h))[i]! =
    if hi : i.val + 1 < K then
      (kappaPrefix n K h)[(⟨i.val + 1, hi⟩ : Fin K)]!
    else
      classOfLetter (iterateU n K)

theorem naiveKappaDynamics_iff (n : ℕ) (K : ℕ) (h : Odd n) (i : Fin K) :
    naiveKappaDynamics n K h i ↔
      (kappaPrefix (collatzU n h) K (collatzU_preserves_odd n h))[i]! =
        if hi : i.val + 1 < K then
          (kappaPrefix n K h)[(⟨i.val + 1, hi⟩ : Fin K)]!
        else
          classOfLetter (iterateU n K) :=
  Iff.rfl

theorem naiveKappaDynamics_holds (n : ℕ) (K : ℕ) (h : Odd n) (i : Fin K) :
    naiveKappaDynamics n K h i :=
  kappaPrefix_get_shift n K h i

end CollatzEabc
