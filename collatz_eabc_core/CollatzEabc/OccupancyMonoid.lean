/-
  CollatzEabc.OccupancyMonoid — Besetzungs-Monoid für EABC-Streaming-Kompression.

  **Ebene 1 (Theorem, GREEN):** abstrakter Zustand Z = (O, T, n) mit
    Z₀ = (∅, ⊤, 0),  Z₁ ⊕ Z₂ = (O₁ ∪ O₂, T₁ ∧ T₂, n₁ + n₂),
  `(Z, ⊕, Z₀)` ist kommutatives Monoid; Block-Fusion ist ordnungsunabhängig.

  **Ebene 2 (Struktur/Brücke):** O = besetzte Kanäle, T = Erstbesetzungs-Trigger,
  n = Ereigniszähler; F : BlockData → M Monoid-Homomorphismus (Kommentar).

  **Ebene 3 (EABC):** n = (N_plus, N_minus, …) — Verknüpfung zu PatternCount folgt.

  Referenz: `eabc_occupancy_tree.py` (`merge_state`), `collatz_eabc_zirkulationshypothese.md`.
  **Kein Collatz-Beweis.**
-/

import Mathlib.Data.Finset.Basic
import Mathlib.Data.List.Perm.Basic
import Mathlib.Algebra.Group.Defs

namespace CollatzEabc

variable {α β : Type*} [DecidableEq α] [LinearOrder β]

/-!
### Ebene 1 — abstrakter Besetzungszustand Z = (O, T, n)
-/

/-- Abstrakter Streaming-Zustand Z = (O, T, n).

* `occupied` — Menge besetzter Kanäle O
* `trigger` — Erstbesetzung T(c); `none` kodiert ⊤ (noch kein Trigger)
* `count` — Ereigniszähler n(c), außerhalb O implizit 0 -/
structure OccupancyState (α β : Type*) [DecidableEq α] where
  occupied : Finset α
  trigger : α → Option β
  count : α → ℕ

namespace OccupancyState

/-- Neutrales Element Z₀ = (∅, ⊤, 0). -/
def identity : OccupancyState α β where
  occupied := ∅
  trigger := fun _ => none
  count := fun _ => 0

/-- Punktweises Meet auf Triggern: min bei beiden definiert, sonst der vorhandene Wert. -/
def triggerMeet (t1 t2 : α → Option β) (c : α) : Option β :=
  match t1 c, t2 c with
  | none, ot => ot
  | ot, none => ot
  | some x, some y => some (min x y)

private theorem triggerMeet_none_left (t : α → Option β) (c : α) :
    triggerMeet (fun _ => none) t c = t c := by
  simp [triggerMeet]

private theorem triggerMeet_none_right (t : α → Option β) (c : α) :
    triggerMeet t (fun _ => none) c = t c := by
  simp only [triggerMeet]
  cases t c <;> rfl

private theorem triggerMeet_comm (t1 t2 : α → Option β) (c : α) :
    triggerMeet t1 t2 c = triggerMeet t2 t1 c := by
  simp only [triggerMeet]
  cases t1 c <;> cases t2 c <;> simp [min_comm]

private theorem triggerMeet_assoc (t1 t2 t3 : α → Option β) (c : α) :
    triggerMeet (triggerMeet t1 t2) t3 c = triggerMeet t1 (triggerMeet t2 t3) c := by
  simp only [triggerMeet]
  rcases t1 c with _ | x <;> rcases t2 c with _ | y <;> rcases t3 c with _ | z <;>
    simp [min_assoc]

/-- Monoid-Merge: Z₁ ⊕ Z₂ = (O₁ ∪ O₂, T₁ ∧ T₂, n₁ + n₂). -/
def merge (z1 z2 : OccupancyState α β) : OccupancyState α β where
  occupied := z1.occupied ∪ z2.occupied
  trigger c := triggerMeet z1.trigger z2.trigger c
  count c := z1.count c + z2.count c

@[ext]
theorem ext {z1 z2 : OccupancyState α β}
    (hO : z1.occupied = z2.occupied)
    (hT : ∀ c, z1.trigger c = z2.trigger c)
    (hn : ∀ c, z1.count c = z2.count c) : z1 = z2 := by
  rcases z1 with ⟨o1, t1, n1⟩
  rcases z2 with ⟨o2, t2, n2⟩
  simp only at hO hT hn
  subst hO
  congr
  · funext c; exact hT c
  · funext c; exact hn c

theorem merge_comm (z1 z2 : OccupancyState α β) : merge z1 z2 = merge z2 z1 := by
  ext c
  · simp [merge, Finset.union_comm]
  · simp [merge, triggerMeet_comm]
  · simp [merge, Nat.add_comm]

theorem merge_identity_left (z : OccupancyState α β) : merge identity z = z := by
  ext c
  · simp [merge, identity]
  · simp [merge, identity, triggerMeet_none_left]
  · simp [merge, identity]

theorem merge_identity_right (z : OccupancyState α β) : merge z identity = z := by
  ext c
  · simp [merge, identity]
  · simp [merge, identity, triggerMeet_none_right]
  · simp [merge, identity]

theorem merge_assoc (z1 z2 z3 : OccupancyState α β) :
    merge (merge z1 z2) z3 = merge z1 (merge z2 z3) := by
  ext c
  · simp [merge, Finset.union_assoc]
  · simp [merge, triggerMeet_assoc]
  · simp [merge, Nat.add_assoc]

theorem merge_comm_middle (x y z : OccupancyState α β) :
    merge x (merge y z) = merge y (merge x z) := by
  rw [← merge_assoc, merge_comm x y, merge_assoc]

instance : CommMonoid (OccupancyState α β) where
  mul := merge
  one := identity
  mul_assoc := merge_assoc
  one_mul := merge_identity_left
  mul_one := merge_identity_right
  mul_comm := merge_comm

/-!
### Folgen-Merge und ordnungsunabhängige Block-Fusion
-/

/-- Lineare Reduktion ⊕_{i=1}^k Z_i über eine Liste (linksfaltend). -/
def foldMerge (zs : List (OccupancyState α β)) : OccupancyState α β :=
  zs.foldl merge identity

theorem foldMerge_nil : foldMerge ([] : List (OccupancyState α β)) = identity := by
  rfl

private theorem foldl_merge_left (z : OccupancyState α β) (zs : List (OccupancyState α β)) :
    List.foldl merge z zs = merge z (List.foldl merge identity zs) := by
  induction zs generalizing z with
  | nil => simp [merge_identity_right]
  | cons w ws ih =>
    calc List.foldl merge z (w :: ws)
        = List.foldl merge (merge z w) ws := by simp [List.foldl_cons]
      _ = merge (merge z w) (List.foldl merge identity ws) := by rw [ih (merge z w)]
      _ = merge z (merge w (List.foldl merge identity ws)) := by rw [merge_assoc]
      _ = merge z (List.foldl merge identity (w :: ws)) := by
          congr 1
          simp [List.foldl_cons, merge_identity_left, ih w]

theorem foldMerge_cons (z : OccupancyState α β) (zs : List (OccupancyState α β)) :
    foldMerge (z :: zs) = merge z (foldMerge zs) := by
  simp only [foldMerge, List.foldl_cons, merge_identity_left]
  exact foldl_merge_left z zs

theorem foldMerge_append (l1 l2 : List (OccupancyState α β)) :
    foldMerge (l1 ++ l2) = merge (foldMerge l1) (foldMerge l2) := by
  induction l1 with
  | nil =>
    simp [foldMerge, merge_identity_left]
  | cons z zs ih =>
    simp only [List.cons_append, foldMerge_cons, ih, merge_assoc]

private theorem foldMerge_swap (x y : OccupancyState α β) (l : List (OccupancyState α β)) :
    foldMerge (x :: y :: l) = foldMerge (y :: x :: l) := by
  simp only [foldMerge_cons]
  exact merge_comm_middle x y (foldMerge l)

/-- Permutation der Blockliste ändert den Endzustand nicht (Kommutativität). -/
theorem foldMerge_perm {l l' : List (OccupancyState α β)} (h : l.Perm l') :
    foldMerge l = foldMerge l' := by
  induction h with
  | nil => rfl
  | cons x h ih =>
    simp only [foldMerge_cons, ih]
  | swap x y l =>
    exact (foldMerge_swap x y l).symm
  | trans h₁ h₂ ih₁ ih₂ =>
    rw [ih₁, ih₂]

/-- **Satz (Streaming-Kompression).**

Es existiert ein kommutatives Monoid M = (Z, ⊕, Z₀) (hier `OccupancyState` mit `merge`/`identity`),
sodass der Endzustand Z_final = ⊕_{i=1}^k Z_i wohldefiniert ist und von der Block-Fusionsreihenfolge
unabhängig bleibt. -/
theorem streaming_compression_monoid :
    (∃ Z0 : OccupancyState α β, Z0 = identity) ∧
    (∀ z1 z2 z3 : OccupancyState α β,
      merge (merge z1 z2) z3 = merge z1 (merge z2 z3)) ∧
    (∀ z1 z2 : OccupancyState α β, merge z1 z2 = merge z2 z1) ∧
    (∀ z : OccupancyState α β, merge identity z = z) ∧
    (∀ z : OccupancyState α β, merge z identity = z) :=
  And.intro ⟨identity, rfl⟩
    (And.intro merge_assoc
      (And.intro merge_comm (And.intro merge_identity_left merge_identity_right)))

theorem streaming_compression_fold_independent (l l' : List (OccupancyState α β))
    (h : l.Perm l') : foldMerge l = foldMerge l' :=
  foldMerge_perm h

/-- Baum- vs. lineare Fusion: jede Permutation der Blockliste liefert denselben Endzustand. -/
theorem fold_merge_independent (l l' : List (OccupancyState α β)) (h : l.Perm l') :
    foldMerge l = foldMerge l' :=
  foldMerge_perm h

/-!
### Ebene 2 — BlockScan: Streaming-Faktorisierung F(P₁ ⊔ P₂) = F(P₁) ⊕ F(P₂)
-/

/-- Block-Scan: Liste lokaler Zustände → Endzustand via ⊕ (Monoid-Reduktion). -/
def blockScan (blocks : List (OccupancyState α β)) : OccupancyState α β :=
  foldMerge blocks

/-- **Satz (Streaming-Faktorisierung).** Disjunkte Blockpartition liefert Monoid-Homomorphismus:
`blockScan (P₁ ++ P₂) = blockScan P₁ ⊕ blockScan P₂`. -/
theorem blockScan_append (p1 p2 : List (OccupancyState α β)) :
    blockScan (p1 ++ p2) = merge (blockScan p1) (blockScan p2) :=
  foldMerge_append p1 p2

/-- Reihenfolge der Blöcke ist irrelevant (`foldMerge_perm`). -/
theorem blockScan_perm {p p' : List (OccupancyState α β)} (h : p.Perm p') :
    blockScan p = blockScan p' :=
  foldMerge_perm h

/-!
Auswertungsfunktoren (Ebene 2, Struktur — nicht Teil des Monoid-Beweises):
  f_D(Z) = D_E,  f_Q(Z) = Q_E,  f_W(Z) = W_E.
Verknüpfung mit Holonomie- und Skalierungsschicht folgt separat (kein PatternCount hier).

### Ebene 3 — EABC-Zählvektor n = (N_plus, N_minus, …)

Spezialisierung des Zählers n auf EABC-Faser-Komponenten — eigener Branch/PR, nicht h12.
-/

end OccupancyState

end CollatzEabc
