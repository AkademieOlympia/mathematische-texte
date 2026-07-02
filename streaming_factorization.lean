import Mathlib

/-!
# Streaming-Faktorisierungssatz (Lean 4 / Mathlib)

Formalisierung des Streaming-Faktorisierungssatzes für abstrakte Monoide
sowie der EABC-Instanz mit komponentenweiser Addition auf `ℤ × ℤ`.

## Hauptresultate

| Name                                | Aussage                                               |
|-------------------------------------|-------------------------------------------------------|
| `streamEval_nil`                    | `streamEval φ [] = 1`                                 |
| `StreamingFactorization`            | `streamEval φ (s₁ ++ s₂) = streamEval φ s₁ * streamEval φ s₂` |
| `streamEvalAdd_nil`                 | `streamEvalAdd φ [] = 0`                              |
| `AdditiveStreamingFactorization`    | `streamEvalAdd φ (s₁ ++ s₂) = streamEvalAdd φ s₁ + streamEvalAdd φ s₂` |
| `eabc_D_additive`                   | `D (s₁ + s₂) = D s₁ + D s₂`                          |
| `eabc_Q_additive`                   | `Q (s₁ + s₂) = Q s₁ + Q s₂`                          |
| `eabc_D_stream_additive`            | D-Additivität über EABC-Streams                       |
| `eabc_Q_stream_additive`            | Q-Additivität über EABC-Streams                       |
| `W_E_not_additive`                  | `W_E = D/Q` ist im Allgemeinen **nicht** additiv      |

## Konzeptionelle Struktur

Ein *Stream* `s : List α` wird in Blöcke `B₁, …, Bₖ` zerlegt.
Eine Blockbewertung `φ : α → M` in ein Monoid `M` ordnet jedem Block
einen Zustandsbeitrag zu. Die Streaming-Auswertung

    F(s) := φ(B₁) ⊗ φ(B₂) ⊗ … ⊗ φ(Bₖ)

faktorisiert über Stromkonkatenation, weil das Monoid assoziativ ist.
-/

-- =============================================================================
-- §1  Abstrakte Streaming-Faktorisierung (multiplikatives Monoid)
-- =============================================================================

section StreamingFactorization

variable {α : Type*} {M : Type*} [Monoid M]

/-- `streamEval φ s` ist das Monoidprodukt der mit `φ` bewerteten Blöcke:
    `F(s) = ∏ i, φ(s[i])`. -/
def streamEval (φ : α → M) (s : List α) : M :=
  (s.map φ).prod

/-- **F [] = 1**  –  der leere Stream liefert das neutrale Element Z₀. -/
@[simp]
theorem streamEval_nil (φ : α → M) : streamEval φ [] = 1 := by
  simp [streamEval]

/-- Einelementiger Stream wertet zum Block-Beitrag aus. -/
@[simp]
theorem streamEval_singleton (φ : α → M) (b : α) : streamEval φ [b] = φ b := by
  simp [streamEval]

/-- **Streaming-Faktorisierungssatz** (multiplikative Form):
    `F(s₁ ‖ s₂) = F(s₁) ⊗ F(s₂)`.

    *Beweis:* Nutze `List.map_append` (Distributivität von `map` über `++`)
    und `List.prod_append` (Produktivität von `prod` über `++`). -/
theorem StreamingFactorization (φ : α → M) (s₁ s₂ : List α) :
    streamEval φ (s₁ ++ s₂) = streamEval φ s₁ * streamEval φ s₂ := by
  simp [streamEval, List.map_append, List.prod_append]

/-- Korollar: volle Assoziativität über drei Teilstreams. -/
theorem StreamingFactorization₃ (φ : α → M) (s₁ s₂ s₃ : List α) :
    streamEval φ (s₁ ++ s₂ ++ s₃) =
      streamEval φ s₁ * streamEval φ s₂ * streamEval φ s₃ := by
  simp [StreamingFactorization, mul_assoc]

end StreamingFactorization

-- =============================================================================
-- §2  Additiver Streaming-Faktorisierungssatz
-- =============================================================================

section AdditiveStreamingFactorization

variable {α : Type*} {M : Type*} [AddMonoid M]

/-- `streamEvalAdd φ s` ist die additive Monoidfaltung der bewerteten Blöcke:
    `F(s) = ∑ i, φ(s[i])`. -/
def streamEvalAdd (φ : α → M) (s : List α) : M :=
  (s.map φ).sum

/-- **F [] = 0**  –  der leere Stream liefert den Nullzustand. -/
@[simp]
theorem streamEvalAdd_nil (φ : α → M) : streamEvalAdd φ [] = 0 := by
  simp [streamEvalAdd]

/-- **Additiver Streaming-Faktorisierungssatz**:
    `F(s₁ ‖ s₂) = F(s₁) + F(s₂)`.

    *Beweis:* Nutze `List.map_append` und `List.sum_append`. -/
theorem AdditiveStreamingFactorization (φ : α → M) (s₁ s₂ : List α) :
    streamEvalAdd φ (s₁ ++ s₂) = streamEvalAdd φ s₁ + streamEvalAdd φ s₂ := by
  simp [streamEvalAdd, List.map_append, List.sum_append]

end AdditiveStreamingFactorization

-- =============================================================================
-- §3  EABC-Instanz: Zustandsraum `ℤ × ℤ`
-- =============================================================================

section EABC

/-- Der EABC-Zustandsraum: `(N₊, N₋) ∈ ℤ × ℤ`.
    Die komponentenweise Addition `(a, b) ⊕ (c, d) = (a+c, b+d)` ist die
    Standard-`AddCommMonoid`-Struktur auf `ℤ × ℤ`, die Lean/Mathlib automatisch
    zur Verfügung stellt. -/
abbrev EABCState := ℤ × ℤ

/-- **Signed Imbalance**: `D(N₊, N₋) := N₊ − N₋`. -/
def D (s : EABCState) : ℤ := s.1 - s.2

/-- **Gesamtdichte**: `Q(N₊, N₋) := N₊ + N₋`. -/
def Q (s : EABCState) : ℤ := s.1 + s.2

/-- **`D` ist additiv unter ⊕**: `D(s₁ ⊕ s₂) = D(s₁) + D(s₂)`.

    *Beweis:* Destrukturierung der Paare, dann `ring` über ℤ. -/
theorem eabc_D_additive (s₁ s₂ : EABCState) :
    D (s₁ + s₂) = D s₁ + D s₂ := by
  obtain ⟨a₁, b₁⟩ := s₁
  obtain ⟨a₂, b₂⟩ := s₂
  simp [D]
  ring

/-- **`Q` ist additiv unter ⊕**: `Q(s₁ ⊕ s₂) = Q(s₁) + Q(s₂)`.

    *Beweis:* Destrukturierung der Paare, dann `ring` über ℤ. -/
theorem eabc_Q_additive (s₁ s₂ : EABCState) :
    Q (s₁ + s₂) = Q s₁ + Q s₂ := by
  obtain ⟨a₁, b₁⟩ := s₁
  obtain ⟨a₂, b₂⟩ := s₂
  simp [Q]
  ring

/-- EABC-Stream-Auswertung über eine Blockbewertungsfunktion `φ`:
    Akkumuliere die Zustandsbeiträge aller Blöcke. -/
def eabcEval {α : Type*} (φ : α → EABCState) (s : List α) : EABCState :=
  streamEvalAdd φ s

/-- **`D`-Faktorisierung über EABC-Streams**:
    `D(F(s₁ ‖ s₂)) = D(F(s₁)) + D(F(s₂))`.

    Folgt aus `AdditiveStreamingFactorization` + `eabc_D_additive`. -/
theorem eabc_D_stream_additive {α : Type*} (φ : α → EABCState) (s₁ s₂ : List α) :
    D (eabcEval φ (s₁ ++ s₂)) = D (eabcEval φ s₁) + D (eabcEval φ s₂) := by
  unfold eabcEval
  rw [AdditiveStreamingFactorization]
  exact eabc_D_additive _ _

/-- **`Q`-Faktorisierung über EABC-Streams**:
    `Q(F(s₁ ‖ s₂)) = Q(F(s₁)) + Q(F(s₂))`.

    Folgt aus `AdditiveStreamingFactorization` + `eabc_Q_additive`. -/
theorem eabc_Q_stream_additive {α : Type*} (φ : α → EABCState) (s₁ s₂ : List α) :
    Q (eabcEval φ (s₁ ++ s₂)) = Q (eabcEval φ s₁) + Q (eabcEval φ s₂) := by
  unfold eabcEval
  rw [AdditiveStreamingFactorization]
  exact eabc_Q_additive _ _

-- =============================================================================
-- §4  Nicht-Additivität von W_E = D / Q
-- =============================================================================

/-- Die **Effizienz-Observable** `W_E := D / Q`.
    Definiert über `ℚ` (rationale Zahlen), damit exakte Division möglich ist;
    `noncomputable` wegen der rationalen Division. -/
noncomputable def W_E (s : EABCState) : ℚ :=
  (D s : ℚ) / (Q s : ℚ)

/-- **`W_E` ist im Allgemeinen nicht additiv.**

    **Gegenbeispiel**: `s₁ = (2, 0)`, `s₂ = (1, 2)`:

    | Observable | `s₁`  | `s₂`  | `s₁ ⊕ s₂ = (3,2)` |
    |-----------|-------|-------|-------------------|
    | D         |  2    | −1    |  1                |
    | Q         |  2    |  3    |  5                |
    | W_E       |  1    | −1/3  | **1/5**           |

    `W_E(s₁) + W_E(s₂) = 1 − 1/3 = 2/3 ≠ 1/5 = W_E(s₁ ⊕ s₂)`.

    **Konsequenz**: `W_E` darf **nicht** inkrementell über Teilstreams summiert
    werden; sie muss nach der vollständigen Aggregation von `D` und `Q`
    neu berechnet werden. -/
theorem W_E_not_additive :
    ∃ (s₁ s₂ : EABCState), W_E (s₁ + s₂) ≠ W_E s₁ + W_E s₂ := by
  refine ⟨(2, 0), (1, 2), ?_⟩
  simp [W_E, D, Q]
  norm_num

end EABC
