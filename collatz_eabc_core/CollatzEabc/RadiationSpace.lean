/-
  CollatzEabc.RadiationSpace — interpretativer Scaffold (rote Schicht).

  **Kein Physikanspruch.** Keine Theoreme über Strahlung, Rydberg-Niveaus oder
  Elektronenschalen. Dient als zukünftige Schnittstelle zwischen metrischen Defekten
  (Schicht A: Prime-Defekte, $D_A$) und strahlungsartigen Lesarten (Schicht C/rot).

  Referenz:
    collatz_eabc_epistemik_schichten.md
    collatz_eabc_epistemik_physik.md

  Verboten: `theorem radiation_space_exists` o.ä. mit physikalischer Bedeutung.
-/

import CollatzEabc.PrefProjection

namespace CollatzEabc

/-!
### Minimaler Zustandsstub (Verknüpfung zu EABC)

`EabcState` koppelt einen EABC-Buchstaben an einen Fensterindex — genug Struktur
für eine spätere Abbildung `epsilon`, ohne Physik zu behaupten.
-/

/-- Diskreter EABC-Zustand auf dem Kreisgraphen (Stub für Schicht A↔rot). -/
structure EabcState where
  /-- Aktuelle EABC-Klasse (E=0, A=1, B=2, C=3). -/
  letter : EabcLetter
  /-- Index eines Gleitfensters / Schritts entlang der Primfolge (rein kombinatorisch). -/
  windowIdx : ℕ

/-- Alias für den interpretativen Scaffold-Namen in der Dokumentation. -/
abbrev State := EabcState

/-!
### RadiationSpace — interpretativer Scaffold (rote Schicht)

Keine Existenzbehauptung, keine Rydberg-Identität, kein Strahlungsraum-Theorem.
`epsilon` ist eine **Typ-signierte Platzhalterabbildung** für künftige Schnittstellen.

Kommentar (nicht formalisiert): Rydberg-Lesart En ~ 1/n^2 als ikonische
Skalierung für Schalenindex `windowIdx` — Interpretation, nicht Konsequenz.
-/

/--
Interpretative scaffold. No physical meaning asserted.
Serves as future interface between metric defects and radiative structures.
-/
class RadiationSpace where
  /-- Trägermenge des noch offenen Strahlungs-Schnitttraums. -/
  carrier : Type
  /-- Platzhalter: effektives epsilon-Feld auf EABC-Zuständen (keine Physik-Semantik). -/
  epsilon : State → carrier

end CollatzEabc
