# EABC Oktaeder-Umgebung — numerische Testdesign

**Status:** Experiment (Modellabbildung)  
**Branch:** `collatz/eabc-euklidische-hebung`, PR #59  
**Code:** `collatz_eabc_oktaeder_flux_test.py` → `collatz_eabc_oktaeder_flux_test.json`  
**Kanonische Theorie:** `collatz_eabc_diskrete_geometrie.md` (Φ_E, ⟨ω_E,h⟩)

**Epistemisches Label:** Die **Oktaeder-Umgebung** ist eine explizite **Modellabbildung** — kein Beweis der
EABC-Vermutung. Sie testet, ob die 1D-C4-Vorhersagen unter einer 8-fach erweiterten Geometrie **konsistent**
bleiben oder **widersprechen**.

---

## 1. Was „Oktaeder-Umgebung“ hier bedeutet

Im Projekt existieren drei verwandte, aber getrennte Begriffe:

| Begriff | Objekt | Rolle in EABC |
|---------|--------|---------------|
| **Oktonion-Schale** Σ_n^(8) | S^7 ⊂ ℝ^8, r_8(n) Gitterpunkte | 8D-Normniveau; Hurwitz-Kette (geplant) |
| **Oktagon-Darstellung** | 8 Punkte + Zentrum in 3D | Prim-Vierlinge, halbe Schalen (`eabc_octagonal_shells_test.py`) |
| **Oktaeder-Umgebung** (dieses Dokument) | Regulärer Oktaeder O_6 ⊂ ℝ^3 | **Dual zum Würfel**, equatorialer C4-Zyklus + polare 8D-Lift-Achse |

**Kanonische Wahl für dieses Experiment:** Der reguläre Oktaeder trägt den C4-Kreis als **Äquator** in der
xy-Ebene; die beiden Pole ±e₃ kodieren die **8-fache Oktonion-Lift-Richtung** (ABCEA → +z, CEABC → −z).

$$\text{Würfel (8 Ecken)} \;\longleftrightarrow\; \text{Oktaeder (6 Ecken, 8 Dreiecksflächen)} \;\supset\; C_4 \cong S^1.$$

---

## 2. Explizite Abbildung C4 → Oktaeder

### 2.1 Knoten

| EABC-Klasse | Oktaeder-Ecke | Koordinaten |
|-------------|---------------|-------------|
| E | +e₁ | (1,0,0) |
| A | +e₂ | (0,1,0) |
| B | −e₁ | (−1,0,0) |
| C | −e₂ | (0,−1,0) |
| P⁺ | +e₃ | (0,0,1) — ABCEA-Pol |
| P⁻ | −e₃ | (0,0,−1) — CEABC-Pol |

### 2.2 Kanten (12 Stück)

**Äquator** (identisch mit C4): EA, AB, BC, CE ∈ E⁺; EC, CB, BA, AE ∈ E⁻.

**Polare Diagonalen** (8-fach Lift): für jedes X ∈ {E,A,B,C} die Kanten X–P⁺ und X–P⁻.

### 2.3 8D-Schalen-Gewicht r_8(n)

Oktonionische Normschale | Appel n^(8): Gewicht w(n) = r_8(n) mit r_8 = Koeffizient von q^n in θ₃(q)^8
(Jacobi; OEIS A000118). Berechnung: Faltung r₄∗r₄ mit r₄(n) = 8 Σ_{d|n, 4∤d} d.

Jedes 5-Fenster-Holonomie-Ereignis bei Prim-Obergrenze p trägt Vorzeichen ω ∈ {±1} und Gewicht w(p).

### 2.4 8→4-Projektion (Koordinatenpaare)

Oktonion-Basis e₁…e₈ wird paarweise auf EABC projiziert:

$$(e_1,e_2)\mapsto E,\;(e_3,e_4)\mapsto A,\;(e_5,e_6)\mapsto B,\;(e_7,e_8)\mapsto C.$$

Polarfluss: ABCEA verteilt +1 auf P⁺-Diagonalen, CEABC auf P⁻-Diagonalen (zusätzlich zum äquatorialen ω_E).

---

## 3. Zu testende Observablen

| Observable | 1D (C4) | Oktaeder-Umgebung |
|------------|---------|-------------------|
| Zirkulation | C_E = N₊ − N₋ | C_oct (äquatorial, ungewichtet) |
| Flussdichte | Φ_E = C_E/S_E | Φ_oct,eq (soll = Φ_E) |
| Schalen-gewichtet | — | Φ_oct,shell = Σ ω·r₈(p) / Σ r₈(p) |
| Harmonische Paarung | ⟨ω_E, h⟩ | ⟨ω_oct, h_eq⟩ auf 4 äquatorialen Kanten |
| Magnetischer Laplace | L_mag auf 4 Knoten | L_mag auf 6 Oktaeder-Knoten |
| Orientierung | ABCEA vs CEABC | P⁺ vs P⁻ Pol-Fluss |

**Zentrale Frage:** Unterstützt Φ_oct,shell ≠ 0 (bzw. ⟨ω,h⟩ ≠ 0) die Vermutung Φ_E ≠ 0, oder kollabiert
das Signal unter Schalen-Gewichtung?

---

## 4. Erwartungen und Falsifikation

| Ergebnis | Lesart |
|----------|--------|
| Φ_oct,eq = Φ_E exakt | Abbildung konsistent (Sanity) |
| Φ_oct,shell ≠ 0 mit gleichem Vorzeichen wie Φ_E | **Unterstützt** EABC-Vermutung in 8D-Lift |
| Φ_oct,shell → 0 während Φ_E ≠ 0 | Schalen-Gewichtung **dämpft** Signal — epistemisch offen |
| ⟨ω,h⟩_oct = ⟨ω,h⟩_C4 auf Äquator | Harmonischer Anteil stabil unter Lift |
| L_mag near-zero-Modus auf O_6 | Analogie zur Übergangsraum-Vermutung |

**Falsifikation:** Φ_oct,eq ≠ Φ_E (Implementierungsfehler) oder ⟨ω,h⟩_oct = 0 bei Φ_E ≠ 0 ohne erklärbare
Projektionsursache.

---

## 5. Ausführung

```bash
python3 collatz_eabc_oktaeder_flux_test.py --max-p 1000000
pytest tests/test_eabc_oktaeder_flux.py -q
```

**JSON:** `collatz_eabc_oktaeder_flux_test.json`

**Querverweise:** `collatz_eabc_hodge_eabc.py`, `collatz_eabc_holonomie_fehlerterm.py`,
`collatz_eabc_diskrete_geometrie.md` §2–§4.
