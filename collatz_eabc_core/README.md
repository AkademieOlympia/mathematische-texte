# collatz_eabc_core

Sorry-freies Lean-4-Paket für die **EABC-Strukturtheorie** der odd-to-odd Collatz-Dynamik (Mathlib 4.29).

**Kein Collatz-Beweis.** Die Uniformitätsvermutung (Stufe E) ist in `CollatzEabc.Open` nur als `Prop` dokumentiert — ohne `sorry`, ohne Beweis.

## Scope (0 `sorry`)

| Modul | Inhalt |
|-------|--------|
| `CollatzEabc.Density` | C-Ketten-Dichte, Tail-Reihen, Mischschranken |
| `CollatzEabc.Uniformity` | LTE-Reset, Valuations-Lemmas |
| `CollatzEabc.Z2Attraktor` | Stufen A–D: `ExceptionSetApprox`, Monotonie, `ExceptionSet`, U-Invarianz |
| `CollatzEabc.Open` | Stufe E: `collatzUniformityConjecture` (offen) |

Forschungsdateien mit `sorry` (Stufe E Strategien) bleiben im Repo-Root: `collatz_z2_attraktor.lean`, `collatz_uniformity_e.lean`.

## Build

```bash
cd collatz_eabc_core
lake update
lake build
```

Sorry-Check:

```bash
grep -r "sorry" CollatzEabc/   # sollte leer sein
```

## Quellen

Kopiert/adaptiert aus Repo-Root:

- `collatz_density_appendix.lean`
- `collatz_uniformity.lean`
- `collatz_z2_attraktor.lean` (nur Stufen A–D)

Offene Punkte: `collatz_offene_punkte.md`
