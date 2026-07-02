# Analyse: EABC-Vierer vs. echte Primzahl-Quadrupletts

## Definitionen

- **Echte Primzahl-Quadrupletts:** Vierer der Form `(n, n+2, n+6, n+8)` mit vier Primzahlen.
- **Allgemeine EABC-Vierer:** Vierer der Form `(n, n+2, n+6, n+8)` fuer beliebiges `n` im Auswertebereich; Primzahlstatus beliebig.

## Setup

- Bereich: `n in [5, 200000]`
- Schrittweite: `2`
- Gepruefte Startwerte: `99998` (Intervalllaenge in n: `199996`)
- Zentrale Groessen: `Delta_mat = n - cbrt((n+2)(n+6)(n+8))`, `Delta_vac = (n+4) - cbrt((n+1)(n+4)(n+7))`.

## Kernbefunde

- Echte Primzahl-Quadrupletts im Bereich: **52**.
- Prime-Count-Verteilung (0..4): `0`: 43666 (43.67%), `1`: 42330 (42.33%), `2`: 12464 (12.46%), `3`: 1486 (1.49%), `4`: 52 (0.05%).
- `Delta_mat`: mean=-5.333177, median=-5.333302, p05=-5.333317, p95=-5.333023.
- `Delta_vac`: mean=0.000152, median=0.000030, p05=0.000016, p95=0.000300.
- Dominante Signaturen: `CCCC`=43666 (43.67%), `CPCC`=10612 (10.61%), `CCPC`=10604 (10.60%), `CCCP`=10561 (10.56%), `PCCC`=10553 (10.55%), `PCPC`=3508 (3.51%), `CPCP`=3491 (3.49%), `PCCP`=1391 (1.39%), `PPCC`=1364 (1.36%), `CCPP`=1363 (1.36%).
- Segmentierung zwischen aufeinanderfolgenden echten Quadrupletts: 51 Intervalle ausgewertet.
- Auffaellige Intervalle (Count- oder Prime3+-Outlier): **6**.

## Special Appearance Patterns (datenbasiert)

- **dominante_signatur**: CCCC (43666/99998 (43.67%))
- **haeufigster_prime_count**: 0 (43666/99998 (43.67%))
- **prime_reich_ge3**: Prime-Count >= 3 (1538/99998 (1.54%))
- **linker_rand_PP**: Signatur beginnt mit PP (2159/99998 (2.16%))
- **rechter_rand_PP**: Signatur endet mit PP (2158/99998 (2.16%))
- **intervall_max_allgemeine_vierer**: [171161, 187631] (general_count=8234, z=3.32)
- **intervall_max_prime3plus_rate**: [11, 101] (prime3plus_rate=22.73%, n=10)

## Top-Intervalle nach Anzahl allgemeiner Vierer

- [171161, 187631]: general_count=8234, avg_prime_count=0.667, dominant_signature=CCCC.
- [122201, 135461]: general_count=6629, avg_prime_count=0.678, dominant_signature=CCCC.
- [144161, 157271]: general_count=6554, avg_prime_count=0.677, dominant_signature=CCCC.
- [88811, 97841]: general_count=4514, avg_prime_count=0.708, dominant_signature=CCCC.
- [34841, 43781]: general_count=4469, avg_prime_count=0.752, dominant_signature=CCCC.

## Prime-reiche Intervalle (Prime-Count >= 3)

- [11, 101]: prime3plus_rate=22.73%, prime3plus_count=10, general_count=44.
- [101, 191]: prime3plus_rate=9.09%, prime3plus_count=4, general_count=44.
- [821, 1481]: prime3plus_rate=6.38%, prime3plus_count=21, general_count=329.
- [191, 821]: prime3plus_rate=6.05%, prime3plus_count=19, general_count=314.
- [1481, 1871]: prime3plus_rate=4.64%, prime3plus_count=9, general_count=194.

## Reproduzierbarkeit

- Skript: `python/analysis/eabc_quadruplet_patterns.py`
- Kommando: `python3 python/analysis/eabc_quadruplet_patterns.py --n-min 5 --n-max 200000 --step 2`

Hinweis: Die Auswertung ist rein arithmetisch-deskriptiv und erhebt keine physikalischen Claims.
