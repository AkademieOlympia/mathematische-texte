# B-infinity Skalen-Audit - v13

## Berichtete Skalenfolge

| N | berichtete Größenordnung p_N | B |
|---:|---:|---:|
| 10^4 | 10^5 | 5.3150 % |
| 10^5 | 10^6 | 5.0865 % |
| 10^6 | 10^7 | 4.9810 % |
| 10^7 | 10^8 | 4.8978 % |

## Berichtete Regressionsmodelle

LOS-Form:

B_LOS(x) = 0.037806 + 0.0687 * ln(ln x)/ln x

Berichtete extrapolierte Konstante: 3.7806 %.

Log-Form:

B_log(x) = 0.042683 + 0.1121/ln x

Berichtete extrapolierte Konstante: 4.2683 %.

## Audit-Grenze

Diese beiden Werte sind modellabhängige Extrapolationen aus vier Skalenpunkten. Die Spanne [3.7806 %, 4.2683 %] ist weder ein Konfidenzintervall noch eine bewiesene Schranke für einen wahren Grenzwert. Aus den im Manuskript gerundeten B-Werten allein lassen sich die angegebenen Regressionskoeffizienten nicht exakt reproduzieren, solange die präzisen x-Werte und ungerundeten Messdaten nicht mitgeliefert werden.

Der robuste Befund lautet nur: Die vier getesteten Punkte zeigen Drift; 4.9 % ist auf dieser Folge kein beobachteter Fixpunkt.

## Kosmologische Referenz

Planck 2018 berichtet in base-LambdaCDM ungefähr Omega_b h^2 = 0.0224 und H0 = 67.4 km/s/Mpc, entsprechend Omega_b ungefähr 0.0493. Der Wert B(10^7)=0.048978 ist daher numerisch nahe, aber kein exakter Vier-Nachkommastellen-Treffer auf den kanonischen Planck-2018-Basiswert.

## Claim-Wall

- [E] Vier-Punkt-Skalenfolge und beobachtete Drift.
- [FIT] 3.7806 % und 4.2683 %: berichtete modellabhängige Extrapolationen.
- [H] B_infinity = Omega_b: offene physikalische Brückenhypothese.
