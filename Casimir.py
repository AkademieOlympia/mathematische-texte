# SageMath 10.5: Thermodynamische Stabilität des BQL
from sage.all import RR, k_B, hbar, c, pi

T = 293 # Temperatur in Kelvin (Bamberg)
d_nm = 100
d = d_nm * 1e-9

# Thermische Energie pro Volumen (ca. k_B * T)
E_therm = k_B * T / d**3

# Casimir-Repulsion (geheilt)
E_casimir = (pi**2 * hbar * c) / (720 * d**3)

verhaeltnis = E_casimir / E_therm
print(f"Stabilitätsfaktor (Casimir/Thermik): {verhaeltnis.n()}")
# Ein Wert > 1 bedeutet, dass die Quantenordnung über das thermische Rauschen dominiert.