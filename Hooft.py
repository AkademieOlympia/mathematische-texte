from decimal import Decimal, getcontext

getcontext().prec = 36

# Numerische Zeugen aus der Simulation
alpha = Decimal('1') / Decimal('137.035999177')
M_X = Decimal('2.0e16') 
M_Planck = Decimal('1.2209e19')

# Berechnung der Monopolmasse
M_M = M_X / alpha
M_M_Planck = M_M / M_Planck  # Das Verhältnis 0.224484...

# Definition der EABC-Gravitationsstrukturkonstante
alpha_G_EABC = M_M_Planck ** 2

print(f"Verhältnis M_M / M_Planck:               {M_M_Planck:.6f}")
print(f"EABC-Gravitationsstrukturkonstante alpha_G: {alpha_G_EABC:.6f}")
print(f"Kehrwert (1 / alpha_G):                   {1 / alpha_G_EABC:.2f}")
