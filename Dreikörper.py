from sage.all import *

# Aktivierung des differentialgeometrischen Pakets von Sage 10.5
M = Manifold(4, 'M', structure='Lorentzian')
X = M.chart('t r th:theta ph:phi')
t, r, th, ph = X[:]

# Definition des arithmetischen Radius als symbolische Funktion/Variable
r_A = var('r_A')
assume(r_A > 0)

# Definition des metrischen Tensors g
g = M.metric('g')
g[0,0] = 1 - r_A/r
g[1,1] = -1 / (1 - r_A/r)
g[2,2] = -r**2
g[3,3] = -r**2 * sin(th)**2

# Berechnung der Christoffel-Symbole (Scheinkräfte der Zahlgeometrie)
print("=== Christoffel-Symbole (Auswahl) ===")
print(g.christoffel_symbols_display())

# Berechnung des Riemannschen Krümmungstensors
Riem = g.riemann()
print("\n=== Riemann-Krümmungstensor (Komponenten) ===")
Riem.display_comp()

# Überprüfung der Einstein-Gleichungen im Vakuum (außerhalb des Kerns)
Ric = g.ricci()
print("\n=== Ricci-Tensor (Sollte außerhalb von r_A verschwinden) ===")
Ric.display()