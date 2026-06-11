# Zusammenfassung: Character-Theoretic Diagonalization (Final mit HL-Analyse)

## Kernbefund

Die Besetzungsasymmetrie der sechs mod-420-Halbkanäle
$\mathcal{R}=\{11,101,191,221,311,401\}$ bei Primzahlvierlingen $(p,p+2,p+6,p+8)$
lässt sich in drei epistemisch getrennte Ebenen zerlegen:

1. **Level 1 (exakte Algebra):** Der kubische Dirichlet-Charakter $\chi_3$ mod $7$
   fällt punktweise mit dem $k{=}2$-DFT-Modus auf $\mathcal{R}$ zusammen
   ($\chi_3|_{\mathcal{R}}=\psi_2$). Die Bias-Zerlegung in Charakterkoordinaten
   ist exakt; $\chi_3$ dominiert algebraisch, nicht asymptotisch.

2. **Level 2 (numerisch, bedingt auf $H$):** Unter der Multinomial-Nullhypothese
   skaliert $|\langle\delta_N,\chi_3\rangle|\sim 1/(6\sqrt{N})$ mit $z\approx 1$
   bis $N_{\max}=10^9$ — reine Fluktuation, kein Chebyshev-Drift.

3. **Level 3 (HL-Absicherung):** Die Hardy–Littlewood-Singulärreihe für
   $(0,2,6,8)$ liefert **identische Koeffizienten** auf allen sechs Kanälen
   (mod $2,3,4,5,7$). Asymptotisch folgt Gleichverteilung $1/6$.
   Die beobachtete Asymmetrie ist statistisch konsistent
   ($\chi^2\approx 4{,}34$, $p\approx 0{,}50$, $\mathrm{KL}\approx 4{,}55\times 10^{-4}$).
   Ein analytischer Beweis der bedingten Äquidistribution bleibt offen (H6).

## Schlussfolgerung

Die $\chi_3$-Dominanz ist **nicht dichtetheoretisch**, sondern durch die
Reduktion des Restklassenraums **algebraisch erzwungen**. HL unterscheidet
die Kanäle nicht; H2 (mod-7-HL-Dichten) und H3 (mod-4/mod-5-HL-Dichten)
sind als Erklärungsmechanismen widerlegt.

## Dateien

| Datei | Inhalt |
|-------|--------|
| `eabc_paper_final.tex` / `.pdf` | Hauptdokument (9 Seiten) |
| `eabc_hl_hypotheses.tex` | Standalone HL-Hypothesen H1–H6 |
| `eabc_hl_coefficient_hypotheses.py` | Numerische HL-Berechnung |
