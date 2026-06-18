# EABC Brachistochrone / Fermat-Prinzip (mit chiraler Doppelbrechung)

**Status:** Modellabbildung — Variationsprinzip, kein Physikanspruch  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Modellabbildung | Analogie | Experiment

**Querverweise:**
- `collatz_eabc_chirale_polarisation.md` — $T_R$, $T_L$, Birefringenz-Analogie
- `collatz_eabc_epistemik_physik.md` — Wegfunktion ohne Zeitdilatation
- `collatz_eabc_holonomie_stufen.md` — drei Stufen, kein SRT
- `collatz_eabc_kritische_abbildung.md` — ABCEA-Trajektorie, Holonomie-Sensor
- `collatz_eabc_brachistochrone.py` — Numerik

---

## 1. Fermat-Prinzip (euklidisch)

$$T = \int \frac{\mathrm{d}s}{v(x,\gamma)} \;\approx\; \sum_j \frac{\Delta s_j}{v_j}.$$

**Kein SRT-Anspruch:** $T$ ist euklidische Laufzeit bei festgelegten $v_j$ — kein Lorentz-Faktor.

Fünf Potenzialfamilien (Modellwahl):
1. $V(x) = \ln x$
2. $V(\gamma) = \sum_n 1/((\gamma-\gamma_n)^2+\varepsilon)$ (ζ-Summe)
3. $V_E = D_E$ (Chiralität)
4. Krümmungshemmung
5. Informationspotential $V = -\ln P$

**Implementierung:** `velocity_from_potential`, `compare_paths_for_potential`.

---

## 2. Chirale Doppelbrechung: $T_R$ und $T_L$

Erweiterung für **zwei Polarisationskanäle** (vgl. `collatz_eabc_chirale_polarisation.md`):

$$T_R = \int \frac{\mathrm{d}s}{v_R},\qquad T_L = \int \frac{\mathrm{d}s}{v_L}.$$

**Analogie:** optisch aktive Medien — unterschiedliche Phasengeschwindigkeit für rechts- ($R$) und linkszirkulare ($L$) Polarisation; Berry-Phase-Systeme mit Kanalabhängigkeit.

**Modell (Brachistochrone):**
$$v_R = v_0 + \alpha\, V_E,\qquad v_L = v_0 - \alpha\, V_E,$$
mit $V_E = D_E$ (globaler Chiralitätsfluss) oder lokalem $D_E$-Proxy entlang der Trajektorie.

**Observable:**
$$\Delta T = T_R - T_L \quad\text{(Birefringenz-Analogie)}.$$

**Label:** Doppelkanal = **Modellabbildung**; optische Aktivität = **Analogie**.

---

## 3. ABCEA-Trajektorie als Testpfad

Der Holonomie-Sensor liefert Knoten $(x_j, \gamma_j)$ entlang ABCEA mit Lücken $(2,4,2,4)$ mod $12$.

Vergleich:
- gerader Polygonzug vs. Gradient-gestörter Pfad
- Halbkreis-Kette in der komplexen Ebene

**Forschungsfrage:** Bleibt der optimale Pfad gerade oder biegt er zu hoher Primdichte / Holonomie?

---

## 4. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $T$ | `travel_time_integral` | `collatz_eabc_brachistochrone` |
| $T_R$, $T_L$ | `travel_time_birefringent` | `collatz_eabc_brachistochrone` |
| $v_R$, $v_L$ | `birefringent_velocity_pair` | `collatz_eabc_brachistochrone` |
| $\Delta T$ | `delta_T` | `travel_time_birefringent` |

---

*Brachistochrone auf dem EABC-Holonomie-Sensor; chirale Erweiterung $T_R/T_L$ als Doppelbrechungs-Analogie — kein relativistisches Eigenzeit-Modell.*
