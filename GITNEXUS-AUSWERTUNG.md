# GitNexus-Auswertung: Mathematische Texte

**Erstellt:** 2. März 2026  
**Index:** Commit `28f2f50` | 387 Dateien | 3.455 Symbole | 6.586 Kanten

---

## 1. Übersicht nach Dateityp

| Typ | Anzahl | Anteil |
|-----|--------|--------|
| **Python** (.py) | 79 | 20 % |
| **LaTeX** (.tex) | 59 | 15 % |
| **Jupyter** (.ipynb) | 34 | 9 % |
| **C++** (.cpp) | 32 | 8 % |
| **JSON** | 16 | 4 % |
| **Markdown** (.md) | 7 | 2 % |
| **Sonstige** (aux, log, out, toc, fls, …) | 160 | 41 % |

---

## 2. Dateien mit den meisten Code-Symbolen (Top 40)

Die symbolreichsten Dateien – Funktionen, Klassen, Methoden, Interfaces:

| Datei | Symbole |
|-------|---------|
| ring_bench_e8_ecm_v8.cpp | 84 |
| Schütte_V0.0.14.cpp | 77 |
| ring_bench_e8_ecm_v7.cpp | 76 |
| ring_bench_e8_ecm_v6.cpp | 76 |
| Schütte_V0.0.11.cpp | 75 |
| Schütte_V0.0.13.cpp | 75 |
| Schütte_V0.0.12.cpp | 75 |
| Schütte_V0.0.10.cpp | 73 |
| Schütte_V0.0.9.cpp | 72 |
| Schütte_V0.0.8.cpp | 70 |
| Schütte_V0.0.6.cpp | 70 |
| Schütte_V0.0.7.cpp | 70 |
| Schütte_V0.0.4.cpp | 70 |
| Schütte_V0.0.5.cpp | 70 |
| Schütte_V0.0.3.cpp | 70 |
| Rabin_Hoff.cpp | 68 |
| bm_qhe_model_v4_4.py | 63 |
| ring_prefilter_fast_mr_toggle.py | 63 |
| Genesis/ring_prefilter_fast_mr_toggle.cpp | 63 |
| bm_qhe_model_v4_3.py | 62 |
| E8_only_V0.0.16.cpp | 56 |
| E8_only_V0.0.15.cpp | 56 |
| ring_prefilter_128_pref_vs_rho_fast128.cpp | 56 |
| ring_prefilter_128_compare_mr_rho.cpp | 54 |
| E8_only_V0.0.14.cpp | 54 |
| ring_prefilter_128_compare_mr.cpp | 53 |
| ring_bench_e8_compare_v5b.cpp | 52 |
| bm_qhe_model_v4.py | 51 |
| Genesis/AntiMiller_Rabin.py | 51 |
| Schütte_V0.0.1.cpp | 48 |
| Schütte_V0.0.2.cpp | 48 |
| bm_qhe_model_v3.py | 45 |
| ring_bench_dataset_rho_boostfixed.cpp | 43 |
| ring_bench_e8_compare_v3.cpp | 41 |
| AlpaFit.py | 41 |
| bm_fib_ico_alphafit_scaler_v3.py | 29 |
| bm_fib_ico_alphafit_scaler.py | 26 |
| ring_bench_dataset_rho_cppint.cpp | 23 |
| Untitled-1.cpp | 21 |

---

## 3. Funktionale Bereiche (Clusters)

Top 20 Module nach Symbolanzahl und Kohäsion (Leiden-Algorithmus):

| Modul | Symbole | Kohäsion |
|-------|---------|----------|
| Genesis | 58 | 60 % |
| Cluster_373 | 49 | 84 % |
| Cluster_371 | 20 | 77 % |
| Cluster_379 | 19 | 76 % |
| Cluster_923 | 17 | 91 % |
| Cluster_385 | 12 | 74 % |
| Cluster_394 | 12 | 86 % |
| Cluster_401 | 11 | 74 % |
| Cluster_841 | 11 | 100 % |
| Cluster_919 | 11 | 96 % |
| Cluster_1023 | 11 | 84 % |
| Cluster_9 | 10 | 92 % |
| Cluster_395 | 10 | 80 % |
| Cluster_396 | 10 | 78 % |
| Cluster_835 | 10 | 95 % |
| Cluster_4 | 9 | 91 % |
| Cluster_931 | 9 | 86 % |
| Cluster_12 | 8 | 80 % |
| Cluster_370 | 8 | 61 % |
| Cluster_393 | 8 | 74 % |

**Genesis** (größter Cluster): Enthält u.a. `AntiMiller_Rabin.py`, `ring_prefilter_fast_mr_toggle.py`, `bm_qhe_model_v4_3.py` – Fokus auf Ring-Filter, Pollard-Rho, QHE-Modelle.

---

## 4. Ausführungsflüsse (Processes)

271 Prozesse erkannt. Beispiele für cross-community Flows:

- `Main → Jacobi` (7 Schritte)
- `Main → Witness` (7 Schritte)
- `Main → Isqrt` (7 Schritte)
- `Chern_projector_marker_stereo → Edges` (7 Schritte)
- `Estimate_qhe_nu_and_chern → Edges` (7 Schritte)
- `Detect_edge_modes → Edges` (7 Schritte)
- `Lz_diagnostics → Edges` (7 Schritte)
- `Lz_shell_splitting → Edges` (7 Schritte)

---

## 5. Graph-Statistik

| Metrik | Wert |
|--------|------|
| Dateien | 387 |
| Symbole (Nodes) | 3.455 |
| Kanten (Relations) | 6.586 |
| Clusters | 1.101 |
| Prozesse | 271 |

---

## 6. Nützliche GitNexus-Befehle

```bash
# Suche im Knowledge Graph
gitnexus query <suchbegriff>

# 360°-Ansicht eines Symbols
gitnexus context <symbol>

# Blast-Radius: Was bricht bei Änderung?
gitnexus impact <symbol>

# Rohdaten-Abfrage
gitnexus cypher "MATCH (f:File) RETURN f.filePath LIMIT 20"
```

---

*Erzeugt mit GitNexus Knowledge Graph aus dem Index vom 2. März 2026.*
