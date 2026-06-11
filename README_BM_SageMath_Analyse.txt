BM-nahe SageMath-Analyse für DNA.gb
===================================

Dateien
-------
- dna_analysis_bm.sage
- DNA.gb

Ziel
----
Diese Variante ergänzt die bisherige Standardanalyse um eine BM-nahe,
diskret-algebraische Sicht:

1. Basen werden quaternionisch/komplex inspiriert kodiert:
   A ->  1
   T -> -1
   C ->  i
   G -> -i

2. Die Sequenz wird auf modulo-12-Residuenklassen projiziert.

3. Für Fenstergrößen 3, 4 und 12 werden spektrale Signaturen berechnet.

4. ORF6 und ORF7a werden als getrennte Segmente verglichen.

Ausführen
---------
Im Verzeichnis mit DNA.gb:

    sage dna_analysis_bm.sage

Ausgabe
-------
Es wird das Verzeichnis

    dna_bm_output/

angelegt, mit u.a.:

- summary.txt
- base_counts.tsv
- top_codons.tsv
- top_4mers.tsv
- mod12_profile.tsv
- complex_mod12.tsv
- complex_mod12_dft.tsv
- cds_bm_signatures.tsv
- cds_comparison.tsv   (falls ORF6 und ORF7a vorhanden sind)

sowie PNG-Grafiken:

- top_codons.png
- top_4mers.png
- mod12_profile.png
- complex_mod12_abs.png
- complex_mod12_dft.png
- mod12_ORF6.png
- mod12_ORF7a.png

Interpretationshinweis
----------------------
Die BM-nahe Kodierung ist eine mathematische Modellierung zur Struktursondierung,
keine etablierte Standardmethode der Bioinformatik. Sie eignet sich für:

- Periodizitäts- und Rastertests
- Segmentvergleiche
- symbolische/algebraische Fingerabdrücke
- mod-12-Profile und Spektren

Naheliegende Erweiterungen
--------------------------
- de-Bruijn-Graph auf Basis der BM-Kodierung
- Vergleich mehrerer Genome oder Isolate
- Korrelationsanalyse zwischen Codonraster und mod-12-Profil
- Quaternionen/Oktonionen mit 4- bzw. 8-dimensionaler Einbettung
- Kopplung an BM-Achsenklassen (z.B. e,a,b,c-Schemata)
