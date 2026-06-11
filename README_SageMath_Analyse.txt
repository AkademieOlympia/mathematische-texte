SageMath-Analyse für DNA.gb
==========================

Dateien:
- dna_analysis.sage   -> SageMath-Skript
- DNA.gb              -> ursprüngliche GenBank-Datei

Verwendung:
1. Beide Dateien in dasselbe Verzeichnis legen.
2. Im Terminal ausführen:

   sage dna_analysis.sage

Falls die GenBank-Datei anders heißt:
- Im Skript unten bei main(path='DNA.gb') den Dateinamen anpassen.

Was das Skript berechnet:
- Metadaten aus der GenBank-Datei
- Basenzählung und GC-Gehalt
- Shannon-Entropie
- Basenverteilung nach Codonposition 1/2/3
- häufigste Codons und 4-mere
- Übersetzung der CDS-Features
- Konsistenzcheck gegen annotierte Proteinsequenzen
- einfache Übergangsmatrix A,C,G,T -> A,C,G,T
- Vektor der Basenzahlen als algebraisches Objekt
- einfache de-Bruijn-artige Kantenliste
