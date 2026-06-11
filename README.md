# Mathematische Texte

Sammlung mathematischer Manuskripte, LaTeX-Quellen und begleitender Skripte (EABC, Geometrie, Zahlentheorie).

## Inhalt

- LaTeX- und PDF-Manuskripte im Repo-Root
- Python-/SageMath-Hilfsskripte für numerische Experimente
- Lean-Teilprojekte (z. B. `ptolemaeus-lean/`) — Build-Artefakte (`.lake/`) werden nicht versioniert

## Große lokale Dateien

Installations-Images (`.dmg`), Prime-Caches (`.npy`) und sehr große CSV-Dateien bleiben absichtlich lokal und sind in `.gitignore` eingetragen.

## PDF bauen

```bash
pdflatex <datei>.tex
```

## Verwandte Repositories

- [eabc-renorm](https://github.com/AkademieOlympia/eabc-renorm) — formaler Lean-Kern der EABC-Renormierung
