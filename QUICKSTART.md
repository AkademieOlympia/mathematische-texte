# FEniCSx Quick Start

## ⚠️ Wichtig: Richtige Python-Version verwenden!

FEniCSx ist nur in der Conda-Umgebung installiert (Python 3.12), **nicht** in pyenv (Python 3.13).

## ✅ Richtig ausführen:

### Methode 1: Skript direkt ausführen (empfohlen)
```bash
./fenicsx_example.py
```
Der Shebang (`#!/opt/miniconda3/envs/fenicsx-env/bin/python`) sorgt dafür, dass das richtige Python verwendet wird.

### Methode 2: Mit Wrapper-Skript
```bash
./run_fenicsx.sh fenicsx_example.py
```

### Methode 3: Vollständiger Python-Pfad
```bash
/opt/miniconda3/envs/fenicsx-env/bin/python fenicsx_example.py
```

## ❌ Falsch:

```bash
# NICHT so - verwendet pyenv Python 3.13:
python fenicsx_example.py
/Users/thomashoffbauer/.pyenv/versions/3.13.11/bin/python fenicsx_example.py
```

## Verfügbare Skripte:

- `fenicsx_example.py` - Einfaches Beispiel
- `fenicsx_standalone.py` - Vollständiges Beispiel mit allen Schritten
- `fenicsx_interaktiv.py` - Für interaktive Nutzung

## Testen:

```bash
# Prüfen Sie, welches Python verwendet wird:
head -1 fenicsx_example.py
# Sollte zeigen: #!/opt/miniconda3/envs/fenicsx-env/bin/python

# Führen Sie es aus:
./fenicsx_example.py
```

## Hilfe bei Problemen:

Falls `ModuleNotFoundError: No module named 'dolfinx'`:
- Sie verwenden das falsche Python
- Verwenden Sie eine der Methoden oben (1, 2 oder 3)
