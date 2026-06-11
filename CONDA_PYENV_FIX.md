# Lösung: Conda-Umgebung mit pyenv

## Problem
Auch nach `conda activate fenicsx-env` verwendet `python` noch pyenv (Python 3.13.11) statt Conda (Python 3.12.8).

## Lösung 1: Vollständigen Python-Pfad verwenden (empfohlen)

```bash
# Direkt ausführen:
/opt/miniconda3/envs/fenicsx-env/bin/python fenicsx_example.py

# Oder mit dem Wrapper-Skript:
./run_fenicsx.sh fenicsx_example.py
```

## Lösung 2: Conda richtig initialisieren

Fügen Sie diese Zeile zu Ihrer `~/.zshrc` hinzu:

```bash
# Conda initialisieren
eval "$(/opt/miniconda3/bin/conda shell.zsh hook)"
```

Dann:
```bash
source ~/.zshrc
conda activate fenicsx-env
which python  # Sollte jetzt /opt/miniconda3/envs/fenicsx-env/bin/python zeigen
```

## Lösung 3: Alias erstellen

Fügen Sie zu `~/.zshrc` hinzu:

```bash
alias python-fenicsx="/opt/miniconda3/envs/fenicsx-env/bin/python"
```

Dann:
```bash
source ~/.zshrc
python-fenicsx fenicsx_example.py
```

## Überprüfung

```bash
# Welches Python wird verwendet?
which python
python --version

# Sollte zeigen:
# /opt/miniconda3/envs/fenicsx-env/bin/python
# Python 3.12.8
```

## Für Jupyter Notebooks

Das Notebook sollte automatisch den richtigen Kernel verwenden, wenn Sie "Python 3.12 (fenicsx-env)" auswählen.
