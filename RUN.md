# Rundweg.py starten

Das Skript nutzt **SageMath** (`sage.all`: exakte Matrizen über `QQ`). Es läuft **nicht** mit normalem `python3` / pyenv.

## Empfohlener Start

```bash
cd "/Users/thomashoffbauer/Desktop/Mathematische Texte"
./run_rundweg.sh
```

Alternativ (gleichwertig):

```bash
sage -python "/Users/thomashoffbauer/Desktop/Mathematische Texte/Rundweg.py"
```

## Sage auf diesem Mac

- Befehl: `sage` (z. B. über Antigravity: `~/.antigravity/antigravity/bin/sage`)
- Version prüfen: `sage --version`

Falls `sage` fehlt: SageMath installieren (z. B. [sagemath.org/download](https://www.sagemath.org/download.html) oder Homebrew `brew install sagemath`) und sicherstellen, dass `sage` im `PATH` liegt.

## Warum nicht `python3 Rundweg.py`?

`ModuleNotFoundError: No module named 'sage'` — das Sage-Modul ist nur im Interpreter `sage -python` verfügbar, nicht in pyenv-Standard-Python.

## Erwartete Ausgabe (Kurz)

Nach dem Start erscheinen u. a. die Matrix **Π_Γ**, der Block `=== ARITHMETISCHE STRUKTUR VON PI_GAMMA ÜBER Q ===` (λ=0, dim 1, Basis `(0,1,0,-1)`; λ=1, dim 3) und die Demo mit lokalem bzw. projiziertem Zustand. Tabellen für das Paper: `PAPER_Q_STRUCTURE.md`.
