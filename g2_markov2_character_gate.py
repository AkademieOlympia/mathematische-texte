#!/usr/bin/env python3
"""Markov-2 Character-Spectrum Gate: K̂^{(2)} auf V₄².

Siehe: scripts/g2_markov2_character_gate.py (kanonische Kopie im Repo).
"""

from __future__ import annotations

import runpy
from pathlib import Path

# Delegation an Repo-Skript, falls vorhanden
_repo = Path.home() / "Projects" / "Kepler-Hurrwitz" / "scripts" / "g2_markov2_character_gate.py"
if _repo.exists():
  runpy.run_path(str(_repo), run_name="__main__")
else:
  raise SystemExit("Repo-Skript g2_markov2_character_gate.py nicht gefunden.")
