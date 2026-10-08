#!/usr/bin/env python3
"""Launcher → Kepler-Hurrwitz/scripts/g2_arith_susy.py"""
from __future__ import annotations

import runpy
from pathlib import Path

_repo = Path.home() / "Projects" / "Kepler-Hurrwitz" / "scripts" / "g2_arith_susy.py"
if not _repo.exists():
    raise SystemExit(f"Nicht gefunden: {_repo}")
runpy.run_path(str(_repo), run_name="__main__")
