#!/usr/bin/env python3
"""Launcher → Kepler-Hurrwitz/scripts/g2_audit_dual_attribution.py"""
from __future__ import annotations

import runpy
from pathlib import Path

_repo = Path.home() / "Projects" / "Kepler-Hurrwitz" / "scripts" / "g2_audit_dual_attribution.py"
if not _repo.exists():
    raise SystemExit(f"Nicht gefunden: {_repo}")
runpy.run_path(str(_repo), run_name="__main__")
