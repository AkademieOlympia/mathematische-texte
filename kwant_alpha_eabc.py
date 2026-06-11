#!/usr/bin/env python3
"""Startet den gleichen Code wie »Kwant Alpha.py« (Pfad ohne Leerzeichen)."""

from pathlib import Path
import runpy


if __name__ == "__main__":
    src = Path(__file__).resolve().parent / "Kwant Alpha.py"
    runpy.run_path(str(src), run_name="__main__")
