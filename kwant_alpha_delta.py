#!/usr/bin/env python3
"""Startet »Kwant Alpha Delta.py« ohne Leerzeichen im Aufruf."""

from pathlib import Path
import runpy


if __name__ == "__main__":
    src = Path(__file__).resolve().parent / "Kwant Alpha Delta.py"
    runpy.run_path(str(src), run_name="__main__")
