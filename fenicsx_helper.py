#!/opt/miniconda3/envs/fenicsx-env/bin/python
"""
Helper-Skript für FEniCSx - wird von SageMath aus aufgerufen
"""
import sys
import json
import os

# MPI-Umgebungsvariablen
os.environ['OMPI_COMM_WORLD_SIZE'] = '1'
os.environ['OMPI_COMM_WORLD_RANK'] = '0'

from mpi4py import MPI
from dolfinx import mesh, fem
import numpy as np

# Lese Code von stdin
code = sys.stdin.read()

# Erfasse print-Ausgaben
import io
from contextlib import redirect_stdout

output_buffer = io.StringIO()

try:
    # Führe Code aus und erfasse Ausgaben
    with redirect_stdout(output_buffer):
        exec(code, globals())
    
    output = output_buffer.getvalue()
    result = {"status": "success", "output": output.strip() if output else "Code erfolgreich ausgeführt"}
except Exception as e:
    import traceback
    result = {"status": "error", "error": str(e), "traceback": traceback.format_exc()[:1000]}

# Ausgabe als JSON (vor MPI-Finalisierung)
print("FENICSX_RESULT:" + json.dumps(result))
sys.stdout.flush()
