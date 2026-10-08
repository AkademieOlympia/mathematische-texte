#!/usr/bin/env bash
# Lokaler G_PC^(V4) Prep-Lauf: m=48 Zellen, gleiche W×U-Topologie wie CoCalc-Full.
# N={200,400,800}, n_iter=40 (Override: N_ITER=... ./run_gp_c_local_grid.sh)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="${ROOT}/gp_c_runs/local"
STAMP="$(date +%Y%m%d_%H%M%S)"
PY="${PY:-python3}"
SEED="${SEED:-2026}"
N_ITER="${N_ITER:-40}"

mkdir -p "${OUT}"

JSON="${OUT}/triple_local_${STAMP}.json"
CSV="${OUT}/triple_local_${STAMP}.csv"
LOG="${OUT}/triple_local_${STAMP}.log"
LATEST_JSON="${OUT}/triple_local_latest.json"
LATEST_CSV="${OUT}/triple_local_latest.csv"
LATEST_LOG="${OUT}/triple_local_latest.log"

echo "=== G_PC local-grid ==="
echo "  python : ${PY}"
echo "  seed   : ${SEED}"
echo "  n_iter : ${N_ITER}"
echo "  json   : ${JSON}"
echo "  csv    : ${CSV}"
echo "  log    : ${LOG}"

"${PY}" "${ROOT}/Topologische Tests.py" \
  --local-grid \
  --seed "${SEED}" \
  --n-iter "${N_ITER}" \
  --bonferroni \
  --export "${JSON}" \
  --export-csv "${CSV}" \
  2>&1 | tee "${LOG}"

cp -f "${JSON}" "${LATEST_JSON}"
cp -f "${CSV}" "${LATEST_CSV}"
cp -f "${LOG}" "${LATEST_LOG}"

echo
echo "Fertig. Latest-Symlinks:"
echo "  ${LATEST_JSON}"
echo "  ${LATEST_CSV}"
echo "  ${LATEST_LOG}"
