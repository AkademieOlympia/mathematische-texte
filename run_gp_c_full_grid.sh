#!/usr/bin/env bash
# CoCalc / rechenstarke Umgebung: voller Grid m=48, N∈{500,1000,5000}, n_iter=1000
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="${ROOT}/gp_c_runs/full"
STAMP="$(date +%Y%m%d_%H%M%S)"
PY="${PY:-python3}"
SEED="${SEED:-2026}"
N_ITER="${N_ITER:-1000}"
ZEROS_JSON="${ZEROS_JSON:-}"

mkdir -p "${OUT}"

JSON="${OUT}/triple_full_${STAMP}.json"
CSV="${OUT}/triple_full_${STAMP}.csv"
LOG="${OUT}/triple_full_${STAMP}.log"

ARGS=(--full-grid --seed "${SEED}" --n-iter "${N_ITER}" --bonferroni
      --export "${JSON}" --export-csv "${CSV}")
if [[ -n "${ZEROS_JSON}" ]]; then
  ARGS+=(--zeros-json "${ZEROS_JSON}")
fi

echo "=== G_PC full-grid (CoCalc) ==="
echo "  n_iter : ${N_ITER}"
echo "  zeros  : ${ZEROS_JSON:-synthetic}"
echo "  out    : ${OUT}"

"${PY}" "${ROOT}/Topologische Tests.py" "${ARGS[@]}" 2>&1 | tee "${LOG}"

cp -f "${JSON}" "${OUT}/triple_full_latest.json"
cp -f "${CSV}" "${OUT}/triple_full_latest.csv"
cp -f "${LOG}" "${OUT}/triple_full_latest.log"
echo "Fertig: ${OUT}/triple_full_latest.{json,csv,log}"
