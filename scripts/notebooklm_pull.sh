#!/usr/bin/env bash
# Holt NotebookLM-Artefakte nach NotebookLM/ im Repo-Root.
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck source=notebooklm_common.sh
source "$(dirname "$0")/notebooklm_common.sh"

OUT_DIR="NotebookLM"
PULL_SLIDES=true
PULL_REPORTS=true
PULL_ALL=false
DRY_RUN=false
SLIDE_FORMAT="pdf"
PDF_REPORT=false

usage() {
  cat <<'EOF'
Usage: bash scripts/notebooklm_pull.sh [OPTIONS]

Lädt Artefakte nach NotebookLM/ (Präfix NotebookLM_*).

Optionen: --slides-only, --reports-only, --all, --pptx, --pdf-report, --dry-run, -h
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --slides-only) PULL_REPORTS=false ;;
    --reports-only) PULL_SLIDES=false ;;
    --all) PULL_ALL=true ;;
    --pptx) SLIDE_FORMAT="pptx" ;;
    --pdf-report) PDF_REPORT=true ;;
    --dry-run) DRY_RUN=true ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unbekannte Option: $1" >&2; usage; exit 1 ;;
  esac
  shift
done

require_nlm
NOTEBOOK_ID=$(resolve_notebook_id "$NOTEBOOK_TITLE")
if [[ -z "${NOTEBOOK_ID:-}" ]]; then
  echo "Notebook nicht gefunden: $NOTEBOOK_TITLE — zuerst: bash scripts/notebooklm_sync.sh" >&2
  exit 3
fi
"$NLM" use "$NOTEBOOK_ID" >/dev/null
mkdir -p "$OUT_DIR"

pull_kind() {
  local label=$1 kind=$2
  shift 2
  local ok=0
  if $DRY_RUN; then
    "$NLM" download "$kind" --dry-run "$@" && ok=1
  else
    "$NLM" download "$kind" --force "$@" && ok=1
  fi
  if [[ "$ok" -eq 0 ]]; then
    echo "  (kein $label — erzeugen: $NLM generate $kind)" >&2
  fi
}

echo "==> Notebook: $NOTEBOOK_TITLE ($NOTEBOOK_ID)"

if $PULL_ALL; then
  TARGET="$OUT_DIR/artefakte"
  mkdir -p "$TARGET"
  $PULL_SLIDES && pull_kind "slide-deck" slide-deck --all "$TARGET" --format "$SLIDE_FORMAT" -n "$NOTEBOOK_ID"
  $PULL_REPORTS && pull_kind "report" report --all "$TARGET" -n "$NOTEBOOK_ID"
else
  $PULL_SLIDES && pull_kind "slide-deck" slide-deck "$OUT_DIR/NotebookLM_slides.$SLIDE_FORMAT" \
    --format "$SLIDE_FORMAT" -n "$NOTEBOOK_ID"
  if $PULL_REPORTS; then
    pull_kind "report" report "$OUT_DIR/NotebookLM_report.md" -n "$NOTEBOOK_ID"
    if $PDF_REPORT && [[ -f "$OUT_DIR/NotebookLM_report.md" ]] && command -v pandoc >/dev/null; then
      pandoc "$OUT_DIR/NotebookLM_report.md" -o "$OUT_DIR/NotebookLM_report.pdf"
    fi
  fi
fi

echo "==> Fertig: $OUT_DIR/"
