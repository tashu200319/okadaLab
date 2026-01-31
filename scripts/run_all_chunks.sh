#!/bin/bash
# data/chunks/chunk_1.csv から順番に1つずつ解析を実行するランナー（標準版）
# - 生成物（ログ含む）は output/ 配下へ
# - chunksの置き場は data/chunks を優先し、無ければ従来の chunks/ にフォールバック
#
# 例:
#   ./scripts/run_all_chunks.sh
#   ./scripts/run_all_chunks.sh --start 1 --end 11
#
# ※ macOS標準のbash(3.2系)でも動くように mapfile は使わない

set -euo pipefail
cd "$(dirname "$0")/.."

SEQ_RATIO=20
MAX_PDBS=50
WORKERS=""
BATCH_SIZE=""
NO_PARALLEL=0
NO_HEATMAP=0
SKIP_PROCESSED=1

START_NUM=""
END_NUM=""

# 標準の置き場（存在する方を自動選択）
DEFAULT_CHUNKS_DIR="data/chunks"
FALLBACK_CHUNKS_DIR="chunks"

CHUNKS_DIR=""
OUTPUT_DIR="output"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/run_all_chunks.sh [options]

Options:
  --seq-ratio N        Sequence ratio percentage (default: 20)
  --max-pdbs N         Maximum PDB entries per ID (default: 50)
  --workers N          Number of parallel workers (default: main.py default)
  --batch-size N       Batch write size (default: main.py default)
  --no-parallel        Disable parallel processing
  --no-heatmap         Disable heatmap generation
  --no-skip            Reprocess all IDs (do not skip processed)
  --start N            Start chunk number (e.g. 1)
  --end N              End chunk number (e.g. 11)
  --chunks-dir PATH    Chunks directory (overrides auto-detect)
  --output-dir PATH    Output directory passed to main.py (default: output)

Examples:
  ./scripts/run_all_chunks.sh
  ./scripts/run_all_chunks.sh --start 1 --end 3 --seq-ratio 20 --max-pdbs 50
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --seq-ratio) SEQ_RATIO="$2"; shift 2 ;;
    --max-pdbs) MAX_PDBS="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --batch-size) BATCH_SIZE="$2"; shift 2 ;;
    --no-parallel) NO_PARALLEL=1; shift ;;
    --no-heatmap) NO_HEATMAP=1; shift ;;
    --no-skip) SKIP_PROCESSED=0; shift ;;
    --start) START_NUM="$2"; shift 2 ;;
    --end) END_NUM="$2"; shift 2 ;;
    --chunks-dir) CHUNKS_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${CHUNKS_DIR}" ]]; then
  if [[ -d "${DEFAULT_CHUNKS_DIR}" ]]; then
    CHUNKS_DIR="${DEFAULT_CHUNKS_DIR}"
  else
    CHUNKS_DIR="${FALLBACK_CHUNKS_DIR}"
  fi
fi

# ログ保存先（生成物は output/ に寄せる）
TS="$(date +%Y%m%d_%H%M%S)"
LOG_ROOT="${OUTPUT_DIR%/}/logs/run_all_chunks_${TS}"
mkdir -p "$LOG_ROOT"

echo "=== run_all_chunks.sh start: $(date) ==="
echo "chunks_dir: $CHUNKS_DIR"
echo "output_dir: $OUTPUT_DIR"
echo "logs: $LOG_ROOT"
echo "params: seq_ratio=$SEQ_RATIO, max_pdbs=$MAX_PDBS, skip_processed=$SKIP_PROCESSED"
echo ""

# chunk_*.csv を数字順に並べる（chunk_10がchunk_2より先になる問題を避ける）
CHUNK_FILES=()
while IFS= read -r line; do
  [[ -n "$line" ]] && CHUNK_FILES+=("$line")
done < <(ls -1 "${CHUNKS_DIR%/}"/chunk_*.csv 2>/dev/null | sed -E 's/.*chunk_([0-9]+)\.csv/\1\t&/' | sort -n | cut -f2-)

if [[ ${#CHUNK_FILES[@]} -eq 0 ]]; then
  echo "❌ No chunk files found in: $CHUNKS_DIR (expected chunk_*.csv)" >&2
  exit 1
fi

run_one() {
  local chunk_file="$1"
  local chunk_base
  chunk_base="$(basename "$chunk_file")"
  local chunk_num
  chunk_num="$(echo "$chunk_base" | sed -E 's/^chunk_([0-9]+)\.csv$/\1/')"

  # start/end フィルタ
  if [[ -n "$START_NUM" && "$chunk_num" -lt "$START_NUM" ]]; then
    return 0
  fi
  if [[ -n "$END_NUM" && "$chunk_num" -gt "$END_NUM" ]]; then
    return 0
  fi

  local log_file="${LOG_ROOT}/chunk_${chunk_num}.log"
  echo "=== Chunk ${chunk_num} start: $(date) ===" | tee "$log_file"
  echo "file: $chunk_file" | tee -a "$log_file"
  echo "log:  $log_file" | tee -a "$log_file"
  echo "" | tee -a "$log_file"

  # main.pyへ渡す引数を組み立て
  local args=( -u main.py
    --file "$chunk_file"
    --seq-ratio "$SEQ_RATIO"
    --max-pdbs "$MAX_PDBS"
    --output-dir "$OUTPUT_DIR"
  )

  if [[ -n "$WORKERS" ]]; then
    args+=( --workers "$WORKERS" )
  fi
  if [[ -n "$BATCH_SIZE" ]]; then
    args+=( --batch-size "$BATCH_SIZE" )
  fi
  if [[ "$NO_PARALLEL" -eq 1 ]]; then
    args+=( --no-parallel )
  fi
  if [[ "$NO_HEATMAP" -eq 1 ]]; then
    args+=( --no-heatmap )
  fi
  if [[ "$SKIP_PROCESSED" -eq 0 ]]; then
    args+=( --no-skip )
  fi

  python "${args[@]}" 2>&1 | tee -a "$log_file"

  echo "" | tee -a "$log_file"
  echo "=== Chunk ${chunk_num} done: $(date) ===" | tee -a "$log_file"
}

for f in "${CHUNK_FILES[@]}"; do
  run_one "$f"
done

echo ""
echo "=== run_all_chunks.sh done: $(date) ==="
echo "logs: $LOG_ROOT"

