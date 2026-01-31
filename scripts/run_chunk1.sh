#!/bin/bash
# chunk_1 の解析を実行し、ログを chunk1_full.log に保存
# python -u でアンバッファ出力 → tee に即反映

set -e
# scripts/ 配下に置く前提。常にリポジトリルートで実行する
cd "$(dirname "$0")/.."

LOG="output/logs/chunk1_full.log"
mkdir -p "output/logs"
echo "=== Chunk 1 解析開始 $(date) ===" | tee "$LOG"
echo "ログ: $LOG (tail -f $LOG で追跡)" | tee -a "$LOG"
echo "" | tee -a "$LOG"

python -u main.py \
  --file data/chunks/chunk_1.csv \
  --max-pdbs 50 \
  --seq-ratio 20 \
  --output-dir output \
  2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== 完了 $(date) ===" | tee -a "$LOG"
