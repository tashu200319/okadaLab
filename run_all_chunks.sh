#!/bin/bash
# 互換ラッパー: ルート直下からでも実行できるように scripts/ 版へ転送
set -euo pipefail
cd "$(dirname "$0")"
exec ./scripts/run_all_chunks.sh "$@"
