#!/bin/bash
# 解析の実行状況を確認するスクリプト

set -e
cd "$(dirname "$0")/../.."

echo "=== 解析プロセス確認 ==="
echo ""

# 1. 実行中のプロセスを確認
echo "📊 実行中のプロセス:"
ps aux | grep "python.*main.py" | grep -v grep || echo "  ❌ 実行中のプロセスはありません"
echo ""

# 2. 最新のログを確認
LATEST_LOG=$(ls -1t output/logs/*.log 2>/dev/null | head -1 || true)
if [ -n "$LATEST_LOG" ] && [ -f "$LATEST_LOG" ]; then
    echo "📄 最新のログ（最後の10行）: $LATEST_LOG"
    tail -10 "$LATEST_LOG"
    echo ""
    
    # ログの最終更新時刻
    echo "🕐 ログの最終更新時刻:"
    stat -f "%Sm" -t "%Y-%m-%d %H:%M:%S" "$LATEST_LOG" 2>/dev/null || \
    ls -l "$LATEST_LOG" | awk '{print $6, $7, $8}'
    echo ""
fi

# 3. 出力ファイルの状態を確認
if [ -f "output/summaries/summary.csv" ]; then
    echo "📊 処理済みID数:"
    # ヘッダー行を除いた行数をカウント
    tail -n +2 output/summaries/summary.csv 2>/dev/null | wc -l | xargs echo "  "
    echo ""
    
    echo "📄 最新の処理結果（最後の3行）:"
    tail -3 output/summaries/summary.csv 2>/dev/null || echo "  ファイルが空です"
    echo ""
fi

# 4. 失敗IDの確認
if [ -f "output/summaries/failed_ids.csv" ]; then
    echo "❌ 失敗ID数:"
    tail -n +2 output/summaries/failed_ids.csv 2>/dev/null | wc -l | xargs echo "  "
    echo ""
fi

# 5. CPU/メモリ使用率（プロセスが見つかった場合）
PID=$(ps aux | grep "python.*main.py" | grep -v grep | awk '{print $2}' | head -1)
if [ ! -z "$PID" ]; then
    echo "💻 リソース使用状況:"
    ps -p $PID -o %cpu,%mem,etime,command 2>/dev/null || echo "  プロセス情報を取得できませんでした"
    echo ""
fi

echo "=== 確認完了 ==="
