#!/bin/bash
# 停止している解析プロセスを終了するスクリプト

echo "=== 停止しているプロセスを確認 ==="
echo ""

# CPU使用率が0%のプロセスを検索
ps aux | grep "python.*main.py" | grep -v grep | while read line; do
    cpu=$(echo $line | awk '{print $3}')
    pid=$(echo $line | awk '{print $2}')
    
    # CPU使用率が0.0%の場合
    if [ "$cpu" = "0.0" ]; then
        echo "⚠️  停止しているプロセスを発見:"
        echo "  PID: $pid, CPU: ${cpu}%"
        echo "  コマンド: $(echo $line | awk '{for(i=11;i<=NF;i++) printf "%s ", $i; print ""}')"
        echo ""
        echo "このプロセスを終了しますか？ (y/n)"
        read -r response
        if [ "$response" = "y" ] || [ "$response" = "Y" ]; then
            kill $pid
            echo "✅ プロセス $pid を終了しました"
        else
            echo "⏸️  プロセスの終了をスキップしました"
        fi
        echo ""
    fi
done

echo "=== 確認完了 ==="
