#!/bin/bash
# main.pyを実行しているすべてのプロセスを強制終了するスクリプト

echo "=== main.pyプロセスを検索中 ==="
echo ""

# main.pyを実行しているすべてのプロセスを検索
PIDS=$(ps aux | grep "python.*main.py" | grep -v grep | awk '{print $2}')

if [ -z "$PIDS" ]; then
    echo "✅ main.pyを実行しているプロセスはありません"
    exit 0
fi

echo "⚠️  以下のプロセスが見つかりました:"
ps aux | grep "python.*main.py" | grep -v grep | while read line; do
    pid=$(echo $line | awk '{print $2}')
    cpu=$(echo $line | awk '{print $3}')
    mem=$(echo $line | awk '{print $4}')
    cmd=$(echo $line | awk '{for(i=11;i<=NF;i++) printf "%s ", $i; print ""}')
    echo "  PID: $pid, CPU: ${cpu}%, MEM: ${mem}%"
    echo "  CMD: ${cmd:0:80}..."
done

echo ""
echo "これらのプロセスをすべて強制終了しますか？ (y/n)"
read -r response

if [ "$response" = "y" ] || [ "$response" = "Y" ]; then
    for pid in $PIDS; do
        echo "🛑 プロセス $pid を終了中..."
        kill -9 $pid 2>/dev/null
        if [ $? -eq 0 ]; then
            echo "  ✅ プロセス $pid を終了しました"
        else
            echo "  ⚠️  プロセス $pid の終了に失敗しました（既に終了している可能性があります）"
        fi
    done
    echo ""
    echo "✅ すべてのプロセスを終了しました"
else
    echo "⏸️  プロセスの終了をキャンセルしました"
fi

echo ""
echo "=== 確認: 残っているプロセス ==="
REMAINING=$(ps aux | grep "python.*main.py" | grep -v grep)
if [ -z "$REMAINING" ]; then
    echo "✅ main.pyを実行しているプロセスはありません"
else
    echo "⚠️  以下のプロセスがまだ実行中です:"
    echo "$REMAINING"
fi
