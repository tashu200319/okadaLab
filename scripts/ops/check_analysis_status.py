#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
解析の実行状況を確認するPythonスクリプト
"""

import os
import subprocess
import time
from pathlib import Path
from datetime import datetime

def check_process():
    """実行中のプロセスを確認"""
    try:
        result = subprocess.run(
            ["ps", "aux"],
            capture_output=True,
            text=True
        )
        lines = result.stdout.split('\n')
        main_processes = [line for line in lines if 'python' in line and 'main.py' in line and 'grep' not in line]
        
        if main_processes:
            print("✅ 実行中のプロセス:")
            for proc in main_processes:
                parts = proc.split()
                if len(parts) >= 11:
                    pid = parts[1]
                    cpu = float(parts[2])
                    mem = float(parts[3])
                    cmd = ' '.join(parts[10:])
                    status = "⚠️  停止中（CPU 0%）" if cpu == 0.0 else "✅ 実行中"
                    print(f"  PID: {pid}, CPU: {cpu}%, MEM: {mem}%, {status}")
                    print(f"  CMD: {cmd[:80]}")
                    
                    # CPU使用率が0%の場合は警告
                    if cpu == 0.0:
                        print(f"  ⚠️  警告: CPU使用率が0%です。プロセスが停止している可能性があります。")
                        print(f"  終了する場合: kill {pid}")
            return True
        else:
            print("❌ 実行中のプロセスはありません")
            return False
    except Exception as e:
        print(f"⚠️  プロセス確認エラー: {e}")
        return False

def check_log_file(log_file="chunk1_full.log"):
    """ログファイルを確認"""
    if not os.path.exists(log_file):
        print(f"⚠️  ログファイルが見つかりません: {log_file}")
        return
    
    print(f"\n📄 ログファイル: {log_file}")
    
    # ファイルサイズと更新時刻
    stat = os.stat(log_file)
    size = stat.st_size
    mtime = datetime.fromtimestamp(stat.st_mtime)
    now = datetime.now()
    age = (now - mtime).total_seconds()
    
    print(f"  サイズ: {size:,} bytes")
    print(f"  最終更新: {mtime.strftime('%Y-%m-%d %H:%M:%S')} ({int(age)}秒前)")
    
    if age > 300:  # 5分以上更新されていない
        print("  ⚠️  警告: 5分以上更新されていません（解析が停止している可能性があります）")
    else:
        print("  ✅ 最近更新されています")
    
    # 最後の10行を表示
    try:
        with open(log_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            if lines:
                print(f"\n  最新のログ（最後の10行）:")
                for line in lines[-10:]:
                    print(f"    {line.rstrip()}")
    except Exception as e:
        print(f"  ⚠️  ログ読み込みエラー: {e}")

def check_output_files():
    """出力ファイルを確認"""
    summary_file = "output/summaries/summary.csv"
    failed_file = "output/summaries/failed_ids.csv"
    
    print("\n📊 出力ファイル:")
    
    # summary.csv
    if os.path.exists(summary_file):
        try:
            with open(summary_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                count = len(lines) - 1 if len(lines) > 1 else 0  # ヘッダーを除く
                print(f"  ✅ summary.csv: {count}件の処理済みID")
                
                if count > 0 and len(lines) > 1:
                    print(f"    最新のID: {lines[-1].split(',')[0] if ',' in lines[-1] else 'N/A'}")
        except Exception as e:
            print(f"  ⚠️  summary.csv読み込みエラー: {e}")
    else:
        print(f"  ⚠️  summary.csvが見つかりません")
    
    # failed_ids.csv
    if os.path.exists(failed_file):
        try:
            with open(failed_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                count = len(lines) - 1 if len(lines) > 1 else 0
                print(f"  ❌ failed_ids.csv: {count}件の失敗ID")
        except Exception as e:
            print(f"  ⚠️  failed_ids.csv読み込みエラー: {e}")
    else:
        print(f"  ℹ️  failed_ids.csvはまだ作成されていません")

def main():
    print("=" * 60)
    print("🔍 解析状況確認")
    print("=" * 60)
    print()
    
    # プロセス確認
    is_running = check_process()
    
    # ログファイル確認
    check_log_file()
    
    # 出力ファイル確認
    check_output_files()
    
    print("\n" + "=" * 60)
    if is_running:
        print("✅ 解析は実行中です")
    else:
        print("❌ 解析は実行されていません")
    print("=" * 60)

if __name__ == "__main__":
    main()
