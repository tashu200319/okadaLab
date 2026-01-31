#!/usr/bin/env python3
"""
failed_ids.csvからIDを抽出して再解析を実行するスクリプト
既にsummary.csvに書き込まれているIDは除外（バッチ書き込み済み）
"""

import pandas as pd
import os
import subprocess
import sys
from pathlib import Path

def extract_ids_from_summary(summary_file, seq_ratio=20):
    """summary.csvからIDを抽出（既に書き込まれているID）"""
    if not os.path.exists(summary_file):
        print(f"⚠️  {summary_file} not found")
        return set()
    
    try:
        df = pd.read_csv(summary_file)
        # seq_ratioが一致するIDを抽出
        if 'seq_ratio' in df.columns:
            df_filtered = df[df['seq_ratio'] == seq_ratio]
            ids = set(df_filtered['uniprotid'].unique().tolist())
        else:
            ids = set(df['uniprotid'].unique().tolist())
        print(f"✅ Extracted {len(ids)} IDs from {summary_file} (already written)")
        return ids
    except Exception as e:
        print(f"❌ Error reading {summary_file}: {e}")
        return set()

def extract_ids_from_failed(failed_file, seq_ratio=20):
    """failed_ids.csvからIDを抽出"""
    if not os.path.exists(failed_file):
        print(f"⚠️  {failed_file} not found")
        return set()
    
    try:
        df = pd.read_csv(failed_file)
        # seq_ratioが一致するIDを抽出
        if 'seq_ratio' in df.columns:
            df_filtered = df[df['seq_ratio'] == seq_ratio]
            ids = set(df_filtered['uniprotid'].unique().tolist())
        else:
            ids = set(df['uniprotid'].unique().tolist())
        print(f"✅ Extracted {len(ids)} IDs from {failed_file}")
        return ids
    except Exception as e:
        print(f"❌ Error reading {failed_file}: {e}")
        return set()

def main():
    base_dir = Path(__file__).resolve().parents[1]

    import argparse
    parser = argparse.ArgumentParser(description="Re-analyze IDs from failed_ids.csv excluding already-written summary IDs")
    parser.add_argument("--summary-file", default=str(base_dir / "output" / "summaries" / "summary.csv"))
    parser.add_argument("--failed-file", default=str(base_dir / "output" / "summaries" / "failed_ids.csv"))
    parser.add_argument("--seq-ratio", type=float, default=20.0)
    parser.add_argument("--max-pdbs", type=int, default=50)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--no-heatmap", action="store_true", help="Disable heatmap generation")
    parser.add_argument("--output-dir", default=str(base_dir / "output"), help="Output directory (passed to main.py)")

    args = parser.parse_args()

    summary_file = args.summary_file
    failed_file = args.failed_file
    seq_ratio = args.seq_ratio
    max_pdbs = args.max_pdbs
    workers = args.workers
    batch_size = args.batch_size
    output_dir = args.output_dir
    
    # IDを抽出
    print("=" * 80)
    print("📋 Extracting IDs from failed_ids.csv for re-analysis")
    print("=" * 80)
    
    # 既にsummary.csvに書き込まれているID（除外対象）
    written_ids = extract_ids_from_summary(summary_file, seq_ratio)
    
    # failed_ids.csvからIDを抽出
    failed_ids = extract_ids_from_failed(failed_file, seq_ratio)
    
    # 再解析対象のIDを決定
    # failed_idsから、既に書き込まれているIDを除外
    target_ids = failed_ids - written_ids
    
    if not target_ids:
        print(f"\n✅ All IDs already written to summary.csv or no IDs to re-analyze")
        print(f"   - Failed IDs: {len(failed_ids)}")
        print(f"   - Already written: {len(written_ids)}")
        return
    
    print(f"\n📊 Total IDs to re-analyze: {len(target_ids)}")
    print(f"   - From failed_ids.csv: {len(failed_ids)}")
    print(f"   - Already written (excluded): {len(written_ids)}")
    print(f"   - Target for re-analysis: {len(target_ids)}")
    
    # IDをファイルに書き出し
    ids_path = Path(output_dir) / "reanalyze" / "reanalyze_ids.txt"
    ids_path.parent.mkdir(parents=True, exist_ok=True)
    with open(ids_path, 'w') as f:
        for uid in sorted(target_ids):
            f.write(f"{uid}\n")
    
    print(f"\n💾 IDs saved to {ids_path}")
    
    # main.pyを実行（--no-skipは不要、既に除外済み）
    print("\n" + "=" * 80)
    print("🚀 Starting re-analysis")
    print("=" * 80)
    
    cmd = [
        sys.executable, "-u", "main.py",
        "--file", str(ids_path),
        "--seq-ratio", str(seq_ratio),
        "--max-pdbs", str(max_pdbs),
        "--workers", str(workers),
        "--batch-size", str(batch_size),
        "--output-dir", str(output_dir),
        "--no-heatmap"
        # --no-skipは不要（既に除外済み）
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print()
    
    # 実行
    subprocess.run(cmd, cwd=str(base_dir))

if __name__ == "__main__":
    main()
