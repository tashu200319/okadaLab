#!/usr/bin/env python3
"""
summary.csvを参照して失敗IDを抽出し、除外リストに追加するスクリプト
"""

import sys
import os
import argparse
from pathlib import Path
import pandas as pd

# プロジェクトルート
BASE_DIR = Path(__file__).resolve().parents[1]

# 親ディレクトリをパスに追加
sys.path.insert(0, str(BASE_DIR))

from core.config import Config


def extract_failed_ids(input_file: str, summary_file: str, seq_ratio: float) -> set:
    """
    入力ファイルとsummary.csvを比較して失敗IDを抽出
    
    Parameters
    ----------
    input_file : str
        入力CSVファイルのパス
    summary_file : str
        summary.csvのパス
    seq_ratio : float
        seq_ratio値
    
    Returns
    -------
    set
        失敗したUniProt IDのセット
    """
    # 入力ファイルからIDを読み込み
    if not os.path.exists(input_file):
        print(f"❌ 入力ファイルが見つかりません: {input_file}")
        return set()
    
    input_df = pd.read_csv(input_file)
    
    # uniprotid列を探す
    if 'uniprotid' in input_df.columns:
        input_ids = set(input_df['uniprotid'].dropna().tolist())
    elif len(input_df.columns) > 0:
        # 最初の列をuniprotidとして扱う
        input_ids = set(input_df.iloc[:, 0].dropna().astype(str).tolist())
        # カンマ区切りの場合は最初の部分だけを取得
        input_ids = {uid.split(',')[0].strip() for uid in input_ids if uid.strip()}
    else:
        print(f"❌ 入力ファイルにIDが見つかりません: {input_file}")
        return set()
    
    print(f"📋 入力ファイルから {len(input_ids)} 個のIDを読み込みました")
    
    # summary.csvから処理済みIDを読み込み
    if not os.path.exists(summary_file):
        print(f"⚠️  summary.csvが見つかりません: {summary_file}")
        print(f"   すべての入力IDが失敗IDとして扱われます")
        return input_ids
    
    try:
        summary_df = pd.read_csv(summary_file)
        
        # seq_ratioでフィルタリング
        if 'seq_ratio' in summary_df.columns:
            processed_df = summary_df[summary_df['seq_ratio'] == seq_ratio]
            processed_ids = set(processed_df['uniprotid'].dropna().tolist())
        else:
            # seq_ratio列がない場合はすべてを処理済みとして扱う
            processed_ids = set(summary_df['uniprotid'].dropna().tolist())
        
        print(f"✅ summary.csvから {len(processed_ids)} 個の処理済みIDを読み込みました (seq_ratio={seq_ratio})")
        
    except Exception as e:
        print(f"⚠️  summary.csvの読み込みエラー: {e}")
        return input_ids
    
    # 失敗ID = 入力ID - 処理済みID
    failed_ids = input_ids - processed_ids
    
    return failed_ids


def main():
    parser = argparse.ArgumentParser(
        description='summary.csvを参照して失敗IDを抽出し、除外リストに追加'
    )
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='入力CSVファイルのパス（例: data/chunks/chunk_1.csv）'
    )
    parser.add_argument(
        '--summary', '-s',
        type=str,
        default='output/summaries/summary.csv',
        help='summary.csvのパス（デフォルト: output/summaries/summary.csv）'
    )
    parser.add_argument(
        '--seq-ratio', '-r',
        type=float,
        default=20.0,
        help='seq_ratio値（デフォルト: 20.0）'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='実際には追加せず、結果のみ表示'
    )
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("📊 失敗ID抽出スクリプト")
    print("=" * 80)
    print(f"入力ファイル: {args.input}")
    print(f"summary.csv: {args.summary}")
    print(f"seq_ratio: {args.seq_ratio}")
    print("=" * 80)
    
    # 失敗IDを抽出
    failed_ids = extract_failed_ids(args.input, args.summary, args.seq_ratio)
    
    if not failed_ids:
        print("\n✅ 失敗IDはありませんでした")
        return
    
    print(f"\n📋 失敗ID: {len(failed_ids)} 個")
    print("\n失敗ID一覧:")
    for fid in sorted(failed_ids):
        print(f"  {fid}")
    
    if args.dry_run:
        print("\n🔍 ドライラン: 除外リストには追加しませんでした")
        return
    
    # 除外リストに追加
    config = Config()
    added_count = config.add_to_excluded_ids(failed_ids)
    
    print("\n" + "=" * 80)
    if added_count > 0:
        print(f"✅ {added_count} 個の失敗IDを除外リストに追加しました")
        print(f"   ファイル: {config.OUTPUT_DIR}/excluded_ids.txt")
    else:
        print("ℹ️  追加されたIDはありませんでした（既に除外リストに含まれています）")
    print("=" * 80)


if __name__ == "__main__":
    main()
