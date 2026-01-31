#!/usr/bin/env python3
"""
除外リストから特定のIDを削除するスクリプト
"""

import sys
import os
import argparse
from pathlib import Path

# プロジェクトルート
BASE_DIR = Path(__file__).resolve().parents[1]

# 親ディレクトリをパスに追加
sys.path.insert(0, str(BASE_DIR))

from core.config import Config


def remove_from_excluded(ids_to_remove: set, input_file: str = None) -> int:
    """
    除外リストからIDを削除
    
    Parameters
    ----------
    ids_to_remove : set
        削除するUniProt IDのセット
    input_file : str, optional
        入力ファイルから失敗IDを読み込む場合
    
    Returns
    -------
    int
        削除されたID数
    """
    config = Config()
    excluded_file = Path(config.OUTPUT_DIR) / "excluded_ids.txt"
    
    if not excluded_file.exists():
        print(f"⚠️  除外リストファイルが見つかりません: {excluded_file}")
        return 0
    
    # 入力ファイルからIDを読み込む場合
    if input_file:
        import pandas as pd
        if os.path.exists(input_file):
            input_df = pd.read_csv(input_file)
            if 'uniprotid' in input_df.columns:
                input_ids = set(input_df['uniprotid'].dropna().tolist())
            else:
                input_ids = set(input_df.iloc[:, 0].dropna().astype(str).tolist())
                input_ids = {uid.split(',')[0].strip() for uid in input_ids if uid.strip()}
            
            # summary.csvから処理済みIDを取得
            summary_file = Path(config.OUTPUT_DIR) / "summaries" / "summary.csv"
            if summary_file.exists():
                summary_df = pd.read_csv(summary_file)
                if 'seq_ratio' in summary_df.columns:
                    processed_df = summary_df[summary_df['seq_ratio'] == 20.0]
                    processed_ids = set(processed_df['uniprotid'].dropna().tolist())
                else:
                    processed_ids = set(summary_df['uniprotid'].dropna().tolist())
                
                # 失敗ID = 入力ID - 処理済みID
                failed_ids = input_ids - processed_ids
                ids_to_remove = failed_ids
                print(f"📋 入力ファイルから {len(failed_ids)} 個の失敗IDを抽出しました")
    
    if not ids_to_remove:
        print("ℹ️  削除するIDがありません")
        return 0
    
    # 現在の除外リストを読み込み
    current_ids = config.load_excluded_ids()
    
    # 削除対象のIDが除外リストに含まれているか確認
    ids_to_remove = ids_to_remove & current_ids
    
    if not ids_to_remove:
        print("ℹ️  指定されたIDは除外リストに含まれていません")
        return 0
    
    # 除外リストから削除
    remaining_ids = current_ids - ids_to_remove
    
    # ファイルを書き直し
    try:
        with open(excluded_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # コメント行と削除対象のIDを除外
        new_lines = []
        for line in lines:
            line_stripped = line.strip()
            # コメント行は保持
            if line_stripped.startswith('#'):
                new_lines.append(line)
            # 空行は保持
            elif not line_stripped:
                new_lines.append(line)
            else:
                # ID行をチェック
                line_id = line_stripped.split()[0] if line_stripped.split() else ""
                if line_id not in ids_to_remove:
                    new_lines.append(line)
        
        # ファイルを書き込み
        with open(excluded_file, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        
        print(f"✅ {len(ids_to_remove)} 個のIDを除外リストから削除しました")
        print(f"   削除されたID: {', '.join(sorted(ids_to_remove))}")
        print(f"   残りの除外ID数: {len(remaining_ids)}")
        
        return len(ids_to_remove)
        
    except Exception as e:
        print(f"❌ エラー: {e}")
        return 0


def keep_only_chunk1_ids():
    """
    chunk_1の除外IDのみを保持し、それ以外を削除
    """
    import pandas as pd
    
    config = Config()
    excluded_file = Path(config.OUTPUT_DIR) / "excluded_ids.txt"
    
    if not excluded_file.exists():
        print(f"⚠️  除外リストファイルが見つかりません: {excluded_file}")
        return 0
    
    # chunk_1.csvから失敗IDを取得（標準: data/chunks/。無ければ旧: chunks/）
    chunk1_file = BASE_DIR / "data" / "chunks" / "chunk_1.csv"
    if not chunk1_file.exists():
        chunk1_file = BASE_DIR / "chunks" / "chunk_1.csv"
    if not chunk1_file.exists():
        print(f"❌ chunk_1.csvが見つかりません: {chunk1_file}")
        return 0
    
    input_df = pd.read_csv(chunk1_file)
    if 'uniprotid' in input_df.columns:
        input_ids = set(input_df['uniprotid'].dropna().tolist())
    else:
        input_ids = set(input_df.iloc[:, 0].dropna().astype(str).tolist())
        input_ids = {uid.split(',')[0].strip() for uid in input_ids if uid.strip()}
    
    # summary.csvから処理済みIDを取得
    summary_file = Path(config.OUTPUT_DIR) / "summaries" / "summary.csv"
    if not summary_file.exists():
        print(f"⚠️  summary.csvが見つかりません: {summary_file}")
        return 0
    
    summary_df = pd.read_csv(summary_file)
    if 'seq_ratio' in summary_df.columns:
        processed_df = summary_df[summary_df['seq_ratio'] == 20.0]
        processed_ids = set(processed_df['uniprotid'].dropna().tolist())
    else:
        processed_ids = set(summary_df['uniprotid'].dropna().tolist())
    
    # chunk_1の失敗ID = 入力ID - 処理済みID
    chunk1_failed_ids = input_ids - processed_ids
    
    print(f"📋 chunk_1の失敗ID: {len(chunk1_failed_ids)} 個")
    
    # 現在の除外リストを読み込み
    current_ids = config.load_excluded_ids()
    
    # chunk_1の失敗ID以外を削除
    ids_to_remove = current_ids - chunk1_failed_ids
    
    if not ids_to_remove:
        print("ℹ️  削除するIDはありません（既にchunk_1のIDのみです）")
        return 0
    
    print(f"🗑️  削除するID: {len(ids_to_remove)} 個")
    print(f"   保持するID（chunk_1）: {len(chunk1_failed_ids)} 個")
    
    # ファイルを書き直し
    try:
        with open(excluded_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # コメント行とchunk_1のIDのみを保持
        new_lines = []
        for line in lines:
            line_stripped = line.strip()
            # コメント行は保持
            if line_stripped.startswith('#'):
                new_lines.append(line)
            # 空行は保持
            elif not line_stripped:
                new_lines.append(line)
            else:
                # ID行をチェック（chunk_1のIDのみ保持）
                line_id = line_stripped.split()[0] if line_stripped.split() else ""
                if line_id in chunk1_failed_ids:
                    new_lines.append(line)
        
        # ファイルを書き込み
        with open(excluded_file, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        
        print(f"✅ {len(ids_to_remove)} 個のIDを除外リストから削除しました")
        print(f"   保持されたID（chunk_1）: {len(chunk1_failed_ids)} 個")
        print(f"   削除されたID: {', '.join(sorted(ids_to_remove))}")
        
        return len(ids_to_remove)
        
    except Exception as e:
        print(f"❌ エラー: {e}")
        return 0


def main():
    parser = argparse.ArgumentParser(
        description='除外リストからIDを削除'
    )
    parser.add_argument(
        '--input', '-i',
        type=str,
        help='入力CSVファイルのパス（このファイルの失敗IDを除外リストから削除）'
    )
    parser.add_argument(
        '--ids',
        nargs='+',
        help='削除するUniProt ID（スペース区切り）'
    )
    parser.add_argument(
        '--file',
        type=str,
        help='削除するIDを1行1つずつ書いたテキストファイル'
    )
    parser.add_argument(
        '--keep-chunk1',
        action='store_true',
        help='chunk_1の除外IDのみを保持し、それ以外を削除'
    )
    
    args = parser.parse_args()
    
    # chunk_1のみ保持する場合
    if args.keep_chunk1:
        keep_only_chunk1_ids()
        return
    
    ids_to_remove = set()
    
    if args.input:
        # 入力ファイルから失敗IDを抽出して削除
        remove_from_excluded(set(), args.input)
        return
    elif args.ids:
        ids_to_remove = set(args.ids)
    elif args.file:
        if os.path.exists(args.file):
            with open(args.file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        ids_to_remove.add(line)
        else:
            print(f"❌ ファイルが見つかりません: {args.file}")
            return
    else:
        parser.print_help()
        return
    
    if ids_to_remove:
        remove_from_excluded(ids_to_remove)


if __name__ == "__main__":
    main()
