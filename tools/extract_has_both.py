#!/usr/bin/env python3
"""
manual_xray_em_check.csv から has_both=true のIDを抽出
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd
from pathlib import Path


def extract_has_both_ids(output_file: str = "has_both_ids.txt",
                         show_all: bool = True):
    """
    has_both列がtrueのIDを抽出
    
    Parameters
    ----------
    output_file : str
        出力テキストファイル名
    show_all : bool
        全IDを表示するか（Falseなら最初の20件のみ）
    """
    
    try:
        # 親ディレクトリのoutputフォルダを参照
        base_dir = Path(__file__).resolve().parent.parent
        csv_file = base_dir / "output" / "manual_xray_em_check.csv"
        
        print(f"📂 Base directory: {base_dir}")
        print(f"📄 CSV file: {csv_file}")
        
        if not csv_file.exists():
            print(f"❌ Error: File not found at {csv_file}")
            return None
        
        # CSVを読み込み
        df = pd.read_csv(csv_file)
        
        print(f"✓ Loaded: {csv_file.name}")
        print(f"Total rows: {len(df):,}")
        print(f"Columns: {', '.join(df.columns)}")
        
        # has_both列の確認
        if 'has_both' not in df.columns:
            print(f"❌ Error: 'has_both' column not found")
            return None
        
        # has_both=Trueの行を抽出（ブール値として直接比較）
        has_both_df = df[df['has_both'] == True]
        
        print(f"\n✓ Found {len(has_both_df):,} entries with has_both=True")
        print(f"  ({len(has_both_df)/len(df)*100:.1f}% of total)")
        
        if len(has_both_df) == 0:
            print("⚠️  No entries found with has_both=True")
            return None
        
        # uniprotid列を使用
        ids = has_both_df['uniprotid'].tolist()
        
        # 統計情報表示
        print(f"\n📊 Statistics:")
        print(f"  X-ray count range: {has_both_df['xray_count'].min()}-{has_both_df['xray_count'].max()}")
        print(f"  EM count range: {has_both_df['em_count'].min()}-{has_both_df['em_count'].max()}")
        print(f"  Avg X-ray: {has_both_df['xray_count'].mean():.1f}")
        print(f"  Avg EM: {has_both_df['em_count'].mean():.1f}")
        
        # テキストファイルに保存（outputディレクトリに）
        output_path = base_dir / "output" / output_file
        with open(output_path, 'w') as f:
            for uid in ids:
                f.write(f"{uid}\n")
        
        print(f"\n💾 Saved to: {output_path}")
        print(f"Total IDs: {len(ids):,}")
        
        # 全件表示 or プレビュー
        if show_all:
            print(f"\n📋 All {len(ids)} IDs:")
            # 20件ごとに区切って表示
            for i in range(0, len(ids), 20):
                batch = ids[i:i+20]
                print(f"\n  [{i+1}-{min(i+20, len(ids))}]")
                for j, uid in enumerate(batch, start=i+1):
                    # X-ray/EM情報も表示
                    row = has_both_df[has_both_df['uniprotid'] == uid].iloc[0]
                    print(f"    {j:4d}. {uid:12s} (X-ray:{row['xray_count']:3d}, EM:{row['em_count']:3d})")
        else:
            print(f"\n📋 Preview (first 20):")
            for i, uid in enumerate(ids[:20], 1):
                row = has_both_df[has_both_df['uniprotid'] == uid].iloc[0]
                print(f"  {i:2d}. {uid:12s} (X-ray:{row['xray_count']:3d}, EM:{row['em_count']:3d})")
            
            if len(ids) > 20:
                print(f"\n  ... and {len(ids) - 20:,} more IDs")
        
        # Pythonリスト形式（1行80文字以内で整形）
        list_output_path = base_dir / "output" / output_file.replace('.txt', '_list.txt')
        with open(list_output_path, 'w') as f:
            f.write('# has_both=True UniProt IDs\n')
            f.write(f'# Total: {len(ids)} IDs\n')
            f.write('# Generated from: manual_xray_em_check.csv\n\n')
            f.write('uniprot_ids = [\n')
            
            # 1行に複数ID（見やすく整形）
            line = '    '
            for i, uid in enumerate(ids):
                line += f'"{uid}", '
                # 1行が長くなったら改行
                if len(line) > 70 or i == len(ids) - 1:
                    if i == len(ids) - 1:
                        line = line.rstrip(', ') + '\n'
                    f.write(line + '\n')
                    line = '    '
            
            f.write(']\n')
        
        print(f"\n🐍 Python list format: {list_output_path}")
        print(f"   (Copy-paste ready for main.py)")
        
        # CSV形式でも保存（詳細情報付き）
        csv_output_path = base_dir / "output" / output_file.replace('.txt', '_details.csv')
        has_both_df[['uniprotid', 'xray_count', 'em_count', 'xray_pdbs', 'em_pdbs']].to_csv(
            csv_output_path, index=False
        )
        print(f"\n📊 Detailed CSV: {csv_output_path}")
        
        return ids
        
    except FileNotFoundError:
        print(f"❌ Error: File '{csv_file}' not found")
        print("Please check the file path")
        return None
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """メイン処理"""
    print("=" * 80)
    print("🔍 Extracting IDs with has_both=True from manual_xray_em_check.csv")
    print("=" * 80)
    print()
    
    # 実行
    ids = extract_has_both_ids(
        output_file="has_both_ids.txt",
        show_all=True  # 全IDを表示
    )
    
    if ids:
        print("\n" + "=" * 80)
        print(f"✅ Extraction completed! {len(ids):,} IDs extracted")
        print("=" * 80)
        print("\n📁 Output files:")
        print("  1. has_both_ids.txt          - Simple list (one ID per line)")
        print("  2. has_both_ids_list.txt     - Python list format (copy-paste ready)")
        print("  3. has_both_ids_details.csv  - Full details (X-ray/EM counts + PDB IDs)")
        print("\n💡 Tip: Copy from has_both_ids_list.txt to use in main.py")
    else:
        print("\n❌ No IDs extracted")


if __name__ == "__main__":
    main()