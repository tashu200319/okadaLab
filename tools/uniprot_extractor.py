#!/usr/bin/env python3
"""
UniProtID抽出ツール（複数ファイル対応）
CSVファイルから polymer entityData 列の UniProt ID を抽出
"""

import sys
import os
# 親ディレクトリをパスに追加（okadaLabプロジェクトと連携する場合）
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import re
import csv
from pathlib import Path
import glob
import argparse

def extract_uniprot_ids(text):
    """
    テキストからUniProtIDを抽出
    
    Parameters:
    -----------
    text : str
        入力テキスト
    
    Returns:
    --------
    set : UniProtIDのセット
    """
    # UniProtIDの正規表現パターン
    uniprot_pattern = r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\b'
    
    matches = re.findall(uniprot_pattern, text)
    uniprot_ids = set()
    
    for match in matches:
        # タプルの場合は最初の要素を取得
        id_value = match[0] if isinstance(match, tuple) else match
        # カンマがある場合は最初のIDのみ
        id_value = id_value.split(',')[0].strip()
        
        # 数字だけのIDやRefSeq IDを除外
        if not re.match(r'^\d+$', id_value) and not id_value.startswith('NR_'):
            uniprot_ids.add(id_value)
    
    return uniprot_ids


def process_csv_file(csv_path, column_index=2):
    """
    CSVファイルから特定の列のUniProtIDを抽出
    
    Parameters:
    -----------
    csv_path : str
        CSVファイルのパス
    column_index : int
        列のインデックス（0始まり）デフォルトは2（C列）
    
    Returns:
    --------
    tuple : (抽出されたUniProtIDのセット, 処理行数)
    """
    all_ids = set()
    line_count = 0
    
    with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        
        # ヘッダーをスキップ
        header = next(reader, None)
        
        # 各行を処理
        for row in reader:
            line_count += 1
            if len(row) > column_index:
                cell_data = row[column_index]
                ids = extract_uniprot_ids(cell_data)
                all_ids.update(ids)
    
    return all_ids, line_count


def process_multiple_files(file_pattern, column_index=2):
    """
    複数のCSVファイルからUniProtIDを抽出
    
    Parameters:
    -----------
    file_pattern : str
        ファイルパターン（例: "*.csv" または "rcsb_*.csv"）
    column_index : int
        列のインデックス
    
    Returns:
    --------
    tuple : (ファイル名をキーとした結果の辞書, 全ファイルからの統合されたUniProtIDのセット)
    """
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"❌ パターン '{file_pattern}' に一致するファイルが見つかりません")
        return {}, set()
    
    all_ids = set()
    results = {}
    
    print(f"\n🔍 {len(files)}個のファイルを処理します...\n")
    
    for i, file_path in enumerate(sorted(files), 1):
        print(f"[{i}/{len(files)}] 処理中: {Path(file_path).name}")
        
        ids, line_count = process_csv_file(file_path, column_index)
        results[Path(file_path).name] = {
            'ids': ids,
            'count': len(ids),
            'lines': line_count
        }
        all_ids.update(ids)
        
        print(f"  ✓ {len(ids)}個のID抽出 ({line_count}行処理)")
    
    return results, all_ids


def save_results(uniprot_ids, output_path='uniprot_ids_unique.txt'):
    """
    結果をファイルに保存（全IDを出力）
    
    Parameters:
    -----------
    uniprot_ids : set or list
        UniProtIDのセットまたはリスト
    output_path : str
        出力ファイルのパス
    """
    sorted_ids = sorted(uniprot_ids)
    
    # 全IDをファイルに保存
    with open(output_path, 'w', encoding='utf-8') as f:
        for id_value in sorted_ids:
            f.write(id_value + '\n')
    
    print(f"\n✅ {len(sorted_ids)}個のユニークなUniProtIDを抽出しました")
    print(f"✅ 全IDを '{output_path}' に保存しました")
    
    # 画面には最初の20件のみ表示（ファイルには全部入ってる）
    print(f"\n📋 最初の20件（ファイルには全{len(sorted_ids)}件保存済み）:")
    for i, id_value in enumerate(sorted_ids[:20], 1):
        print(f"  {i:2d}. {id_value}")
    
    if len(sorted_ids) > 20:
        print(f"  ... 他 {len(sorted_ids) - 20} 件（全て {output_path} に保存済み）")


def save_detailed_report(results, all_ids, output_path='uniprot_extraction_report.txt'):
    """
    詳細レポートを保存（統計情報 + 全ID一覧）
    
    Parameters:
    -----------
    results : dict
        ファイルごとの結果
    all_ids : set
        全てのUniProtID
    output_path : str
        レポートファイルのパス
    """
    sorted_ids = sorted(all_ids)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("UniProtID 抽出レポート\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"処理ファイル数: {len(results)}\n")
        f.write(f"統合後のユニークID数: {len(all_ids)}\n\n")
        
        f.write("ファイル別詳細:\n")
        f.write("-" * 70 + "\n")
        
        for filename, data in sorted(results.items()):
            f.write(f"\n{filename}\n")
            f.write(f"  - 処理行数: {data['lines']:,}\n")
            f.write(f"  - 抽出ID数: {data['count']:,}\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("抽出された全UniProtID一覧:\n")
        f.write("=" * 70 + "\n\n")
        
        # 全IDを番号付きで出力
        for i, uid in enumerate(sorted_ids, 1):
            f.write(f"{i:4d}. {uid}\n")
        
        f.write("\n" + "=" * 70 + "\n")
    
    print(f"✅ 詳細レポート（全{len(sorted_ids)}件のID含む）を '{output_path}' に保存しました")


def main():
    """
    メイン処理
    """
    parser = argparse.ArgumentParser(
        description='Extract UniProt IDs from CSV files (polymer entityData column)'
    )
    parser.add_argument(
        '--pattern',
        help='File pattern (e.g., "rcsb_*.csv" or "../研究室/rcsb_*.csv")'
    )
    parser.add_argument(
        '--file',
        help='Single CSV file path'
    )
    parser.add_argument(
        '--output',
        default='uniprot_ids_unique.txt',
        help='Output file name (default: uniprot_ids_unique.txt)'
    )
    parser.add_argument(
        '--column',
        type=int,
        default=2,
        help='Column index to extract from (0-based, default: 2 for C column)'
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("UniProtID 抽出ツール（複数ファイル対応）")
    print("=" * 70)
    
    # コマンドライン引数がある場合
    if args.pattern:
        print(f"\n🔍 パターン '{args.pattern}' で処理します\n")
        results, all_ids = process_multiple_files(args.pattern, column_index=args.column)
        
        if not all_ids:
            print("❌ UniProtIDが見つかりませんでした")
            return
        
        save_results(all_ids, args.output)
        report_name = args.output.replace('.txt', '_report.txt')
        save_detailed_report(results, all_ids, report_name)
        
    elif args.file:
        print(f"\n🔍 '{args.file}' を処理中...\n")
        
        if not Path(args.file).exists():
            print(f"❌ エラー: ファイル '{args.file}' が見つかりません")
            return
        
        uniprot_ids, line_count = process_csv_file(args.file, column_index=args.column)
        
        if not uniprot_ids:
            print("❌ UniProtIDが見つかりませんでした")
            return
        
        save_results(uniprot_ids, args.output)
        
    else:
        # 対話式
        print("\n処理方法を選択してください:")
        print("1. 単一ファイルを処理")
        print("2. 複数ファイルをまとめて処理（パターン指定）")
        print("3. カレントディレクトリの全CSVファイルを処理")
        
        choice = input("\n選択 [1/2/3]: ").strip()
        
        if choice == "1":
            # 単一ファイル処理
            file_path = input("\nファイルパス: ").strip()
            
            if not Path(file_path).exists():
                print(f"❌ エラー: ファイル '{file_path}' が見つかりません")
                return
            
            print(f"\n🔍 '{file_path}' を処理中...")
            uniprot_ids, line_count = process_csv_file(file_path, column_index=args.column)
            
            if not uniprot_ids:
                print("❌ UniProtIDが見つかりませんでした")
                return
            
            output_path = input("\n出力ファイル名（Enter で 'uniprot_ids_unique.txt'）: ").strip()
            if not output_path:
                output_path = 'uniprot_ids_unique.txt'
            
            save_results(uniprot_ids, output_path)
        
        elif choice == "2":
            # パターン指定で複数ファイル処理
            pattern = input("\nファイルパターン（例: rcsb_*.csv）: ").strip()
            
            results, all_ids = process_multiple_files(pattern, column_index=args.column)
            
            if not all_ids:
                print("❌ UniProtIDが見つかりませんでした")
                return
            
            output_path = input("\n出力ファイル名（Enter で 'uniprot_ids_all_unique.txt'）: ").strip()
            if not output_path:
                output_path = 'uniprot_ids_all_unique.txt'
            
            save_results(all_ids, output_path)
            save_detailed_report(results, all_ids)
        
        elif choice == "3":
            # カレントディレクトリの全CSVファイル
            print("\n現在のディレクトリの全CSVファイルを処理します")
            
            results, all_ids = process_multiple_files("*.csv", column_index=args.column)
            
            if not all_ids:
                print("❌ UniProtIDが見つかりませんでした")
                return
            
            save_results(all_ids, 'uniprot_ids_all_unique.txt')
            save_detailed_report(results, all_ids)
        
        else:
            print("❌ 無効な選択です")
            return
    
    print("\n" + "=" * 70)
    print("処理完了！")
    print("=" * 70)


if __name__ == "__main__":
    main()