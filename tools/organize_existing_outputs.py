#!/usr/bin/env python3
"""
既存のscore_detailsとscatterファイルを整理するスクリプト
- AND検索の有無で分類
- X-ray AND EM の両方を持つかどうかで分類
"""

import sys
import os
import shutil
from pathlib import Path
from typing import List, Set

# プロジェクトルート
BASE_DIR = Path(__file__).resolve().parents[1]

# 親ディレクトリをパスに追加
sys.path.insert(0, str(BASE_DIR))

from core.uniprot_handler import UniprotData


def check_has_both_methods(uniprotid: str) -> bool:
    """
    UniProt IDがX-rayとEMの両方を持っているかチェック
    
    Returns
    -------
    bool
        両方持っていればTrue
    """
    try:
        unidata = UniprotData(uniprotid)
        pdbdata = unidata.getpdbdata({"X-ray", "EM"})
        
        xray_count = sum(1 for col in pdbdata.columns 
                        if pdbdata.at['method', col] == 'X-ray')
        em_count = sum(1 for col in pdbdata.columns 
                      if pdbdata.at['method', col] == 'EM')
        
        return xray_count > 0 and em_count > 0
    
    except Exception as e:
        print(f"    ⚠️  Error checking {uniprotid}: {e}")
        return False


def extract_uniprotid_from_filename(filename: str) -> str:
    """
    ファイル名からUniProt IDを抽出
    例: score_details_P01308_20.csv → P01308
    例: score_details_P01308_20_scatter.png → P01308
    """
    basename = Path(filename).stem
    # score_details_ または scatter を除去
    parts = basename.replace('score_details_', '').replace('_scatter', '').split('_')
    return parts[0]


def organize_score_details(
    input_dir: Path = None,
    create_subdirs: bool = True,
    dry_run: bool = False
):
    """
    score_detailsファイルを整理
    
    構造:
    output/score_details/
    ├── with_and_search/        # AND検索済み（X-ray AND EMの両方持つ）
    ├── without_and_search/     # AND検索なし（どちらか片方のみ）
    └── [既存ファイル]
    
    Parameters
    ----------
    input_dir : Path
        score_detailsディレクトリ
    create_subdirs : bool
        サブディレクトリを作成するか
    dry_run : bool
        実際には移動せず、処理内容のみ表示
    """
    if input_dir is None:
        input_dir = BASE_DIR / "output" / "score_details"
    
    input_dir = Path(input_dir)
    
    if not input_dir.exists():
        print(f"❌ Directory not found: {input_dir}")
        return
    
    print("=" * 80)
    print("Organizing score_details files")
    print("=" * 80)
    print(f"Input: {input_dir}")
    print(f"Dry run: {dry_run}")
    print("=" * 80)
    
    # サブディレクトリ作成
    with_and_dir = input_dir / "with_and_search"
    without_and_dir = input_dir / "without_and_search"
    
    if create_subdirs and not dry_run:
        with_and_dir.mkdir(exist_ok=True)
        without_and_dir.mkdir(exist_ok=True)
        print(f"📁 Created: {with_and_dir}")
        print(f"📁 Created: {without_and_dir}")
    
    # score_details_*.csvファイルを取得
    csv_files = list(input_dir.glob("score_details_*.csv"))
    
    if not csv_files:
        print("⚠️  No score_details files found")
        return
    
    print(f"\nFound {len(csv_files)} files")
    print("-" * 80)
    
    with_and_count = 0
    without_and_count = 0
    error_count = 0
    
    for i, csv_file in enumerate(csv_files, 1):
        uniprotid = extract_uniprotid_from_filename(csv_file.name)
        print(f"({i}/{len(csv_files)}) {csv_file.name} [{uniprotid}]", end=" ")
        
        try:
            has_both = check_has_both_methods(uniprotid)
            
            if has_both:
                dest_dir = with_and_dir
                category = "WITH AND"
                with_and_count += 1
            else:
                dest_dir = without_and_dir
                category = "WITHOUT AND"
                without_and_count += 1
            
            dest_file = dest_dir / csv_file.name
            
            print(f"→ {category}")
            
            if not dry_run:
                shutil.move(str(csv_file), str(dest_file))
        
        except Exception as e:
            print(f"→ ERROR: {e}")
            error_count += 1
    
    # サマリー
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"✓ WITH AND search (both X-ray & EM):    {with_and_count}")
    print(f"✓ WITHOUT AND search (X-ray or EM):     {without_and_count}")
    print(f"✗ Errors:                                {error_count}")
    print("=" * 80)
    
    if dry_run:
        print("\n⚠️  DRY RUN: No files were actually moved")
    else:
        print(f"\n✅ Files organized in: {input_dir}")


def organize_scatter(
    input_dir: Path = None,
    create_subdirs: bool = True,
    dry_run: bool = False
):
    """
    scatterファイルを整理
    
    構造:
    output/score_details/scatter/
    ├── with_and_search/        # X-ray AND EMの両方持つ
    ├── without_and_search/     # どちらか片方のみ
    └── [既存ファイル]
    
    Parameters
    ----------
    input_dir : Path
        scatterディレクトリ
    create_subdirs : bool
        サブディレクトリを作成するか
    dry_run : bool
        実際には移動せず、処理内容のみ表示
    """
    if input_dir is None:
        input_dir = BASE_DIR / "output" / "score_details" / "scatter"
    
    input_dir = Path(input_dir)
    
    if not input_dir.exists():
        print(f"❌ Directory not found: {input_dir}")
        return
    
    print("=" * 80)
    print("Organizing scatter plot files")
    print("=" * 80)
    print(f"Input: {input_dir}")
    print(f"Dry run: {dry_run}")
    print("=" * 80)
    
    # サブディレクトリ作成
    with_and_dir = input_dir / "with_and_search"
    without_and_dir = input_dir / "without_and_search"
    
    if create_subdirs and not dry_run:
        with_and_dir.mkdir(exist_ok=True)
        without_and_dir.mkdir(exist_ok=True)
        print(f"📁 Created: {with_and_dir}")
        print(f"📁 Created: {without_and_dir}")
    
    # scatter画像ファイルを取得
    png_files = list(input_dir.glob("score_details_*_scatter.png"))
    
    if not png_files:
        print("⚠️  No scatter plot files found")
        return
    
    print(f"\nFound {len(png_files)} files")
    print("-" * 80)
    
    with_and_count = 0
    without_and_count = 0
    error_count = 0
    
    for i, png_file in enumerate(png_files, 1):
        uniprotid = extract_uniprotid_from_filename(png_file.name)
        print(f"({i}/{len(png_files)}) {png_file.name} [{uniprotid}]", end=" ")
        
        try:
            has_both = check_has_both_methods(uniprotid)
            
            if has_both:
                dest_dir = with_and_dir
                category = "WITH AND"
                with_and_count += 1
            else:
                dest_dir = without_and_dir
                category = "WITHOUT AND"
                without_and_count += 1
            
            dest_file = dest_dir / png_file.name
            
            print(f"→ {category}")
            
            if not dry_run:
                shutil.move(str(png_file), str(dest_file))
        
        except Exception as e:
            print(f"→ ERROR: {e}")
            error_count += 1
    
    # サマリー
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"✓ WITH AND search (both X-ray & EM):    {with_and_count}")
    print(f"✓ WITHOUT AND search (X-ray or EM):     {without_and_count}")
    print(f"✗ Errors:                                {error_count}")
    print("=" * 80)
    
    if dry_run:
        print("\n⚠️  DRY RUN: No files were actually moved")
    else:
        print(f"\n✅ Files organized in: {input_dir}")


def main():
    """メイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Organize existing score_details and scatter files'
    )
    parser.add_argument(
        '--score-details-dir',
        type=str,
        default=None,
        help='score_details directory (default: ./output/score_details/)'
    )
    parser.add_argument(
        '--scatter-dir',
        type=str,
        default=None,
        help='scatter directory (default: ./output/score_details/scatter/)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without actually moving files'
    )
    parser.add_argument(
        '--skip-score-details',
        action='store_true',
        help='Skip organizing score_details files'
    )
    parser.add_argument(
        '--skip-scatter',
        action='store_true',
        help='Skip organizing scatter files'
    )
    
    args = parser.parse_args()
    
    print("\n📂 File Organization Tool")
    print("=" * 80)
    
    # score_detailsを整理
    if not args.skip_score_details:
        organize_score_details(
            input_dir=Path(args.score_details_dir) if args.score_details_dir else None,
            dry_run=args.dry_run
        )
        print("\n")
    
    # scatterを整理
    if not args.skip_scatter:
        organize_scatter(
            input_dir=Path(args.scatter_dir) if args.scatter_dir else None,
            dry_run=args.dry_run
        )
    
    print("\n✅ Organization complete!")
    
    if args.dry_run:
        print("\n💡 Tip: Remove --dry-run to actually move the files")


if __name__ == "__main__":
    main()