#!/usr/bin/env python3
"""
score_detailsファイルから Distance vs Score (log) の散布図を作成
Y軸を score_log / distance_std から選択可能
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, LogLocator, MaxNLocator
from pathlib import Path
from typing import List, Optional


def plot_score_scatter(csv_file: str, 
                       output_dir: str,
                       y_axis: str = 'score_log',
                       point_size: int = 5,
                       alpha: float = 0.4,
                       dpi: int = 300) -> str:
    """
    1つのscore_detailsファイルから散布図を作成
    
    Parameters
    ----------
    csv_file : str
        score_detailsファイルのパス
    output_dir : str
        出力ディレクトリ
    y_axis : str
        Y軸の種類 ('score_log' or 'distance_std')
    point_size : int
        点のサイズ
    alpha : float
        点の透明度
    dpi : int
        画像解像度
    
    Returns
    -------
    str
        保存した画像のパス
    """
    # データ読み込み
    # 一部の score_details_*.csv が実体として PNG 等のバイナリになっている場合があるため、
    # 先頭シグネチャと読み込みエラーを見て安全にスキップする
    try:
        with open(csv_file, 'rb') as f:
            sig = f.read(8)
        if sig.startswith(b'\x89PNG\r\n\x1a\n'):
            print(f"⚠️  Skipping {csv_file}: looks like PNG (wrong extension)")
            return None
    except Exception as e:
        print(f"⚠️  Skipping {csv_file}: failed to read file header ({e})")
        return None

    try:
        df = pd.read_csv(csv_file)
    except UnicodeDecodeError:
        print(f"⚠️  Skipping {csv_file}: UnicodeDecodeError (not a text CSV)")
        return None
    except pd.errors.ParserError as e:
        print(f"⚠️  Skipping {csv_file}: ParserError ({e})")
        return None
    
    # 必要な列をチェック
    required_cols = ['distance mean']
    if y_axis == 'score_log':
        required_cols.append('score')
    elif y_axis == 'distance_std':
        required_cols.append('distance std')
    
    if not all(col in df.columns for col in required_cols):
        print(f"⚠️  Skipping {csv_file}: missing required columns")
        return None
    
    # UniProt IDをファイル名から抽出
    basename = os.path.basename(csv_file)
    # 例: score_details_P01308_20.csv → P01308
    uniprotid = basename.replace('score_details_', '').split('_')[0]
    
    # プロット作成
    fig, ax = plt.subplots(figsize=(8, 8))
    
    x_data = df['distance mean']
    
    # Y軸のデータとラベルを設定
    if y_axis == 'score_log':
        # scoreが0以下の値を除外（logが計算できないため）
        valid_mask = df['score'] > 0
        x_data = x_data[valid_mask]
        y_data = df.loc[valid_mask, 'score']  # 対数スケール表示するので生の値を使用
        y_label = 'Score'
        title_suffix = 'Distance vs Score (log scale)'
    else:  # distance_std
        y_data = df['distance std']
        y_label = 'Stddev (Å)'
        title_suffix = 'Distance Mean vs Std'
    
    # 散布図プロット
    ax.scatter(
        x_data, 
        y_data,
        s=point_size,
        alpha=alpha,
        color='steelblue'
    )
    
    # score_logの場合はY軸を対数スケールに設定し、目盛りを細かく
    if y_axis == 'score_log':
        ax.set_yscale('log')
        # 対数軸の副目盛りを1~9の位置に表示（目盛り線のみ）
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1, 10), numticks=100))
        # 副目盛りのラベルは2, 4, 6, 8のみ表示
        def log_formatter(x, pos):
            if x < 10:  # 10未満は表示しない
                return ''
            # x = 20, 40, 60, 80 → 2, 4, 6, 8
            # x = 200, 400, 600, 800 → 2, 4, 6, 8
            # x = 2000, 4000, 6000, 8000 → 2, 4, 6, 8
            
            # xを文字列にして最初の桁を取得
            x_str = f'{int(x)}'
            first_digit = int(x_str[0])
            
            # 2, 4, 6, 8のみ表示
            if first_digit in [2, 4, 6, 8] and len(x_str) > 1:
                return f'{first_digit}'
            return ''
        ax.yaxis.set_minor_formatter(plt.FuncFormatter(log_formatter))
        ax.grid(True, which='major', alpha=0.3, linestyle='-', linewidth=0.5)
        ax.grid(True, which='minor', alpha=0.1, linestyle='-', linewidth=0.3)
    else:
        ax.grid(alpha=0.3)
    
    ax.set_xlabel('Distance (Å)', fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(f'{uniprotid} - {title_suffix}\n({len(x_data)} residue pairs)', 
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    # 保存先を決定
    output_file = os.path.join(
        output_dir,
        basename.replace('.csv', f'_{y_axis}_scatter.png')
    )
    
    # ★ 既存ファイルがあればスキップ
    if os.path.exists(output_file):
        print(f"↩️  Skip (already exists): {os.path.basename(output_file)}")
        plt.close()
        return None
    
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file



def plot_all_score_scatter(input_dir: str,
                           output_dir: str,
                           y_axis: str = 'score_log',
                           point_size: int = 5,
                           alpha: float = 0.4,
                           dpi: int = 300,
                           max_files: int = None) -> List[str]:
    """
    全てのscore_detailsファイルから散布図を一括生成
    
    Parameters
    ----------
    input_dir : str
        score_detailsファイルがあるディレクトリ
    output_dir : str
        出力ディレクトリ
    y_axis : str
        Y軸の種類 ('score_log' or 'distance_std')
    point_size : int
        点のサイズ
    alpha : float
        点の透明度
    dpi : int
        画像解像度
    max_files : int, optional
        処理する最大ファイル数（テスト用）
    
    Returns
    -------
    list
        生成された画像ファイルのパスリスト
    """
    # 出力ディレクトリ作成
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📁 Created output directory: {output_dir}")
    
    # score_detailsファイルを検索
    pattern = os.path.join(input_dir, "score_details_*_20.csv")
    csv_files = sorted(glob.glob(pattern))
    
    if not csv_files:
        print(f"❌ No score_details files found: {pattern}")
        return []
    
    # max_files制限
    if max_files:
        csv_files = csv_files[:max_files]
    
    print("=" * 80)
    print(f"Creating scatter plots for {len(csv_files)} files")
    print("=" * 80)
    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")
    print(f"Y-axis: {y_axis}")
    print(f"Settings: point_size={point_size}, alpha={alpha}, dpi={dpi}")
    print("=" * 80)
    
    output_files = []
    
    for i, csv_file in enumerate(csv_files, 1):
        basename = os.path.basename(csv_file)
        print(f"({i}/{len(csv_files)}) Processing {basename}...", end=" ")
        
        try:
            output_file = plot_score_scatter(
                csv_file,
                output_dir=output_dir,
                y_axis=y_axis,
                point_size=point_size,
                alpha=alpha,
                dpi=dpi
            )
            
            if output_file:
                output_files.append(output_file)
                print(f"✓ Saved to {os.path.basename(output_file)}")
            else:
                print("✗ Skipped")
                
        except Exception as e:
            print(f"✗ Error: {e}")
            continue
    
    # サマリー
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"✓ Successfully created: {len(output_files)} plots")
    print(f"📂 Saved to: {output_dir}")
    print("=" * 80)
    
    return output_files

def reorganize_existing_plots(scatter_dir: str = None):
    """既存の散布図を整理（新旧両方の命名規則に対応）"""
    import shutil
    import glob

    # デフォルトはリポジトリ相対（絶対パスを避けて移動耐性を上げる）
    if scatter_dir is None:
        base_dir = Path(__file__).resolve().parents[1]
        scatter_dir = str(base_dir / "output" / "score_details" / "scatter")
    
    # with_and_search内にサブディレクトリ作成
    with_and_search_dir = os.path.join(scatter_dir, 'with_and_search')
    score_log_dir = os.path.join(with_and_search_dir, 'with_score_log')
    distance_std_dir = os.path.join(with_and_search_dir, 'with_distance_std')
    
    os.makedirs(score_log_dir, exist_ok=True)
    os.makedirs(distance_std_dir, exist_ok=True)
    
    moved_score = 0
    moved_std = 0
    
    # 1. with_and_search/ 直下のファイルを移動
    # 新しい命名規則: *_score_log_scatter.png
    for file in glob.glob(os.path.join(with_and_search_dir, '*_score_log_scatter.png')):
        dest = os.path.join(score_log_dir, os.path.basename(file))
        if os.path.exists(dest):
            os.remove(dest)
        shutil.move(file, score_log_dir)
        moved_score += 1
    
    # 新しい命名規則: *_distance_std_scatter.png
    for file in glob.glob(os.path.join(with_and_search_dir, '*_distance_std_scatter.png')):
        dest = os.path.join(distance_std_dir, os.path.basename(file))
        if os.path.exists(dest):
            os.remove(dest)
        shutil.move(file, distance_std_dir)
        moved_std += 1
    
    # 古い命名規則: score_details_*_scatter.png（suffix無し）
    # これらは score_log として扱う
    for file in glob.glob(os.path.join(with_and_search_dir, 'score_details_*_scatter.png')):
        # 新しい命名規則のファイルは除外（すでに処理済み）
        basename = os.path.basename(file)
        if '_score_log_scatter.png' in basename or '_distance_std_scatter.png' in basename:
            continue
        
        # 古い形式のファイルを score_log に移動
        dest = os.path.join(score_log_dir, basename)
        if os.path.exists(dest):
            os.remove(dest)
        shutil.move(file, score_log_dir)
        moved_score += 1
        print(f"  📦 Moved old format: {basename}")
    
    # 2. with_score_log/ フォルダの中身を統合
    old_score_log_dir = os.path.join(scatter_dir, 'with_score_log')
    if os.path.exists(old_score_log_dir):
        # score_log ファイルを移動
        for file in glob.glob(os.path.join(old_score_log_dir, '*_score_log_scatter.png')):
            dest = os.path.join(score_log_dir, os.path.basename(file))
            if os.path.exists(dest):
                os.remove(dest)
            shutil.move(file, score_log_dir)
            moved_score += 1
        
        # 古い命名規則のファイルも移動
        for file in glob.glob(os.path.join(old_score_log_dir, 'score_details_*_scatter.png')):
            basename = os.path.basename(file)
            if '_score_log_scatter.png' in basename or '_distance_std_scatter.png' in basename:
                continue
            dest = os.path.join(score_log_dir, basename)
            if os.path.exists(dest):
                os.remove(dest)
            shutil.move(file, score_log_dir)
            moved_score += 1
        
        # distance_std ファイルを移動
        for file in glob.glob(os.path.join(old_score_log_dir, '*_distance_std_scatter.png')):
            dest = os.path.join(distance_std_dir, os.path.basename(file))
            if os.path.exists(dest):
                os.remove(dest)
            shutil.move(file, distance_std_dir)
            moved_std += 1
        
        # 空になったwith_score_log/フォルダを削除
        try:
            os.rmdir(old_score_log_dir)
            print(f"🗑️  Removed empty directory: {old_score_log_dir}")
        except OSError as e:
            remaining = os.listdir(old_score_log_dir)
            print(f"⚠️  Directory not empty ({len(remaining)} files remaining), keeping: {old_score_log_dir}")
    
    print(f"\n✅ Reorganized: {moved_score} score_log, {moved_std} distance_std plots")
    print(f"📂 New structure:")
    print(f"   {score_log_dir}")
    print(f"   {distance_std_dir}")

def main():
    """メイン処理"""
    import argparse
    import os

    # デフォルトはリポジトリ相対（絶対パスを避けて移動耐性を上げる）
    base_dir = Path(__file__).resolve().parents[1]
    default_output_dir = base_dir / "output"

    parser = argparse.ArgumentParser(
        description="Generate scatter plots with score_log or distance_std from score_details files"
    )
    parser.add_argument(
        "--input-dir",
        default=str(default_output_dir / "score_details" / "with_and_search"),
        help="Input directory containing score_details CSV files",
    )
    parser.add_argument(
        "--output-root",
        default=str(default_output_dir / "score_details" / "scatter"),
        help="Root output directory (default: .../scatter)",
    )
    parser.add_argument(
        "--y-axis",
        choices=["score_log", "distance_std"],
        default="score_log",
        help="Y-axis type: score_log or distance_std (default: score_log)",
    )
    parser.add_argument(
        "--both",
        action="store_true",
        help="Generate both score_log and distance_std plots",
    )
    parser.add_argument("--point-size", type=int, default=5)
    parser.add_argument("--alpha", type=float, default=0.4)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--max-files", type=int, default=None)

    args = parser.parse_args()

    # --bothが指定された場合は両方のプロットを生成
    if args.both:
        y_axes = ["score_log", "distance_std"]
    else:
        y_axes = [args.y_axis]

    all_output_files = []
    for y_axis in y_axes:
        # y_axis に応じて出力先を自動で決める
        if y_axis == "score_log":
            output_dir = os.path.join(args.output_root, "with_score_log")
        else:
            output_dir = os.path.join(args.output_root, "with_distance_std")

        output_files = plot_all_score_scatter(
            input_dir=args.input_dir,
            output_dir=output_dir,
            y_axis=y_axis,
            point_size=args.point_size,
            alpha=args.alpha,
            dpi=args.dpi,
            max_files=args.max_files,
        )
        all_output_files.extend(output_files)
    
    if args.both:
        print("\n" + "=" * 80)
        print("Total Summary")
        print("=" * 80)
        print(f"✓ Successfully created: {len(all_output_files)} plots (both types)")
        print("=" * 80)

if __name__ == "__main__":
    main()
