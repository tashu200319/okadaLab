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
    df = pd.read_csv(csv_file)
    
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


def main():
    """メイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Generate scatter plots with score_log or distance_std from score_details files'
    )
    parser.add_argument(
        '--input-dir',
        default='/Users/tashiroshuya/Desktop/okadaLab/output/score_details/with_and_search',
        help='Input directory containing score_details CSV files'
    )
    parser.add_argument(
        '--output-dir',
        default='/Users/tashiroshuya/Desktop/okadaLab/output/score_details/scatter/with_score_log',
        help='Output directory for scatter plots'
    )
    parser.add_argument(
        '--y-axis',
        choices=['score_log', 'distance_std'],
        default='score_log',
        help='Y-axis type: score_log or distance_std (default: score_log)'
    )
    parser.add_argument(
        '--point-size',
        type=int,
        default=5,
        help='Point size (default: 5)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.4,
        help='Point transparency (default: 0.4)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Image resolution (default: 300)'
    )
    parser.add_argument(
        '--max-files',
        type=int,
        default=None,
        help='Maximum number of files to process (for testing)'
    )
    
    args = parser.parse_args()
    
    plot_all_score_scatter(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        y_axis=args.y_axis,
        point_size=args.point_size,
        alpha=args.alpha,
        dpi=args.dpi,
        max_files=args.max_files
    )


if __name__ == "__main__":
    main()