"""
解析結果の可視化モジュール
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def generate_heatmap(score: pd.DataFrame, output_path: str = None):
    """
    スコアのヒートマップを生成
    
    Parameters
    ----------
    score : pd.DataFrame
        スコアデータ（残基ペア × スコア）
    output_path : str, optional
        保存先パス。Noneの場合は表示のみ
    """
    # 最終残基ペアから行列サイズを決定
    n0, n1 = score.iloc[-1, 0].split(', ')
    df1 = pd.DataFrame(np.zeros((int(n1), int(n1))))
    df1[:] = np.nan
    
    def Q(x, df):
        """スコアを行列に配置"""
        x00, x01 = x[0].split(', ')
        df.loc[int(x00) - 1, int(x01) - 1] = x[4]  # score列
    
    score.apply(Q, df=df1, axis=1)
    
    # ヒートマップ描画
    plt.figure(figsize=(12, 10))
    sns.heatmap(df1, cmap='viridis', cbar_kws={'label': 'Score'})
    plt.xlabel('Residue Number')
    plt.ylabel('Residue Number')
    plt.title('DSA Score Heatmap')
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved to {output_path}")
    else:
        plt.show()
    
    plt.close()
    
    return df1


def plot_distance_distribution(distance: pd.DataFrame, output_path: str = None):
    """
    距離分布のヒストグラムを作成
    
    Parameters
    ----------
    distance : pd.DataFrame
        距離データ
    output_path : str, optional
        保存先パス
    """
    # 距離データの抽出（最初の2列はID情報）
    dist_values = distance.iloc[:, 2:].values.flatten()
    dist_values = dist_values[~np.isnan(dist_values)]
    
    plt.figure(figsize=(10, 6))
    plt.hist(dist_values, bins=50, edgecolor='black', alpha=0.7)
    plt.xlabel('Distance (Å)')
    plt.ylabel('Frequency')
    plt.title('Distance Distribution')
    plt.grid(True, alpha=0.3)
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Distance distribution saved to {output_path}")
    else:
        plt.show()
    
    plt.close()


def plot_cis_analysis(cis_dist: pd.DataFrame, cis_threshold: float = 3.3,
                     output_path: str = None):
    """
    cis結合の解析結果を可視化
    
    Parameters
    ----------
    cis_dist : pd.DataFrame
        cis距離データ
    cis_threshold : float
        cis判定閾値
    output_path : str, optional
        保存先パス
    """
    if len(cis_dist) == 0:
        print("No cis pairs detected")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 距離の分布
    dist_cols = cis_dist.columns[2:]  # 最初の2列はID情報
    dist_data = []
    for col in dist_cols:
        if col not in ['distance mean', 'distance std', 'score', 
                       'cis_cnt', 'trans_cnt']:
            dist_data.extend(cis_dist[col].dropna().values)
    
    axes[0, 0].hist(dist_data, bins=30, edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(cis_threshold, color='r', linestyle='--', 
                       label=f'Threshold={cis_threshold}Å')
    axes[0, 0].set_xlabel('Distance (Å)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Cis Distance Distribution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 平均距離
    if 'distance mean' in cis_dist.columns:
        axes[0, 1].hist(cis_dist['distance mean'].dropna(), bins=30, 
                       edgecolor='black', alpha=0.7, color='green')
        axes[0, 1].set_xlabel('Mean Distance (Å)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Mean Distance per Residue Pair')
        axes[0, 1].grid(True, alpha=0.3)
    
    # 標準偏差
    if 'distance std' in cis_dist.columns:
        axes[1, 0].hist(cis_dist['distance std'].dropna(), bins=30, 
                       edgecolor='black', alpha=0.7, color='orange')
        axes[1, 0].set_xlabel('Std Distance (Å)')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Std Distance per Residue Pair')
        axes[1, 0].grid(True, alpha=0.3)
    
    # スコア分布
    if 'score' in cis_dist.columns:
        axes[1, 1].hist(cis_dist['score'].dropna(), bins=30, 
                       edgecolor='black', alpha=0.7, color='purple')
        axes[1, 1].set_xlabel('Score')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].set_title('Score Distribution')
        axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Cis analysis plot saved to {output_path}")
    else:
        plt.show()
    
    plt.close()


def plot_summary_comparison(summary_df: pd.DataFrame, output_path: str = None):
    """
    複数タンパク質のサマリー比較
    
    Parameters
    ----------
    summary_df : pd.DataFrame
        サマリーデータ（複数のUniProt IDを含む）
    output_path : str, optional
        保存先パス
    """
    if len(summary_df) < 2:
        print("Need at least 2 entries for comparison")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    
    metrics = [
        ('Entries', 'Number of PDB Entries'),
        ('Chains', 'Number of Chains'),
        ('Length', 'Sequence Length'),
        ('Resolution', 'Resolution (Å)'),
        ('UMF', 'UMF Score'),
        ('cis/Length(%)', 'Cis Ratio (%)')
    ]
    
    for idx, (metric, title) in enumerate(metrics):
        row = idx // 3
        col = idx % 3
        
        if metric in summary_df.columns:
            axes[row, col].bar(range(len(summary_df)), 
                              summary_df[metric].astype(float))
            axes[row, col].set_xlabel('Protein Index')
            axes[row, col].set_ylabel(metric)
            axes[row, col].set_title(title)
            axes[row, col].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Summary comparison saved to {output_path}")
    else:
        plt.show()
    
    plt.close()