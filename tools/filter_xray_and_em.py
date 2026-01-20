#!/usr/bin/env python3
"""
既存の解析結果から X-ray AND EM の両方を持つUniProt IDを抽出
"""

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
from core.uniprot_handler import UniprotData
from typing import List, Dict

# このスクリプト tools/filter_xray_and_em.py から見たプロジェクトルート
BASE_DIR = Path(__file__).resolve().parents[1]



def check_uniprot_has_both_methods(uniprotid: str) -> Dict[str, List[str]]:
    """
    UniProt IDがX-rayとEMの両方のPDB構造を持っているかチェック
    
    Returns
    -------
    dict
        {'xray': [pdb_ids], 'em': [pdb_ids], 'has_both': bool}
    """
    try:
        unidata = UniprotData(uniprotid)
        
        # 全PDBデータを取得
        pdbdata = unidata.getpdbdata({"X-ray", "EM"})
        
        # 手法ごとに分類
        xray_pdbs = []
        em_pdbs = []
        
        for pdb_id in pdbdata.columns:
            method = pdbdata.at['method', pdb_id]
            if method == 'X-ray':
                xray_pdbs.append(pdb_id)
            elif method == 'EM':
                em_pdbs.append(pdb_id)
        
        return {
            'uniprotid': uniprotid,
            'xray': xray_pdbs,
            'em': em_pdbs,
            'xray_count': len(xray_pdbs),
            'em_count': len(em_pdbs),
            'has_both': len(xray_pdbs) > 0 and len(em_pdbs) > 0
        }
    
    except Exception as e:
        print(f"Error checking {uniprotid}: {e}")
        return {
            'uniprotid': uniprotid,
            'xray': [],
            'em': [],
            'xray_count': 0,
            'em_count': 0,
            'has_both': False
        }


def filter_from_summary(summary_file: str = None,
                        output_dir: str = None):
    """
    summary.csv から X-ray AND EM の両方を持つエントリをフィルタリング
    """
    # デフォルトパスを BASE_DIR から組み立て
    if summary_file is None:
        summary_file = BASE_DIR / "output" / "summaries" / "summary.csv"
    else:
        summary_file = Path(summary_file)

    if output_dir is None:
        output_dir = BASE_DIR / "output2"
    else:
        output_dir = Path(output_dir)

    # output2 ディレクトリを作成
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"📁 Created directory: {output_dir}")

    if not summary_file.exists():
        print(f"❌ Error: {summary_file} not found")
        return

    print("=" * 80)
    print(f"Filtering {summary_file}")
    print("=" * 80)
    
    # summary.csv を読み込み
    df = pd.read_csv(summary_file)
    
    # ユニークなUniProt IDリストを取得
    unique_uniprots = df['uniprotid'].unique()
    print(f"\nTotal UniProt IDs in summary: {len(unique_uniprots)}")
    
    # 各UniProt IDをチェック
    results = []
    both_methods_ids = []
    
    for idx, uniprotid in enumerate(unique_uniprots, 1):
        print(f"({idx}/{len(unique_uniprots)}) Checking {uniprotid}...", end=" ")
        
        info = check_uniprot_has_both_methods(uniprotid)
        results.append(info)
        
        if info['has_both']:
            print(f"✅ X-ray: {info['xray_count']}, EM: {info['em_count']}")
            both_methods_ids.append(uniprotid)
        else:
            print(f"❌ X-ray: {info['xray_count']}, EM: {info['em_count']}")
    
    # フィルタリング
    filtered_df = df[df['uniprotid'].isin(both_methods_ids)]
    
    # output2 ディレクトリに保存
    output_file = output_dir / "summary_xray_and_em.csv"
    filtered_df.to_csv(output_file, index=False)

    
    # 統計情報
    print("\n" + "=" * 80)
    print("Results:")
    print("=" * 80)
    print(f"Total UniProts checked: {len(unique_uniprots)}")
    print(f"UniProts with BOTH X-ray AND EM: {len(both_methods_ids)}")
    print(f"Percentage: {len(both_methods_ids)/len(unique_uniprots)*100:.1f}%")
    print(f"\nFiltered summary saved to: {output_file}")
    print("=" * 80)
    
    # 詳細レポート保存
    report_file = os.path.join(output_dir, "xray_em_report.csv")
    results_df = pd.DataFrame(results)
    results_df.to_csv(report_file, index=False)
    print(f"Detailed report saved to: {report_file}")
    
    return both_methods_ids


def filter_from_links(links_file: str = "./output/links/uniprot_pdb_links.csv",
                     output_dir: str = "./output2/"):
    """
    uniprot_pdb_links.csv から X-ray AND EM の両方を持つエントリをフィルタリング
    """
    import os
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if not pd.io.common.file_exists(links_file):
        print(f"❌ Error: {links_file} not found")
        return
    
    print("\n" + "=" * 80)
    print("Filtering PDB links for UniProts with BOTH methods")
    print("=" * 80)
    
    # links.csv を読み込み
    df = pd.read_csv(links_file)
    
    # ユニークなUniProt IDを取得
    unique_uniprots = df['uniprotid'].unique()
    
    # X-ray AND EM を持つUniProt IDを特定
    both_methods_ids = []
    
    for uniprotid in unique_uniprots:
        info = check_uniprot_has_both_methods(uniprotid)
        if info['has_both']:
            both_methods_ids.append(uniprotid)
    
    # フィルタリング
    filtered_df = df[df['uniprotid'].isin(both_methods_ids)]
    
    # output2 ディレクトリに保存
    output_file = os.path.join(output_dir, "uniprot_pdb_links_xray_and_em.csv")
    filtered_df.to_csv(output_file, index=False)
    
    print(f"Filtered {len(filtered_df)} PDB links")
    print(f"Saved to: {output_file}")


def generate_comparison_report(summary_file: str = None,
                               output_dir: str = None):

    if summary_file is None:
        summary_file = BASE_DIR / "output2" / "summary_xray_and_em.csv"  # ←ここへ変更
    else:
        summary_file = Path(summary_file)

    if output_dir is None:
        output_dir = BASE_DIR / "output2"
    else:
        output_dir = Path(output_dir)

    df = pd.read_csv(summary_file)
    unique_uniprots = df['uniprotid'].unique()
    
    xray_only = []
    em_only = []
    both = []
    
    for uniprotid in unique_uniprots:
        info = check_uniprot_has_both_methods(uniprotid)
        
        if info['xray_count'] > 0 and info['em_count'] > 0:
            both.append(info)
        elif info['xray_count'] > 0:
            xray_only.append(info)
        elif info['em_count'] > 0:
            em_only.append(info)
    
    print(f"\nX-ray only: {len(xray_only)}")
    print(f"EM only: {len(em_only)}")
    print(f"Both X-ray AND EM: {len(both)}")
    
    # 詳細保存
    report = {
        'category': [],
        'uniprotid': [],
        'xray_count': [],
        'em_count': []
    }
    
    for info in xray_only:
        report['category'].append('X-ray only')
        report['uniprotid'].append(info['uniprotid'])
        report['xray_count'].append(info['xray_count'])
        report['em_count'].append(0)
    
    for info in em_only:
        report['category'].append('EM only')
        report['uniprotid'].append(info['uniprotid'])
        report['xray_count'].append(0)
        report['em_count'].append(info['em_count'])
    
    for info in both:
        report['category'].append('Both')
        report['uniprotid'].append(info['uniprotid'])
        report['xray_count'].append(info['xray_count'])
        report['em_count'].append(info['em_count'])
    
    report_df = pd.DataFrame(report)
    report_file = "./output/summaries/method_comparison.csv"
    report_df.to_csv(report_file, index=False)
    print(f"\nComparison report saved to: {report_file}")
    
    return report_df


if __name__ == "__main__":
    import sys
    
    print("\nX-ray AND EM Filter Tool")
    print("=" * 80)
    print("This tool extracts UniProt IDs that have BOTH X-ray and EM structures")
    print("=" * 80)
    
    # オプション選択
    print("\nOptions:")
    print("1. Filter summary.csv")
    print("2. Filter uniprot_pdb_links.csv")
    print("3. Generate comparison report")
    print("4. All of the above")
    
    choice = input("\nSelect option (1-4): ").strip()
    
    if choice == '1':
        both_ids = filter_from_summary()
    elif choice == '2':
        filter_from_links()
    elif choice == '3':
        generate_comparison_report()
    elif choice == '4':
        both_ids = filter_from_summary()
        filter_from_links()
        generate_comparison_report()
    else:
        print("Invalid choice. Running option 4 (all).")
        both_ids = filter_from_summary()
        filter_from_links()
        generate_comparison_report()
    
    print("\n✅ Done!")
