#!/usr/bin/env python3
"""
未解析UniProt IDの診断と自動解析スクリプト
"""

import os
import pandas as pd
from typing import List, Tuple, Dict
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.config import Config
from core.uniprot_handler import UniprotData

def get_remaining_uniprots(seq_ratio: float = 20) -> List[str]:
    """未解析のUniProt IDリストを取得"""
    unique_file = "./output/links/unique_uniprots.csv"
    summary_file = "./output/summaries/summary.csv"
    
    if not os.path.exists(unique_file):
        print(f"❌ Error: {unique_file} not found")
        return []
    
    unique_df = pd.read_csv(unique_file)
    total_uniprots = set(unique_df['uniprotid'].tolist())
    
    if not os.path.exists(summary_file):
        return list(total_uniprots)
    
    summary_df = pd.read_csv(summary_file)
    completed_df = summary_df[summary_df['seq_ratio'] == seq_ratio]
    completed_uniprots = set(completed_df['uniprotid'].tolist())
    
    remaining = total_uniprots - completed_uniprots
    return list(remaining)


def diagnose_uniprot(uniprotid: str, max_pdbs: int = 50) -> Dict:
    """UniProt IDの解析可能性を診断"""
    config = Config()
    result = {
        'uniprotid': uniprotid,
        'status': 'unknown',
        'reason': '',
        'pdb_count': 0,
        'can_analyze': False
    }
    
    try:
        # UniProtデータ取得
        unidata = UniprotData(uniprotid)
        
        # 基本情報取得
        try:
            fullname = unidata.get_fullname()
            organism = unidata.get_organism()
            result['fullname'] = fullname
            result['organism'] = organism
        except Exception as e:
            result['status'] = 'invalid_uniprot'
            result['reason'] = f"Cannot fetch UniProt data: {str(e)}"
            return result
        
        # PDBリスト取得
        try:
            if config.USE_AND_SEARCH and config.USE_XRAY and config.USE_EM:
                xray_pdblist = unidata.pdblist({"X-ray"})
                em_pdblist = unidata.pdblist({"EM"})
                pdblist = list(set(xray_pdblist) & set(em_pdblist))
            else:
                pdblist = unidata.pdblist(config.METHODS_SELECTED)
            
            # max_pdbs制限を適用
            if max_pdbs and len(pdblist) > max_pdbs:
                pdblist = pdblist[:max_pdbs]
            
            result['pdb_count'] = len(pdblist)
            
            if len(pdblist) < config.PDB_THRESHOLD:
                result['status'] = 'insufficient_pdbs'
                result['reason'] = f"Only {len(pdblist)} PDBs (need {config.PDB_THRESHOLD})"
                return result
            
        except Exception as e:
            result['status'] = 'pdb_fetch_error'
            result['reason'] = f"Cannot fetch PDB list: {str(e)}"
            return result
        
        # 配列データ取得テスト
        try:
            fasta = unidata.fasta()
            if not fasta or len(fasta) == 0:
                result['status'] = 'no_sequence'
                result['reason'] = "No sequence data available"
                return result
        except Exception as e:
            result['status'] = 'sequence_error'
            result['reason'] = f"Cannot fetch sequence: {str(e)}"
            return result
        
        # 解析可能と判定
        result['status'] = 'analyzable'
        result['reason'] = f"{len(pdblist)} PDBs available"
        result['can_analyze'] = True
        
    except Exception as e:
        result['status'] = 'error'
        result['reason'] = str(e)
    
    return result


def diagnose_all_remaining(seq_ratio: float = 20, max_pdbs: int = 50):
    """全ての未解析IDを診断"""
    print("=" * 80)
    print("Diagnosing Remaining UniProt IDs")
    print("=" * 80)
    
    remaining = get_remaining_uniprots(seq_ratio)
    
    if not remaining:
        print("\n🎉 All UniProt IDs have been analyzed!")
        return [], []
    
    print(f"\n📋 Found {len(remaining)} remaining UniProt IDs")
    print(f"Parameters: seq_ratio={seq_ratio}%, max_pdbs={max_pdbs}\n")
    
    analyzable = []
    unanalyzable = []
    
    for i, uniprotid in enumerate(remaining, 1):
        print(f"({i}/{len(remaining)}) Diagnosing {uniprotid}...", end=" ")
        
        diagnosis = diagnose_uniprot(uniprotid, max_pdbs)
        
        if diagnosis['can_analyze']:
            print(f"✅ {diagnosis['reason']}")
            analyzable.append(diagnosis)
        else:
            print(f"❌ {diagnosis['status']}: {diagnosis['reason']}")
            unanalyzable.append(diagnosis)
    
    # サマリー表示
    print(f"\n{'=' * 80}")
    print("Diagnosis Summary")
    print(f"{'=' * 80}")
    print(f"✅ Analyzable: {len(analyzable)}")
    print(f"❌ Cannot analyze: {len(unanalyzable)}")
    
    if unanalyzable:
        print(f"\n{'=' * 80}")
        print("Reasons for Failed Analysis:")
        print(f"{'=' * 80}")
        
        # 理由ごとに集計
        reasons = {}
        for item in unanalyzable:
            status = item['status']
            if status not in reasons:
                reasons[status] = []
            reasons[status].append(item['uniprotid'])
        
        for status, ids in reasons.items():
            print(f"\n{status} ({len(ids)} IDs):")
            for uid in ids[:5]:  # 最初の5つを表示
                detail = next(d for d in unanalyzable if d['uniprotid'] == uid)
                print(f"  - {uid}: {detail['reason']}")
            if len(ids) > 5:
                print(f"  ... and {len(ids) - 5} more")
    
    return analyzable, unanalyzable


def save_diagnosis_report(analyzable: List[Dict], unanalyzable: List[Dict], 
                          output_file: str = "./output/diagnosis_report.csv"):
    """診断結果をCSVに保存"""
    all_results = analyzable + unanalyzable
    
    if not all_results:
        return
    
    df = pd.DataFrame(all_results)
    df = df.sort_values('can_analyze', ascending=False)
    df.to_csv(output_file, index=False)
    
    print(f"\n📊 Diagnosis report saved to: {output_file}")


def create_analyzable_list(analyzable: List[Dict], 
                           output_file: str = "./output/analyzable_uniprots.txt"):
    """解析可能なUniProt IDをテキストファイルに保存"""
    if not analyzable:
        return
    
    uniprotids = [item['uniprotid'] for item in analyzable]
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(uniprotids))
    
    print(f"📝 Analyzable UniProt IDs saved to: {output_file}")
    print(f"\n💡 You can edit main.py and use:")
    print(f"   uniprot_ids = open('{output_file}').read().strip().split('\\n')")


def main():
    """メイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Diagnose remaining UniProt IDs')
    parser.add_argument('--seq-ratio', type=float, default=20,
                       help='Sequence ratio used in main.py (default: 20)')
    parser.add_argument('--max-pdbs', type=int, default=50,
                       help='Maximum number of PDBs to process (default: 50)')
    parser.add_argument('--analyze', action='store_true',
                       help='Automatically run analysis on analyzable IDs')
    
    args = parser.parse_args()
    
    # 診断実行
    analyzable, unanalyzable = diagnose_all_remaining(args.seq_ratio, args.max_pdbs)
    
    if not analyzable and not unanalyzable:
        return
    
    # レポート保存
    save_diagnosis_report(analyzable, unanalyzable)
    create_analyzable_list(analyzable)
    
    # 自動解析オプション
    if args.analyze and analyzable:
        print(f"\n{'=' * 80}")
        print("Starting Automatic Analysis")
        print(f"{'=' * 80}")
        
        response = input(f"\n⚠️  Analyze {len(analyzable)} UniProt IDs? (y/n): ")
        
        if response.lower() == 'y':
            # main.pyのrun_analysisを呼び出す代わりに、
            # ユーザーにmain.pyを編集して実行してもらう
            print("\n📝 To analyze these IDs:")
            print("1. Edit main.py")
            print("2. Set: uniprot_ids = open('./output/analyzable_uniprots.txt').read().strip().split('\\n')")
            print("3. Run: python main.py")
        else:
            print("\n⏸️  Analysis skipped")
    
    print(f"\n{'=' * 80}")


if __name__ == "__main__":
    main()