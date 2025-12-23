#!/usr/bin/env python3
"""
チェーン数不足の疑いがあるIDを詳細調査
実際にprep()を実行してチェーン数を確認
"""

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.config import Config
from core.uniprot_handler import UniprotData
from core.structure_analyzer import CifData

def investigate_chains(uniprotid: str, max_pdbs: int = 50, seq_ratio: float = 20):
    """実際にチェーン数を調査"""
    print(f"\n{'='*80}")
    print(f"Investigating {uniprotid}")
    print(f"{'='*80}")
    
    config = Config()
    
    try:
        unidata = UniprotData(uniprotid)
        fullname = unidata.get_fullname()
        organism = unidata.get_organism()
        
        print(f"Protein: {fullname}")
        print(f"Organism: {organism}")
        
        # PDBリスト取得
        if config.USE_AND_SEARCH and config.USE_XRAY and config.USE_EM:
            xray_pdblist = unidata.pdblist({"X-ray"})
            em_pdblist = unidata.pdblist({"EM"})
            pdblist = list(set(xray_pdblist) & set(em_pdblist))
        else:
            pdblist = unidata.pdblist(config.METHODS_SELECTED)
            
        if max_pdbs and len(pdblist) > max_pdbs:
            pdblist = pdblist[:max_pdbs]
        
        print(f"\nPDB count: {len(pdblist)}")
        print(f"PDBs: {pdblist[:10]}..." if len(pdblist) > 10 else f"PDBs: {pdblist}")
        
        # 各PDBのチェーン情報を調査
        valid_chains = []
        chain_details = []
        
        for pdbid in pdblist[:10]:  # 最初の10個を詳細調査
            try:
                cifdata = CifData(pdbid)
                mut_judge = cifdata.mutationjudge(unidata.get_id(), pdbid)
                
                # 正常またはsubstitutionのみカウント
                if mut_judge in ['normal', 'substitution']:
                    # 位置情報チェック
                    beg, end = unidata.position(pdbid)
                    if beg is not None and end is not None:
                        valid_chains.append(pdbid)
                        chain_details.append({
                            'pdbid': pdbid,
                            'mutation': mut_judge,
                            'position': f"{beg}-{end}",
                            'length': end - beg + 1
                        })
                        print(f"  ✓ {pdbid}: {mut_judge}, pos {beg}-{end}")
                    else:
                        print(f"  ✗ {pdbid}: No position data")
                else:
                    print(f"  ✗ {pdbid}: {mut_judge} (excluded)")
                    
            except Exception as e:
                print(f"  ✗ {pdbid}: Error - {str(e)}")
        
        print(f"\n{'='*40}")
        print(f"Valid chains found: {len(valid_chains)}")
        print(f"Chain threshold: {config.CHAIN_THRESHOLD}")
        
        if len(valid_chains) >= config.CHAIN_THRESHOLD:
            print(f"✅ SUFFICIENT - Should be analyzable!")
            verdict = "sufficient"
        else:
            print(f"❌ INSUFFICIENT - Need {config.CHAIN_THRESHOLD - len(valid_chains)} more")
            verdict = "insufficient"
        
        return {
            'uniprotid': uniprotid,
            'fullname': fullname,
            'pdb_count': len(pdblist),
            'valid_chains': len(valid_chains),
            'threshold': config.CHAIN_THRESHOLD,
            'verdict': verdict,
            'chain_details': chain_details
        }
        
    except Exception as e:
        print(f"❌ Investigation failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def investigate_multiple(uniprotids: list, max_pdbs: int = 50):
    """複数IDを一括調査"""
    results = []
    
    print("="*80)
    print(f"Chain Count Investigation")
    print(f"Investigating {len(uniprotids)} UniProt IDs")
    print("="*80)
    
    for uid in uniprotids:
        result = investigate_chains(uid, max_pdbs)
        if result:
            results.append(result)
    
    # サマリー表示
    print("\n" + "="*80)
    print("Investigation Summary")
    print("="*80)
    
    sufficient = [r for r in results if r['verdict'] == 'sufficient']
    insufficient = [r for r in results if r['verdict'] == 'insufficient']
    
    print(f"\n✅ Sufficient chains: {len(sufficient)} IDs")
    for r in sufficient:
        print(f"   {r['uniprotid']}: {r['valid_chains']} chains")
    
    print(f"\n❌ Insufficient chains: {len(insufficient)} IDs")
    for r in insufficient:
        print(f"   {r['uniprotid']}: {r['valid_chains']}/{r['threshold']} chains")
    
    # CSV保存
    summary_df = pd.DataFrame([{
        'uniprotid': r['uniprotid'],
        'fullname': r['fullname'],
        'pdb_count': r['pdb_count'],
        'valid_chains': r['valid_chains'],
        'threshold': r['threshold'],
        'verdict': r['verdict']
    } for r in results])
    
    output_file = "./output/chain_investigation.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\n📊 Results saved to: {output_file}")
    
    return results


def main():
    """メイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Investigate chain counts')
    parser.add_argument('--file', default='./output/score_details_diagnosis.csv',
                       help='Diagnosis CSV file')
    parser.add_argument('--max-pdbs', type=int, default=50,
                       help='Maximum PDBs to process')
    parser.add_argument('--uniprotid', help='Investigate single UniProt ID')
    
    args = parser.parse_args()
    
    if args.uniprotid:
        # 単一ID調査
        investigate_chains(args.uniprotid, args.max_pdbs)
    else:
        # CSVから possibly_insufficient_chains を抽出
        if not os.path.exists(args.file):
            print(f"❌ File not found: {args.file}")
            return
        
        df = pd.read_csv(args.file)
        target_ids = df[df['reason'] == 'possibly_insufficient_chains']['uniprotid'].tolist()
        
        if not target_ids:
            print("✅ No IDs with possibly_insufficient_chains found")
            return
        
        print(f"Found {len(target_ids)} IDs to investigate:")
        print(target_ids)
        
        investigate_multiple(target_ids, args.max_pdbs)


if __name__ == "__main__":
    import os
    main()