#!/usr/bin/env python3
"""
scoredetailsの解析状況をチェックし、未解析の原因を詳細診断するスクリプト
main.pyの処理ロジックに基づいて原因を特定
"""

import os
import glob
import pandas as pd
from typing import List, Dict, Tuple

# main.pyからインポート（必要に応じて調整）
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from core.config import Config
    from core.uniprot_handler import UniprotData
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False
    print("⚠️  Warning: config.py or uniprot_handler.py not found")
    print("   File existence check only mode")


def get_all_uniprots(unique_file: str = "./output/links/unique_uniprots.csv") -> List[str]:
    """全UniProt IDリストを取得"""
    if not os.path.exists(unique_file):
        print(f"❌ Error: {unique_file} not found")
        return []
    
    unique_df = pd.read_csv(unique_file)
    return unique_df['uniprotid'].tolist()


def diagnose_not_analyzed_reason(uniprotid: str, seq_ratio: float = 20, 
                                 max_pdbs: int = 50) -> Dict:
    """
    未解析の原因を詳細診断（main.pyの処理フローに基づく）
    """
    diagnosis = {
        'uniprotid': uniprotid,
        'reason': 'unknown',
        'details': '',
        'can_retry': False,
        'suggestion': ''
    }
    
    if not IMPORTS_AVAILABLE:
        diagnosis['reason'] = 'cannot_diagnose'
        diagnosis['details'] = 'Import modules not available'
        return diagnosis
    
    try:
        config = Config()
        
        # ステップ1: UniProtデータ取得
        try:
            unidata = UniprotData(uniprotid)
        except Exception as e:
            diagnosis['reason'] = 'uniprot_fetch_error'
            diagnosis['details'] = f"Failed to fetch UniProt data: {str(e)}"
            diagnosis['can_retry'] = 'Connection' in str(e) or 'aborted' in str(e)
            diagnosis['suggestion'] = 'Retry after some time' if diagnosis['can_retry'] else 'Check UniProt ID validity'
            return diagnosis
        
        # ステップ2: 基本情報取得
        try:
            fullname = unidata.get_fullname()
            organism = unidata.get_organism()
            diagnosis['fullname'] = fullname
            diagnosis['organism'] = organism
        except Exception as e:
            diagnosis['reason'] = 'metadata_error'
            diagnosis['details'] = f"Cannot get protein info: {str(e)}"
            diagnosis['can_retry'] = False
            return diagnosis
        
        # ステップ3: PDBリスト取得とカウント
        try:
            if config.USE_AND_SEARCH and config.USE_XRAY and config.USE_EM:
                xray_pdblist = unidata.pdblist({"X-ray"})
                em_pdblist = unidata.pdblist({"EM"})
                pdblist = list(set(xray_pdblist) & set(em_pdblist))
            else:
                pdblist = unidata.pdblist(config.METHODS_SELECTED)
                
            original_pdb_count = len(pdblist)
            
            # max_pdbs制限を適用
            if max_pdbs and len(pdblist) > max_pdbs:
                pdblist = pdblist[:max_pdbs]
            
            pdb_count = len(pdblist)
            diagnosis['pdb_count'] = pdb_count
            diagnosis['original_pdb_count'] = original_pdb_count
            
            # PDB_THRESHOLDチェック
            if pdb_count < config.PDB_THRESHOLD:
                diagnosis['reason'] = 'insufficient_pdbs'
                diagnosis['details'] = f"Only {pdb_count} PDBs (need {config.PDB_THRESHOLD})"
                diagnosis['can_retry'] = False
                diagnosis['suggestion'] = f'Need at least {config.PDB_THRESHOLD} PDB entries'
                return diagnosis
                
        except Exception as e:
            diagnosis['reason'] = 'pdb_list_error'
            diagnosis['details'] = f"Cannot fetch PDB list: {str(e)}"
            diagnosis['can_retry'] = True
            diagnosis['suggestion'] = 'Retry or check network connection'
            return diagnosis
        
        # ステップ4: prep()相当の処理（簡易版）
        try:
            # 配列データ取得
            fasta = unidata.fasta()
            if not fasta or len(fasta) == 0:
                diagnosis['reason'] = 'no_sequence'
                diagnosis['details'] = 'No sequence data available'
                diagnosis['can_retry'] = False
                return diagnosis
            
            diagnosis['sequence_length'] = len(fasta)
            
            # PDBごとの位置情報チェック（簡易）
            valid_pdbs = 0
            invalid_pdbs = []
            
            for pdbid in pdblist[:5]:  # 最初の5つだけチェック
                try:
                    beg, end = unidata.position(pdbid)
                    if beg is not None and end is not None:
                        valid_pdbs += 1
                    else:
                        invalid_pdbs.append(pdbid)
                except:
                    invalid_pdbs.append(pdbid)
            
            if valid_pdbs == 0:
                diagnosis['reason'] = 'no_valid_pdb_positions'
                diagnosis['details'] = f'No valid position data in checked PDBs: {invalid_pdbs}'
                diagnosis['can_retry'] = False
                diagnosis['suggestion'] = 'PDB entries may lack proper mapping'
                return diagnosis
            
        except Exception as e:
            diagnosis['reason'] = 'sequence_processing_error'
            diagnosis['details'] = f"Error in sequence processing: {str(e)}"
            diagnosis['can_retry'] = False
            return diagnosis
        
        # ステップ5: チェーン数推定（実際の処理なしで推定）
        # 実際に処理しないと正確にはわからないが、PDB数から推定
        estimated_chains = valid_pdbs * 2  # 平均的に各PDBから2チェーン取れると仮定
        
        if estimated_chains < config.CHAIN_THRESHOLD:
            diagnosis['reason'] = 'possibly_insufficient_chains'
            diagnosis['details'] = f'Estimated {estimated_chains} chains (need {config.CHAIN_THRESHOLD})'
            diagnosis['can_retry'] = False
            diagnosis['suggestion'] = 'May fail at chain count threshold'
            return diagnosis
        
        # ここまで来たら解析可能のはず
        diagnosis['reason'] = 'should_be_analyzable'
        diagnosis['details'] = f'{pdb_count} PDBs, {len(fasta)} residues'
        diagnosis['can_retry'] = True
        diagnosis['suggestion'] = 'Should be analyzable. May have failed during processing.'
        
    except Exception as e:
        diagnosis['reason'] = 'diagnosis_error'
        diagnosis['details'] = f"Error during diagnosis: {str(e)}"
        diagnosis['can_retry'] = False
    
    return diagnosis


def check_scoredetails_exists(uniprotid: str, seq_ratio: float = 20, 
                              diagnose_reason: bool = False, max_pdbs: int = 50) -> Dict:
    """特定のUniProt IDのscoredetailsファイルが存在するかチェック"""
    result = {
        'uniprotid': uniprotid,
        'seq_ratio': seq_ratio,
        'status': 'not_found',
        'file_path': None,
        'file_size': 0,
        'row_count': 0,
        'has_data': False,
        'diagnosis': None
    }
    
    # ファイルパスを構築
    file_path = f"./output/score_details/score_details_{uniprotid}_{int(seq_ratio)}.csv"
    result['file_path'] = file_path
    
    # ファイルの存在チェック
    if not os.path.exists(file_path):
        result['status'] = 'not_analyzed'
        
        # 詳細診断を実行
        if diagnose_reason and IMPORTS_AVAILABLE:
            result['diagnosis'] = diagnose_not_analyzed_reason(uniprotid, seq_ratio, max_pdbs)
        
        return result
    
    # ファイルサイズチェック
    file_size = os.path.getsize(file_path)
    result['file_size'] = file_size
    
    if file_size == 0:
        result['status'] = 'empty_file'
        return result
    
    # CSVの内容チェック
    try:
        df = pd.read_csv(file_path)
        result['row_count'] = len(df)
        
        if len(df) == 0:
            result['status'] = 'empty_data'
            return result
        
        # データが存在する
        result['status'] = 'completed'
        result['has_data'] = True
        
    except Exception as e:
        result['status'] = 'read_error'
        result['error'] = str(e)
        return result
    
    return result


def check_all_scoredetails(seq_ratio: float = 20, verbose: bool = True, 
                           diagnose: bool = False, max_pdbs: int = 50):
    """全UniProt IDのscoredetails解析状況をチェック"""
    print("=" * 80)
    print(f"Checking scoredetails Analysis Status (seq_ratio={seq_ratio}%)")
    print("=" * 80)
    
    all_uniprots = get_all_uniprots()
    
    if not all_uniprots:
        return None, None
    
    print(f"\n📋 Total UniProt IDs: {len(all_uniprots)}")
    if diagnose and IMPORTS_AVAILABLE:
        print("🔍 Detailed diagnosis mode: ON")
    elif diagnose and not IMPORTS_AVAILABLE:
        print("⚠️  Detailed diagnosis mode: Limited (missing imports)")
    
    results = []
    status_counts = {}
    reason_counts = {}
    
    for i, uniprotid in enumerate(all_uniprots, 1):
        if verbose:
            print(f"({i}/{len(all_uniprots)}) Checking {uniprotid}...", end=" ")
        
        result = check_scoredetails_exists(uniprotid, seq_ratio, diagnose, max_pdbs)
        results.append(result)
        
        status = result['status']
        status_counts[status] = status_counts.get(status, 0) + 1
        
        # 未解析の場合、原因をカウント
        if status == 'not_analyzed' and result.get('diagnosis'):
            reason = result['diagnosis']['reason']
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
        
        if verbose:
            if result['status'] == 'completed':
                print(f"✅ ({result['row_count']} rows)")
            elif result['status'] == 'not_analyzed':
                if diagnose and result.get('diagnosis'):
                    reason = result['diagnosis']['reason']
                    print(f"❌ {reason}")
                else:
                    print(f"❌ Not analyzed")
            elif result['status'] == 'empty_file':
                print(f"⚠️  Empty file")
            elif result['status'] == 'empty_data':
                print(f"⚠️  No data")
            else:
                print(f"❌ {result['status']}")
    
    # サマリー表示
    print(f"\n{'=' * 80}")
    print("Summary")
    print(f"{'=' * 80}")
    
    total = len(all_uniprots)
    for status, count in sorted(status_counts.items()):
        percentage = count / total * 100
        emoji = "✅" if status == "completed" else "❌"
        print(f"{emoji} {status:20s}: {count:4d} ({percentage:5.1f}%)")
    
    print(f"{'=' * 80}")
    print(f"Total:                   {total:4d}")
    
    # 詳細診断結果を表示
    if diagnose and reason_counts:
        print_reason_analysis(results, reason_counts)
    
    return results, status_counts


def print_reason_analysis(results: List[Dict], reason_counts: Dict):
    """未解析の原因別サマリーを表示"""
    print(f"\n{'=' * 80}")
    print("Not Analyzed - Reason Breakdown")
    print(f"{'=' * 80}")
    
    # 原因ごとにカウント
    for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
        print(f"\n❌ {reason}: {count} IDs")
        
        # 各原因の例を表示
        examples = [r for r in results 
                   if r['status'] == 'not_analyzed' 
                   and r.get('diagnosis', {}).get('reason') == reason][:3]
        
        for ex in examples:
            diag = ex['diagnosis']
            print(f"   • {ex['uniprotid']}: {diag.get('details', 'No details')}")
            if diag.get('suggestion'):
                print(f"     → {diag['suggestion']}")
    
    # リトライ可能なものをカウント
    can_retry = [r for r in results 
                 if r['status'] == 'not_analyzed' 
                 and r.get('diagnosis', {}).get('can_retry')]
    
    if can_retry:
        print(f"\n{'=' * 80}")
        print(f"🔄 Can Retry: {len(can_retry)} IDs")
        print("   These may succeed on retry (network errors, temporary issues)")
        for r in can_retry[:5]:
            print(f"   • {r['uniprotid']}")


def save_diagnosis_report(results: List[Dict], 
                         output_file: str = "./output/score_details_diagnosis.csv"):
    """診断結果を詳細CSVに保存"""
    if not results:
        return
    
    report_data = []
    for r in results:
        row = {
            'uniprotid': r['uniprotid'],
            'status': r['status'],
            'file_path': r['file_path'],
            'file_size': r['file_size'],
            'row_count': r['row_count'],
            'has_data': r['has_data']
        }
        
        if r.get('diagnosis'):
            diag = r['diagnosis']
            row['reason'] = diag.get('reason', '')
            row['details'] = diag.get('details', '')
            row['can_retry'] = diag.get('can_retry', False)
            row['suggestion'] = diag.get('suggestion', '')
            row['pdb_count'] = diag.get('pdb_count', '')
            row['fullname'] = diag.get('fullname', '')
            row['organism'] = diag.get('organism', '')
        
        report_data.append(row)
    
    df = pd.DataFrame(report_data)
    df = df.sort_values(['status', 'uniprotid'])
    df.to_csv(output_file, index=False)
    
    print(f"\n📊 Diagnosis report saved to: {output_file}")


def create_remaining_list(results: List[Dict], 
                         output_file: str = "./output/remaining_score_details.txt"):
    """未解析のUniProt IDリストをテキストファイルに保存"""
    remaining = [r['uniprotid'] for r in results if not r['has_data']]
    
    if not remaining:
        print("\n🎉 All UniProt IDs have been analyzed!")
        return
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(remaining))
    
    print(f"📝 Remaining UniProt IDs saved to: {output_file}")
    print(f"   Total remaining: {len(remaining)}")


def create_retryable_list(results: List[Dict],
                         output_file: str = "./output/retryable_score_details.txt"):
    """リトライ可能なIDリストを保存"""
    retryable = [r['uniprotid'] for r in results 
                 if r['status'] == 'not_analyzed' 
                 and r.get('diagnosis', {}).get('can_retry')]
    
    if not retryable:
        return
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(retryable))
    
    print(f"🔄 Retryable UniProt IDs saved to: {output_file}")
    print(f"   Total retryable: {len(retryable)}")


def main():
    """メイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Check scoredetails analysis status with detailed diagnosis'
    )
    parser.add_argument(
        '--seq-ratio', 
        type=float, 
        default=20,
        help='Sequence ratio to check (default: 20)'
    )
    parser.add_argument(
        '--max-pdbs',
        type=int,
        default=50,
        help='Maximum PDBs used in analysis (default: 50)'
    )
    parser.add_argument(
        '--quiet', 
        action='store_true',
        help='Show summary only (no per-ID output)'
    )
    parser.add_argument(
        '--diagnose', 
        action='store_true',
        help='Perform detailed diagnosis for not_analyzed IDs'
    )
    parser.add_argument(
        '--save-report', 
        action='store_true',
        dest='save_report',
        help='Save detailed report to CSV'
    )
    
    args = parser.parse_args()
    
    # 解析状況チェック
    results, status_counts = check_all_scoredetails(
        seq_ratio=args.seq_ratio,
        verbose=not args.quiet,
        diagnose=args.diagnose,
        max_pdbs=args.max_pdbs
    )
    
    if not results:
        return
    
    # 未解析リスト作成
    create_remaining_list(results)
    
    # リトライ可能リスト作成
    if args.diagnose:
        create_retryable_list(results)
    
    # レポート保存
    if args.save_report or args.diagnose:
        save_diagnosis_report(results)
    
    print(f"\n{'=' * 80}\n")


if __name__ == "__main__":
    main()