#!/usr/bin/env python3
"""
RefSeq ID → UniProt ID 変換スクリプト
変換後は main.py または ./scripts/run_all_chunks.sh で解析できます
"""

import requests
import time
import pandas as pd
from typing import List, Dict, Tuple

def refseq_to_uniprot(refseq_ids: List[str]) -> Dict[str, str]:
    """
    RefSeq IDをUniProt IDに変換
    
    Parameters
    ----------
    refseq_ids : list
        RefSeq IDのリスト (NM_, NR_, NP_ など)
    
    Returns
    -------
    dict
        {RefSeq ID: UniProt ID} のマッピング辞書
    """
    if not refseq_ids:
        return {}
    
    print(f"🔄 Converting {len(refseq_ids)} RefSeq IDs to UniProt...")
    
    # UniProt ID Mapping API を使用
    url = "https://rest.uniprot.org/idmapping/run"
    
    payload = {
        "from": "RefSeq_Protein",  # NP_, XP_ など
        "to": "UniProtKB",
        "ids": ",".join(refseq_ids)
    }
    
    try:
        # ジョブ開始
        response = requests.post(url, data=payload)
        response.raise_for_status()
        job_data = response.json()
        job_id = job_data.get("jobId")
        
        if not job_id:
            print("❌ Failed to start conversion job")
            return {}
        
        print(f"  Job ID: {job_id}")
        print(f"  Waiting for results...", end="")
        
        # 結果を待つ（最大60秒）
        result_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        
        for i in range(60):
            time.sleep(1)
            status_response = requests.get(result_url)
            status_data = status_response.json()
            
            if "results" in status_data or "failedIds" in status_data:
                print(" Done!")
                break
            print(".", end="", flush=True)
        else:
            print(" Timeout!")
            return {}
        
        # 結果取得
        results_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
        results_response = requests.get(results_url)
        results_data = results_response.json()
        
        # マッピング辞書作成
        mapping = {}
        for result in results_data.get("results", []):
            refseq_id = result.get("from")
            uniprot_id = result.get("to", {}).get("primaryAccession")
            if refseq_id and uniprot_id:
                mapping[refseq_id] = uniprot_id
        
        return mapping
        
    except Exception as e:
        print(f"\n❌ Error during conversion: {str(e)}")
        return {}


def try_alternative_mapping(refseq_ids: List[str]) -> Dict[str, str]:
    """
    代替方法: NM_/NR_ を含むRefSeqの場合、Gene情報経由で変換を試みる
    """
    print(f"\n🔄 Trying alternative mapping for NM_/NR_ IDs...")
    
    mapping = {}
    
    for refseq_id in refseq_ids:
        try:
            # UniProt検索APIを使用
            search_url = f"https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": f"xref:refseq-{refseq_id}",
                "format": "json",
                "size": 1
            }
            
            response = requests.get(search_url, params=params)
            data = response.json()
            
            if data.get("results"):
                uniprot_id = data["results"][0].get("primaryAccession")
                if uniprot_id:
                    mapping[refseq_id] = uniprot_id
                    print(f"  ✓ {refseq_id} → {uniprot_id}")
            
            time.sleep(0.3)  # API制限回避
            
        except Exception as e:
            print(f"  ✗ {refseq_id}: {str(e)}")
    
    return mapping


def convert_and_save(refseq_list: List[str], output_file: str = "./output/refseq_to_uniprot_mapping.csv"):
    """
    RefSeqリストを変換して結果を保存
    """
    print("=" * 80)
    print("RefSeq → UniProt Conversion")
    print("=" * 80)
    
    if not refseq_list:
        print("❌ No RefSeq IDs provided")
        return [], []
    
    # RefSeqのタイプを判別
    np_ids = [rid for rid in refseq_list if rid.startswith("NP_") or rid.startswith("XP_")]
    nm_nr_ids = [rid for rid in refseq_list if rid.startswith("NM_") or rid.startswith("NR_")]
    other_ids = [rid for rid in refseq_list if rid not in np_ids and rid not in nm_nr_ids]
    
    print(f"\n📋 Input RefSeq IDs:")
    print(f"  NP_/XP_ (Protein): {len(np_ids)}")
    print(f"  NM_/NR_ (RNA): {len(nm_nr_ids)}")
    print(f"  Other: {len(other_ids)}")
    
    # NP_/XP_ (タンパク質) を変換
    mapping = {}
    if np_ids:
        print(f"\n🧬 Converting protein RefSeq IDs...")
        mapping.update(refseq_to_uniprot(np_ids))
    
    # NM_/NR_ (RNA) を代替方法で変換
    if nm_nr_ids:
        print(f"\n🧬 Converting RNA RefSeq IDs (may take longer)...")
        mapping.update(try_alternative_mapping(nm_nr_ids))
    
    # 結果サマリー
    print(f"\n{'=' * 80}")
    print("Conversion Results")
    print(f"{'=' * 80}")
    print(f"✅ Successfully converted: {len(mapping)}/{len(refseq_list)}")
    print(f"❌ Failed to convert: {len(refseq_list) - len(mapping)}")
    
    # 変換されたUniProt IDリスト
    uniprot_ids = list(mapping.values())
    failed_ids = [rid for rid in refseq_list if rid not in mapping]
    
    # 結果を表示
    if mapping:
        print(f"\n✅ Converted mappings:")
        for refseq_id, uniprot_id in list(mapping.items())[:10]:
            print(f"  {refseq_id} → {uniprot_id}")
        if len(mapping) > 10:
            print(f"  ... and {len(mapping) - 10} more")
    
    if failed_ids:
        print(f"\n❌ Failed to convert:")
        for rid in failed_ids[:10]:
            print(f"  {rid}")
        if len(failed_ids) > 10:
            print(f"  ... and {len(failed_ids) - 10} more")
    
    # CSV保存
    if mapping:
        df = pd.DataFrame([
            {"RefSeq_ID": k, "UniProt_ID": v} 
            for k, v in mapping.items()
        ])
        df.to_csv(output_file, index=False)
        print(f"\n💾 Mapping saved to: {output_file}")
    
    # UniProt IDリストを保存
    uniprot_list_file = "./output/converted_uniprots.txt"
    if uniprot_ids:
        with open(uniprot_list_file, 'w') as f:
            f.write('\n'.join(uniprot_ids))
        print(f"💾 UniProt IDs saved to: {uniprot_list_file}")
    
    return uniprot_ids, failed_ids


def main():
    """メイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Convert RefSeq IDs to UniProt IDs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert RefSeq IDs from file
  python refseq_to_uniprot_converter.py --input refseq_list.txt
  
  # Convert specific IDs
  python refseq_to_uniprot_converter.py --ids NM_001178430 NR_132154
  
  # After conversion, analyze with:
  ./scripts/run_all_chunks.sh
        """
    )
    
    parser.add_argument('--input', '-i', type=str,
                       help='Input file with RefSeq IDs (one per line)')
    parser.add_argument('--ids', nargs='+',
                       help='RefSeq IDs to convert (space-separated)')
    parser.add_argument('--output', '-o', type=str,
                       default='./output/refseq_to_uniprot_mapping.csv',
                       help='Output CSV file (default: ./output/refseq_to_uniprot_mapping.csv)')
    
    args = parser.parse_args()
    
    # RefSeq IDリストを取得
    refseq_list = []
    
    if args.input:
        try:
            with open(args.input, 'r') as f:
                refseq_list = [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"❌ Error reading input file: {e}")
            return
    elif args.ids:
        refseq_list = args.ids
    else:
        # デフォルト: あなたが提供したRefSeq IDリスト
        refseq_list = [
            "NM_001178430",
            "NR_132154",
            "NR_132184",
            "NR_132206"
        ]
        print("💡 Using default RefSeq IDs (specify --input or --ids to use custom list)")
    
    # 変換実行
    uniprot_ids, failed_ids = convert_and_save(refseq_list, args.output)
    
    # 次のステップを提示
    if uniprot_ids:
        print(f"\n{'=' * 80}")
        print("Next Steps")
        print(f"{'=' * 80}")
        print("1. Analyze converted UniProt IDs:")
        print("   ./scripts/run_all_chunks.sh")
        print()
        print("2. Or directly analyze in main.py:")
        print("   uniprot_ids = open('./output/converted_uniprots.txt').read().strip().split('\\n')")


if __name__ == "__main__":
    main()