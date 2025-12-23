#!/usr/bin/env python3
"""
UniProt IDチェッカー - 高速化版
X-rayとEMの両方を持つかチェック + 並列処理 + スキップ機能 + 最適化
"""

import sys
import os
import time
import threading

# 親ディレクトリをパスに追加
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd
from core.uniprot_handler import UniprotData
from core.config import Config
from typing import List, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed

# グローバル変数でリクエスト制御
request_lock = threading.Lock()
last_request_time = 0
MIN_REQUEST_INTERVAL = 0.05  # 50ms


def rate_limited_request():
    """レート制限を考慮したリクエスト制御"""
    global last_request_time
    with request_lock:
        current_time = time.time()
        time_since_last = current_time - last_request_time
        if time_since_last < MIN_REQUEST_INTERVAL:
            time.sleep(MIN_REQUEST_INTERVAL - time_since_last)
        last_request_time = time.time()


def load_processed_ids(output_file: str) -> Dict[str, Dict]:
    """既に処理済みのUniProt IDとその結果を読み込む"""
    if not os.path.exists(output_file):
        return {}
    
    try:
        df = pd.read_csv(output_file)
        processed = {}
        
        for _, row in df.iterrows():
            uid = row['uniprotid']
            processed[uid] = {
                'uniprotid': uid,
                'status': row.get('status', 'OK'),
                'has_both': row.get('has_both', False),
                'xray_count': row.get('xray_count', 0),
                'em_count': row.get('em_count', 0),
                'xray': row.get('xray_pdbs', '').split(', ') if pd.notna(row.get('xray_pdbs')) else [],
                'em': row.get('em_pdbs', '').split(', ') if pd.notna(row.get('em_pdbs')) else []
            }
            processed[uid]['xray'] = [x for x in processed[uid]['xray'] if x]
            processed[uid]['em'] = [x for x in processed[uid]['em'] if x]
        
        return processed
    except Exception as e:
        print(f"⚠️  Warning: Could not load existing results: {e}")
        return {}


def check_uniprot_optimized(uniprotid: str, retry_count: int = 2) -> Dict:
    """
    最適化されたUniProtチェック（リトライ機能付き）
    """
    for attempt in range(retry_count):
        try:
            rate_limited_request()  # レート制限
            
            unidata = UniprotData(uniprotid)
            pdbdata = unidata.getpdbdata({"X-ray", "EM"})
            
            if pdbdata is None or pdbdata.empty:
                return {
                    'uniprotid': uniprotid,
                    'status': 'OK',
                    'xray': [],
                    'em': [],
                    'has_both': False,
                }

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
                'status': 'OK',
                'xray': xray_pdbs,
                'em': em_pdbs,
                'has_both': bool(xray_pdbs and em_pdbs),
            }

        except Exception as e:
            if attempt < retry_count - 1:
                time.sleep(0.5)  # リトライ前に待機
                continue
            else:
                return {
                    'uniprotid': uniprotid,
                    'status': 'ERROR_FETCH',
                    'xray': [],
                    'em': [],
                    'has_both': False,
                }


def save_intermediate_results(new_results: List[Dict], output_file: str):
    """中間結果を保存（追記式・重複排除）"""
    if not new_results:
        return
        
    new_df = pd.DataFrame([{
        'uniprotid': r['uniprotid'],
        'status': r.get('status', 'OK'),
        'has_both': r['has_both'],
        'xray_count': r['xray_count'],
        'em_count': r['em_count'],
        'xray_pdbs': ', '.join(r['xray'][:10]),
        'em_pdbs': ', '.join(r['em'][:10]),
    } for r in new_results])

    if os.path.exists(output_file):
        old_df = pd.read_csv(output_file)
        combined = pd.concat([old_df, new_df], ignore_index=True)
        combined = combined.drop_duplicates(subset="uniprotid", keep="last")
        combined.to_csv(output_file, index=False)
    else:
        # outputディレクトリがなければ作成
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        new_df.to_csv(output_file, index=False)


def check_manual_list_optimized(
    uniprot_ids: List[str],
    save_report: bool = True,
    output_file: str = "./output/manual_xray_em_check.csv",
    max_workers: int = 20,
    skip_processed: bool = True,
    progress_interval: int = 100
) -> Dict:
    """
    最適化版のチェック処理
    
    Parameters
    ----------
    uniprot_ids : list
        チェックするUniProt IDのリスト
    save_report : bool
        結果をCSVに保存するか
    output_file : str
        出力ファイル名
    max_workers : int
        並列処理のスレッド数（デフォルト: 20）
    skip_processed : bool
        既に処理済みのIDをスキップするか
    progress_interval : int
        何件ごとに進捗を保存するか（デフォルト: 100件）
    """
    config = Config()
    use_and_search = (config.SEARCH_MODE == "AND")
    
    # 既存結果の読み込み
    processed_results = {}
    skipped_ids = []
    
    if skip_processed:
        processed_results = load_processed_ids(output_file)
        if processed_results:
            print(f"📋 Found {len(processed_results)} already processed IDs")
            skipped_ids = [uid for uid in uniprot_ids if uid in processed_results]
            if skipped_ids:
                print(f"⏭️  Skipping {len(skipped_ids)} already processed IDs")
    
    ids_to_process = [uid for uid in uniprot_ids if uid not in processed_results]
    
    print("=" * 80)
    print(f"🚀 OPTIMIZED: Using parallel processing with {max_workers} workers")
    print("=" * 80)
    print(f"Total IDs: {len(uniprot_ids)}")
    print(f"Already processed: {len(skipped_ids)}")
    print(f"New IDs to check: {len(ids_to_process)}")
    print(f"💾 Auto-save every {progress_interval} IDs")
    print()
    
    results = []
    both_ids = []
    xray_only_ids = []
    em_only_ids = []
    neither_ids = []
    error_ids = []
    
    # スキップ分を追加
    for uid in skipped_ids:
        info = processed_results[uid]
        results.append(info)
        
        xray_count = info.get('xray_count', len(info.get('xray', [])))
        em_count = info.get('em_count', len(info.get('em', [])))
        
        if info.get('has_both'):
            both_ids.append(uid)
        elif xray_count > 0:
            xray_only_ids.append(uid)
        elif em_count > 0:
            em_only_ids.append(uid)
        else:
            neither_ids.append(uid)
    
    # 並列処理
    if ids_to_process:
        print(f"\n{'─' * 80}")
        print("Processing new IDs:")
        print(f"{'─' * 80}\n")
        
        total = len(ids_to_process)
        completed = 0
        start_time = time.time()
        temp_results = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_id = {
                executor.submit(check_uniprot_optimized, uid): uid
                for uid in ids_to_process
            }
            
            for future in as_completed(future_to_id):
                uid = future_to_id[future]
                completed += 1

                try:
                    info = future.result()
                except Exception as e:
                    info = {
                        "uniprotid": uid,
                        "status": "ERROR_THREAD",
                        "xray": [],
                        "em": [],
                        "has_both": False,
                    }

                status = info.get("status", "OK")
                xray = info.get("xray") or []
                em = info.get("em") or []
                xray_count = len(xray)
                em_count = len(em)
                has_both = (xray_count > 0 and em_count > 0)

                normalized = {
                    "uniprotid": info.get("uniprotid", uid),
                    "status": status,
                    "xray": xray,
                    "em": em,
                    "has_both": has_both,
                    "xray_count": xray_count,
                    "em_count": em_count,
                }
                results.append(normalized)
                temp_results.append(normalized)

                # 進捗表示（100件ごと）
                if completed % 100 == 0 or completed == total:
                    elapsed = time.time() - start_time
                    rate = completed / elapsed if elapsed > 0 else 0
                    remaining = (total - completed) / rate if rate > 0 else 0
                    print(f"📊 Progress: {completed}/{total} ({completed/total*100:.1f}%) | "
                          f"Speed: {rate:.1f} IDs/sec | ETA: {remaining/60:.1f} min")

                # 分類
                if status != "OK":
                    error_ids.append(uid)
                elif has_both:
                    both_ids.append(uid)
                elif xray_count > 0:
                    xray_only_ids.append(uid)
                elif em_count > 0:
                    em_only_ids.append(uid)
                else:
                    neither_ids.append(uid)

                # 定期保存（データロス防止）
                if save_report and len(temp_results) >= progress_interval:
                    save_intermediate_results(temp_results, output_file)
                    temp_results = []

        # 最後の残りを保存
        if save_report and temp_results:
            save_intermediate_results(temp_results, output_file)

    # サマリー表示
    print(f"\n{'=' * 80}")
    print("Summary")
    print(f"{'=' * 80}")
    print(f"✅ Has BOTH X-ray AND EM:  {len(both_ids):4d} ({len(both_ids)/len(uniprot_ids)*100:5.1f}%)")
    print(f"📊 X-ray only:            {len(xray_only_ids):4d} ({len(xray_only_ids)/len(uniprot_ids)*100:5.1f}%)")
    print(f"🔬 EM only:               {len(em_only_ids):4d} ({len(em_only_ids)/len(uniprot_ids)*100:5.1f}%)")
    print(f"❌ Neither:               {len(neither_ids):4d} ({len(neither_ids)/len(uniprot_ids)*100:5.1f}%)")
    if error_ids:
        print(f"⚠️  Errors:                {len(error_ids):4d} ({len(error_ids)/len(uniprot_ids)*100:5.1f}%)")
    print(f"{'=' * 80}")
    
    if save_report:
        print(f"\n💾 Final results saved to {output_file}")
        
        # BOTH持ちのIDをテキストファイルにも保存
        if both_ids:
            both_file = output_file.replace('.csv', '_both_only.txt')
            with open(both_file, 'w') as f:
                f.write('\n'.join(sorted(both_ids)))
            print(f"💾 IDs with BOTH saved to: {both_file}")
    
    return {
        'both': both_ids,
        'xray_only': xray_only_ids,
        'em_only': em_only_ids,
        'neither': neither_ids,
        'error': error_ids,
        'results': results
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Optimized UniProt X-ray/EM checker')
    parser.add_argument('--file', required=True, help='File containing UniProt IDs')
    parser.add_argument('--output', default='./output/manual_xray_em_check.csv')
    parser.add_argument('--workers', type=int, default=20, help='Parallel workers (default: 20)')
    parser.add_argument('--save-interval', type=int, default=100, help='Auto-save every N IDs (default: 100)')
    parser.add_argument('--no-skip', action='store_true', help='Re-check all IDs')
    
    args = parser.parse_args()
    
    try:
        with open(args.file, 'r') as f:
            uniprot_ids = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"❌ Error: File not found: {args.file}")
        return
    
    if not uniprot_ids:
        print("❌ Error: No UniProt IDs found")
        return
    
    check_manual_list_optimized(
        uniprot_ids,
        output_file=args.output,
        max_workers=args.workers,
        skip_processed=not args.no_skip,
        progress_interval=args.save_interval
    )


if __name__ == "__main__":
    main()