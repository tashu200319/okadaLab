#!/usr/bin/env python3
"""
DSA (Distance Structure Analysis) 高速化版
並列処理 + PDBキャッシュ + バッチ書き込み対応
PDBダウンロードの競合対策追加
失敗ID自動記録機能追加
"""

import os
import re
import pandas as pd
import csv
import shutil
import datetime
import pytz
from typing import List, Tuple, Optional, Dict, Set
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count, Manager
import traceback
import time

from core.config import Config
from core.uniprot_handler import UniprotData
from core.structure_analyzer import CifData
from core.sequence_processor import (
    convert_three, trim_sequence, trim2_sequence, 
    sort_sequence, getcoord
)
from core.distance_calculator import getdistance2, getscore
from core.report_generator import generate_log_content, export_to_csv
from core.visualization import generate_heatmap

from tools.failed_id_manager import FailedIDManager, classify_error_type


BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = BASE_DIR / "output"
ARCHIVE_DIR = OUTPUT_DIR / "archive"
TRIM_DIR = ARCHIVE_DIR / "trimsequence"
BACKUP_DIR = ARCHIVE_DIR / "summary_backup"

# ===== グローバルPDBロック管理 =====
# プロセス間で共有するロック辞書
pdb_locks = None

def cleanup_pdb_files_after_analysis(uniprotid: str, pdb_used: List[str]):
    """
    解析完了後に使用したPDBファイルを削除
    
    Parameters
    ----------
    uniprotid : str
        処理したUniProt ID
    pdb_used : list
        使用したPDB IDのリスト
    """
    import os
    from pathlib import Path
    
    deleted_count = 0
    
    for pdbid in pdb_used:
        # pdb_files/ ディレクトリ
        pdb_file = Path(f"pdb_files/{pdbid.lower()}.cif")
        if pdb_file.exists():
            try:
                pdb_file.unlink()
                deleted_count += 1
            except Exception as e:
                print(f"  ⚠️ Failed to delete {pdb_file}: {e}")
        
        # atom_coord/ ディレクトリ（こちらは残す場合はコメントアウト）
        # atom_csv = Path(f"atom_coord/{pdbid.upper()}.csv")
        # if atom_csv.exists():
        #     try:
        #         atom_csv.unlink()
        #     except Exception:
        #         pass
    
    if deleted_count > 0:
        print(f"  🧹 Cleaned {deleted_count} PDB files for {uniprotid}")

def init_worker(lock_dict):
    """ワーカープロセス初期化時にロック辞書を設定"""
    global pdb_locks
    pdb_locks = lock_dict

def acquire_pdb_lock(pdbid: str, timeout: float = 30.0) -> bool:
    """
    PDBファイルのダウンロードロックを取得
    
    Parameters
    ----------
    pdbid : str
        PDB ID
    timeout : float
        タイムアウト時間（秒）
    
    Returns
    -------
    bool
        ロック取得成功ならTrue
    """
    global pdb_locks
    if pdb_locks is None:
        return True  # ロック辞書がない場合はスキップ
    
    pdbid_upper = pdbid.upper()
    start_time = time.time()
    
    # ロックが解放されるまで待機
    while True:
        # アトミックな操作でロックを試行
        if pdb_locks.get(pdbid_upper, False) == False:
            pdb_locks[pdbid_upper] = True
            return True
        
        # タイムアウトチェック
        if time.time() - start_time > timeout:
            print(f"⚠️  Lock timeout for {pdbid}")
            return False
        
        time.sleep(0.1)  # 100ms待機

def release_pdb_lock(pdbid: str):
    """PDBファイルのダウンロードロックを解放"""
    global pdb_locks
    if pdb_locks is None:
        return
    
    pdbid_upper = pdbid.upper()
    pdb_locks[pdbid_upper] = False



# ===== PDBキャッシュ管理（ロック対応版） =====
def check_pdb_cache(pdbid: str) -> bool:
    """PDBファイルがキャッシュ済みかチェック"""
    pdb_path = Path("pdb_files") / f"{pdbid.lower()}.cif"
    atom_path = Path("atom_coord") / f"{pdbid.upper()}.csv"
    return pdb_path.exists() or atom_path.exists()

def safe_download_pdb(pdbid: str, verbose: bool = False) -> bool:
    """
    ロック機構付きPDBダウンロード
    
    Returns
    -------
    bool
        成功ならTrue
    """
    # すでにキャッシュにあればスキップ
    if check_pdb_cache(pdbid):
        if verbose:
            print(f"  💾 Cache hit: {pdbid}")
        return True
    
    # ロック取得
    if not acquire_pdb_lock(pdbid):
        print(f"  ⚠️  Failed to acquire lock: {pdbid}")
        return False
    
    try:
        # ダブルチェック（ロック取得中に他のプロセスがDLした可能性）
        if check_pdb_cache(pdbid):
            if verbose:
                print(f"  💾 Cache hit (after lock): {pdbid}")
            return True
        
        # 実際のダウンロード処理
        cifdata = CifData(pdbid)
        
        if verbose:
            print(f"  ⬇️  Downloaded: {pdbid}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Download failed: {pdbid} - {e}")
        return False
        
    finally:
        # 必ずロックを解放
        release_pdb_lock(pdbid)

def get_cache_stats(pdblist: List[str]) -> Dict[str, int]:
    """キャッシュヒット率を取得"""
    cached = sum(1 for pdb in pdblist if check_pdb_cache(pdb))
    return {
        'total': len(pdblist),
        'cached': cached,
        'to_download': len(pdblist) - cached,
        'hit_rate': round(cached / len(pdblist) * 100, 1) if pdblist else 0
    }

# ===== 既存関数（変更なし） =====
def setup_archive_dirs():
    """archive用ディレクトリを生成"""
    TRIM_DIR.mkdir(parents=True, exist_ok=True)
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)

def auto_archive_output():
    """output直下のCSV整理"""
    moved = 0
    for f in OUTPUT_DIR.glob("trimsequence_*.csv"):
        shutil.move(str(f), TRIM_DIR / f.name)
        moved += 1
    for f in OUTPUT_DIR.glob("summary_backup_*.csv"):
        shutil.move(str(f), BACKUP_DIR / f.name)
        moved += 1
    if moved > 0:
        print(f"\n✨ Auto-archive completed. {moved} files organized.\n")

def load_processed_uniprots(summary_file: str, seq_ratio: float) -> set:
    """処理済みUniProt IDを取得"""
    if not os.path.exists(summary_file):
        return set()
    try:
        df = pd.read_csv(summary_file)
        processed = df[df['seq_ratio'] == seq_ratio]['uniprotid'].tolist()
        print(f"📋 Found {len(processed)} already processed UniProt IDs")
        return set(processed)
    except Exception as e:
        print(f"⚠️  Warning: {e}")
        return set()

def clean_pdb_files_selective(keep_pdblist: List[str], verbose: bool = True):
    """指定されたPDBリストに含まれないファイルを削除"""
    keep_pdbs = set([pdb.upper() for pdb in keep_pdblist])
    deleted_count = 0
    
    if os.path.exists('pdb_files'):
        for filename in os.listdir('pdb_files'):
            if filename.startswith('.'):
                continue
            pdb_id = filename.split('.')[0].upper()
            if pdb_id not in keep_pdbs:
                filepath = os.path.join('pdb_files', filename)
                os.remove(filepath)
                deleted_count += 1
    
    if os.path.exists('atom_coord'):
        for filename in os.listdir('atom_coord'):
            if filename.startswith('.'):
                continue
            pdb_id = filename.split('.')[0].upper()
            if pdb_id not in keep_pdbs:
                filepath = os.path.join('atom_coord', filename)
                os.remove(filepath)
    
    if verbose and deleted_count > 0:
        print(f"🧹 Cleaned up {deleted_count} old PDB files")

def count_pdb(uniprotid: str, methods: Optional[set] = None, 
              negative_pdbid: str = "", max_pdbs: int = None) -> bool:
    """PDB数をカウント"""
    cfg = Config()
    if methods is None:
        methods = cfg.METHODS_SELECTED
    
    unidata = UniprotData(uniprotid)
    methods_to_search = []
    if cfg.USE_XRAY:
        methods_to_search.append("X-ray")
    if cfg.USE_NMR:
        methods_to_search.append("NMR")
    if cfg.USE_EM:
        methods_to_search.append("EM")

    if cfg.SEARCH_MODE == "AND" and len(methods_to_search) > 1:
        pdb_sets = []
        for method in methods_to_search:
            pdb_sets.append(set(unidata.pdblist({method})))
        pdblist = list(set.intersection(*pdb_sets))
    else:
        pdblist = unidata.pdblist(cfg.METHODS_SELECTED)
    
    if negative_pdbid:
        negative_list = re.split(r'[,\s]+', negative_pdbid.strip())
        negative_list_upper = [neg.upper() for neg in negative_list]
        pdblist = [item for item in pdblist 
                   if item.upper() not in negative_list_upper]
    
    if max_pdbs is not None and len(pdblist) > max_pdbs:
        pdblist = pdblist[:max_pdbs]
    
    return len(pdblist) >= cfg.PDB_THRESHOLD

def prep(uniprotid: str, methods: Optional[set] = None, 
         negative_pdbid: str = "", max_pdbs: int = None,
         verbose: bool = True) -> Tuple:
    """データ準備"""
    cfg = Config()
    if methods is None:
        methods = cfg.METHODS_SELECTED
    
    unidata = UniprotData(uniprotid)
    uniprotids = unidata.get_id()
    id_str = str(uniprotids)
    fasta = unidata.fasta()
    sequence = convert_three(fasta)
    seqdata = pd.DataFrame(sequence, columns=[id_str])
    len_seqdata = len(seqdata)
    
    methods_to_search = []
    if cfg.USE_XRAY:
        methods_to_search.append("X-ray")
    if cfg.USE_NMR:
        methods_to_search.append("NMR")
    if cfg.USE_EM:
        methods_to_search.append("EM")

    if cfg.SEARCH_MODE == "AND" and len(methods_to_search) > 1:
        pdb_sets = []
        for method in methods_to_search:
            pdb_sets.append(set(unidata.pdblist({method})))
        pdblist = list(set.intersection(*pdb_sets))
    else:
        pdblist = unidata.pdblist(cfg.METHODS_SELECTED)
    
    if negative_pdbid:
        negative_list = re.split(r'[,\s]+', negative_pdbid.strip())
        negative_list_upper = [neg.upper() for neg in negative_list]
        pdblist = [item for item in pdblist 
                   if item.upper() not in negative_list_upper]
    
    if max_pdbs is not None and len(pdblist) > max_pdbs:
        if verbose:
            print(f"  Limiting to first {max_pdbs} PDB entries (out of {len(pdblist)} available)")
        pdblist = pdblist[:max_pdbs]
    
    # キャッシュ統計表示
    if verbose:
        cache_stats = get_cache_stats(pdblist)
        print(f"  📦 Cache: {cache_stats['cached']}/{cache_stats['total']} "
              f"({cache_stats['hit_rate']}% hit rate)")
        print(f"  📥 Will download: {cache_stats['to_download']} PDB entries")
        if pdblist:
            pdb_list_str = ", ".join(pdblist[:10])  # 最初の10個を表示
            if len(pdblist) > 10:
                pdb_list_str += f", ... (+{len(pdblist)-10} more)"
            print(f"  📋 PDB list: {pdb_list_str}")
        print(f"  Processing {len(pdblist)} PDB entries ...")
    
    nor_pdblist = []
    sub_pdblist = []
    chi_pdblist = []
    din_pdblist = []
    skipped_count = 0
    
    for n, pdbid in enumerate(pdblist):
        try:
            cifdata = CifData(pdbid)
            mut_judge = cifdata.mutationjudge(uniprotids, pdbid)
            
            if verbose:
                cache_marker = "💾" if check_pdb_cache(pdbid) else "⬇️"
                print(f" {cache_marker} ({n+1}/{len(pdblist)}) {pdbid} {mut_judge}")
            
            if mut_judge == 'normal':
                nor_pdblist.append(pdbid)
            elif mut_judge == 'substitution':
                sub_pdblist.append(pdbid)
            elif mut_judge == 'chimera':
                chi_pdblist.append(pdbid)
            elif mut_judge == 'delins':
                din_pdblist.append(pdbid)
            else:
                if verbose:
                    print(f"  Skipping {pdbid}: {mut_judge}")
                skipped_count += 1
                continue
            
            beg, end = unidata.position(pdbid)
            
            if beg is None or end is None:
                if verbose:
                    print(f"  Skipping {pdbid}: no valid position data")
                skipped_count += 1
                if pdbid in nor_pdblist:
                    nor_pdblist.remove(pdbid)
                elif pdbid in sub_pdblist:
                    sub_pdblist.remove(pdbid)
                elif pdbid in chi_pdblist:
                    chi_pdblist.remove(pdbid)
                elif pdbid in din_pdblist:
                    din_pdblist.remove(pdbid)
                continue
            
            df_beg = pd.DataFrame(index=list(range(beg-1)))
            df_end = pd.DataFrame(index=list(range(len_seqdata - end)))
            seq = cifdata.getsequence(uniprotids)
            seq = pd.concat([df_beg, seq, df_end])
            seq.reset_index(inplace=True, drop=True)
            seqdata = pd.concat([seqdata, seq], axis=1)
            
        except Exception as e:
            if verbose:
                print(f"  Error processing {pdbid}: {str(e)}")
            skipped_count += 1
            continue
    
    all_pdblist = [nor_pdblist, sub_pdblist, chi_pdblist, din_pdblist]
    
    if verbose:
        total = sum(len(lst) for lst in all_pdblist)
        print(f" ✓ Data Preparation: {total}/{len(pdblist)} processed")
        if skipped_count > 0:
            print(f" Skipped: {skipped_count} entries")
    
    return seqdata, all_pdblist

def run_DSA(uniprotid: str, seqdata, export: bool, seqtype: str,
            methods: Optional[set] = None, seq_ratio: float = 80,
            cis_threshold: float = 3.3, dirpath: str = "./output/",
            verbose: bool = True):
    """DSA解析を実行"""
    cfg = Config()
    if methods is None:
        methods = cfg.METHODS_SELECTED
    
    unidata = UniprotData(uniprotid)
    uniprotids = unidata.get_id()
    str_ids = str(uniprotids)
    fasta = unidata.fasta()
    sequence = convert_three(fasta)
    
    unidata.getpdbdata(methods)
    
    if seqdata is None or len(seqdata) == 0:
        print(f"Error: seqdata is empty for {uniprotid}")
        return None, "", None, None
    
    trimsequence = sort_sequence(str_ids, seqdata, seq_ratio)
    if trimsequence is None or len(trimsequence) == 0:
        print(f"Error: trimsequence is empty for {uniprotid}")
        return None, "", None, None
    
    trimsequence.to_csv(os.path.join(dirpath, f"trimsequence_{uniprotid}.csv"), 
                        index=False)
    trimseqcol = trimsequence.columns.values[1:]
    
    if len(trimseqcol) <= cfg.CHAIN_THRESHOLD - 1:
        print(f"Error: Not enough chains for {uniprotid}")
        return None, "", None, None
    
    atomcoord = getcoord(trimsequence, uniprotid)
    if atomcoord is None or len(atomcoord) == 0:
        print(f"Error: atomcoord is empty for {uniprotid}")
        return None, "", None, None
    
    distance = getdistance2(atomcoord)
    if distance is None or len(distance) == 0:
        print(f"Error: distance is empty for {uniprotid}")
        return None, "", None, None
    
    score = getscore(distance, 0)
    
    from itertools import combinations
    residue_pairs = list(combinations(atomcoord.index, 2))
    residue_num1_list = [pair[0] + 1 for pair in residue_pairs]
    residue_num2_list = [pair[1] + 1 for pair in residue_pairs]
    residue_num_df = pd.DataFrame({
        'residue_num1': residue_num1_list,
        'residue_num2': residue_num2_list
    })
    distance_cols = distance.columns[2:]
    distance_data_df = distance[distance_cols].copy()
    merged_df = pd.concat([residue_num_df, distance_data_df], axis=1)
    # merged_df.to_csv(os.path.join(dirpath, f"distance_{uniprotid}.csv"), 
    #                  index=False, header=False)
    
    cis_index = []
    for col in distance.columns.values.tolist()[2:]:
        try:
            tmp = distance.query(f'`{col}`<=@cis_threshold').index.to_list()
        except (SyntaxError, pd.errors.ParserError):
            tmp = distance[distance[col] <= cis_threshold].index.tolist()
        cis_index.extend(tmp)
    
    if not cis_index:
        cis_info = [[0, 0, 0, 0, 0]]
    else:
        cis_index = sorted(set(cis_index))
        cis_dist = distance.iloc[cis_index, :]
        cis_dist_mean = cis_dist.iloc[:, 2:].mean(axis=None)
        cis_dist_std = cis_dist.iloc[:, 2:].stack().std(ddof=0)
        cis_score_mean = (cis_dist.iloc[:, 2:].mean(axis=1) / 
                          (cis_dist.iloc[:, 2:].std(axis=1).replace(0, 0.0001))).mean()
        cis_num = len(cis_dist)
        mix = 0
        cis_info = [[cis_dist_mean, cis_dist_std, cis_score_mean, cis_num, mix]]
    
    summary_df = generate_log_content(unidata.pdbdata, len(sequence), 
                                      trimsequence, score, cis_info)
    
    header_line = " ".join(summary_df.columns.astype(str).tolist())
    value_line = " ".join(str(v) for v in summary_df.iloc[0].tolist())
    log_text = header_line + "\n" + value_line
    
    return score, log_text, summary_df, distance

def save_score_details(uniprotid: str, score: pd.DataFrame, distance: pd.DataFrame,
                       seq_ratio: float, dirpath: str, 
                       existing_details: pd.DataFrame = None):
    """score detailsファイルを保存"""
    first_col = distance.iloc[:, 0]
    residue_nums = first_col.str.split(', ', expand=True)
    residue_num1 = residue_nums[0].astype(int)
    residue_num2 = residue_nums[1].astype(int)
    
    distance_values = distance.iloc[:, 2:]
    distance_min = distance_values.min(axis=1)
    distance_max = distance_values.max(axis=1)
    
    new_details = pd.DataFrame({
        'uniprotid': uniprotid,
        'residue_num1': residue_num1.values,
        'residue_num2': residue_num2.values,
        'distance mean': score['distance mean'].values,
        'distance std': score['distance std'].values,
        'distance min': distance_min.values,
        'distance max': distance_max.values,
        'score': score['score'].values
    })
    
    score_details_dir = os.path.join(dirpath, "score_details")
    details_file = os.path.join(
        score_details_dir, 
        f"score_details_{uniprotid}_{int(seq_ratio)}.csv"
    )
    new_details.to_csv(details_file, index=False)
    
    return None

def save_summary_statistics(uniprotid: str, fullName: str, organism: str,
                           score: pd.DataFrame, summary_df: pd.DataFrame,
                           seq_ratio: float, dirpath: str,
                           existing_stats: pd.DataFrame = None):
    """統計サマリーを保存"""
    mean_distance = score['distance mean'].mean()
    std_distance = score['distance std'].mean()
    mean_score = score['score'].mean()
    umf = float(summary_df.iloc[0]['UMF'])
    
    new_stats = pd.DataFrame({
        'uniprotid': [uniprotid],
        'seq_ratio': [seq_ratio],
        'fullName': [fullName],
        'organism': [organism],
        'Entries': [int(summary_df.iloc[0]['Entries'])],
        'Chains': [int(summary_df.iloc[0]['Chains'])],
        'Length': [int(summary_df.iloc[0]['Length'])],
        'mean_distance': [round(mean_distance, 3)],
        'mean_std': [round(std_distance, 3)],
        'UMF': [round(umf, 3)],
        'mean_score': [round(mean_score, 3)]
    })
    
    stats_file = os.path.join(dirpath, "summaries", "summary_statistics.csv")
    
    if not os.path.exists(stats_file):
        new_stats.to_csv(stats_file, index=False, mode='w')
    else:
        new_stats.to_csv(stats_file, index=False, mode='a', header=False)
    
    if existing_stats is not None and len(existing_stats) > 0:
        combined_stats = pd.concat([existing_stats, new_stats], ignore_index=True)
    else:
        combined_stats = new_stats
    
    return combined_stats

def save_uniprot_pdb_links(uniprotid: str, pdblist: List[str], 
                           seq_ratio: float, dirpath: str,
                           existing_links: pd.DataFrame = None):
    """UniProt-PDBリンクを保存"""
    new_links = pd.DataFrame({
        'uniprotid': [uniprotid] * len(pdblist),
        'pdbid': pdblist,
        'seq_ratio': [seq_ratio] * len(pdblist)
    })
    
    links_file = os.path.join(dirpath, "links", "uniprot_pdb_links.csv")
    
    if not os.path.exists(links_file):
        new_links.to_csv(links_file, index=False, mode='w')
    else:
        new_links.to_csv(links_file, index=False, mode='a', header=False)
    
    if existing_links is not None and len(existing_links) > 0:
        combined_links = pd.concat([existing_links, new_links], ignore_index=True)
    else:
        combined_links = new_links
    
    return combined_links

def setup_archive_dirs(output_dir: str) -> Tuple[str, str]:
    """archiveディレクトリを作成"""
    archive_dir = os.path.join(output_dir, "archive")
    trim_dir = os.path.join(archive_dir, "trimsequence")
    backup_dir = os.path.join(archive_dir, "summary_backup")

    for d in (archive_dir, trim_dir, backup_dir):
        if not os.path.exists(d):
            os.makedirs(d)

    return trim_dir, backup_dir

def archive_trim_and_summary(
    base_output_dir: str,
    trim_dir: str,
    backup_dir: str,
    verbose: bool = True,
) -> None:
    """trimsequence/summary_backupをarchiveへ移動"""
    moved = 0

    for name in os.listdir(base_output_dir):
        path = os.path.join(base_output_dir, name)
        if not os.path.isfile(path):
            continue

        if name.startswith("trimsequence_") and name.endswith(".csv"):
            dst = os.path.join(trim_dir, name)
            shutil.move(path, dst)
            moved += 1

        elif name.startswith("summary_backup_") and name.endswith(".csv"):
            dst = os.path.join(backup_dir, name)
            shutil.move(path, dst)
            moved += 1

    if verbose and moved > 0:
        print(f"🗂 Archived {moved} file(s).\n")

# ===== バッチ書き込み用クラス =====
class BatchWriter:
    """CSVバッチ書き込みマネージャー"""
    def __init__(self, batch_size: int = 50):
        self.batch_size = batch_size
        self.summary_buffer = []
        self.stats_buffer = []
        self.links_buffer = []
        
    def add_summary(self, row_dict: dict):
        self.summary_buffer.append(row_dict)
        
    def add_stats(self, stats_df: pd.DataFrame):
        self.stats_buffer.append(stats_df)
        
    def add_links(self, links_df: pd.DataFrame):
        self.links_buffer.append(links_df)
    
    def should_flush(self) -> bool:
        return len(self.summary_buffer) >= self.batch_size
    
    def flush(self, dirpath: str, fieldnames: List[str]):
        """バッファを一括書き込み"""
        if not self.summary_buffer:
            return
        
        # summary.csv書き込み
        filename = os.path.join(dirpath, "summaries", "summary.csv")
        with open(filename, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            for row in self.summary_buffer:
                writer.writerow(row)
        
        # summary_statistics.csv書き込み
        if self.stats_buffer:
            stats_file = os.path.join(dirpath, "summaries", "summary_statistics.csv")
            combined_stats = pd.concat(self.stats_buffer, ignore_index=True)
            if not os.path.exists(stats_file):
                combined_stats.to_csv(stats_file, index=False, mode='w')
            else:
                combined_stats.to_csv(stats_file, index=False, mode='a', header=False)
        
        # uniprot_pdb_links.csv書き込み
        if self.links_buffer:
            links_file = os.path.join(dirpath, "links", "uniprot_pdb_links.csv")
            combined_links = pd.concat(self.links_buffer, ignore_index=True)
            if not os.path.exists(links_file):
                combined_links.to_csv(links_file, index=False, mode='w')
            else:
                combined_links.to_csv(links_file, index=False, mode='a', header=False)
        
        count = len(self.summary_buffer)
        self.summary_buffer.clear()
        self.stats_buffer.clear()
        self.links_buffer.clear()
        
        print(f"💾 Batch written: {count} entries")
def cleanup_batch_pdb_files(keep_latest: int = 0):
    """
    pdb_filesディレクトリ全体をクリーンアップ
    
    Parameters
    ----------
    keep_latest : int
        最新N件のファイルを残す（0=全削除）
    """
    import os
    from pathlib import Path
    
    pdb_dir = Path("pdb_files")
    if not pdb_dir.exists():
        return
    
    # 全PDBファイルを取得
    pdb_files = list(pdb_dir.glob("*.cif"))
    
    if keep_latest > 0:
        # 更新日時でソートして古いものから削除
        pdb_files.sort(key=lambda x: x.stat().st_mtime, reverse=True)
        to_delete = pdb_files[keep_latest:]
    else:
        to_delete = pdb_files
    
    deleted_count = 0
    for pdb_file in to_delete:
        try:
            pdb_file.unlink()
            deleted_count += 1
        except Exception as e:
            print(f"  ⚠️ Failed to delete {pdb_file}: {e}")
    
    if deleted_count > 0:
        print(f"  🧹 Batch cleanup: Deleted {deleted_count} PDB files")
        
        # ディスク容量を表示
        import shutil
        disk_usage = shutil.disk_usage("/")
        free_gb = disk_usage.free / (1024**3)
        print(f"  💾 Free disk space: {free_gb:.1f} GB")

# ===== 並列処理用ワーカー関数 =====
def process_single_uniprot(args: tuple) -> dict:
    """
    単一のUniProt IDを処理（並列化用）
    
    Returns
    -------
    dict: 処理結果 {'success': bool, 'uniprotid': str, 'data': dict or 'error': str, 'error_type': str}
    """
    uniprotid, seq_ratio, max_pdbs, config_dict, clean_old_pdbs = args
    
    try:
        # Configを再構築（プロセス間で共有できないため）
        config = Config()
        dirpath = config.OUTPUT_DIR
        
        print(f"📄 [{uniprotid}] Starting...")
        
        unidata = UniprotData(uniprotid)
        fullName = unidata.get_fullname()
        organism = unidata.get_organism()
        
        # PDBリスト取得
        methods_to_search = []
        if config.USE_XRAY:
            methods_to_search.append("X-ray")
        if config.USE_NMR:
            methods_to_search.append("NMR")
        if config.USE_EM:
            methods_to_search.append("EM")

        if config.SEARCH_MODE == "AND" and len(methods_to_search) > 1:
            pdb_sets = []
            for method in methods_to_search:
                pdb_sets.append(set(unidata.pdblist({method})))
            current_pdblist = list(set.intersection(*pdb_sets))
        else:
            current_pdblist = unidata.pdblist(config.METHODS_SELECTED)
        
        # PDB数チェック
        if not count_pdb(uniprotid, methods=config.METHODS_SELECTED,
                       negative_pdbid=config.NEGATIVE_PDBID,
                       max_pdbs=max_pdbs):
            return {
                'success': False,
                'uniprotid': uniprotid,
                'error': 'Less than threshold PDB entries',
                'error_type': 'PDB_THRESHOLD'
            }
        
        # データ準備
        seqdata, all_pdblist = prep(uniprotid, 
                                   methods=config.METHODS_SELECTED,
                                   negative_pdbid=config.NEGATIVE_PDBID,
                                   max_pdbs=max_pdbs,
                                   verbose=False)  # 並列時はログ抑制
        
        seqdata1 = seqdata.filter(like=uniprotid)
        
        seqtype = 'nor+sub'
        pdbtuple = tuple(all_pdblist[0] + all_pdblist[1])
        pdb_used = all_pdblist[0] + all_pdblist[1]
        
        seqdata2 = seqdata.loc[:, seqdata.columns.str.startswith(pdbtuple)]
        norsub_seqdata = pd.concat([seqdata1, seqdata2], axis=1)
        
        # DSA実行
        sc_all, log_all, df_all, dist_all = run_DSA(
            uniprotid, norsub_seqdata, config.EXPORT, seqtype,
            methods=config.METHODS_SELECTED,
            seq_ratio=seq_ratio,
            cis_threshold=config.CIS_THRESHOLD,
            dirpath=dirpath,
            verbose=False
        )
        
        if df_all is None or len(df_all) == 0:
            return {
                'success': False,
                'uniprotid': uniprotid,
                'error': 'df_all is empty',
                'error_type': 'EMPTY_SUMMARY'
            }
        
        # score details保存
        save_score_details(uniprotid, sc_all, dist_all, seq_ratio, dirpath, None)
        
        # 統計データ作成
        mean_distance = sc_all['distance mean'].mean()
        std_distance = sc_all['distance std'].mean()
        mean_score = sc_all['score'].mean()
        umf = float(df_all.iloc[0]['UMF'])
        
        stats_df = pd.DataFrame({
            'uniprotid': [uniprotid],
            'seq_ratio': [seq_ratio],
            'fullName': [fullName],
            'organism': [organism],
            'Entries': [int(df_all.iloc[0]['Entries'])],
            'Chains': [int(df_all.iloc[0]['Chains'])],
            'Length': [int(df_all.iloc[0]['Length'])],
            'mean_distance': [round(mean_distance, 3)],
            'mean_std': [round(std_distance, 3)],
            'UMF': [round(umf, 3)],
            'mean_score': [round(mean_score, 3)]
        })
        
        # リンクデータ作成
        links_df = pd.DataFrame({
            'uniprotid': [uniprotid] * len(pdb_used),
            'pdbid': pdb_used,
            'seq_ratio': [seq_ratio] * len(pdb_used)
        })
        
        # summaryデータ作成
        row_dict = {
            'uniprotid': uniprotid,
            'seq_ratio': seq_ratio,
            'fullName': fullName,
            'organism': organism
        }
        for col in df_all.columns:
            if col in ['Entries', 'Chains', 'Length', 'Length(%)', 'Resolution', 
                      'UMF', 'cis/Length(%)', 'mean_cisDist', 'std_cisDist', 
                      'mean_cisScore', 'cis', 'mix']:
                row_dict[col] = df_all.iloc[0][col]
        
        print(f"✅ [{uniprotid}] Completed")
        
        # 🔴 追加：使用済みPDBファイルを削除（一時的に無効化）
        # cleanup_pdb_files_after_analysis(uniprotid, pdb_used)
        
        return {
            'success': True,
            'uniprotid': uniprotid,
            'data': {
                'row_dict': row_dict,
                'stats_df': stats_df,
                'links_df': links_df
            }
        }
        
    except Exception as e:
        # 🔴 追加：失敗時もクリーンアップ（一時的に無効化）
        # cleanup_pdb_files_after_analysis(uniprotid, [])
        
        error_msg = str(e)
        return {
            'success': False,
            'uniprotid': uniprotid,
            'error': error_msg,
            'error_type': classify_error_type(error_msg),
            'traceback': traceback.format_exc()
        }
    

# ===== メイン処理（並列版） =====
def main():
    """メイン処理（並列化対応 + 失敗ID自動記録）"""
    
    import argparse
    
    # ===== コマンドライン引数パーサー追加 =====
    parser = argparse.ArgumentParser(description='DSA Analysis - High Performance Mode')
    parser.add_argument('--ids', nargs='+', help='Specific UniProt IDs to analyze (e.g., --ids P01308 P00789)')
    parser.add_argument('--file', help='File containing UniProt IDs (one per line)')  
    parser.add_argument('--seq-ratio', type=float, default=20, help='Sequence ratio percentage (default: 20)')
    parser.add_argument('--max-pdbs', type=int, default=50, help='Maximum PDB entries per ID (default: 50)')
    parser.add_argument('--workers', type=int, default=7, help='Number of parallel workers (default: 7)')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch write size (default: 50)')
    parser.add_argument('--test-mode', action='store_true', help='Enable test mode')
    parser.add_argument('--test-count', type=int, default=433, help='Number of IDs for test mode (default: 433)')
    parser.add_argument('--no-parallel', action='store_true', help='Disable parallel processing (for debugging)')
    parser.add_argument('--skip-processed', action='store_true', default=True, help='Skip already processed IDs (default: True)')
    parser.add_argument('--no-skip', dest='skip_processed', action='store_false', help='Reprocess all IDs')
    
    args = parser.parse_args()
    
    # ===== ユーザー設定（引数から取得） =====
    seq_ratio = args.seq_ratio
    max_pdbs = args.max_pdbs
    
    # 並列処理設定
    ENABLE_PARALLEL = not args.no_parallel
    MAX_WORKERS = args.workers
    BATCH_SIZE = args.batch_size
    
    # テストモード設定
    TEST_MODE = args.test_mode
    TEST_COUNT = args.test_count
    
    # ===== Config初期化とID読み込み =====
    config = Config()
    
    # コマンドライン引数でIDが指定されていればそれを使用
    if args.ids:
        uniprot_ids = args.ids
        print(f"\n{'='*80}")
        print(f"🎯 SPECIFIC ID MODE")
        print(f"{'='*80}")
        print(f"Analyzing {len(uniprot_ids)} specific ID(s):")
        for uid in uniprot_ids:
            print(f"  • {uid}")
        print(f"{'='*80}\n")
    
    # ファイルが指定されている場合
    elif args.file:
        try:
            with open(args.file, 'r') as f:
                uniprot_ids = [line.strip() for line in f if line.strip()]
                print(f"\n{'='*80}")
                print(f"📄 FILE MODE")
                print(f"{'='*80}")
                print(f"Loaded {len(uniprot_ids)} IDs from {args.file}")
                print(f"{'='*80}\n")
        except FileNotFoundError:
            print(f"❌ Error: File not found: {args.file}")
            return
    else:
        # UniProt IDリストを読み込み
        uniprot_ids = config.load_uniprot_ids()
        
        # テストモード：指定数だけに制限
        if TEST_MODE and len(uniprot_ids) > TEST_COUNT:
            print(f"\n{'='*80}")
            print(f"🧪 TEST MODE ENABLED")
            print(f"{'='*80}")
            print(f"Original IDs: {len(uniprot_ids)}")
            print(f"Test IDs: {TEST_COUNT}")
            print(f"{'='*80}\n")
            uniprot_ids = uniprot_ids[:TEST_COUNT]
    
    use_pdb_search_results = False
    skip_processed = args.skip_processed
    clean_old_pdbs = True
    
    # 🔵 テストモード：指定数だけに制限
    if TEST_MODE and len(uniprot_ids) > TEST_COUNT:
        print(f"\n{'='*80}")
        print(f"🧪 TEST MODE ENABLED")
        print(f"{'='*80}")
        print(f"Original IDs: {len(uniprot_ids)}")
        print(f"Test IDs: {TEST_COUNT}")
        print(f"{'='*80}\n")
        uniprot_ids = uniprot_ids[:TEST_COUNT]
    
    use_pdb_search_results = False
    skip_processed = True
    clean_old_pdbs = True
    
    dirpath = config.OUTPUT_DIR
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    
    trim_archive_dir, backup_archive_dir = setup_archive_dirs(dirpath)
    archive_trim_and_summary(dirpath, trim_archive_dir, backup_archive_dir, 
                            verbose=config.VERBOSE)
    
    summaries_dir = os.path.join(dirpath, "summaries")
    score_details_dir = os.path.join(dirpath, "score_details")
    links_dir = os.path.join(dirpath, "links")
    
    for subdir in [summaries_dir, score_details_dir, links_dir]:
        if not os.path.exists(subdir):
            os.makedirs(subdir)
    
    # UniProt ID読み込み（既存ロジック - フォールバック）
    if not uniprot_ids and use_pdb_search_results:
        from pathlib import Path
        csv_file = Path("./output/links/unique_uniprots.csv")
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            uniprot_ids = df['uniprotid'].tolist()
            # 🔵 テストモード適用
            if TEST_MODE and len(uniprot_ids) > TEST_COUNT:
                uniprot_ids = uniprot_ids[:TEST_COUNT]
    
    # ===== 除外IDリスト読み込み =====
    EXCLUDED_IDS = config.load_excluded_ids()
    
    if EXCLUDED_IDS:
        uniprot_ids = [uid for uid in uniprot_ids if uid not in EXCLUDED_IDS]
    
    # ===== 失敗ID管理の初期化 =====
    failed_manager = FailedIDManager()
    MAX_RETRIES = 3
    
    # 処理済みIDと失敗IDを両方除外
    summary_file = os.path.join(summaries_dir, "summary.csv")
    processed_ids = load_processed_uniprots(summary_file, seq_ratio) if skip_processed else set()

    # 失敗IDもスキップ対象に追加
    failed_ids = failed_manager.get_failed_ids(seq_ratio, max_retries=MAX_RETRIES)
    skip_ids = processed_ids | failed_ids

    # 🔴 修正：--idsで指定された場合はスキップ処理をしない
    # さらに--fileで指定された場合もスキップしない（テスト用）
    if skip_processed and skip_ids and not args.ids and not args.file:  # ← ここに `and not args.file` を追加
        uniprot_ids = [uid for uid in uniprot_ids if uid not in skip_ids]
        if processed_ids or failed_ids:
            print(f"🚫 Skipping: {len(processed_ids)} processed + {len(failed_ids)} failed IDs")
    
    if not uniprot_ids:
        print("\n✅ All UniProt IDs already processed!")
        return
    
    print("\n" + "=" * 80)
    print(f"🚀 DSA Analysis - High Performance Mode")
    if TEST_MODE:
        print(f"🧪 TEST MODE: Analyzing {len(uniprot_ids)} IDs")
    print("=" * 80)
    print(f"📊 Total IDs: {len(uniprot_ids)}")
    print(f"⚙️  Parallel: {ENABLE_PARALLEL} (Workers: {MAX_WORKERS if ENABLE_PARALLEL else 1})")
    print(f"💾 Batch size: {BATCH_SIZE}")
    print(f"📦 PDB cache: Enabled with lock protection")
    print(f"Parameters: seq_ratio={seq_ratio}%, max_pdbs={max_pdbs}")
    print("=" * 80 + "\n")
    
    # ... 以下は既存のコードと同じ（変更なし）
    # サマリーファイル設定
    filename = os.path.join(summaries_dir, "summary.csv")
    fieldnames = [
        'uniprotid', 'seq_ratio', 'fullName', 'organism',
        'Entries', 'Chains', 'Length', 'Length(%)', 'Resolution', 'UMF',
        'cis/Length(%)', 'mean_cisDist', 'std_cisDist', 
        'mean_cisScore', 'cis', 'mix'
    ]
    
    # 既存ファイルバックアップ
    if os.path.exists(filename):
        jst = pytz.timezone('Asia/Tokyo')
        timestamp = datetime.datetime.now(jst).strftime('%Y%m%d_%H%M%S')
        backup_filename = os.path.join(dirpath, f"summary_backup_{timestamp}.csv")
        shutil.copy2(filename, backup_filename)
        print(f'📋 Backup: {backup_filename}')
    else:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
    
    # バッチライター初期化
    batch_writer = BatchWriter(batch_size=BATCH_SIZE)
    
    # 並列処理実行
    if ENABLE_PARALLEL:
        print(f"🔄 Processing {len(uniprot_ids)} IDs with {MAX_WORKERS} workers...\n")
        
        # マネージャーでプロセス間共有ロック辞書を作成
        with Manager() as manager:
            lock_dict = manager.dict()
            
            config_dict = {}
            task_args = [
                (uid, seq_ratio, max_pdbs, config_dict, clean_old_pdbs)
                for uid in uniprot_ids
            ]
            
            success_count = 0
            error_count = 0
            
            # initializerでロック辞書を渡す
            with ProcessPoolExecutor(max_workers=MAX_WORKERS, 
                                   initializer=init_worker, 
                                   initargs=(lock_dict,)) as executor:
                futures = {
                    executor.submit(process_single_uniprot, args): args[0]
                    for args in task_args
                }
                
                completed_batch = 0  # 🔴 追加

                for i, future in enumerate(as_completed(futures), 1):
                    uniprotid = futures[future]
                    try:
                        result = future.result()
                        
                        if result['success']:
                            batch_writer.add_summary(result['data']['row_dict'])
                            batch_writer.add_stats(result['data']['stats_df'])
                            batch_writer.add_links(result['data']['links_df'])
                            success_count += 1
                            
                            if batch_writer.should_flush():
                                batch_writer.flush(dirpath, fieldnames)
                            
                        else:
                            # 失敗を記録
                            error_type = result.get('error_type', classify_error_type(result.get('error', 'Unknown')))
                            failed_manager.record_failure(
                                uniprotid=result['uniprotid'],
                                seq_ratio=seq_ratio,
                                error_type=error_type,
                                error_message=result.get('error', 'Unknown error')
                            )
                            print(f"⚠️  [{uniprotid}] {error_type}: {result.get('error', 'Unknown')}")
                            error_count += 1
                        
                        # 進捗表示
                        print(f"Progress: {i}/{len(uniprot_ids)} "
                              f"(✓{success_count} ✗{error_count})")
                        
                        # 🔴 追加：7件ごとにバッチクリーンアップ（一時的に無効化）
                        # if i % 7 == 0:  # 7件ごとに
                        #     cleanup_batch_pdb_files()
                        #     print(f"🧹 Batch cleanup completed at {i}/{len(uniprot_ids)}")
                        
                    except Exception as e:
                        print(f"❌ [{uniprotid}] Fatal error: {e}")
                        error_count += 1
            
            # 残りのバッファを書き込み
            batch_writer.flush(dirpath, fieldnames)
        
    else:
        # シーケンシャル処理（デバッグ用）
        print(f"Processing {len(uniprot_ids)} IDs sequentially...\n")
        
        success_count = 0
        error_count = 0
        
        for i, uniprotid in enumerate(uniprot_ids, 1):
            config_dict = {}
            args = (uniprotid, seq_ratio, max_pdbs, config_dict, clean_old_pdbs)
            result = process_single_uniprot(args)
            
            if result['success']:
                batch_writer.add_summary(result['data']['row_dict'])
                batch_writer.add_stats(result['data']['stats_df'])
                batch_writer.add_links(result['data']['links_df'])
                success_count += 1
                
                if batch_writer.should_flush():
                    batch_writer.flush(dirpath, fieldnames)
            else:
                # 失敗を記録
                error_type = result.get('error_type', classify_error_type(result.get('error', 'Unknown')))
                failed_manager.record_failure(
                    uniprotid=result['uniprotid'],
                    seq_ratio=seq_ratio,
                    error_type=error_type,
                    error_message=result.get('error', 'Unknown error')
                )
                error_count += 1
            
            print(f"Progress: {i}/{len(uniprot_ids)} (✓{success_count} ✗{error_count})")
        
        batch_writer.flush(dirpath, fieldnames)
    
    # 失敗記録を保存
    failed_manager.save()
    
    # 失敗IDを除外リストに自動追加
    if config.VERBOSE:
        print("\n" + "=" * 80)
        print("📝 Updating excluded IDs list...")
        print("=" * 80)
    
    added_count = config.update_excluded_from_failed(
        seq_ratio=seq_ratio,
        max_retries=MAX_RETRIES
    )
    
    if added_count > 0 and config.VERBOSE:
        print(f"✅ Added {added_count} failed IDs to excluded list")
        print(f"   File: {config.OUTPUT_DIR}/excluded_ids.txt")
    
    # 失敗統計を表示
    stats = failed_manager.get_statistics(seq_ratio)
    if stats['total'] > 0:
        print("\n" + "=" * 80)
        print("📊 Failed ID Statistics")
        print("=" * 80)
        print(f"Total failed: {stats['total']}")
        print("\nBy error type:")
        for error_type, count in sorted(stats['by_error_type'].items(), 
                                       key=lambda x: x[1], reverse=True):
            print(f"  {error_type}: {count}")
        print("\nRetry distribution:")
        for retry_count, count in sorted(stats['retry_distribution'].items()):
            print(f"  {retry_count} attempt(s): {count} IDs")
        print("=" * 80)
    
    # 最終整理
    archive_trim_and_summary(dirpath, trim_archive_dir, backup_archive_dir, 
                            verbose=config.VERBOSE)
    
    print("\n" + "=" * 80)
    print("✅ All processing completed!")
    if TEST_MODE:
        print(f"🧪 TEST MODE: Completed {len(uniprot_ids)} IDs")
    print(f"📊 Final: ✓{success_count} successful, ✗{error_count} failed")
    if (success_count + error_count) > 0:
        print(f"   Success rate: {success_count/(success_count+error_count)*100:.1f}%")
    print("=" * 80)

if __name__ == "__main__":
    config = Config()
    trim_dir, backup_dir = setup_archive_dirs(config.OUTPUT_DIR)
    main()
    auto_archive_output()


