#!/usr/bin/env python3
"""
DSA (Distance Structure Analysis) メインスクリプト
UniProtIDから構造データを解析し、距離行列とスコアを計算します
統計ファイル出力機能付き
"""

import os
import re
import pandas as pd
import csv
import shutil
import datetime
import pytz
from typing import List, Tuple, Optional

from config import Config
from uniprot_handler import UniprotData
from structure_analyzer import CifData
from sequence_processor import (
    convert_three, trim_sequence, trim2_sequence, 
    sort_sequence, getcoord
)
from distance_calculator import getdistance2, getscore
from report_generator import generate_log_content, export_to_csv
from visualization import generate_heatmap


def clean_pdb_files_selective(keep_pdblist: List[str], verbose: bool = True):
    """
    指定されたPDBリストに含まれないファイルを削除
    
    Parameters
    ----------
    keep_pdblist : list
        保持するPDB IDのリスト
    verbose : bool
        削除ログを表示するか
    """
    keep_pdbs = set([pdb.upper() for pdb in keep_pdblist])
    deleted_count = 0
    
    # pdb_filesディレクトリのクリーンアップ
    if os.path.exists('pdb_files'):
        for filename in os.listdir('pdb_files'):
            if filename.startswith('.'):
                continue
            # ファイル名からPDB IDを抽出（例: 1abc.cif → 1ABC）
            pdb_id = filename.split('.')[0].upper()
            if pdb_id not in keep_pdbs:
                filepath = os.path.join('pdb_files', filename)
                os.remove(filepath)
                deleted_count += 1
                if verbose:
                    print(f"  Deleted: {filename}")
    
    # atom_coordディレクトリのクリーンアップ
    if os.path.exists('atom_coord'):
        for filename in os.listdir('atom_coord'):
            if filename.startswith('.'):
                continue
            # ファイル名からPDB IDを抽出（例: 1ABC.csv → 1ABC）
            pdb_id = filename.split('.')[0].upper()
            if pdb_id not in keep_pdbs:
                filepath = os.path.join('atom_coord', filename)
                os.remove(filepath)
                if verbose:
                    print(f"  Deleted: {filename}")
    
    if verbose and deleted_count > 0:
        print(f"Cleaned up {deleted_count} old PDB files")
    elif verbose:
        print("No old PDB files to clean")


def count_pdb(uniprotid: str, methods: Optional[set] = None, 
              negative_pdbid: str = "", max_pdbs: int = None) -> bool:
    """選択した構造決定法に限定してPDB数をカウント"""
    cfg = Config()
    
    if methods is None:
        methods = cfg.METHODS_SELECTED
    
    unidata = UniprotData(uniprotid)
    pdblist = unidata.pdblist(methods)
    
    if negative_pdbid:
        negative_list = re.split(r'[,\s]+', negative_pdbid.strip())
        negative_list_upper = [neg.upper() for neg in negative_list]
        pdblist = [item for item in pdblist 
                   if item.upper() not in negative_list_upper]
    
    # PDB数制限
    if max_pdbs is not None and len(pdblist) > max_pdbs:
        pdblist = pdblist[:max_pdbs]
    
    return len(pdblist) >= cfg.PDB_THRESHOLD


def prep(uniprotid: str, methods: Optional[set] = None, 
         negative_pdbid: str = "", max_pdbs: int = None,
         verbose: bool = True) -> Tuple:
    """データ準備: PDBファイルをダウンロードし配列情報を整理"""
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
    pdblist = unidata.pdblist(methods)
    
    if negative_pdbid:
        negative_list = re.split(r'[,\s]+', negative_pdbid.strip())
        negative_list_upper = [neg.upper() for neg in negative_list]
        pdblist = [item for item in pdblist 
                   if item.upper() not in negative_list_upper]
    
    # PDB数制限（ここで制限をかける）
    if max_pdbs is not None and len(pdblist) > max_pdbs:
        if verbose:
            print(f"  Limiting to first {max_pdbs} PDB entries (out of {len(pdblist)} available)")
        pdblist = pdblist[:max_pdbs]
    
    if verbose:
        print(f"  Processing {len(pdblist)} PDB entries ...")
    
    nor_pdblist = []
    sub_pdblist = []
    chi_pdblist = []
    din_pdblist = []
    
    for n, pdbid in enumerate(pdblist):
        cifdata = CifData(pdbid)
        mut_judge = cifdata.mutationjudge(uniprotids, pdbid)
        
        if verbose:
            print(f" ({n+1}/{len(pdblist)}) judge: {pdbid} {mut_judge}")
        
        if mut_judge == 'normal':
            nor_pdblist.append(pdbid)
        elif mut_judge == 'substitution':
            sub_pdblist.append(pdbid)
        elif mut_judge == 'chimera':
            chi_pdblist.append(pdbid)
        elif mut_judge == 'delins':
            din_pdblist.append(pdbid)
        else:
            continue
        
        beg, end = unidata.position(pdbid)
        df_beg = pd.DataFrame(index=list(range(beg-1)))
        df_end = pd.DataFrame(index=list(range(len_seqdata - end)))
        seq = cifdata.getsequence(uniprotids)
        seq = pd.concat([df_beg, seq, df_end])
        seq.reset_index(inplace=True, drop=True)
        seqdata = pd.concat([seqdata, seq], axis=1)
    
    all_pdblist = [nor_pdblist, sub_pdblist, chi_pdblist, din_pdblist]
    
    if verbose:
        total = sum(len(lst) for lst in all_pdblist)
        print(f" Data Preparation Finished: {total}/{len(pdblist)} PDB entries, "
              f"{len(seqdata.columns)-1} chains as {uniprotid}")
        print(f" (Normal: {len(nor_pdblist)}, Substitution: {len(sub_pdblist)}, "
              f"Chimera: {len(chi_pdblist)}, DelIns: {len(din_pdblist)})")
    
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
    
    # pdbdataを取得（summary_df生成に必要）
    if methods is None:
        methods = cfg.METHODS_SELECTED
    unidata.getpdbdata(methods)
    
    if seqdata is None or len(seqdata) == 0:
        print(f"Error: seqdata is empty for {uniprotid}")
        return None, "", None
    
    trimsequence = sort_sequence(str_ids, seqdata, seq_ratio)
    if trimsequence is None or len(trimsequence) == 0:
        print(f"Error: trimsequence is empty for {uniprotid}")
        return None, "", None
    
    trimsequence.to_csv(os.path.join(dirpath, f"trimsequence_{uniprotid}.csv"), 
                        index=False)
    trimseqcol = trimsequence.columns.values[1:]
    
    # CHAIN_THRESHOLDチェック
    if len(trimseqcol) <= cfg.CHAIN_THRESHOLD - 1:
        print(f"Error: Not enough chains for {uniprotid} (found {len(trimseqcol)}, need {cfg.CHAIN_THRESHOLD})")
        return None, "", None
    
    atomcoord = getcoord(trimsequence)
    if atomcoord is None or len(atomcoord) == 0:
        print(f"Error: atomcoord is empty for {uniprotid}")
        return None, "", None
    
    distance = getdistance2(atomcoord)
    if distance is None or len(distance) == 0:
        print(f"Error: distance is empty for {uniprotid}")
        return None, "", None
    
    score = getscore(distance, 0)
    
    # 距離データのCSV出力
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
    merged_df.to_csv(os.path.join(dirpath, f"distance_{uniprotid}.csv"), 
                     index=False, header=False)
    
    # cis解析
    cis_index = []
    for col in distance.columns.values.tolist()[2:]:
        tmp = distance.query(f'`{col}`<=@cis_threshold').index.to_list()
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
    
    # サマリーDataFrame生成
    summary_df = generate_log_content(unidata.pdbdata, len(sequence), 
                                      trimsequence, score, cis_info)
    
    # ログテキスト生成
    header_line = " ".join(summary_df.columns.astype(str).tolist())
    value_line = " ".join(str(v) for v in summary_df.iloc[0].tolist())
    log_text = header_line + "\n" + value_line
    
    return score, log_text, summary_df


def save_score_details(uniprotid: str, score: pd.DataFrame, seq_ratio: float, 
                       dirpath: str, existing_details: pd.DataFrame = None):
    """アプローチ1: 残基ペアごとの詳細データを保存"""
    # 新しいエントリを作成
    new_details = score[['residue pair', 'distance mean', 'distance std', 'score']].copy()
    new_details.insert(0, 'uniprotid', uniprotid)
    new_details.insert(1, 'seq_ratio', seq_ratio)
    
    # 既存データと結合
    if existing_details is not None and len(existing_details) > 0:
        # 同じuniprotidとseq_ratioのデータを削除
        existing_details = existing_details[
            ~((existing_details['uniprotid'] == uniprotid) & 
              (existing_details['seq_ratio'] == seq_ratio))
        ]
        combined_details = pd.concat([existing_details, new_details], ignore_index=True)
    else:
        combined_details = new_details
    
    # ファイル保存
    details_file = os.path.join(dirpath, "score_details.csv")
    combined_details.to_csv(details_file, index=False)
    
    return combined_details


def save_summary_statistics(uniprotid: str, fullName: str, organism: str,
                           score: pd.DataFrame, summary_df: pd.DataFrame,
                           seq_ratio: float, dirpath: str,
                           existing_stats: pd.DataFrame = None):
    """アプローチ2: タンパク質全体の統計サマリーを保存"""
    # 統計値を計算
    mean_distance = score['distance mean'].mean()
    std_distance = score['distance std'].mean()
    mean_score = score['score'].mean()
    
    # UMFは既存のsummary_dfから取得
    umf = float(summary_df.iloc[0]['UMF'])
    
    # 新しいエントリを作成
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
    
    # 既存データと結合
    if existing_stats is not None and len(existing_stats) > 0:
        # 同じuniprotidとseq_ratioのデータを削除
        existing_stats = existing_stats[
            ~((existing_stats['uniprotid'] == uniprotid) & 
              (existing_stats['seq_ratio'] == seq_ratio))
        ]
        combined_stats = pd.concat([existing_stats, new_stats], ignore_index=True)
    else:
        combined_stats = new_stats
    
    # ファイル保存
    stats_file = os.path.join(dirpath, "summary_statistics.csv")
    combined_stats.to_csv(stats_file, index=False)
    
    return combined_stats


def save_uniprot_pdb_links(uniprotid: str, pdblist: List[str], 
                           seq_ratio: float, dirpath: str,
                           existing_links: pd.DataFrame = None):
    """UniProt IDとPDB IDのリンクファイルを作成
    
    Parameters
    ----------
    uniprotid : str
        UniProt ID
    pdblist : list
        PDB IDのリスト
    seq_ratio : float
        使用した配列比率
    dirpath : str
        出力ディレクトリ
    existing_links : pd.DataFrame
        既存のリンクデータ
    
    Returns
    -------
    pd.DataFrame
        更新されたリンクデータ
    """
    # 新しいリンクデータを作成
    new_links = pd.DataFrame({
        'uniprotid': [uniprotid] * len(pdblist),
        'pdbid': pdblist,
        'seq_ratio': [seq_ratio] * len(pdblist)
    })
    
    # 既存データと結合
    if existing_links is not None and len(existing_links) > 0:
        # 同じuniprotidとseq_ratioのデータを削除
        existing_links = existing_links[
            ~((existing_links['uniprotid'] == uniprotid) & 
              (existing_links['seq_ratio'] == seq_ratio))
        ]
        combined_links = pd.concat([existing_links, new_links], ignore_index=True)
    else:
        combined_links = new_links
    
    # ファイル保存
    links_file = os.path.join(dirpath, "uniprot_pdb_links.csv")
    combined_links.to_csv(links_file, index=False)
    
    return combined_links


def main():
    """メイン処理"""
    # === ユーザー設定可能なパラメータ ===
    seq_ratio = 20  # 解析に使用する配列長の割合(%)
    max_pdbs = 50   # 処理するPDB数の上限（Noneで無制限）
    uniprot_ids = ["P25156"]  # UniProt IDリスト
    clean_old_pdbs = True  # 前のUniProt IDのPDBファイルを削除するか
    # ====================================
    
    config = Config()
    
    # 出力ディレクトリ設定
    dirpath = config.OUTPUT_DIR
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    
    print(f"Output directory: {os.path.abspath(dirpath)}")
    print(f"Parameters: seq_ratio={seq_ratio}%, max_pdbs={max_pdbs if max_pdbs else 'unlimited'}")
    print(f"Clean old PDBs: {clean_old_pdbs}")
    print("=" * 80)
    
    # サマリーファイル設定
    filename = os.path.join(dirpath, "summary.csv")
    fieldnames = [
        'uniprotid', 'seq_ratio', 'fullName', 'organism',
        'Entries', 'Chains', 'Length', 'Length(%)', 'Resolution', 'UMF',
        'cis/Length(%)', 'mean_cisDist', 'std_cisDist', 
        'mean_cisScore', 'cis', 'mix'
    ]
    
    # 既存データ読み込み
    existing_data = []
    if os.path.exists(filename):
        jst = pytz.timezone('Asia/Tokyo')
        timestamp = datetime.datetime.now(jst).strftime('%Y%m%d_%H%M%S')
        backup_filename = os.path.join(dirpath, f"summary_backup_{timestamp}.csv")
        shutil.copy2(filename, backup_filename)
        print(f'Backup created: {backup_filename}')
        
        with open(filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                existing_data.append(row)
    else:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
    
    # 統計ファイルの既存データ読み込み
    details_file = os.path.join(dirpath, "score_details.csv")
    stats_file = os.path.join(dirpath, "summary_statistics.csv")
    links_file = os.path.join(dirpath, "uniprot_pdb_links.csv")
    
    existing_details = None
    existing_stats = None
    existing_links = None
    
    if os.path.exists(details_file):
        existing_details = pd.read_csv(details_file)
        print(f'Loaded existing score_details.csv')
    
    if os.path.exists(stats_file):
        existing_stats = pd.read_csv(stats_file)
        print(f'Loaded existing summary_statistics.csv')
    
    if os.path.exists(links_file):
        existing_links = pd.read_csv(links_file)
        print(f'Loaded existing uniprot_pdb_links.csv')
    
    # 各IDを処理
    for uniprotid in uniprot_ids:
        try:
            print("=" * 80)
            print(f"Processing {uniprotid} ...")
            
            unidata = UniprotData(uniprotid)
            fullName = unidata.get_fullname()
            organism = unidata.get_organism()
            print(f"{fullName} from {organism}")
            
            # 現在のUniProt IDに必要なPDBリストを取得
            unidata_temp = UniprotData(uniprotid)
            current_pdblist = unidata_temp.pdblist(config.METHODS_SELECTED)
            
            # 不要なPDBファイルを削除
            if clean_old_pdbs:
                print(f"\nCleaning old PDB files...")
                clean_pdb_files_selective(current_pdblist, verbose=config.VERBOSE)
            
            # PDB数チェック（max_pdbsを渡す）
            if not count_pdb(uniprotid, methods=config.METHODS_SELECTED,
                           negative_pdbid=config.NEGATIVE_PDBID,
                           max_pdbs=max_pdbs):
                print("Less than threshold PDB entries")
                continue
            
            # データ準備（max_pdbsを渡す）
            seqdata, all_pdblist = prep(uniprotid, 
                                       methods=config.METHODS_SELECTED,
                                       negative_pdbid=config.NEGATIVE_PDBID,
                                       max_pdbs=max_pdbs,
                                       verbose=config.VERBOSE)
            
            seqdata1 = seqdata.filter(like=uniprotid)
            
            # Normal + Substitution解析
            seqtype = 'nor+sub'
            pdbtuple = tuple(all_pdblist[0] + all_pdblist[1])
            pdb_used = all_pdblist[0] + all_pdblist[1]  # リスト形式で保持
            print(f"\n### normal & mutant ###")
            print(f"PDB: {pdbtuple}")
            print(f"{len(pdbtuple)} entries were processed")
            
            seqdata2 = seqdata.loc[:, seqdata.columns.str.startswith(pdbtuple)]
            norsub_seqdata = pd.concat([seqdata1, seqdata2], axis=1)
            
            # DSA実行（seq_ratioを渡す）
            sc_all, log_all, df_all = run_DSA(
                uniprotid, norsub_seqdata, config.EXPORT, seqtype,
                methods=config.METHODS_SELECTED,
                seq_ratio=seq_ratio,  # 設定したseq_ratioを使用
                cis_threshold=config.CIS_THRESHOLD,
                dirpath=dirpath,
                verbose=config.VERBOSE
            )
            
            if df_all is None or len(df_all) == 0:
                print(f"Error: df_all is empty for {uniprotid}")
                continue
            
            # アプローチ1: 残基ペアごとの詳細データを保存
            print(f"\nSaving score details...")
            existing_details = save_score_details(
                uniprotid, sc_all, seq_ratio, dirpath, existing_details
            )
            print(f"  -> score_details.csv updated")
            
            # アプローチ2: タンパク質全体の統計サマリーを保存
            print(f"Saving summary statistics...")
            existing_stats = save_summary_statistics(
                uniprotid, fullName, organism, sc_all, df_all, 
                seq_ratio, dirpath, existing_stats
            )
            print(f"  -> summary_statistics.csv updated")
            
            # UniProt-PDBリンクを保存
            print(f"Saving UniProt-PDB links...")
            existing_links = save_uniprot_pdb_links(
                uniprotid, pdb_used, seq_ratio, dirpath, existing_links
            )
            print(f"  -> uniprot_pdb_links.csv updated ({len(pdb_used)} PDB entries)")
            
            # エントリ作成(既存のsummary.csv用)
            row0 = df_all.iloc[0]
            new_entry = {
                'uniprotid': uniprotid,
                'seq_ratio': float(seq_ratio),  # 設定したseq_ratioを使用
                'fullName': fullName,
                'organism': organism,
                'Entries': int(row0['Entries']),
                'Chains': int(row0['Chains']),
                'Length': int(row0['Length']),
                'Length(%)': float(row0['Length(%)']),
                'Resolution': float(row0['Resolution']),
                'UMF': float(row0['UMF']),
                'cis/Length(%)': float(row0['cis/Length(%)']),
                'mean_cisDist': float(row0['mean_cisDist']),
                'std_cisDist': float(row0['std_cisDist']),
                'mean_cisScore': float(row0['mean_cisScore']),
                'cis': int(row0['cis']),
                'mix': int(row0['mix']),
            }
            
            # テキスト出力
            txtfilepath = os.path.join(dirpath, 
                                      f"{uniprotid}_{seq_ratio}_summary.txt")
            with open(txtfilepath, 'w') as f:
                print(f"### {uniprotid} (seq_ratio={seq_ratio}) Summary ###", 
                      file=f)
                print(fullName, file=f)
                print(organism, file=f)
                print(df_all.to_string(index=False), file=f)
            
            # 既存データ更新
            found = False
            for i, existing_entry in enumerate(existing_data):
                if (existing_entry['uniprotid'] == new_entry['uniprotid'] and
                    float(existing_entry['seq_ratio']) == new_entry['seq_ratio']):
                    if config.OVERWRITE:
                        existing_data[i] = {k: str(v) for k, v in new_entry.items()}
                    found = True
                    break
            if not found:
                existing_data.append({k: str(v) for k, v in new_entry.items()})
            
            print(f"Processing {uniprotid} Finished")
            
        except Exception as e:
            print(f"Error processing {uniprotid}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # CSV書き出し(既存のsummary.csv)
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in existing_data:
            writer.writerow(row)
    
    print(f"\n{'=' * 80}")
    print(f"Updated files:")
    print(f"  1. {os.path.abspath(filename)} (既存のサマリー)")
    print(f"  2. {os.path.abspath(details_file)} (残基ペアごとの詳細)")
    print(f"  3. {os.path.abspath(stats_file)} (タンパク質ごとの統計)")
    print(f"  4. {os.path.abspath(links_file)} (UniProt-PDBリンク)")
    print(f"{'=' * 80}")
    print("Job Completed")


if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    main()