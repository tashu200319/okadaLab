"""
配列処理モジュール (型不一致修正版)
"""

import os
import pandas as pd
import gzip
from mimetypes import guess_type
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from typing import Dict
from core.structure_analyzer import downloadpdb


def _open_cif(pdbid: str):
    """CIFファイルを開く(gzip対応)"""
    file = pdbid.lower() + ".cif"
    ciffile = "pdb_files/" + file
    
    if guess_type(file)[1] == "gzip":
        return gzip.open(ciffile, mode='rt')
    else:
        return open(ciffile)


def convert_three(sequence: str) -> list:
    """1文字アミノ酸コードを3文字コードに変換"""
    dic = {
        "A": "ALA", "B": "D|N", "C": "CYS", "D": "ASP", 
        "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", 
        "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET",
        "N": "ASN", "O": "HYP", "P": "PRO", "Q": "GLN", 
        "R": "ARG", "S": "SER", "T": "THR", "U": "SEC", 
        "V": "VAL", "W": "TRP", "X": "any", "Y": "TYR", 
        "Z": "E|Q"
    }
    return [dic[char] for char in sequence]


def trim_sequence(sequencedata: pd.DataFrame, seq_ratio: float = 80) -> pd.DataFrame:
    """座標データが存在しているアミノ酸残基の割合がseq_ratio以下であれば削除"""
    sequencedata.dropna(subset=sequencedata.columns[0], inplace=True)
    seqlen = len(sequencedata)
    
    delchain = [
        chain for chain, item in sequencedata.items() 
        if 100 - (item.isnull().sum() / seqlen * 100) < seq_ratio
    ]
    
    seqdata = sequencedata.drop(columns=delchain)
    seqdata.dropna(inplace=True)
    
    return seqdata


def trim2_sequence(sequencedata: pd.DataFrame, seq_ratio: float = 80) -> pd.DataFrame:
    """seq_id重複の場合、最初だけを残して残りは削除"""
    seq = sequencedata.iloc[:, 1:].map(
        lambda x: int(x.split(', ')[1]) if isinstance(x, str) else x
    )
    
    duplicate_indices = set()
    
    for column in seq.columns:
        duplicates = seq[column].duplicated(keep='first')
        duplicate_indices.update(seq[duplicates].index)
    
    duplicate_indices = sorted(list(duplicate_indices))
    trim2_seq = sequencedata.drop(index=duplicate_indices)
    
    return trim2_seq


def _diff(uniprotid: str, df1: pd.Series, df2: pd.Series, shift: int = 0) -> int:
    """配列の一致度を計算"""
    diff = pd.concat([df1, df2.shift(shift)], axis=1)
    diff.dropna(inplace=True)
    diff.drop_duplicates(subset=uniprotid, ignore_index=True, inplace=True)
    return (diff.iloc[:, 0] == diff.iloc[:, 1]).sum()


def sort_sequence(uniprotid: str, sequencedata: pd.DataFrame, seq_ratio: float) -> pd.DataFrame:
    """配列のソートとアライメント調整"""
    seq = sequencedata.map(lambda x: x.split(', ')[0] if isinstance(x, str) else x)
    
    trimdata = trim_sequence(seq, seq_ratio)
    trimdata.drop_duplicates(subset=uniprotid, ignore_index=True, inplace=True)
    trimdata.reset_index(inplace=True, drop=True)
    trimdata = trimdata.T
    
    columns = trimdata.columns
    IDs = []
    
    for col in columns:
        diff = trimdata[trimdata[col] != trimdata.at[uniprotid, col]].index
        if len(diff) != 0:
            IDs.extend(diff)
            trimdata.drop(diff, inplace=True)
    
    uniseq = seq[uniprotid]
    
    for ID in IDs:
        difseq = seq[ID]
        unique = _diff(uniprotid, uniseq, difseq)
        
        if unique > 10:
            continue
        
        num = 1
        unique = 0
        
        while unique < 10 and num < 100:
            unique = _diff(uniprotid, uniseq, difseq, num)
            num = (-num) + 1 if num < 0 else -num
        
        if unique > 10:
            diff = sequencedata[ID].shift((-num) + 1 if num > 0 else -num)
            loc = sequencedata.columns.get_loc(ID)
            sequencedata.drop(ID, axis=1, inplace=True)
            sequencedata.insert(loc, ID, diff)
        else:
            sequencedata.drop(ID, axis=1, inplace=True)
    
    sorted_seqdata = trim_sequence(sequencedata, seq_ratio)
    uniq_sorted_seqdata = trim2_sequence(sorted_seqdata)
    
    return uniq_sorted_seqdata


def getcoord(
    trimsequence: pd.DataFrame,
    uniprotid: str,
    *,
    verbose: bool = False,
    logger=None,
) -> pd.DataFrame:
    """原子座標を取得(PDBファイルから直接) - 型不一致を修正"""
    atomcoord = pd.DataFrame(trimsequence.iloc[:, 0])
    atomindex = atomcoord.index.tolist()
    trimseq = trimsequence.iloc[:, 1:].map(
        lambda x: int(x.split(', ')[1]) if isinstance(x, str) else x
    )
    columns = trimseq.columns.tolist()
    
    pdbids = {}
    for col in columns:
        pdbid, strand_id = col.split(' ')
        pdbids.setdefault(pdbid, []).append(strand_id)
    
    for pdbid, chain_id in pdbids.items():
        if not downloadpdb(pdbid):
            continue
        
        try:
            with _open_cif(pdbid) as handle:
                mmcifdict = MMCIF2Dict(handle)
        except Exception as e:
            if logger:
                logger.warning(f"Error reading CIF file for {pdbid}: {e}")
            elif verbose:
                print(f"Error reading CIF file for {pdbid}: {e}")
            continue

        
        try:
            model_num = mmcifdict.get("_atom_site.pdbx_PDB_model_num", [])
            asym_id = mmcifdict.get("_atom_site.auth_asym_id", [])
            comp_id = mmcifdict.get("_atom_site.auth_comp_id", [])
            seq_id = mmcifdict.get("_atom_site.auth_seq_id", [])
            atom_id = mmcifdict.get("_atom_site.auth_atom_id", [])
            cartn_x = mmcifdict.get("_atom_site.Cartn_x", [])
            cartn_y = mmcifdict.get("_atom_site.Cartn_y", [])
            cartn_z = mmcifdict.get("_atom_site.Cartn_z", [])
            alt_id = mmcifdict.get("_atom_site.label_alt_id", [])
            group_PDB = mmcifdict.get("_atom_site.group_PDB", [])
            
            lengths = [len(model_num), len(asym_id), len(comp_id), len(seq_id),
                      len(atom_id), len(cartn_x), len(cartn_y), len(cartn_z),
                      len(alt_id), len(group_PDB)]
            if len(set(lengths)) > 1:
                min_length = min(lengths)
                model_num = model_num[:min_length]
                asym_id = asym_id[:min_length]
                comp_id = comp_id[:min_length]
                seq_id = seq_id[:min_length]
                atom_id = atom_id[:min_length]
                cartn_x = cartn_x[:min_length]
                cartn_y = cartn_y[:min_length]
                cartn_z = cartn_z[:min_length]
                alt_id = alt_id[:min_length]
                group_PDB = group_PDB[:min_length]
            
            struct = pd.DataFrame({
                "model_num": model_num,
                "asym_id": asym_id,
                "comp_id": comp_id,
                "seq_id": seq_id,
                "atom_id": atom_id,
                "Cartn_x": cartn_x,
                "Cartn_y": cartn_y,
                "Cartn_z": cartn_z,
                "alt_id": alt_id,
                "group_PDB": group_PDB
            })
            
            struct["Cartn_x"] = pd.to_numeric(struct["Cartn_x"], errors='coerce')
            struct["Cartn_y"] = pd.to_numeric(struct["Cartn_y"], errors='coerce')
            struct["Cartn_z"] = pd.to_numeric(struct["Cartn_z"], errors='coerce')
            
            struct = struct[
                (struct['atom_id'] == 'CA') & 
                (struct['group_PDB'] == 'ATOM')
            ]
            
            struct['original_index'] = struct.index
            alt_id_dot = struct[struct['alt_id'].str.contains(r'\.', na=False)]
            alt_id_not_dot = struct[~struct['alt_id'].str.contains(r'\.', na=False)]
            alt_id_not_dot_unique = alt_id_not_dot.drop_duplicates(
                subset=['seq_id', 'atom_id']
            )
            struct = pd.concat([alt_id_dot, alt_id_not_dot_unique])
            struct = struct.sort_values('original_index')
            struct = struct.drop(columns=['original_index', 'alt_id', 'group_PDB', 'model_num', 'atom_id'])
            
        except Exception as e:
            if logger:
                logger.warning(f"Error extracting atom coordinates for {pdbid}: {e}")
            elif verbose:
                print(f"Error extracting atom coordinates for {pdbid}: {e}")
            continue

        
        for chain in chain_id:
            seq_num = trimseq[pdbid + ' ' + chain]
            
            chaindata = struct[struct["asym_id"] == chain].copy()
            
            # ====== デバッグログ（verbose制御） ======
            if verbose:
                print(f"  [{pdbid} {chain}] chaindata before filter: {len(chaindata)} rows")
                if len(chaindata) > 0:
                    print(f"    seq_id range: {chaindata['seq_id'].min()} - {chaindata['seq_id'].max()}")
                    print(f"    seq_id sample: {chaindata['seq_id'].head(3).tolist()}")
                    print(f"    seq_id dtype: {chaindata['seq_id'].dtype}")

            # 🔴 修正: seq_idを整数型に統一
            try:
                chaindata['seq_id'] = pd.to_numeric(chaindata['seq_id'], errors='coerce')
                chaindata = chaindata.dropna(subset=['seq_id'])
                chaindata['seq_id'] = chaindata['seq_id'].astype(int)
            except Exception as e:
                if logger:
                    logger.warning(f"[{pdbid} {chain}] Failed to convert seq_id to int: {e}")
                elif verbose:
                    print(f"    ⚠️  Failed to convert seq_id to int: {e}")
                continue

            # 🔴 修正: 整数インデックスを設定
            chaindata = chaindata.set_index('seq_id')

            # 重複インデックスがあれば最初だけ残す
            if chaindata.index.duplicated().any():
                chaindata = chaindata[~chaindata.index.duplicated(keep='first')]

            # ====== フィルタ対象 seq_num の確認ログ ======
            if verbose:
                if len(seq_num) > 0:
                    t = type(seq_num.iloc[0])
                    sample = seq_num.head(3).tolist()
                else:
                    t = 'empty'
                    sample = []
                print(f"    seq_num to filter (type: {t}): {sample}")
                print(
                    f"    chaindata.index (type: {chaindata.index.dtype}): "
                    f"{chaindata.index[:3].tolist() if len(chaindata) > 0 else []}"
                    )

            # 🔴 修正: seq_numも整数に統一してからフィルタリング
            seq_num_int = seq_num.astype(int)

            coord = chaindata[['comp_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']].reindex(
                seq_num_int.tolist()
            )

            # ====== フィルタ後の結果 ======
            if verbose:
                print(f"    coord after filter: {len(coord)} rows")

            if len(coord) == 0 or coord.isna().all().all():
                if logger:
                    logger.debug(f"[{pdbid} {chain}] No matching coordinates found")
                elif verbose:
                    print(f"    ⚠️  No matching coordinates found!")
                continue

            
            # インデックスをatomindexに戻す
            coord.index = seq_num.index
            coord = pd.concat([seq_num, coord], axis=1)
            coord.drop(columns=pdbid + ' ' + chain, inplace=True)
            coord.rename(columns={'comp_id': pdbid + ' ' + chain}, inplace=True)
            coord.index = atomindex
            atomcoord = pd.concat([atomcoord, coord], axis=1)
    
    if len(atomcoord.columns) <= 1:
        if logger:
            logger.warning(f"No coordinate columns added for {uniprotid}")
        elif verbose:
            print(f"Warning: No coordinate columns added for {uniprotid}")
        return atomcoord
    
    if verbose:
        print(f"Debug: Before dropna - {len(atomcoord)} rows, {len(atomcoord.columns)} columns")
    
    coord_cols = atomcoord.columns[1:]
    atomcoord = atomcoord.dropna(subset=coord_cols)
    
    if verbose:
        if atomcoord.empty:
            print(f"Debug: atomcoord is empty AFTER dropna for {uniprotid}")
        else:
            print(f"Debug: After dropna - {len(atomcoord)} rows")

    return atomcoord 