"""
配列処理モジュール
"""

import os
import pandas as pd
from typing import Dict


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
            print(f"{ID} is not used due to sequence alignment failure")
            sequencedata.drop(ID, axis=1, inplace=True)
    
    sorted_seqdata = trim_sequence(sequencedata, seq_ratio)
    uniq_sorted_seqdata = trim2_sequence(sorted_seqdata)
    
    return uniq_sorted_seqdata


def getcoord(trimsequence: pd.DataFrame) -> pd.DataFrame:
    """原子座標を取得"""
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
        struct = pd.read_csv(f'atom_coord/{pdbid}.csv')
        struct["asym_id"] = struct["asym_id"].astype(str)
        struct = struct[struct["atom_id"] == "CA"]
        struct.drop(columns=['model_num', 'atom_id'], inplace=True)
        
        for chain in chain_id:
            seq_num = trimseq[pdbid + ' ' + chain]
            seq_num.index = seq_num.tolist()
            
            chaindata = struct[struct["asym_id"] == chain]
            chaindata.index = chaindata["seq_id"].tolist()
            
            if chaindata['seq_id'].duplicated().any():
                chaindata = chaindata.drop_duplicates(subset='seq_id', keep='first')
            
            coord = chaindata[['comp_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']]
            coord = chaindata[['comp_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']].filter(
                items=seq_num.tolist(), axis=0
            )
            
            coord = pd.concat([seq_num, coord], axis=1)
            coord.drop(columns=pdbid + ' ' + chain, inplace=True)
            coord.rename(columns={'comp_id': pdbid + ' ' + chain}, inplace=True)
            coord.index = atomindex
            atomcoord = pd.concat([atomcoord, coord], axis=1)
    
    atomcoord.dropna(inplace=True)
    
    return atomcoord
