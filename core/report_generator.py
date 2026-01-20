"""
解析結果のレポート生成モジュール
"""

import os
import numpy as np
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP
from typing import List, Tuple


def generate_log_content(pdbdata: pd.DataFrame, len_sequence: int,
                        trimsequence: pd.DataFrame, score: pd.DataFrame,
                        cis_info: List[List]) -> pd.DataFrame:
    """解析結果のサマリーを生成"""
    cis_dist_mean, cis_dist_std, cis_score_mean, cis_num, mix = cis_info[0]
    
    cols = trimsequence.columns.values[1:]
    pdbids = [i.split(' ')[0] for i in cols]
    
    reso_list = []
    for pdbid in pdbids:
        reso = pdbdata.at['resolution', pdbid]
        if reso is not None:
            reso = ''.join(char for char in str(reso) 
                          if char.isdigit() or char == '.')
            if reso:
                reso_list.append(float(reso))
    
    if reso_list:
        reso_ave = Decimal(str(np.mean(reso_list))).quantize(
            Decimal('0.01'), rounding=ROUND_HALF_UP
        )
    else:
        reso_ave = Decimal('0.00')
    
    seted = sorted(set(pdbids), key=pdbids.index)
    
    return pd.DataFrame({
        'Entries': [len(seted)],
        'Chains': [len(pdbids)],
        'Length': [len(trimsequence)],
        'Length(%)': [round((len(trimsequence) * 100 / len_sequence), 1)],
        'Resolution': [reso_ave],
        'UMF': [round((score["distance mean"] / score["distance std"]).mean(), 1)],
        'cis/Length(%)': [round((cis_num * 100 / len(trimsequence)), 2)],
        'mean_cisDist': [round(cis_dist_mean, 2)],
        'std_cisDist': [round(cis_dist_std, 2)],
        'mean_cisScore': [round(cis_score_mean, 2)],
        'cis': [cis_num],
        'mix': [mix]
    })


def export_to_csv(uniprotid: str, seq_ratio: float, outputdataname: str,
                 outputdata: pd.DataFrame, seqtype: str, dirpath: str):
    """データをCSVファイルに出力"""
    filepath = os.path.join(
        dirpath, 
        f"{uniprotid}_{str(seq_ratio)}_{outputdataname}_{seqtype}.csv"
    )
    outputdata.to_csv(filepath, index=False)
