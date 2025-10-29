"""
原子間距離計算モジュール
"""

import numpy as np
import pandas as pd
from itertools import combinations
from numba import jit


@jit(nopython=True)
def calculat(atom1: np.ndarray, atom2: np.ndarray) -> float:
    """
    2つの原子間の距離を計算（高精度版）
    
    Parameters
    ----------
    atom1 : np.ndarray
        原子1の座標 [x, y, z]
    atom2 : np.ndarray
        原子2の座標 [x, y, z]
    
    Returns
    -------
    float
        距離（Å）
    """
    xyz = atom1 - atom2
    xyz = np.rint(xyz * 1000)
    dis = np.sqrt(np.sum(xyz**2))
    return dis / 1000


def getdistance2(atomcoord: pd.DataFrame) -> pd.DataFrame:
    """
    全残基ペア間の距離を計算
    
    Parameters
    ----------
    atomcoord : pd.DataFrame
        原子座標データ
    
    Returns
    -------
    pd.DataFrame
        距離行列
    """
    id_col = atomcoord.iloc[:, 0].name
    cols = atomcoord.iloc[:, 1::4].columns.tolist()
    
    distance = pd.DataFrame({
        id_col: [str(int(a)+1) + ", " + str(int(b)+1) 
                for a, b in combinations(map(str, atomcoord.index), 2)],
        "residue pair": [resi0 + ", " + resi1 
                        for resi0, resi1 in combinations(atomcoord[id_col], 2)],
        **{col: np.nan for col in cols}
    })
    
    combination = list(combinations(range(len(atomcoord)), 2))
    
    for i, col in enumerate(cols):
        i = (i * 4) + 2
        atoms = atomcoord.iloc[:, i:i+3].to_numpy()
        distance[col] = [calculat(atoms[n1], atoms[n2]) 
                        for n1, n2 in combination]
    
    return distance


def getscore(distance: pd.DataFrame, ddof: int = 0) -> pd.DataFrame:
    """
    距離の平均・標準偏差・スコアを計算
    
    Parameters
    ----------
    distance : pd.DataFrame
        距離データ
    ddof : int
        標準偏差の自由度補正（0: 母集団, 1: 標本）
    
    Returns
    -------
    pd.DataFrame
        スコアデータ
    """
    dis = distance.iloc[:, 2:]
    means = dis.mean(axis='columns')
    stds = dis.std(axis='columns', ddof=ddof)
    stds = stds.map(lambda x: 0.0001 if x == 0 else x)
    
    column0 = distance.columns[0]
    
    return pd.DataFrame({
        column0: distance[column0],
        "residue pair": distance["residue pair"],
        "distance mean": means,
        "distance std": stds,
        "score": means / stds
    })


def getscore_cis(distance: pd.DataFrame, ddof: int = 0) -> pd.DataFrame:
    """
    cis配列の距離スコアを計算（getscoreと同一機能）
    """
    return getscore(distance, ddof)
