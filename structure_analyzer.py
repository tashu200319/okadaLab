"""
PDB/mmCIF構造ファイルの解析モジュール
"""

import os
import gzip
from mimetypes import guess_type
import pandas as pd
from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from typing import List


pdb_list = PDBList()


def downloadpdb(pdbid: str):
    """
    PDBファイルをダウンロード
    
    Parameters
    ----------
    pdbid : str
        PDB ID (例: 1ABC)
    """
    pdb_list.retrieve_pdb_file(pdbid, pdir="pdb_files/", file_format="mmCif")


def _open(pdbid: str):
    """CIFファイルを開く（gzip対応）"""
    file = pdbid.lower() + ".cif"
    ciffile = "pdb_files/" + file
    
    if guess_type(file)[1] == "gzip":
        return gzip.open(ciffile, mode='rt')
    else:
        return open(ciffile)


class CifData:
    """
    Cifファイルを解析し、配列情報を取得
    
    参考:
    https://mmcif.pdbj.org/docs/pdb_to_pdbx_correspondences.html#DBREF
    https://mmcif.pdbj.org/dictionaries/mmcif_pdbx_v50.dic/Items/index.html
    """
    
    def __init__(self, pdbid: str):
        """
        Parameters
        ----------
        pdbid : str
            PDB ID
        """
        self.pdbid = pdbid
        downloadpdb(self.pdbid)
        
        with _open(self.pdbid) as handle:
            mmcifdict = MMCIF2Dict(handle)
        
        # struct_ref_seqの構築
        self.struct_ref_seq = pd.DataFrame({
            "strand_id": mmcifdict["_struct_ref_seq.pdbx_strand_id"],
            "accession": [i.upper() for i in mmcifdict["_struct_ref_seq.pdbx_db_accession"]],
            "seq_align_beg": mmcifdict["_struct_ref_seq.seq_align_beg"],
            "seq_align_end": mmcifdict["_struct_ref_seq.seq_align_end"]
        })
        
        pdb_strand_id = mmcifdict["_pdbx_poly_seq_scheme.pdb_strand_id"]
        for i, struct_strand_id in enumerate(self.struct_ref_seq["strand_id"]):
            self.struct_ref_seq.at[i, "sort_index"] = pdb_strand_id.index(struct_strand_id)
        
        self.struct_ref_seq.sort_values("sort_index", inplace=True)
        
        # struct_ref_seq_difの構築
        try:
            self.struct_ref_seq_dif = pd.DataFrame({
                "strand_id": mmcifdict["_struct_ref_seq_dif.pdbx_pdb_strand_id"],
                "seq_num": mmcifdict["_struct_ref_seq_dif.pdbx_auth_seq_num"],
                "db_seq_num": mmcifdict["_struct_ref_seq_dif.pdbx_seq_db_seq_num"],
                "details": [i.lower() for i in mmcifdict["_struct_ref_seq_dif.details"]]
            })
            # フィルタリング
            self.struct_ref_seq_dif = self.struct_ref_seq_dif[
                ~self.struct_ref_seq_dif['details'].isin([
                    "expression tag", "linker", "conflict", "microheterogeneity"
                ])
            ]
        except KeyError:
            self.struct_ref_seq_dif = pd.DataFrame({
                "strand_id": [], "seq_num": [], "db_seq_num": []
            })
        
        # Chain情報の構築
        self._build_chain_info(mmcifdict)
        
        # 原子座標の取得
        self._extract_atom_coord(mmcifdict)
    
    def _build_chain_info(self, mmcifdict):
        """Chain情報を構築"""
        self.chain = []
        self.chainid = []
        self.hetero_info = []
        self.ind = -1
        hetero_pdb_seq_num = ""
        
        for pdb_mon_id, pdb_seq_num, hetero, chainid in zip(
            mmcifdict["_pdbx_poly_seq_scheme.pdb_mon_id"],
            mmcifdict["_pdbx_poly_seq_scheme.pdb_seq_num"],
            mmcifdict["_pdbx_poly_seq_scheme.hetero"],
            mmcifdict["_pdbx_poly_seq_scheme.pdb_strand_id"]
        ):
            self.ind += 1
            
            if hetero == "n":
                hetero_pdb_seq_num = ""
                if pdb_mon_id != "?":
                    self.chain.append(pdb_mon_id + ", " + pdb_seq_num)
                    self.chainid.append(chainid)
                else:
                    self.chain.append(None)
                    self.chainid.append(chainid)
            else:
                if pdb_seq_num == hetero_pdb_seq_num:
                    self.hetero_info.append(self.ind)
                    continue
                else:
                    if pdb_mon_id != "?":
                        self.chain.append(pdb_mon_id + ", " + pdb_seq_num)
                        self.chainid.append(chainid)
                        hetero_pdb_seq_num = pdb_seq_num
                    else:
                        self.chain.append(None)
                        self.chainid.append(chainid)
        
        # sort_indexの更新
        for j, strandid in enumerate(self.struct_ref_seq["strand_id"]):
            self.struct_ref_seq.at[j, "sort_index"] = self.chainid.index(strandid)
    
    def _extract_atom_coord(self, mmcifdict):
        """原子座標を抽出してCSVに保存"""
        atom_coord = pd.DataFrame({
            "model_num": mmcifdict["_atom_site.pdbx_PDB_model_num"],
            "asym_id": mmcifdict["_atom_site.auth_asym_id"],
            "comp_id": mmcifdict["_atom_site.auth_comp_id"],
            "seq_id": mmcifdict["_atom_site.auth_seq_id"],
            "atom_id": mmcifdict["_atom_site.auth_atom_id"],
            "Cartn_x": mmcifdict["_atom_site.Cartn_x"],
            "Cartn_y": mmcifdict["_atom_site.Cartn_y"],
            "Cartn_z": mmcifdict["_atom_site.Cartn_z"],
            "alt_id": mmcifdict["_atom_site.label_alt_id"],
            "group_PDB": mmcifdict["_atom_site.group_PDB"],
            "ins_code": mmcifdict["_atom_site.pdbx_PDB_ins_code"]
        })
        
        atom_coord["asym_id"] = atom_coord["asym_id"].astype(str)
        
        # alt_idの処理
        atom_coord['original_index'] = atom_coord.index
        alt_id_dot = atom_coord[atom_coord['alt_id'].str.contains('\\.')]
        alt_id_not_dot = atom_coord[~atom_coord['alt_id'].str.contains('\\.')]
        alt_id_not_dot_unique = alt_id_not_dot.drop_duplicates(
            subset=['seq_id', 'atom_id']
        )
        
        atom_coord = pd.concat([alt_id_dot, alt_id_not_dot_unique])
        atom_coord = atom_coord.sort_values('original_index')
        atom_coord = atom_coord.drop(columns=['original_index'])
        atom_coord = atom_coord[
            (atom_coord['group_PDB'] == 'ATOM')
        ].drop(columns=['alt_id', 'group_PDB'])
        
        if not os.path.exists('atom_coord/'):
            os.makedirs('atom_coord/')
        
        atom_coord.to_csv(f'atom_coord/{self.pdbid}.csv', index=False)
    
    def mutationjudge(self, uniprotids: List[str], pdbid: str) -> str:
        """
        変異の判定
        
        Returns
        -------
        str
            "normal", "substitution", "chimera", "delins", または 
            "UniProt ID mismatch"
        """
        m_pd = self.struct_ref_seq[["strand_id", "accession"]]
        unim_pd = m_pd[m_pd["accession"].isin(uniprotids)]
        
        if unim_pd["accession"].count() == 0:
            return "UniProt ID mismatch"
        
        if unim_pd.duplicated().sum() != 0:
            return "chimera"
        
        m_id = list(unim_pd["strand_id"])
        mdif_pd = self.struct_ref_seq_dif[
            self.struct_ref_seq_dif["strand_id"].isin(m_id)
        ]
        
        if len(mdif_pd) == 0:
            return "normal"
        
        mdif_pd_details = mdif_pd["details"].unique()
        
        if "engineered mutation" in mdif_pd_details:
            return "substitution"
        elif "microheterogeneity" in mdif_pd_details:
            return "normal"
        
        s_list = list(m_pd["strand_id"])
        if len(s_list) != len(set(s_list)):
            return "chimera"
        
        for i in m_id:
            strand_mdif_pd = mdif_pd[mdif_pd["strand_id"] == i]
            seq_num_list = list(strand_mdif_pd["seq_num"])
            db_seq_num_list = list(strand_mdif_pd["db_seq_num"])
            
            if len(seq_num_list) != len(set(seq_num_list)):
                return "delins"
            elif len(db_seq_num_list) != len(set(db_seq_num_list)):
                return "delins"
        
        return "substitution"
    
    def getsequence(self, uniprotids: List[str]) -> pd.DataFrame:
        """
        配列情報を取得
        
        Parameters
        ----------
        uniprotids : list
            UniProt IDのリスト
        
        Returns
        -------
        pd.DataFrame
            配列データ
        """
        firstLoop = True
        struct = self.struct_ref_seq[
            self.struct_ref_seq['accession'].isin(uniprotids)
        ].drop_duplicates(subset=["strand_id"])
        
        for row in struct.itertuples():
            if row.accession in uniprotids:
                sort_index = int(row.sort_index)
                align_beg = sort_index + int(row.seq_align_beg) - 1
                align_end = sort_index + int(row.seq_align_end)
                chain = self.chain[align_beg:align_end]
                
                mutat_info = self.struct_ref_seq_dif[
                    self.struct_ref_seq_dif["strand_id"] == row.strand_id
                ].drop(columns='strand_id')
                
                if len(mutat_info) != 0:
                    # deletion処理
                    deletion = mutat_info[(mutat_info["seq_num"] == '?')].index
                    if len(deletion) != 0:
                        mutat_info.drop(deletion, inplace=True)
                        chain_num = pd.Series(chain).map(
                            lambda x: int(x.split(', ')[1]) if isinstance(x, str) else x
                        ).diff()
                        deletion = chain_num[(chain_num != 1)].dropna()
                        for index, i in zip(deletion.index, deletion):
                            chain[index:index] = [None] * int(i)
                    
                    # insertion処理
                    insertion = mutat_info[(mutat_info["db_seq_num"] == '?')]["seq_num"]
                    if len(insertion) != 0:
                        mutat_info.drop(insertion.index, inplace=True)
                        insertion = insertion.values.tolist()
                        for i in chain:
                            if isinstance(i, str):
                                for n in insertion:
                                    if n == i.split(', ')[1]:
                                        insertion.remove(n)
                                        chain.remove(i)
                    
                    # delins処理
                    dup_mutat = mutat_info[
                        mutat_info.duplicated(subset=["seq_num"], keep=False)
                    ]
                    if len(dup_mutat) != 0:
                        mutat_info.drop(dup_mutat.index, inplace=True)
                        for i in dup_mutat['seq_num'].drop_duplicates():
                            chain_num = pd.Series(chain).map(
                                lambda x: int(x.split(', ')[1]) if isinstance(x, str) else x
                            )
                            num = len(dup_mutat[dup_mutat['seq_num'] == i]) - 1
                            index = chain_num[chain_num == int(i)].index[0] + 1
                            chain[index:index] = [None] * num
                    
                    # ins処理
                    dup_mutat = mutat_info[
                        mutat_info.duplicated(subset=["db_seq_num"], keep=False)
                    ]
                    if len(dup_mutat) != 0:
                        mutat_info.drop(dup_mutat.index, inplace=True)
                        insertion = []
                        for i in dup_mutat['db_seq_num'].drop_duplicates():
                            insertion += dup_mutat[
                                dup_mutat['db_seq_num'] == i
                            ]["seq_num"].reset_index(drop=True).drop([0]).values.tolist()
                        
                        m = 0
                        for i in range(len(chain)):
                            i = chain[i + m]
                            if isinstance(i, str):
                                for n in insertion:
                                    if n == i.split(', ')[1]:
                                        insertion.remove(n)
                                        chain.remove(i)
                                        m -= 1
                
                if firstLoop:
                    firstLoop = False
                    sequence = pd.DataFrame(
                        chain, 
                        columns=[self.pdbid + ' ' + row.strand_id]
                    )
                else:
                    strand = pd.Series(
                        chain, 
                        name=self.pdbid + ' ' + row.strand_id
                    )
                    sequence = pd.concat([sequence, strand], axis=1)
        
        if firstLoop:
            sequence = pd.DataFrame()
        
        return sequence