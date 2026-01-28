"""
PDB/mmCIF構造ファイルの解析モジュール (修正版)
"""

import os
import gzip
import time
import threading
from mimetypes import guess_type
import pandas as pd
import requests
from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from typing import List


_pdb_list = None
def _get_pdb_list():
    global _pdb_list
    if _pdb_list is None:
        _pdb_list = PDBList()
    return _pdb_list

# レート制限用（PDBダウンロード間隔制御）
_pdb_download_lock = threading.Lock()
_last_pdb_download_time = 0.0
_DEFAULT_PDB_DELAY_SEC = 2.0  # デフォルト2秒間隔
_DEFAULT_PDB_TIMEOUT_SEC = 120.0  # タイムアウト（秒）。大きいCIF用
_DEFAULT_PDB_RETRIES = 3  # リトライ回数
_RCSB_CIF_URL = "https://files.rcsb.org/view/{pdbid}.cif"  # HTTPSで安定


def _download_via_requests(pdbid: str, timeout_sec: float, out_path: str) -> bool:
    """requests で RCSB HTTPS から CIF を取得。成功時 True。"""
    url = _RCSB_CIF_URL.format(pdbid=pdbid.lower())
    try:
        r = requests.get(url, timeout=timeout_sec, stream=True)
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        return True
    except Exception:
        return False


def _downloadpdb_with_retry(pdbid: str, cif_file: str, delay_sec: float) -> bool:
    """リトライ・タイムアウト付きで CIF を取得。HTTPS 優先、失敗時 PDBList にフォールバック。"""
    timeout_sec = float(os.environ.get("DSA_PDB_TIMEOUT_SEC", _DEFAULT_PDB_TIMEOUT_SEC))
    max_retries = int(os.environ.get("DSA_PDB_RETRIES", _DEFAULT_PDB_RETRIES))
    os.makedirs("pdb_files", exist_ok=True)

    for attempt in range(max_retries):
        if _download_via_requests(pdbid, timeout_sec, cif_file):
            return True
        backoff = min(5 * (2 ** attempt), 60)
        if attempt < max_retries - 1:
            time.sleep(backoff)

    try:
        _get_pdb_list().retrieve_pdb_file(pdbid, pdir="pdb_files/", file_format="mmCif", overwrite=False)
        return True
    except Exception as e:
        print(f"Desired structure not found or download failed. '{pdbid}': {e}")
        return False


def downloadpdb(pdbid: str, delay_sec: float = None) -> bool:
    """
    PDB CIF をダウンロード（レート制限・リトライ・タイムアウト対策付き）
    
    - RCSB HTTPS (files.rcsb.org) を優先。失敗時は PDBList にフォールバック。
    - リトライ: 環境変数 DSA_PDB_RETRIES またはデフォルト 3 回、指数バックオフ。
    - タイムアウト: 環境変数 DSA_PDB_TIMEOUT_SEC またはデフォルト 120 秒。
    - 間隔: 環境変数 DSA_PDB_DELAY_SEC またはデフォルト 2 秒。
    """
    global _last_pdb_download_time
    
    if delay_sec is None:
        delay_sec = float(os.environ.get("DSA_PDB_DELAY_SEC", _DEFAULT_PDB_DELAY_SEC))
    
    cif_file = f"pdb_files/{pdbid.lower()}.cif"
    cache_hit = os.path.exists(cif_file)
    
    # 空ファイル（0バイト）を検出して再ダウンロード
    if cache_hit:
        try:
            if os.path.getsize(cif_file) == 0:
                # 空ファイルを削除して再ダウンロード
                os.remove(cif_file)
                cache_hit = False
        except OSError:
            # ファイルサイズ取得に失敗した場合は再ダウンロード
            cache_hit = False
    
    if not cache_hit and delay_sec > 0:
        with _pdb_download_lock:
            elapsed = time.time() - _last_pdb_download_time
            if elapsed < delay_sec:
                time.sleep(delay_sec - elapsed)
            _last_pdb_download_time = time.time()
    
    if not cache_hit:
        ok = _downloadpdb_with_retry(pdbid, cif_file, delay_sec)
        if ok and delay_sec > 0:
            with _pdb_download_lock:
                _last_pdb_download_time = time.time()
        return ok
    return True


def _open(pdbid: str):
    """CIFファイルを開く(gzip対応)"""
    file = pdbid.lower() + ".cif"
    ciffile = "pdb_files/" + file
    
    if not os.path.exists(ciffile):
        raise FileNotFoundError(f"No such file or directory: '{ciffile}'")
    
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
    
    def __init__(self, pdbid: str, skip_download: bool = False):
        """
        Parameters
        ----------
        pdbid : str
            PDB ID
        skip_download : bool
            Trueの場合、PDBダウンロードをスキップ（既にダウンロード済みの場合）
        """
        self.pdbid = pdbid
        if not skip_download:
            downloadpdb(self.pdbid)
        
        try:
            with _open(self.pdbid) as handle:
                mmcifdict = MMCIF2Dict(handle)
        except Exception as e:
            print(f"Error reading CIF file for {pdbid}: {e}")
            raise
        
        # struct_ref_seqの構築
        try:
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
        except KeyError as e:
            print(f"Warning: Missing struct_ref_seq data in {pdbid}: {e}")
            self.struct_ref_seq = pd.DataFrame(columns=["strand_id", "accession", "seq_align_beg", "seq_align_end"])
        
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
        
        # 原子座標の取得（atom_coord保存は不要になったためスキップ）
        # self._extract_atom_coord(mmcifdict)
    
    def _build_chain_info(self, mmcifdict):
        """Chain情報を構築"""
        self.chain = []
        self.chainid = []
        self.hetero_info = []
        self.ind = -1
        hetero_pdb_seq_num = ""
        
        try:
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
                if strandid in self.chainid:
                    self.struct_ref_seq.at[j, "sort_index"] = self.chainid.index(strandid)
        except KeyError as e:
            print(f"Warning: Missing chain info in {self.pdbid}: {e}")
    
    def _extract_atom_coord(self, mmcifdict):
        """原子座標を抽出してCSVに保存 (修正版: データ長チェック追加)"""
        try:
            # 各フィールドを取得
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
            ins_code = mmcifdict.get("_atom_site.pdbx_PDB_ins_code", [])
            
            # データ長チェック
            lengths = [len(model_num), len(asym_id), len(comp_id), len(seq_id),
                      len(atom_id), len(cartn_x), len(cartn_y), len(cartn_z),
                      len(alt_id), len(group_PDB), len(ins_code)]
            
            if len(set(lengths)) > 1:
                print(f"Warning: Inconsistent data lengths in {self.pdbid}: {lengths}")
                min_length = min(lengths)
                print(f"  Truncating to minimum length: {min_length}")
                
                # 最小の長さに合わせる
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
                ins_code = ins_code[:min_length]
            
            # DataFrameを作成
            atom_coord = pd.DataFrame({
                "model_num": model_num,
                "asym_id": asym_id,
                "comp_id": comp_id,
                "seq_id": seq_id,
                "atom_id": atom_id,
                "Cartn_x": cartn_x,
                "Cartn_y": cartn_y,
                "Cartn_z": cartn_z,
                "alt_id": alt_id,
                "group_PDB": group_PDB,
                "ins_code": ins_code
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
            print(f"  Saved {len(atom_coord)} atoms to atom_coord/{self.pdbid}.csv")
            
        except Exception as e:
            print(f"Error extracting atom coordinates for {self.pdbid}: {e}")
            import traceback
            traceback.print_exc()
            raise
    
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
                        insertion_list = insertion.values.tolist()
                        
                        # 削除する要素のインデックスを収集
                        indices_to_remove = []
                        for idx, elem in enumerate(chain):
                            if isinstance(elem, str):
                                seq_num = elem.split(', ')[1]
                                if seq_num in insertion_list:
                                    indices_to_remove.append(idx)
                                    try:
                                        insertion_list.remove(seq_num)
                                    except ValueError:
                                        pass
                
                        # 後ろから削除(インデックスのずれを防ぐ)
                        for idx in sorted(indices_to_remove, reverse=True):
                            del chain[idx]
                
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
                        insertion_list = []
                        for i in dup_mutat['db_seq_num'].drop_duplicates():
                            insertion_list += dup_mutat[
                                dup_mutat['db_seq_num'] == i
                            ]["seq_num"].reset_index(drop=True).drop([0]).values.tolist()
                        
                        # 削除する要素のインデックスを収集
                        indices_to_remove = []
                        for idx, elem in enumerate(chain):
                            if isinstance(elem, str):
                                seq_num = elem.split(', ')[1]
                                if seq_num in insertion_list:
                                    indices_to_remove.append(idx)
                                    try:
                                        insertion_list.remove(seq_num)
                                    except ValueError:
                                        pass
                        
                        # 後ろから削除(インデックスのずれを防ぐ)
                        for idx in sorted(indices_to_remove, reverse=True):
                            del chain[idx]
                
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