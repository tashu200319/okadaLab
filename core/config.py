"""
DSA プロジェクト設定ファイル（拡張版）
UniProt ID管理機能を追加
"""

import os
import pandas as pd
from typing import List, Set, Optional
from pathlib import Path


class Config:
    """解析パラメータの設定クラス"""
    
    def __init__(self, **kwargs):
        """
        設定の初期化
        
        Parameters
        ----------
        **kwargs : dict
            カスタム設定（例: SEQ_RATIO=70, USE_EM=True）
        """
        # 構造決定法の選択
        self.USE_XRAY = True
        self.USE_NMR = False
        self.USE_EM = True
        
        # 検索モード追加
        self.SEARCH_MODE = "OR"  # "AND" or "OR"
        
        # 解析閾値
        self.PDB_THRESHOLD = 1          # 最小PDB数
        self.CHAIN_THRESHOLD = 3        # 最小Chain数（標準偏差計算のため）
        
        # 配列解析パラメータ
        self.SEQ_RATIO = 20             # 使用する配列長の割合(%)
        
        # cis解析パラメータ
        self.CIS_THRESHOLD = 3.3        # cisプロリン結合のCA-CA間距離(Å)
        self.PROC_CIS = True            # cis解析を実行するか
        
        # 除外するPDB ID（カンマまたはスペース区切り）
        self.NEGATIVE_PDBID = ""
        
        # 入力設定
        # 入力設定（標準: data/chunks/。無ければ従来の chunks/ にフォールバック）
        if Path("./data/chunks/chunk_1.csv").exists():
            self.INPUT_FILE = "./data/chunks/chunk_1.csv"
        else:
            self.INPUT_FILE = "./chunks/chunk_1.csv"  # デフォルトの入力CSVファイル（旧）
        
        # 出力設定
        self.OUTPUT_DIR = "./output/"
        self.EXPORT = True              # CSV出力するか
        self.HEATMAP = True             # ヒートマップ生成するか
        self.OVERWRITE = True           # データ上書きするか
        self.VERBOSE = True             # 処理状況を表示するか
        
        # カスタム設定で上書き
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: Unknown config parameter '{key}' ignored")
    
    @property
    def METHODS_SELECTED(self):
        """選択された構造決定法を返す"""
        methods = []
        if self.USE_XRAY:
            methods.append("X-ray")
        if self.USE_NMR:
            methods.append("NMR")
        if self.USE_EM:
            methods.append("EM")
        
        # 何も選択されていない場合は全種類
        if not methods:
            methods = ["X-ray", "NMR", "EM"]
        
        return set(methods)
    
    def display(self):
        """現在の設定を表示"""
        print("=" * 60)
        print("DSA Configuration")
        print("=" * 60)
        print(f"構造決定法:")
        print(f"  X-ray: {self.USE_XRAY}")
        print(f"  NMR: {self.USE_NMR}")
        print(f"  EM: {self.USE_EM}")
        print(f"  選択: {self.METHODS_SELECTED}")
        print(f"\n閾値:")
        print(f"  PDB_THRESHOLD: {self.PDB_THRESHOLD}")
        print(f"  CHAIN_THRESHOLD: {self.CHAIN_THRESHOLD}")
        print(f"\n配列解析:")
        print(f"  SEQ_RATIO: {self.SEQ_RATIO}%")
        print(f"\ncis解析:")
        print(f"  CIS_THRESHOLD: {self.CIS_THRESHOLD} Å")
        print(f"  PROC_CIS: {self.PROC_CIS}")
        print(f"\n除外PDB: {self.NEGATIVE_PDBID if self.NEGATIVE_PDBID else 'なし'}")
        print(f"\n出力設定:")
        print(f"  OUTPUT_DIR: {self.OUTPUT_DIR}")
        print(f"  EXPORT: {self.EXPORT}")
        print(f"  HEATMAP: {self.HEATMAP}")
        print(f"  OVERWRITE: {self.OVERWRITE}")
        print(f"  VERBOSE: {self.VERBOSE}")
        print("=" * 60)
    
    # ===== 🆕 UniProt ID管理機能 =====
    
    def load_uniprot_ids(self, 
                        filename: str = "has_both_ids_list.txt",
                        fallback_csv: str = "unique_uniprots.csv") -> List[str]:
        """UniProt IDリストを読み込む（最小修正版）"""
        output_dir = Path(self.OUTPUT_DIR)
        
        # 1. 抽出した精鋭リストを最優先でチェック（追加箇所）
        target_csv = output_dir / "only_true_results.csv"
        if target_csv.exists():
            return self._load_uniprot_from_csv(target_csv)
        
        # 2. 既存の処理（引数エラーを回避するために維持）
        main_file = output_dir / filename
        if main_file.exists():
            ids = self._load_uniprot_from_txt(main_file)
            if ids: return ids
        
        fallback_file = output_dir / "links" / fallback_csv
        if fallback_file.exists():
            ids = self._load_uniprot_from_csv(fallback_file)
            if ids: return ids
            
        return []
    
    def _load_uniprot_from_txt(self, filepath: Path) -> List[str]:
        """
        Pythonリスト形式またはテキスト形式のファイルから読み込み
        
        Format 1 (Pythonリスト):
        ```
        uniprot_ids = [
            "P62258", "Q13144", ...
        ]
        ```
        
        Format 2 (プレーンテキスト):
        ```
        P62258
        Q13144
        ```
        """
        try:
            with open(filepath, 'r') as f:
                content = f.read()
            
            # Format 1: uniprot_ids = [...] の形式
            if 'uniprot_ids = [' in content or 'uniprot_ids=[' in content:
                start = content.find('[')
                end = content.rfind(']')
                if start != -1 and end != -1:
                    list_str = content[start:end+1]
                    
                    import ast
                    ids = ast.literal_eval(list_str)
                    
                    if self.VERBOSE:
                        print(f"📋 Loaded {len(ids)} UniProt IDs from {filepath.name}")
                    return ids
            
            # Format 2: 通常のテキストファイル（1行1ID）
            ids = []
            for line in content.split('\n'):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                for part in line.split(','):
                    uid = part.strip().strip('"').strip("'")
                    if uid and not uid.startswith('#'):
                        ids.append(uid)
            
            if ids and self.VERBOSE:
                print(f"📋 Loaded {len(ids)} UniProt IDs from {filepath.name}")
            
            return ids
                
        except Exception as e:
            if self.VERBOSE:
                print(f"❌ Error loading {filepath}: {e}")
            return []
    
    def _load_uniprot_from_csv(self, filepath: Path) -> List[str]:
        """CSVファイルからUniProt IDを読み込み"""
        try:
            df = pd.read_csv(filepath)
            
            # uniprotid列を探す
            possible_columns = ['uniprotid', 'UniProtID', 'uniprot_id', 'UNIPROTID']
            column_name = None
            
            for col in possible_columns:
                if col in df.columns:
                    column_name = col
                    break
            
            if column_name:
                ids = df[column_name].dropna().unique().tolist()
            elif len(df.columns) == 1:
                ids = df.iloc[:, 0].dropna().unique().tolist()
            else:
                if self.VERBOSE:
                    print(f"⚠️  Warning: Could not find uniprotid column in {filepath}")
                return []
            
            if self.VERBOSE:
                print(f"📋 Loaded {len(ids)} UniProt IDs from {filepath.name}")
            
            return ids
            
        except Exception as e:
            if self.VERBOSE:
                print(f"❌ Error loading {filepath}: {e}")
            return []
    
    def load_excluded_ids(self, 
                         filename: str = "excluded_ids.txt") -> Set[str]:
        """
        除外IDリストを読み込む
        
        Parameters
        ----------
        filename : str
            除外IDファイル名
        
        Returns
        -------
        Set[str]
            除外するUniProt IDのセット
        """
        # config/excluded_ids.txt を優先
        config_dir = Path("./config")
        config_file = config_dir / filename
        if config_file.exists():
            return self._load_excluded_from_file(config_file)
        
        # output/excluded_ids.txt をチェック
        output_file = Path(self.OUTPUT_DIR) / filename
        if output_file.exists():
            return self._load_excluded_from_file(output_file)
        
        # なければ空セットを返す（エラーではない）
        if self.VERBOSE:
            print(f"ℹ️  No excluded IDs file found (optional)")
        
        return set()
    
    def _load_excluded_from_file(self, filepath: Path) -> Set[str]:
        """除外IDファイルを読み込み"""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            ids = set()
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    for uid in line.replace(',', ' ').split():
                        uid = uid.strip().strip('"').strip("'")
                        if uid:
                            ids.add(uid)
            
            if self.VERBOSE:
                print(f"🚫 Loaded {len(ids)} excluded IDs from {filepath.name}")
            
            return ids
            
        except Exception as e:
            if self.VERBOSE:
                print(f"❌ Error loading {filepath}: {e}")
            return set()
    
    def add_to_excluded_ids(self, 
                           new_ids: Set[str],
                           filename: str = "excluded_ids.txt",
                           max_retries: int = 3) -> int:
        """
        除外IDリストに新しいIDを追加
        
        Parameters
        ----------
        new_ids : Set[str]
            追加するUniProt IDのセット
        filename : str
            除外IDファイル名
        max_retries : int
            この回数以上失敗したIDのみ追加
        
        Returns
        -------
        int
            追加されたID数
        """
        if not new_ids:
            return 0
        
        # output/excluded_ids.txt を使用
        output_file = Path(self.OUTPUT_DIR) / filename
        
        # 既存の除外IDを読み込み
        existing_ids = set()
        if output_file.exists():
            existing_ids = self._load_excluded_from_file(output_file)
        
        # 新規IDのみ抽出
        ids_to_add = new_ids - existing_ids
        
        if not ids_to_add:
            if self.VERBOSE:
                print(f"ℹ️  No new IDs to add to excluded list")
            return 0
        
        try:
            # ファイルに追記
            with open(output_file, 'a') as f:
                if not output_file.exists() or output_file.stat().st_size == 0:
                    # 新規ファイルの場合はヘッダーを追加
                    f.write("# Excluded UniProt IDs\n")
                    f.write("# Auto-generated from failed analysis\n\n")
                else:
                    # 既存ファイルに追記する場合は区切り追加
                    f.write(f"\n# Added from analysis (max_retries={max_retries})\n")
                
                for uid in sorted(ids_to_add):
                    f.write(f"{uid}\n")
            
            if self.VERBOSE:
                print(f"➕ Added {len(ids_to_add)} IDs to {output_file}")
            
            return len(ids_to_add)
            
        except Exception as e:
            if self.VERBOSE:
                print(f"❌ Error saving to {output_file}: {e}")
            return 0
    
    def update_excluded_from_failed(self,
                                    failed_ids_file: str = None,
                                    seq_ratio: float = 20,
                                    max_retries: int = 3) -> int:
        """
        failed_ids.csvから除外IDリストを自動更新
        
        Parameters
        ----------
        failed_ids_file : str, optional
            failed_ids.csvのパス（Noneの場合はデフォルト）
        seq_ratio : float
            対象のseq_ratio
        max_retries : int
            この回数以上失敗したIDを追加
        
        Returns
        -------
        int
            追加されたID数
        """
        if failed_ids_file is None:
            failed_ids_file = os.path.join(self.OUTPUT_DIR, "summaries", "failed_ids.csv")
        
        if not os.path.exists(failed_ids_file):
            if self.VERBOSE:
                print(f"ℹ️  No failed_ids.csv found")
            return 0
        
        try:
            # failed_ids.csvを読み込み
            df = pd.read_csv(failed_ids_file)
            
            # 型変換（文字列で読み込まれた場合の対策）
            if 'seq_ratio' in df.columns:
                df['seq_ratio'] = pd.to_numeric(df['seq_ratio'], errors='coerce')
            if 'retry_count' in df.columns:
                df['retry_count'] = pd.to_numeric(df['retry_count'], errors='coerce').fillna(0).astype(int)
            
            # 指定seq_ratioで最大リトライ回数を超えたIDを抽出
            mask = (df['seq_ratio'] == seq_ratio) & (df['retry_count'] >= max_retries)
            failed_ids = set(df[mask]['uniprotid'].tolist())
            
            if not failed_ids:
                if self.VERBOSE:
                    print(f"ℹ️  No IDs with {max_retries}+ failures for seq_ratio={seq_ratio}")
                return 0
            
            # 除外リストに追加
            count = self.add_to_excluded_ids(failed_ids, max_retries=max_retries)
            
            return count
            
        except Exception as e:
            if self.VERBOSE:
                print(f"❌ Error updating from failed_ids.csv: {e}")
            return 0
    
    def get_filtered_ids(self, 
                        uniprot_ids: List[str] = None,
                        load_excluded: bool = True) -> List[str]:
        """
        フィルタリング済みのUniProt IDリストを取得
        
        Parameters
        ----------
        uniprot_ids : List[str], optional
            元のIDリスト（Noneの場合は自動読み込み）
        load_excluded : bool
            除外IDリストを適用するか
        
        Returns
        -------
        List[str]
            フィルタリング後のIDリスト
        """
        # IDリスト読み込み
        if uniprot_ids is None:
            uniprot_ids = self.load_uniprot_ids()
        
        if not uniprot_ids:
            return []
        
        # 除外IDを適用
        if load_excluded:
            excluded_ids = self.load_excluded_ids()
            if excluded_ids:
                original_count = len(uniprot_ids)
                uniprot_ids = [uid for uid in uniprot_ids if uid not in excluded_ids]
                
                if self.VERBOSE:
                    removed = original_count - len(uniprot_ids)
                    if removed > 0:
                        print(f"🚫 Filtered out {removed} excluded IDs")
        
        return uniprot_ids