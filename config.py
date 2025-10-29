"""
DSA プロジェクト設定ファイル
"""

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