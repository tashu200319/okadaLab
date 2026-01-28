#!/usr/bin/env python3
"""
失敗したUniProt IDを自動記録する機能
tools/ ディレクトリ用
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd
import datetime
import pytz
from pathlib import Path
from typing import Set, Optional


class FailedIDManager:
    """失敗したUniProt IDを管理するクラス"""
    
    def __init__(self, failed_ids_file: Optional[str] = None):
        """
        Parameters
        ----------
        failed_ids_file : str, optional
            失敗ID記録ファイルのパス（Noneの場合は自動検出）
        """
        try:
            if failed_ids_file is None:
                # プロジェクトルートを自動検出
                try:
                    current_dir = Path(__file__).resolve().parent.parent
                    failed_ids_file = str(current_dir / "output" / "summaries" / "failed_ids.csv")
                except Exception as e:
                    # フォールバック: 相対パスを使用
                    failed_ids_file = "output/summaries/failed_ids.csv"
            self.failed_ids_file = failed_ids_file
            self.failed_ids_df = self._load_failed_ids()
        except Exception as e:
            # 初期化に失敗した場合は空のDataFrameを使用
            print(f"⚠️  Warning: FailedIDManager initialization failed: {e}", flush=True)
            import traceback
            traceback.print_exc()
            self.failed_ids_file = failed_ids_file or "output/summaries/failed_ids.csv"
            self.failed_ids_df = pd.DataFrame(columns=['uniprotid', 'seq_ratio', 'error_type', 
                                                        'error_message', 'timestamp', 'retry_count'])
    
    def _load_failed_ids(self) -> pd.DataFrame:
        """既存の失敗ID記録を読み込む"""
        if not os.path.exists(self.failed_ids_file):
            return pd.DataFrame(columns=['uniprotid', 'seq_ratio', 'error_type', 
                                        'error_message', 'timestamp', 'retry_count'])
        
        try:
            # ファイルが空でないか確認
            if os.path.getsize(self.failed_ids_file) == 0:
                return pd.DataFrame(columns=['uniprotid', 'seq_ratio', 'error_type', 
                                            'error_message', 'timestamp', 'retry_count'])
            
            df = pd.read_csv(self.failed_ids_file)
            
            # 空のDataFrameの場合は空のDataFrameを返す
            if len(df) == 0:
                return pd.DataFrame(columns=['uniprotid', 'seq_ratio', 'error_type', 
                                            'error_message', 'timestamp', 'retry_count'])
            
            # 必要なカラムが存在するか確認
            required_columns = ['uniprotid', 'seq_ratio', 'error_type', 
                              'error_message', 'timestamp', 'retry_count']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                print(f"⚠️  Warning: Missing columns in failed_ids.csv: {missing_columns}")
                return pd.DataFrame(columns=required_columns)
            
            # seq_ratioをfloat型に変換（文字列で読み込まれた場合の対策）
            if 'seq_ratio' in df.columns:
                df['seq_ratio'] = pd.to_numeric(df['seq_ratio'], errors='coerce')
            # retry_countをint型に変換（文字列で読み込まれた場合の対策）
            if 'retry_count' in df.columns:
                df['retry_count'] = pd.to_numeric(df['retry_count'], errors='coerce').fillna(0).astype(int)
            
            print(f"📋 Loaded {len(df)} failed ID records")
            return df
        except Exception as e:
            print(f"⚠️  Warning: Could not load failed IDs: {e}", flush=True)
            import traceback
            traceback.print_exc()
            return pd.DataFrame(columns=['uniprotid', 'seq_ratio', 'error_type', 
                                        'error_message', 'timestamp', 'retry_count'])
    
    def get_failed_ids(self, seq_ratio: float, max_retries: int = 3) -> Set[str]:
        """
        指定したseq_ratioで失敗したIDのセットを取得
        
        Parameters
        ----------
        seq_ratio : float
            seq_ratio値
        max_retries : int
            最大リトライ回数（一時的なエラーの場合のみ適用）
        
        Returns
        -------
        set
            スキップすべきUniProt IDのセット
        """
        if len(self.failed_ids_df) == 0:
            return set()
        
        # 確実な失敗理由（1回でスキップ）
        CERTAIN_ERROR_TYPES = {
            'EMPTY_TRIMSEQUENCE', 'NOT_ENOUGH_CHAINS', 'EMPTY_SUMMARY',
            'PDB_THRESHOLD', 'EMPTY_SEQDATA', 'EMPTY_ATOMCOORD', 'EMPTY_DISTANCE'
        }
        
        # 一時的なエラー（max_retries回まで再試行）
        TEMPORARY_ERROR_TYPES = {'TIMEOUT', 'DOWNLOAD_ERROR', 'OTHER_ERROR'}
        
        # 指定seq_ratioの失敗IDを取得
        seq_mask = self.failed_ids_df['seq_ratio'] == seq_ratio
        failed_df = self.failed_ids_df[seq_mask]
        
        if len(failed_df) == 0:
            return set()
        
        failed_ids = set()
        
        # 確実なエラー: 1回失敗したら即スキップ
        certain_mask = failed_df['error_type'].isin(CERTAIN_ERROR_TYPES) & \
                      (failed_df['retry_count'] >= 1)
        certain_ids = set(failed_df[certain_mask]['uniprotid'].tolist())
        failed_ids.update(certain_ids)
        
        # 一時的なエラー: max_retries回以上失敗したらスキップ
        temporary_mask = failed_df['error_type'].isin(TEMPORARY_ERROR_TYPES) & \
                        (failed_df['retry_count'] >= max_retries)
        temporary_ids = set(failed_df[temporary_mask]['uniprotid'].tolist())
        failed_ids.update(temporary_ids)
        
        if failed_ids:
            certain_count = len(certain_ids)
            temporary_count = len(temporary_ids)
            if certain_count > 0 and temporary_count > 0:
                print(f"🚫 Skipping {len(failed_ids)} IDs ({certain_count} certain errors, {temporary_count} temporary errors)")
            elif certain_count > 0:
                print(f"🚫 Skipping {len(failed_ids)} IDs ({certain_count} certain errors)")
            else:
                print(f"🚫 Skipping {len(failed_ids)} IDs ({temporary_count} failed {max_retries}+ times)")
        
        return failed_ids
    
    def record_failure(self, uniprotid: str, seq_ratio: float, 
                      error_type: str, error_message: str):
        """
        失敗を記録
        
        Parameters
        ----------
        uniprotid : str
            UniProt ID
        seq_ratio : float
            seq_ratio値
        error_type : str
            エラータイプ（例: "PDB_THRESHOLD", "EMPTY_TRIMSEQUENCE", "NOT_ENOUGH_CHAINS"）
        error_message : str
            詳細エラーメッセージ
        """
        jst = pytz.timezone('Asia/Tokyo')
        timestamp = datetime.datetime.now(jst).strftime('%Y-%m-%d %H:%M:%S')
        
        # 既存の記録があるかチェック
        mask = (self.failed_ids_df['uniprotid'] == uniprotid) & \
               (self.failed_ids_df['seq_ratio'] == seq_ratio) & \
               (self.failed_ids_df['error_type'] == error_type)
        
        if mask.any():
            # retry_countをインクリメント
            self.failed_ids_df.loc[mask, 'retry_count'] += 1
            self.failed_ids_df.loc[mask, 'timestamp'] = timestamp
            self.failed_ids_df.loc[mask, 'error_message'] = error_message
        else:
            # 新規記録
            new_record = pd.DataFrame({
                'uniprotid': [uniprotid],
                'seq_ratio': [seq_ratio],
                'error_type': [error_type],
                'error_message': [error_message],
                'timestamp': [timestamp],
                'retry_count': [1]
            })
            self.failed_ids_df = pd.concat([self.failed_ids_df, new_record], 
                                          ignore_index=True)
    
    def save(self):
        """失敗記録をファイルに保存"""
        if len(self.failed_ids_df) > 0:
            # ディレクトリ作成
            os.makedirs(os.path.dirname(self.failed_ids_file), exist_ok=True)
            
            # retry_countで降順ソート（失敗回数が多い順）
            self.failed_ids_df = self.failed_ids_df.sort_values(
                by=['retry_count', 'uniprotid'], 
                ascending=[False, True]
            )
            
            self.failed_ids_df.to_csv(self.failed_ids_file, index=False)
            print(f"💾 Saved {len(self.failed_ids_df)} failed records to {self.failed_ids_file}")
    
    def get_statistics(self, seq_ratio: float) -> dict:
        """
        失敗統計を取得
        
        Parameters
        ----------
        seq_ratio : float
            seq_ratio値
        
        Returns
        -------
        dict
            統計情報
        """
        if len(self.failed_ids_df) == 0:
            return {'total': 0, 'by_error_type': {}, 'retry_distribution': {}}
        
        mask = self.failed_ids_df['seq_ratio'] == seq_ratio
        filtered_df = self.failed_ids_df[mask]
        
        stats = {
            'total': len(filtered_df),
            'by_error_type': filtered_df['error_type'].value_counts().to_dict() if len(filtered_df) > 0 else {},
            'retry_distribution': filtered_df['retry_count'].value_counts().to_dict() if len(filtered_df) > 0 else {}
        }
        
        return stats
    
    def clear_failed_ids(self, seq_ratio: float = None, uniprotid: str = None):
        """
        失敗記録をクリア
        
        Parameters
        ----------
        seq_ratio : float, optional
            特定のseq_ratioの記録のみクリア
        uniprotid : str, optional
            特定のUniProt IDの記録のみクリア
        """
        if seq_ratio is not None and uniprotid is not None:
            # 特定のseq_ratioとuniprotidの組み合わせのみクリア
            mask = (self.failed_ids_df['seq_ratio'] == seq_ratio) & \
                   (self.failed_ids_df['uniprotid'] == uniprotid)
            self.failed_ids_df = self.failed_ids_df[~mask]
            print(f"🗑️  Cleared record for {uniprotid} (seq_ratio={seq_ratio})")
        elif seq_ratio is not None:
            # 特定のseq_ratioの記録のみクリア
            mask = self.failed_ids_df['seq_ratio'] == seq_ratio
            count = mask.sum()
            self.failed_ids_df = self.failed_ids_df[~mask]
            print(f"🗑️  Cleared {count} records for seq_ratio={seq_ratio}")
        elif uniprotid is not None:
            # 特定のuniprotidの記録のみクリア
            mask = self.failed_ids_df['uniprotid'] == uniprotid
            count = mask.sum()
            self.failed_ids_df = self.failed_ids_df[~mask]
            print(f"🗑️  Cleared {count} records for {uniprotid}")
        else:
            # 全てクリア
            count = len(self.failed_ids_df)
            self.failed_ids_df = pd.DataFrame(columns=['uniprotid', 'seq_ratio', 'error_type', 
                                                       'error_message', 'timestamp', 'retry_count'])
            print(f"🗑️  Cleared all {count} records")


def classify_error_type(error_message: str) -> str:
    """
    エラーメッセージからエラータイプを分類
    
    Parameters
    ----------
    error_message : str
        エラーメッセージ
    
    Returns
    -------
    str
        エラータイプ
    """
    error_lower = error_message.lower()
    
    if 'less than threshold' in error_lower or 'pdb_threshold' in error_lower:
        return 'PDB_THRESHOLD'
    elif 'trimsequence is empty' in error_lower or 'trimsequence' in error_lower:
        return 'EMPTY_TRIMSEQUENCE'
    elif 'not enough chains' in error_lower or 'chain_threshold' in error_lower:
        return 'NOT_ENOUGH_CHAINS'
    elif 'seqdata is empty' in error_lower:
        return 'EMPTY_SEQDATA'
    elif 'atomcoord is empty' in error_lower:
        return 'EMPTY_ATOMCOORD'
    elif 'distance is empty' in error_lower:
        return 'EMPTY_DISTANCE'
    elif 'df_all is empty' in error_lower:
        return 'EMPTY_SUMMARY'
    elif 'timeout' in error_lower:
        return 'TIMEOUT'
    elif 'download' in error_lower or 'fetch' in error_lower:
        return 'DOWNLOAD_ERROR'
    else:
        return 'OTHER_ERROR'


def main():
    """ツールとして実行時のメイン処理"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Failed UniProt ID Manager Tool'
    )
    parser.add_argument(
        '--seq-ratio',
        type=int,
        default=20,
        help='SEQ_RATIO value (default: 20)'
    )
    parser.add_argument(
        '--stats',
        action='store_true',
        help='Show statistics'
    )
    parser.add_argument(
        '--clear',
        action='store_true',
        help='Clear failed records'
    )
    parser.add_argument(
        '--clear-id',
        type=str,
        help='Clear records for specific UniProt ID'
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='List all failed IDs'
    )
    
    args = parser.parse_args()
    
    manager = FailedIDManager()
    
    if args.stats:
        # 統計表示
        stats = manager.get_statistics(args.seq_ratio)
        print("\n" + "=" * 80)
        print(f"📊 Failed ID Statistics (seq_ratio={args.seq_ratio})")
        print("=" * 80)
        print(f"Total failed: {stats['total']}")
        
        if stats['by_error_type']:
            print("\nBy error type:")
            for error_type, count in sorted(stats['by_error_type'].items(), 
                                           key=lambda x: x[1], reverse=True):
                print(f"  {error_type}: {count}")
        
        if stats['retry_distribution']:
            print("\nRetry distribution:")
            for retry_count, count in sorted(stats['retry_distribution'].items()):
                print(f"  {retry_count} attempt(s): {count} IDs")
        print("=" * 80)
    
    elif args.clear:
        # クリア
        if args.clear_id:
            manager.clear_failed_ids(seq_ratio=args.seq_ratio, uniprotid=args.clear_id)
        else:
            confirm = input(f"⚠️  Clear all records for seq_ratio={args.seq_ratio}? (y/n): ")
            if confirm.lower() == 'y':
                manager.clear_failed_ids(seq_ratio=args.seq_ratio)
        manager.save()
    
    elif args.list:
        # リスト表示
        if len(manager.failed_ids_df) == 0:
            print("No failed IDs recorded")
        else:
            mask = manager.failed_ids_df['seq_ratio'] == args.seq_ratio
            filtered_df = manager.failed_ids_df[mask]
            
            if len(filtered_df) == 0:
                print(f"No failed IDs for seq_ratio={args.seq_ratio}")
            else:
                print(f"\n📋 Failed IDs (seq_ratio={args.seq_ratio}):")
                print(filtered_df.to_string(index=False))
    
    else:
        # デフォルト: 統計表示
        stats = manager.get_statistics(args.seq_ratio)
        failed_ids = manager.get_failed_ids(args.seq_ratio, max_retries=3)
        
        print("\n" + "=" * 80)
        print("Failed ID Manager")
        print("=" * 80)
        print(f"Total failed records: {stats['total']}")
        print(f"IDs to skip (3+ failures): {len(failed_ids)}")
        print("\nUse --stats for detailed statistics")
        print("Use --list to see all failed IDs")
        print("Use --clear to remove records")
        print("=" * 80)


if __name__ == "__main__":
    main()