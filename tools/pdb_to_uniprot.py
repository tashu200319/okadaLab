#!/usr/bin/env python3
"""
PDB IDからUniProt IDを検索するモジュール（1対1マッピング版）
"""

import os
import requests
import pandas as pd
from typing import List, Tuple, Optional

def get_primary_uniprot_from_pdb(pdbid: str) -> Tuple[Optional[str], List[str]]:
    """
    PDB IDから最初のUniProt IDを取得（PDBe APIの順序通り）
    
    Parameters
    ----------
    pdbid : str
        PDB ID (例: "1ABC")
    
    Returns
    -------
    tuple
        (primary_uniprot_id, all_uniprot_ids)
        primary_uniprot_id: 最初のUniProt ID（なければNone）
        all_uniprot_ids: 全てのUniProt IDのリスト
    """
    pdbid = pdbid.lower()
    
    # PDBe API を使用
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdbid}"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        uniprot_ids = []
        
        # データ構造: {pdbid: {'UniProt': {uniprot_id: {...}, ...}}}
        if pdbid in data and 'UniProt' in data[pdbid]:
            uniprot_dict = data[pdbid]['UniProt']
            # 辞書のキーの順序を保持（Python 3.7+）
            uniprot_ids = list(uniprot_dict.keys())
        
        # 最初のIDを返す
        primary = uniprot_ids[0] if uniprot_ids else None
        return primary, uniprot_ids
    
    except Exception as e:
        print(f"Error fetching UniProt for {pdbid}: {e}")
        return None, []


def batch_pdb_to_uniprot(pdbids: List[str], output_file: str = "./output/links/pdb_to_uniprot_results.csv") -> pd.DataFrame:
    """
    複数のPDB IDを一括処理してCSVに保存（1対1マッピング版・追記対応）
    """
    results = []
    duplicate_check = {}
    
    print("=" * 80)
    print("PDB to UniProt Reverse Search (1-to-1 Mapping)")
    print("=" * 80)
    
    for i, pdbid in enumerate(pdbids, 1):
        pdbid = pdbid.strip().upper()
        print(f"({i}/{len(pdbids)}) Searching UniProt for {pdbid}...")
        
        primary_uniprot, all_uniprots = get_primary_uniprot_from_pdb(pdbid)
        
        if primary_uniprot:
            if len(all_uniprots) > 1:
                others = ', '.join(all_uniprots[1:])
                print(f"  ✓ Primary: {primary_uniprot} (Others: {others})")
            else:
                print(f"  ✓ Found: {primary_uniprot}")
            
            if primary_uniprot not in duplicate_check:
                duplicate_check[primary_uniprot] = []
            duplicate_check[primary_uniprot].append(pdbid)
            
            results.append({
                'pdbid': pdbid,
                'uniprotid': primary_uniprot,
                'alternative_uniprots': ', '.join(all_uniprots[1:]) if len(all_uniprots) > 1 else ''
            })
        else:
            print(f"  ✗ No UniProt ID found")
    
    # 重複通知
    print(f"\n{'=' * 80}")
    print("Duplicate UniProt IDs Check:")
    print("=" * 80)
    
    has_duplicates = False
    for uniprot_id, pdb_list in duplicate_check.items():
        if len(pdb_list) > 1:
            has_duplicates = True
            print(f"⚠ UniProt {uniprot_id} is mapped to multiple PDBs: {', '.join(pdb_list)}")
    
    if not has_duplicates:
        print("✓ No duplicate UniProt IDs found (all mappings are unique)")
    
    # DataFrameに変換して保存
    if results:
        df = pd.DataFrame(results)
        
        # ===== 追記方式に変更 =====
        if not os.path.exists(output_file):
            # ファイルが存在しない場合はヘッダー付きで作成
            df.to_csv(output_file, index=False, mode='w')
            print(f"\n📝 Created new file: {output_file}")
        else:
            # 既存ファイルに追記（ヘッダーなし）
            df.to_csv(output_file, index=False, mode='a', header=False)
            print(f"\n📝 Appended to existing file: {output_file}")
        # ==========================
        
        print(f"\n{'=' * 80}")
        print(f"Results saved to: {output_file}")
        print(f"New mappings: {len(results)}")
        print(f"Unique UniProts in this batch: {len(duplicate_check)}")
        print(f"{'=' * 80}")
        return df
    else:
        print("\nNo results found.")
        return pd.DataFrame(columns=['pdbid', 'uniprotid', 'alternative_uniprots'])


def get_unique_uniprots(csv_file: str = "./output/links/pdb_to_uniprot_results.csv",
                        output_file: str = "./output/links/unique_uniprots.csv",
                        append_mode: bool = True) -> List[str]:
    """
    CSVからユニークなUniProt IDのリストを取得
    
    Parameters
    ----------
    csv_file : str
        入力CSVファイル
    output_file : str
        ユニークなUniProt IDを保存するCSVファイル
    append_mode : bool
        True=既存データに追加, False=上書き
    
    Returns
    -------
    list
        ユニークなUniProt IDのリスト
    """
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found")
        return []
    
    # 新しいデータを読み込み
    new_df = pd.read_csv(csv_file)
    
    # 既存データがあれば読み込んで結合
    if append_mode and os.path.exists(output_file):
        existing_df = pd.read_csv(output_file)
        print(f"📋 Existing data: {len(existing_df)} entries")
        
        # 既存データと新データを結合
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
        
        # 重複を除いて最初のPDB IDを保持
        unique_df = combined_df.drop_duplicates(subset='uniprotid', keep='first')
        
        added_count = len(unique_df) - len(existing_df)
    else:
        # 重複を除いて最初のPDB IDを保持
        unique_df = new_df.drop_duplicates(subset='uniprotid', keep='first')
        added_count = len(unique_df)
    
    print(f"\n{'=' * 80}")
    print("Unique UniProt IDs Extraction")
    print(f"{'=' * 80}")
    print(f"New PDB-UniProt mappings: {len(new_df)}")
    print(f"Total unique UniProt IDs: {len(unique_df)}")
    print(f"New UniProts added: {added_count}")
    
    # 重複があるUniProt IDを表示（新データ内での重複のみ）
    duplicates = new_df[new_df.duplicated(subset='uniprotid', keep=False)]
    if not duplicates.empty:
        print(f"\n📋 UniProt IDs with multiple PDB entries in new data (keeping first only):")
        for uniprot in duplicates['uniprotid'].unique():
            pdbs = new_df[new_df['uniprotid'] == uniprot]['pdbid'].tolist()
            kept_pdb = pdbs[0]
            removed_pdbs = pdbs[1:]
            print(f"  {uniprot}:")
            print(f"    ✓ Kept: {kept_pdb}")
            print(f"    ✗ Removed: {', '.join(removed_pdbs)}")
    else:
        print("\n✓ No duplicate UniProt IDs found in new data")
    
    # ユニークなリストを保存
    unique_df.to_csv(output_file, index=False)
    print(f"\n{'=' * 80}")
    print(f"Unique mappings saved to: {output_file}")
    print(f"{'=' * 80}")
    
    return unique_df['uniprotid'].tolist()


if __name__ == "__main__":
    # 出力ディレクトリを作成
    output_dir = "./output/links"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 検索したいPDB IDリスト
    test_pdbids = [
        '3j5p', '5lzp', '6a90', '6csx', '6fkf', '6jcw', '6nby', '6rbk',
        '3j7h', '5n8n', '6a91', '6ct0', '6fu8', '6jcz', '6nht', '6rbn',
        '3j9c', '5ni1', '6a95', '6cud', '6g2j', '6jd1', '6nhv',
        '3j9i', '5nik', '6acf', '6cvj', '6giq', '6mb3', '6niy',
        '3j9j', '5nv3', '6ach', '6cvm', '6gyb', '6mdr', '6njo',
        '3j9q', '5oof', '6adm', '6dcq', '6gyn', '6mgv', '6njp',
        '3jak', '5syc', '6adq', '6dde', '6gyo', '6mgw', '6nq1',
        '3jal', '5sye', '6agf', '6ddf', '6gyu', '6mho', '6nq2',
        '3jar', '5syf', '6ap1', '6dg7', '6h3i', '6mhq', '6nr3',
        '3jas', '5syg', '6aui', '6djm', '6h3n', '6mhs', '6nsj',
        '3jat', '5t4d', '6az0', '6djn', '6hbc', '6mhv', '6nsk',
        '3jbs', '5tfy', '6b3j', '6dnf', '6hbu', '6mhy', '6nt3',
        '3jck', '5tj6', '6baj', '6dpu', '6hcy', '6mks', '6nt4',
        '3jcu', '5u6o', '6bcj', '6dpv', '6hiq', '6mlm', '6nt9',
        '3jcz', '5uj9', '6bco', '6dpw', '6hjp', '6mrt', '6nyf',
        '3jd0', '5uja', '6bcq', '6drj', '6hls', '6mru', '6nyj',
        '3jd1', '5vfo', '6bdf', '6drv', '6hn5', '6msm', '6nyn',
        '3jd2', '5vy3', '6bgl', '6dso', '6hu9', '6mwq', '6nyv',
        '3jd4', '5vy4', '6bhu', '6du8', '6hug', '6mzb', '6nzu',
        '4ci0', '5vy5', '6bjc', '6dw1', '6huj', '6mzu', '6nzw',
        '4d3e', '5w0s', '6bly', '6dwb', '6hum', '6mzv', '6nzz',
        '4v1a', '5w3e', '6bmf', '6e0g', '6huo', '6mzx', '6o1n',
        '5a0q', '5w3l', '6bqr', '6e1h', '6hwh', '6mzy', '6o1p',
        '5a1a', '5w3m', '6bqv', '6e1k', '6hzm', '6n06', '6o1u',
        '5a63', '5w3s', '6btm', '6e1m', '6i1y', '6n09', '6o20',
        '5foj', '5w5e', '6c0v', '6e3y', '6i2k', '6n1h', '6o2r',
        '5fij', '5w5f', '6c1d', '6e7b', '6i53', '6n23', '6o7t',
        '5fik', '5w68', '6c24', '6e7p', '6idf', '6n24', '6o7u',
        '5fil', '5w81', '6c26', '6e9z', '6igz', '6n25', '6o81',
        '5fim', '5wc3', '6c6l', '6ebk', '6ihb', '6n26', '6o85',
        '5fin', '5wj9', '6c70', '6ebl', '6ijj', '6n27', '6o9z',
        '5g05', '5wq7', '6c96', '6et5', '6ijo', '6n28', '6oeu',
        '5gaq', '5wq8', '6c9a', '6eti', '6ilk', '6n2d', '6oij',
        '5h1q', '5xb1', '6c9i', '6eu2', '6irs', '6n2y', '6on2',
        '5h3o', '5xnl', '6c9k', '6eu3', '6irt', '6n2z', '6qm7',
        '5i68', '5xnm', '6caj', '6ezj', '6iyc', '6n30', '6qm8',
        '5irx', '5xno', '6cjq', '6ezm', '6j0b', '6n4b', '6qnt',
        '5irz', '5xtb', '6cjt', '6ezn', '6j5t', '6nb3', '6qp6',
        '5is0', '5xte', '6cju', '6f1t', '6j8e', '6nbb', '6qpc',
        '5jx1', '5yi5', '6cmx', '6f1u', '6j8g', '6nbc', '6qpi',
        '5k0z', '5z1w', '6cnm', '6f1v', '6j8h', '6nbd', '6r3q',
        '5k12', '5z96', '6cnn', '6f1y', '6j8i', '6nbf', '6r4p',
        '5kmg', '5zdh', '6co7', '6f1z', '6j8j', '6nbh', '6r7x',
        '5lkh', '5zji', '6coy', '6fi95', '6jb1', '6nbq', '6rao',
        '5lki', '5zx5', '6coz', '6fhl', '6jcv', '6nbx', '6rap'
    ]
    
    # Step 1: PDB → UniProt 検索
    df = batch_pdb_to_uniprot(test_pdbids)
    
    if df is not None and not df.empty:
        # Step 2: ユニークなUniProt IDを抽出
        unique_uniprots = get_unique_uniprots(append_mode=False)  # やり直しなのでFalse
        
        print(f"\n{'=' * 80}")
        print("✅ PDB to UniProt mapping completed!")
        print(f"{'=' * 80}")
        print(f"📁 Output files created:")
        print(f"  1. ./output/links/pdb_to_uniprot_results.csv (all mappings)")
        print(f"  2. ./output/links/unique_uniprots.csv (unique UniProts only)")
        print(f"\n🚀 Next step: Run analysis with ./scripts/run_all_chunks.sh (or main.py --ids/--file)")
        print(f"{'=' * 80}")
    else:
        print("\n❌ No PDB-UniProt mappings found.")