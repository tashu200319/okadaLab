#!/usr/bin/env python3
"""
filtered.csvからX-rayとEM両方を持つもののみ抽出
"""

import pandas as pd
import argparse

def extract_both_only(
    filtered_csv: str = 'output/summaries/filtered.csv',
    check_csv: str = 'output/filtered_xray_em_check.csv',
    output_file: str = 'output/summaries/filtered_both.csv'
):
    """
    filtered.csvからX-rayとEM両方を持つものだけ抽出
    """
    # データ読み込み
    print("📖 Loading data...")
    filtered_df = pd.read_csv(filtered_csv)
    check_df = pd.read_csv(check_csv)
    
    print(f"  filtered.csv: {len(filtered_df)} rows")
    print(f"  check results: {len(check_df)} rows")
    
    # 両方持つIDを抽出
    both_ids = check_df[(check_df['xray_count'] > 0) & (check_df['em_count'] > 0)]['uniprotid'].tolist()
    
    print(f"\n✅ Both X-ray and EM: {len(both_ids)} IDs")
    
    # filtered.csvから抽出
    both_df = filtered_df[filtered_df['uniprotid'].isin(both_ids)]
    
    # 保存
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    both_df.to_csv(output_file, index=False)
    
    print(f"\n💾 Saved: {output_file}")
    print(f"   Rows: {len(both_df)}")
    print(f"   Percentage: {len(both_df)/len(filtered_df)*100:.1f}% of filtered.csv")


def main():
    parser = argparse.ArgumentParser(description='Extract IDs with both X-ray and EM')
    parser.add_argument('--filtered', default='output/summaries/filtered.csv')
    parser.add_argument('--check', default='output/filtered_xray_em_check.csv')
    parser.add_argument('--output', default='output/summaries/filtered_both.csv')
    
    args = parser.parse_args()
    
    extract_both_only(
        filtered_csv=args.filtered,
        check_csv=args.check,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
