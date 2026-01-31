#!/usr/bin/env python3
"""
Minimal evaluation of linearity between
distance mean and distance std per CSV.

Output columns:
- uniprot
- n_pairs
- pearson_r
- r2
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd

def extract_uniprot_id(csv_path: str) -> str:
    base = os.path.basename(csv_path)
    return base.replace("score_details_", "").split("_")[0]

def eval_one_csv(csv_path: str):
    # 一部の score_details_*_20.csv は実体が PNG などのバイナリの可能性があるため、
    # UnicodeDecodeError などの読み込みエラーが発生した場合はスキップする
    try:
        df = pd.read_csv(csv_path)
    except UnicodeDecodeError:
        print(f"Skipping (UnicodeDecodeError): {csv_path}")
        return None

    if "distance mean" not in df.columns or "distance std" not in df.columns:
        return None

    x = pd.to_numeric(df["distance mean"], errors="coerce")
    y = pd.to_numeric(df["distance std"], errors="coerce")

    mask = x.notna() & y.notna()
    x = x[mask].to_numpy()
    y = y[mask].to_numpy()

    n = len(x)
    if n < 10:
        return None

    r = float(np.corrcoef(x, y)[0, 1])
    r2 = r * r

    return {
        "uniprot": extract_uniprot_id(csv_path),
        "n_pairs": n,
        "pearson_r": r,
        "r2": r2,
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "score_details_*_20.csv")))

    rows = []
    for f in files:
        res = eval_one_csv(f)
        if res:
            rows.append(res)

    df = pd.DataFrame(rows).sort_values("pearson_r", ascending=False)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    print(f"Saved: {args.out_csv} ({len(df)} rows)")

if __name__ == "__main__":
    main()
