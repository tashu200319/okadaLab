#!/usr/bin/env python3
"""
manual_xray_em_check.csv から True 行だけ抽出して only_true_results.csv を作る前処理。

標準:
- input : output/manual_xray_em_check.csv
- output: output/only_true_results.csv

実行例:
  python scripts/prep_true_results.py
  python scripts/prep_true_results.py --input output/manual_xray_em_check.csv --output output/only_true_results.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main() -> int:
    root = _repo_root()
    default_input = root / "output" / "manual_xray_em_check.csv"
    default_output = root / "output" / "only_true_results.csv"

    parser = argparse.ArgumentParser(
        description="Filter True rows from manual_xray_em_check.csv to only_true_results.csv"
    )
    parser.add_argument("--input", "-i", default=str(default_input), help="Input CSV path")
    parser.add_argument("--output", "-o", default=str(default_output), help="Output CSV path")
    parser.add_argument(
        "--true-col",
        type=int,
        default=2,
        help="0-based column index that contains True flag (default: 2 = 3rd column)",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise SystemExit(f"❌ Error: input not found: {input_path}")

    df = pd.read_csv(input_path)
    if df.shape[1] <= args.true_col:
        raise SystemExit(
            f"❌ Error: input has only {df.shape[1]} columns, but --true-col={args.true_col}"
        )

    col = df.iloc[:, args.true_col]
    mask = (col == True) | (col.astype(str).str.upper() == "TRUE")
    filtered = df[mask]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(output_path, index=False)

    print("-" * 30)
    print("✅ 抽出完了!")
    print(f"📂 入力ファイル: {input_path}")
    print(f"📄 出力ファイル: {output_path}")
    print(f"🔢 該当件数: {len(filtered)} 件")
    print("-" * 30)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

