#!/usr/bin/env python3
"""
only_true_results.csv から data/chunks/chunk_*.csv を作る chunk 生成スクリプト。

標準:
- input : output/only_true_results.csv
- output: data/chunks/

実行例:
  python scripts/make_chunks.py
  python scripts/make_chunks.py --input output/only_true_results.csv --out-dir data/chunks --chunk-size 200
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _extract_uniprot_ids(df: pd.DataFrame) -> List[str]:
    # 1) uniprotid列があればそれを使う
    for col in ["uniprotid", "UniProtID", "uniprot_id", "UNIPROTID"]:
        if col in df.columns:
            raw = df[col].dropna().astype(str).tolist()
            break
    else:
        # 2) 1列目を使う
        raw = df.iloc[:, 0].dropna().astype(str).tolist()

    ids: List[str] = []
    for rid in raw:
        uid = rid.strip().split(",")[0].strip()
        if not uid:
            continue
        if uid.lower() == "uniprotid":  # ヘッダー混入対策
            continue
        ids.append(uid)
    # 重複削除しつつ順序維持
    seen = set()
    uniq = []
    for uid in ids:
        if uid in seen:
            continue
        seen.add(uid)
        uniq.append(uid)
    return uniq


def main() -> int:
    root = _repo_root()
    default_input = root / "output" / "only_true_results.csv"
    default_out_dir = root / "data" / "chunks"

    parser = argparse.ArgumentParser(description="Split UniProt IDs into chunk CSVs")
    parser.add_argument("--input", "-i", default=str(default_input), help="Input CSV path")
    parser.add_argument(
        "--out-dir",
        "-o",
        default=str(default_out_dir),
        help="Output directory for chunk_*.csv",
    )
    parser.add_argument("--chunk-size", type=int, default=200, help="Rows per chunk (default: 200)")

    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out_dir)

    if not input_path.exists():
        raise SystemExit(f"❌ Error: input not found: {input_path}")

    df = pd.read_csv(input_path)
    ids = _extract_uniprot_ids(df)
    if not ids:
        raise SystemExit(f"❌ Error: No UniProt IDs found in: {input_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    chunk_size = int(args.chunk_size)
    chunk_num = 0
    for i in range(0, len(ids), chunk_size):
        chunk_num += 1
        chunk_ids = ids[i : i + chunk_size]
        out_df = pd.DataFrame({"uniprotid": chunk_ids})
        out_df.to_csv(out_dir / f"chunk_{chunk_num}.csv", index=False)

    print(f"✅ {len(ids)} IDs split into {chunk_num} files in '{out_dir}/' folder.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

