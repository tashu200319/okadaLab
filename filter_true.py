#!/usr/bin/env python3
"""
互換ラッパー（旧: filter_true.py）。

新しい入口は `python scripts/prep_true_results.py` です。
このファイルは過去の手順互換のために残し、scripts版へ委譲します。
"""

from __future__ import annotations

import runpy
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parent
    target = root / "scripts" / "prep_true_results.py"
    if not target.exists():
        raise SystemExit(f"❌ scripts/prep_true_results.py not found: {target}")
    runpy.run_path(str(target), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
