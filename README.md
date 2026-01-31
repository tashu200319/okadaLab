# okadaLab DSA (Distance Structure Analysis)

## リポジトリ構成（来年・別の人向け）

- **`scripts/`**: 実行用スクリプト（この中のものを使う）
- **`data/`**: 入力データ置き場（この中に集約するのが標準）
  - `data/chunks/`: `chunk_1.csv` など
  - `data/raw/`: RCSBの巨大CSVなど元データ
- **`output/`**: 生成物（解析結果・ログなど。git管理しない）

当面は互換のため `chunks/` も読めますが、**新しく置くなら `data/chunks/` が標準**です。

## 連番でchunkを1つずつ実行する（標準）

`data/chunks/chunk_1.csv` → `data/chunks/chunk_2.csv` → … の順に、**1個ずつ** `main.py` を回します。  
ログは “生成物” として **`output/` 配下** に自動で保存されます。

### 基本

```bash
./scripts/run_all_chunks.sh
```

### 範囲指定（例: chunk_1〜chunk_3 だけ）

```bash
./scripts/run_all_chunks.sh --start 1 --end 3
```

### main.py に渡す代表的な引数

```bash
./scripts/run_all_chunks.sh --seq-ratio 20 --max-pdbs 50
./scripts/run_all_chunks.sh --no-parallel --workers 1
./scripts/run_all_chunks.sh --no-heatmap
```

### ログの場所

- `output/logs/run_all_chunks_YYYYMMDD_HHMMSS/chunk_N.log`

別ターミナルで以下が便利です。

```bash
tail -f output/logs/run_all_chunks_*/chunk_1.log
```

## パスの方針（整理しやすくするためのルール）

- **生成物は `output/` に寄せる**（ログや解析結果は `output/` 配下）
- **パスの優先順位は「引数 > Config > デフォルト相対パス」**
  - 例: `main.py --output-dir ./output` が最優先
  - 指定がなければ `core/config.py` の `OUTPUT_DIR`（デフォルト `./output/`）
