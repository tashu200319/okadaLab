# okadaLab DSA (Distance Structure Analysis)

## 卒業研究テーマ（概要）

本プロジェクトは **DSA（Distance Structure Analysis）** を用いて、異なる構造決定手法（**X線結晶構造解析（X-ray）** および **電子顕微鏡法（cryo-EM / EM）**）で得られたタンパク質立体構造の違いを **定量的に検証**することを目的としています。  
特に、タンパク質全体の体積変化（膨張・収縮）を評価する指標として **UMF（UnMorphness Factor）** を算出し、解析します。

## 解析システム構成

- **Main Engine (`main.py`)**: 全体フロー制御、並列処理、出力（summary/score_details/heatmap等）
- **Core Modules (`core/`)**: 設定、データ取得、構造解析、距離計算、レポート生成などの基幹処理
- **Tools (`tools/`)**: 解析対象の選定、結果の再集計、可視化などの補助ツール

## 解析フロー（概略）

- **ステップ1：解析対象の選定**
  - 例: `tools/check_xray_and_em_manual.py` などで X-ray / EM 両方を持つ UniProt を選別
- **ステップ2：データ取得と前処理**
  - UniProt からメタデータ取得: `core/uniprot_handler.py`
  - RCSB PDB から mmCIF ダウンロード: `core/structure_analyzer.py`
- **ステップ3：配列整列と領域特定**
  - 複数構造（Chain）間で共通範囲を同定: `core/sequence_processor.py`
  - `SEQ_RATIO`（デフォルト 20%）に基づき解析範囲を決定: `core/config.py`
- **ステップ4：距離計算とスコアリング**
  - 残基ペア距離（例: CA-CA）を計算: `core/distance_calculator.py`
  - 残基ペアごとの score: \(\mathrm{score} = \frac{\mathrm{distance\ mean}}{\mathrm{distance\ std}}\)
  - UMF（本コードの定義）: 全残基ペアの score を平均した値
  - 出力整形: `core/report_generator.py`
- **ステップ5：可視化と分析**
  - ヒートマップ（残基ペアごとのスコア分布）: `core/visualization.py`
  - 散布図（Distance vs Score / Std）: `tools/plot_score_log.py`

## 主な設定項目（`core/config.py`）

- **`USE_XRAY` / `USE_EM`**: 使用する構造決定手法
- **`SEQ_RATIO`**: 使用する配列長の割合
- **`CHAIN_THRESHOLD`**: 解析に含める最小チェーン数

## 出力成果物（`output/`）

- **`output/summaries/summary.csv`**: UMF、解像度など主要指標の一覧
- **`output/score_details/score_details_<UniProt>_<seq_ratio_int>.csv`**: 残基ペアごとの詳細（平均距離、標準偏差、score等）
- **`output/score_details/scatter/`**: 散布図（PNG。`with_score_log/` や `with_distance_std/` を含む）
- **`output/logs/`**: 実行ログ（chunkごと）

## 可視化：散布図（`tools/plot_score_log.py`）

`output/score_details/with_and_search/` 等の `score_details_*.csv` から散布図を生成します。

- **score_log**: 横軸=平均距離、縦軸=score（対数表示）
- **distance_std**: 横軸=平均距離、縦軸=距離の標準偏差

例:

```bash
python tools/plot_score_log.py --y-axis score_log
python tools/plot_score_log.py --y-axis distance_std
python tools/plot_score_log.py --both
```

## リポジトリ構成（来年・別の人向け）

- **`scripts/`**: 実行用スクリプト（この中のものを使う）
- **`data/`**: 入力データ置き場（この中に集約するのが標準）
  - `data/chunks/`: `chunk_1.csv` など
  - `data/raw/`: RCSBの巨大CSVなど元データ
- **`output/`**: 生成物（解析結果・ログなど。git管理しない）

当面は互換のため `chunks/` も読めますが、**新しく置くなら `data/chunks/` が標準**です。

## 配布（外付けHDDで渡す場合）

このプロジェクトは外付けHDDで配布します。作業するPCに `okadaLab/` フォルダごとコピーして使ってください。

- **正本**: コピーした `okadaLab/` フォルダ
- **生成物**: `output/` に出ます（git管理外）

## セットアップ（推奨: venv）

別PCでも環境差分で壊れないよう、プロジェクト専用のPython環境（venv）を作ることを推奨します。

```bash
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -r requirements.txt
```

## 連番でchunkを1つずつ実行する（標準）

`data/chunks/chunk_1.csv` → `data/chunks/chunk_2.csv` → … の順に、**1個ずつ** `main.py` を回します。  
ログは “生成物” として **`output/` 配下** に自動で保存されます。

## chunk作成（前処理）

手元の結果から `data/chunks/chunk_*.csv` を作り直したい時の標準手順です。

### 1) True行だけ抽出（任意）

`output/manual_xray_em_check.csv` から True 行だけ抜いた `output/only_true_results.csv` を作ります。

```bash
python scripts/prep_true_results.py
```

### 2) chunkを生成

`output/only_true_results.csv` から `data/chunks/chunk_*.csv` を作ります。

```bash
python scripts/make_chunks.py
```

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
