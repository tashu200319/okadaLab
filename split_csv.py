import pandas as pd
import os

# 入力は “生成物” として output/ 配下が標準。
# 互換のため、ルート直下 only_true_results.csv もフォールバックで探す。
input_file = 'output/only_true_results.csv'
if os.path.exists(input_file):
    df = pd.read_csv(input_file)
else:
    input_file = 'only_true_results.csv'
    if os.path.exists(input_file):
        df = pd.read_csv(input_file)
    else:
        print(f"❌ Error: {input_file} not found.")
        raise SystemExit(1)

# chunkの標準出力先は data/chunks/（新しい置き場）
out_dir = "data/chunks"
os.makedirs(out_dir, exist_ok=True)

chunk_size = 200
for i in range(0, len(df), chunk_size):
    chunk_df = df.iloc[i:i + chunk_size]
    chunk_num = i // chunk_size + 1
    chunk_df.to_csv(f"{out_dir}/chunk_{chunk_num}.csv", index=False)
print(f"✅ {len(df)} IDs split into {chunk_num} files in '{out_dir}/' folder.")
