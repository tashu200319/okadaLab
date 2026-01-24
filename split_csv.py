import pandas as pd
import os

input_file = 'only_true_results.csv'
if os.path.exists(input_file):
    df = pd.read_csv(input_file)
    os.makedirs('chunks', exist_ok=True)
    chunk_size = 200
    for i in range(0, len(df), chunk_size):
        chunk_df = df.iloc[i:i + chunk_size]
        chunk_num = i // chunk_size + 1
        chunk_df.to_csv(f'chunks/chunk_{chunk_num}.csv', index=False)
    print(f"✅ {len(df)} IDs split into {chunk_num} files in 'chunks/' folder.")
else:
    print(f"❌ Error: {input_file} not found.")
