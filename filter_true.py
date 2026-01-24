import pandas as pd
import os

# 設定
input_file = './output/manual_xray_em_check.csv'
output_file = './output/only_true_results.csv'

if not os.path.exists(input_file):
    print(f"❌ Error: {input_file} が見つかりません。")
else:
    # CSV読み込み
    df = pd.read_csv(input_file)
    
    # 3列目（インデックス2）が True または "True" のものを抽出
    # bool型と文字列型の両方に対応
    mask = (df.iloc[:, 2] == True) | (df.iloc[:, 2].astype(str).str.upper() == 'TRUE')
    filtered_df = df[mask]
    
    # 保存
    filtered_df.to_csv(output_file, index=False)
    
    print("-" * 30)
    print(f"✅ 抽出完了!")
    print(f"📂 入力ファイル: {input_file}")
    print(f"📄 出力ファイル: {output_file}")
    print(f"🔢 該当件数: {len(filtered_df)} 件")
    print("-" * 30)
