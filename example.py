#!/usr/bin/env python3
"""
DSA解析の実行例
"""

import os
import pandas as pd
from config import Config
from main import count_pdb, prep, run_DSA
from uniprot_handler import UniprotData
from visualization import generate_heatmap, plot_cis_analysis


def example1_single_protein():
    """例1: 単一タンパク質の解析"""
    print("=" * 80)
    print("Example 1: Single Protein Analysis")
    print("=" * 80)
    
    # 設定
    config = Config(
        SEQ_RATIO=80,
        USE_XRAY=True,
        USE_EM=False,
        VERBOSE=True
    )
    
    uniprotid = "P01308"  # インスリン
    
    # UniProtデータ取得
    unidata = UniprotData(uniprotid)
    print(f"\nProtein: {unidata.get_fullname()}")
    print(f"Organism: {unidata.get_organism()}")
    
    # PDB数チェック
    if not count_pdb(uniprotid, methods=config.METHODS_SELECTED):
        print("Not enough PDB entries")
        return
    
    # データ準備
    seqdata, all_pdblist = prep(
        uniprotid,
        methods=config.METHODS_SELECTED,
        verbose=config.VERBOSE
    )
    
    print(f"\nNormal PDBs: {len(all_pdblist[0])}")
    print(f"Substitution PDBs: {len(all_pdblist[1])}")
    
    # 解析実行
    seqdata1 = seqdata.filter(like=uniprotid)
    pdbtuple = tuple(all_pdblist[0] + all_pdblist[1])
    seqdata2 = seqdata.loc[:, seqdata.columns.str.startswith(pdbtuple)]
    norsub_seqdata = pd.concat([seqdata1, seqdata2], axis=1)
    
    score, log_text, summary_df = run_DSA(
        uniprotid,
        norsub_seqdata,
        export=True,
        seqtype='nor+sub',
        methods=config.METHODS_SELECTED,
        seq_ratio=config.SEQ_RATIO,
        cis_threshold=config.CIS_THRESHOLD,
        dirpath=config.OUTPUT_DIR,
        verbose=config.VERBOSE
    )
    
    if summary_df is not None:
        print("\n" + "=" * 80)
        print("Results:")
        print("=" * 80)
        print(summary_df.to_string(index=False))
        
        # ヒートマップ生成
        heatmap_path = os.path.join(
            config.OUTPUT_DIR,
            f"{uniprotid}_{config.SEQ_RATIO}_heatmap.png"
        )
        generate_heatmap(score, heatmap_path)


def example2_multiple_proteins():
    """例2: 複数タンパク質の比較解析"""
    print("\n" + "=" * 80)
    print("Example 2: Multiple Protein Comparison")
    print("=" * 80)
    
    config = Config(SEQ_RATIO=70, VERBOSE=False)
    
    # 複数のUniProt ID
    uniprot_ids = ["P01308", "P00004", "P00720"]
    
    results = []
    
    for uniprotid in uniprot_ids:
        try:
            print(f"\nProcessing {uniprotid}...")
            
            unidata = UniprotData(uniprotid)
            fullname = unidata.get_fullname()
            
            if not count_pdb(uniprotid, methods=config.METHODS_SELECTED):
                print(f"  Skipped (not enough PDBs)")
                continue
            
            seqdata, all_pdblist = prep(
                uniprotid,
                methods=config.METHODS_SELECTED,
                verbose=False
            )
            
            seqdata1 = seqdata.filter(like=uniprotid)
            pdbtuple = tuple(all_pdblist[0] + all_pdblist[1])
            seqdata2 = seqdata.loc[:, seqdata.columns.str.startswith(pdbtuple)]
            norsub_seqdata = pd.concat([seqdata1, seqdata2], axis=1)
            
            score, log_text, summary_df = run_DSA(
                uniprotid,
                norsub_seqdata,
                export=False,
                seqtype='nor+sub',
                methods=config.METHODS_SELECTED,
                seq_ratio=config.SEQ_RATIO,
                verbose=False
            )
            
            if summary_df is not None:
                summary_df['uniprotid'] = uniprotid
                summary_df['fullName'] = fullname
                results.append(summary_df)
                print(f"  Success: {fullname}")
        
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    # 結果の結合と表示
    if results:
        combined = pd.concat(results, ignore_index=True)
        print("\n" + "=" * 80)
        print("Comparison Results:")
        print("=" * 80)
        print(combined[['uniprotid', 'fullName', 'Entries', 'Chains', 
                       'Length', 'UMF', 'cis']].to_string(index=False))


def example3_custom_parameters():
    """例3: カスタムパラメータでの解析"""
    print("\n" + "=" * 80)
    print("Example 3: Custom Parameters")
    print("=" * 80)
    
    # カスタム設定
    config = Config(
        SEQ_RATIO=60,           # 配列使用率を下げる
        CIS_THRESHOLD=3.5,      # cis閾値を変更
        USE_XRAY=True,
        USE_EM=True,            # EMも含める
        CHAIN_THRESHOLD=5,      # より多くのChainを要求
        VERBOSE=True
    )
    
    print(f"\nCustom Configuration:")
    print(f"  SEQ_RATIO: {config.SEQ_RATIO}%")
    print(f"  CIS_THRESHOLD: {config.CIS_THRESHOLD} Å")
    print(f"  Methods: {config.METHODS_SELECTED}")
    print(f"  CHAIN_THRESHOLD: {config.CHAIN_THRESHOLD}")
    
    uniprotid = "P01308"
    
    # 通常通り解析実行
    if count_pdb(uniprotid, methods=config.METHODS_SELECTED):
        seqdata, all_pdblist = prep(
            uniprotid,
            methods=config.METHODS_SELECTED,
            verbose=config.VERBOSE
        )
        
        print(f"\nWith custom parameters, found {len(seqdata.columns)-1} chains")


def example4_visualization_only():
    """例4: 既存データの可視化のみ"""
    print("\n" + "=" * 80)
    print("Example 4: Visualization Only")
    print("=" * 80)
    
    config = Config()
    uniprotid = "P01308"
    
    # 既存のスコアファイルを読み込み（存在する場合）
    score_file = os.path.join(
        config.OUTPUT_DIR,
        f"{uniprotid}_{config.SEQ_RATIO}_score_nor+sub.csv"
    )
    
    if os.path.exists(score_file):
        print(f"\nLoading existing score data from {score_file}")
        score = pd.read_csv(score_file)
        
        # ヒートマップ生成
        heatmap_path = os.path.join(
            config.OUTPUT_DIR,
            f"{uniprotid}_heatmap_regenerated.png"
        )
        generate_heatmap(score, heatmap_path)
        print(f"Heatmap regenerated: {heatmap_path}")
    else:
        print(f"\nScore file not found: {score_file}")
        print("Please run analysis first (Example 1)")


if __name__ == "__main__":
    # 実行したい例を選択
    print("\nDSA Analysis Examples\n")
    print("1. Single protein analysis")
    print("2. Multiple protein comparison")
    print("3. Custom parameters")
    print("4. Visualization only")
    
    choice = input("\nSelect example (1-4, or 'all'): ").strip()
    
    if choice == '1':
        example1_single_protein()
    elif choice == '2':
        example2_multiple_proteins()
    elif choice == '3':
        example3_custom_parameters()
    elif choice == '4':
        example4_visualization_only()
    elif choice.lower() == 'all':
        example1_single_protein()
        example2_multiple_proteins()
        example3_custom_parameters()
        example4_visualization_only()
    else:
        print("Invalid choice. Running Example 1 by default.")
        example1_single_protein()
    
    print("\n" + "=" * 80)
    print("Examples completed!")
    print("=" * 80)