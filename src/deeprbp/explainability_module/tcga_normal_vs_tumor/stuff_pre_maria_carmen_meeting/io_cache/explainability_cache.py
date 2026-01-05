
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/explainability_cache.py

import os
import pandas as pd
from typing import Dict
from argparse import Namespace

from ....explainability_module.main_explainer import (
    run_explainability_pipeline,
    save_results,
)
from ....data_loading.config_loader import ConfigParser
from ....training_module.model import PredictorModel

from ..core.wilcoxon_tests import compute_rbp_transcript_wilcoxon

def load_or_run_explainer(
    condition_label: str,
    cfg: ConfigParser,
    args_cond: Namespace,
    dataset,
    model: PredictorModel,
    getBM,
) -> Dict[str, pd.DataFrame]:
    """
    Compute or load DeepRBP explainability outputs for a single condition.

    If df_scores_TxRBP_per_sample.csv and result_table.csv already exist
    in args_cond.output_dir, they are loaded from disk and the heavy
    DeepLIFT pass is skipped.

    Returns
    -------
    dict
        Keys:
        - 'df_per_sample' : Transcript × (Sample×RBP) matrix
        - 'result_table'  : gene-level summary table
        - 'df_scores_TxRBP', 'df_scores_GxRBP', 'df_scores_HLxRBP'
          (may be None if loaded from cache).
    """
    out_dir = args_cond.output_dir
    os.makedirs(out_dir, exist_ok=True)

    per_sample_path = os.path.join(out_dir, "df_scores_TxRBP_per_sample.csv")
    result_table_path = os.path.join(out_dir, "result_table.csv")

    if os.path.exists(per_sample_path) and os.path.exists(result_table_path):
        print(f"\n[tcga_compare] Found existing explainability outputs for {condition_label}:")
        print(f"  → {out_dir}")
        print("[tcga_compare] Skipping DeepRBP attribution recomputation.")

        df_per_sample = pd.read_csv(per_sample_path, header=[0, 1], index_col=0)
        result_table = pd.read_csv(result_table_path, index_col=0)

        return {
            "df_per_sample": df_per_sample,
            "result_table": result_table,
            "df_scores_TxRBP": None,
            "df_scores_GxRBP": None,
            "df_scores_HLxRBP": None,
        }

    print(f"\n[tcga_compare] No cached explainability outputs for {condition_label}.")
    print(f"[tcga_compare] Running DeepRBP explainability for {condition_label}...")
    res = run_explainability_pipeline(cfg, dataset, model, getBM, analyze_hidden_layer=False)

    save_results(
        path_save_results=out_dir,
        df_scores_TxRBP=res["df_scores_TxRBP"],
        df_scores_GxRBP=res["df_scores_GxRBP"],
        result_table=res["result_table"],
        df_scores_HLxRBP=res["df_scores_HLxRBP"],
        df_per_sample=res["df_per_sample"],
    )
    print(f"[tcga_compare] Saved explainability outputs for {condition_label} to:\n  → {out_dir}")
    return res

def load_or_run_rbp_tx_wilcoxon(
    GxN_both,
    GxT_both,
    stats_dir: str,
    filename: str = "RBPxTranscript_wilcoxon.csv",
    overwrite: bool = False,
    **wilcoxon_kwargs,
) -> pd.DataFrame:
    """
    Heavy step: compute or load Wilcoxon tests per (Transcript, RBP).

    If the CSV exists in stats_dir and overwrite=False, it is loaded
    from disk. Otherwise, compute with compute_rbp_transcript_wilcoxon
    and save.

    Parameters
    ----------
    GxN_both, GxT_both : pd.DataFrame
        Score matrices (Normal vs Tumor) restricted to RBPs expressed in both.
    stats_dir : str
        Folder where the CSV will be stored / loaded.
    filename : str
        Name of the CSV file to use.
    overwrite : bool
        If True, force recomputation even if file exists.
    wilcoxon_kwargs :
        Extra kwargs passed to compute_rbp_transcript_wilcoxon (e.g. verbose=False).

    Returns
    -------
    pd.DataFrame
        The RBPxTranscript Wilcoxon results.
    """
    os.makedirs(stats_dir, exist_ok=True)
    path = os.path.join(stats_dir, filename)

    if (not overwrite) and os.path.exists(path):
        print(f"\n[tcga_compare] Found cached RBPxTranscript Wilcoxon results in:")
        print(f"  → {path}")
        df = pd.read_csv(path)
        print(f"[tcga_compare]  • #rows     : {df.shape[0]}")
        if "RBP" in df.columns:
            print(f"[tcga_compare]  • #RBPs     : {df['RBP'].nunique()}")
        if "Transcript_ID" in df.columns:
            print(f"[tcga_compare]  • #Tx total : {df['Transcript_ID'].nunique()}")
        return df

    print("\n[tcga_compare] Computing RBPxTranscript Wilcoxon tests (no cache found)...")
    df = compute_rbp_transcript_wilcoxon(GxN_both, GxT_both, **wilcoxon_kwargs)
    df.to_csv(path, index=False)
    print(f"[tcga_compare] 💾 Saved RBPxTranscript-level tests to:\n  → {path}")
    print(f"[tcga_compare]  • #rows     : {df.shape[0]}")
    print(f"[tcga_compare]  • #RBPs     : {df['RBP'].nunique()}")
    print(f"[tcga_compare]  • #Tx total : {df['Transcript_ID'].nunique()}")

    return df