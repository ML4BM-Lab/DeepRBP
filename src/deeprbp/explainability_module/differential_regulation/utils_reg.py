# src/deeprbp/explainability_module/differential_regulation/utils_reg.py

import re
import pandas as pd
from typing import Dict
from argparse import Namespace

from ...explainability_module.main_explainer import (
    run_explainability_pipeline,
    save_results,
)
from ...data_loading.config_loader import ConfigParser
from ...training_module.model import PredictorModel

def _load_or_run_explainer(
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


def _slugify_condition(name: str) -> str:
    """
    Turn arbitrary condition labels into safe folder names.
    - Keeps alphanumerics, '_', '-', '.'
    - Collapses whitespace into '_'
    """
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_\-\.]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "Condition"

def _parse_conditions(raw: str) -> List[str]:
    """
    Parse conditions from:
      - "A,B,C"
      - "A;B;C"
      - multiple --conditions flags are NOT supported here; keep it simple.
    """
    if not raw or not raw.strip():
        return []
    # accept comma or semicolon
    parts = re.split(r"[;,]", raw.strip())
    return [p.strip() for p in parts if p.strip()]
