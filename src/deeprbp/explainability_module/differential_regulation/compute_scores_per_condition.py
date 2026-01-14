# src/deeprbp/explainability_module/differential_regulation/compute_scores_per_condition.py

import os
import argparse
from argparse import Namespace
from datetime import datetime
from typing import List

from ...data_loading.config_loader import ConfigParser
from ...training_module.model import PredictorModel
from ...explainability_module.main_explainer import data_preparation

from .utils_reg import (
    _load_or_run_explainer,
    _slugify_condition,
    _parse_conditions,
)

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Compute DeepRBP DeepLIFT explainability scores independently for multiple conditions.")
    p.add_argument("--config_path", required=True, help="YAML configuration for the explainer (loaded by ConfigParser).")
    p.add_argument("--model_ckpt_path", required=True, help="Checkpoint of the trained DeepRBP predictor.")
    p.add_argument("--output_dir", required=True, help="Base output directory. Will create output_dir/explainability_scores/<Condition>/")
    p.add_argument("--select_category", required=True, 
     help="Dataset category / tissue identifier. "
             "This value MUST be present and correctly defined in the input metadata, "
             "otherwise the dataset cannot be loaded.")
    p.add_argument("--conditions", required=True, 
     help="Comma-separated list of condition labels to process "
             "(e.g. 'Control,Treated'). Each condition MUST exist in the metadata.")
    return p.parse_args(args=argv)

def main(argv=None) -> None:
    args = parse_args(argv)

    conditions: List[str] = _parse_conditions(args.conditions)
    if len(conditions) == 0:
        raise ValueError(
            "No conditions provided. Use --conditions 'ConditionA,ConditionB'."
        )

    print("\n" + "=" * 88)
    print("[diff_reg / scores] Step 1: compute DeepRBP DeepLIFT scores per condition")
    print(f"[diff_reg / scores]  • select_category : {args.select_category}")
    print(f"[diff_reg / scores]  • conditions      : {conditions}")
    print(f"[diff_reg / scores]  • output_dir      : {args.output_dir}")
    print("=" * 88 + "\n")

    # ---- Prepare output dirs ----
    scores_root = os.path.join(args.output_dir, "explainability_scores")
    os.makedirs(scores_root, exist_ok=True)

    # ---- Load config ----
    print("[diff_reg / scores] Loading explainer configuration...")
    cfg = ConfigParser(args.config_path)
     
    # ---- Load model ----
    print("[diff_reg / scores] Loading trained DeepRBP predictor...")
    print(f"[diff_reg / scores]  • Config YAML    : {args.config_path}")
    print(f"[diff_reg / scores]  • Model CKPT     : {args.model_ckpt_path}")
    model = PredictorModel.load_from_checkpoint(args.model_ckpt_path)

    # ---- Process each condition independently ----
    for cond in conditions:
        cond_slug = _slugify_condition(cond)
        cond_out_dir = os.path.join(scores_root, cond_slug)
        os.makedirs(cond_out_dir, exist_ok=True)

        print("\n" + "-" * 88)
        print(f"[diff_reg / scores] Processing condition: {cond}")
        print(f"[diff_reg / scores]  • condition_slug : {cond_slug}")
        print(f"[diff_reg / scores]  • output_dir     : {cond_out_dir}")
        print("-" * 88)

        # Update condition selection in config
        print(f"[diff_reg / scores] Updating config: select_condition = [{cond}]")
        cfg.update("select_condition", [cond])

        # Build condition-specific Namespace
        args_cond = Namespace(**vars(args))
        args_cond.output_dir = cond_out_dir

        # ---- Load dataset ----
        print("[diff_reg / scores] Loading dataset for condition...")
        ds_cond, getBM = data_preparation(cfg, args_cond, args.select_category)

        n_samples = None
        n_rbps = None
        try:
            n_samples = int(ds_cond.features["scaled_rbp_df"].shape[0])
            n_rbps = int(ds_cond.features["scaled_rbp_df"].shape[1])
            print(f"[diff_reg / scores]  • samples : {n_samples}")
            print(f"[diff_reg / scores]  • RBPs    : {n_rbps}")
        except Exception:
            print("[diff_reg / scores]  • Dataset size information not available.")

        # ---- Run or load explainability ----
        print("[diff_reg / scores] Running / loading DeepRBP explainability...")
        res = _load_or_run_explainer(
            condition_label=cond,
            cfg=cfg,
            args_cond=args_cond,
            dataset=ds_cond,
            model=model,
            getBM=getBM,
        )

        if isinstance(res, dict) and "df_per_sample" in res:
            try:
                print(
                    "[diff_reg / scores]  • df_per_sample shape:",
                    res["df_per_sample"].shape,
                )
            except Exception:
                pass

    print("\n" + "=" * 88)
    print("[diff_reg / scores] ✅ Done.")
    print(f"[diff_reg / scores] Manifest saved to:\n  → {run_manifest_path}")
    print("=" * 88 + "\n")


if __name__ == "__main__":
    main()