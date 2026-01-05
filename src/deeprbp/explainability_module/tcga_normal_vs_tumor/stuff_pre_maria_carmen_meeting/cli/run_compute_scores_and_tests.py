
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/cli/run_compute_scores_and_tests.py

"""
CLI step 1/3 - Compute DeepRBP explainability scores and per-transcript Wilcoxon tests.

High-level pipeline
-------------------
1) Load configuration + trained DeepRBP predictor.
2) Build NORMAL and TUMOR datasets for the requested TCGA category.
3) Classify RBPs by expression:
     - RBPs expressed in both groups  → used in Normal vs Tumor comparisons.
     - RBPs expressed only in tumor   → used later in tumor-only analyses.
   → Saves:  explainability_scores/qc/expression_threshold.json
             (epsilon, rbps_both, rbps_tumor_only)

4) Run (or reuse from cache) DeepRBP explainability for NORMAL and TUMOR:
     - Per-sample Tx(sample,RBP) score matrices
     - Standard result_table and aggregated scores
   → Saves:  explainability_scores/Normal/
               df_scores_TxRBP_per_sample.csv
               df_scores_TxRBP.csv
               df_scores_GxRBP.csv
               result_table.csv
             explainability_scores/Tumor/
               (same structure)

5) Filter score matrices to RBPs expressed in both groups and generate a QC
   histogram comparing global score distributions (Normal vs Tumor).
   → Saves:  explainability_scores/qc/scores_distribution_both.png

6) Run (or reuse from cache) Mann-Whitney tests per (Transcript, RBP),
   comparing NORMAL vs TUMOR across samples.
   → Saves:  explainability_scores/stats/RBPxTranscript_wilcoxon.csv

            Transcript_ID              RBP  n_normal  n_tumor  median_normal  median_tumor  median_diff   log2_fc   p_value  p_adj_fdr
# 0        ENST00000610897  ENSG00000004478        11       73       0.003816      0.009408     0.005592  1.301710  0.000011   0.000225
# 1        ENST00000615252  ENSG00000004478        11       73      -0.001155     -0.003984    -0.002829 -1.785852  0.000048   0.000537
# 2        ENST00000378609  ENSG00000004478        11       73       0.000237      0.000354     0.000117  0.577753  0.329752   0.420076
# 3        ENST00000461893  ENSG00000004478        11       73       0.000111      0.000558     0.000447  2.331568  0.000222   0.001483
# 4        ENST00000471354  ENSG00000004478        11       73       0.000025      0.000172     0.000147  2.763309  0.000009   0.000198
# ...                  ...              ...       ...      ...            ...           ...          ...       ...       ...        ...
# 6004511  ENST00000504061  ENSG00000280165        11       73      -0.000160      0.000028     0.000188 -2.489005  0.198387   0.277857
# 6004512  ENST00000369476  ENSG00000280165        11       73      -0.000862     -0.001294    -0.000432 -0.585676  0.276918   0.364067
# 6004513  ENST00000482244  ENSG00000280165        11       73       0.001687      0.001726     0.000039  0.032654  0.968270   0.976952
# 6004514  ENST00000476116  ENSG00000280165        11       73       0.001084      0.000427    -0.000657  1.343787  0.071345   0.120425
# 6004515  ENST00000362018  ENSG00000280165        11       73       0.000095      0.000192     0.000096  1.011163  0.049717   0.089837

The outputs of this script are then consumed by:
  - run_normal_vs_tumor_rbp_programs.py  (Section A - RBP programs)
  - run_tumor_only_rbp_targets.py        (Section B - tumor-only targets, WIP)
"""

import os
import json
import argparse
from argparse import Namespace
import numpy as np
import pandas as pd

from ....data_loading.config_loader import ConfigParser
from ....data_preparation.prepare_data import PrepareData
from ....training_module.model import PredictorModel
from ....explainability_module.main_explainer import data_preparation

from ..core.expression_thresholds import classify_rbps_by_expression
from ..core.rbp_score_matrix_ops import filter_score_matrix_by_rbps
from ..io_cache.explainability_cache import load_or_run_explainer, load_or_run_rbp_tx_wilcoxon
from ..plots.plot_distribution import plot_scores_distribution_two_groups

_DEFAULTS = dict(
    config_path="/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml", # si quieres usar el train:"/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat2_luego_borra.yaml"
    model_ckpt_path="/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt",
    scaler_dir="/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data",
    output_dir="/scratch/jsanchoz/DeepRBP/output/results/explain_tcga_compare/LIHC",
    select_category="Liver_Hepatocellular_Carcinoma",
    normal_condition="Solid_Tissue_Normal",
    tumor_condition="Primary_Tumor"
)

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Compute DeepRBP scores for Normal vs Tumor and test per RBP (Mann-Whitney)."
    )
    p.add_argument("--config_path", default=_DEFAULTS["config_path"], help="YAML del explainer (usa select_category, etc.)")
    p.add_argument("--model_ckpt_path", default=_DEFAULTS["model_ckpt_path"], help="Checkpoint del predictor entrenado")
    p.add_argument("--scaler_dir", default=_DEFAULTS["scaler_dir"], help="Directorio con scaler/sigma del predictor")
    p.add_argument("--output_dir", default=_DEFAULTS["output_dir"], help="Directorio base de salida")
    p.add_argument("--select_category", default=_DEFAULTS["select_category"], help="Categoría (ej: Liver_Hepatocellular_Carcinoma)")
    p.add_argument("--normal_condition", default=_DEFAULTS["normal_condition"])
    p.add_argument("--tumor_condition",  default=_DEFAULTS["tumor_condition"])
    p.add_argument("--rbps_focus", default=None, help="Comma-separated list of RBP IDs or symbols to always prioritize (e.g. 'SLU7,SRSF3,SRSF1').",)
    return p.parse_args(args=argv)


# ---------- main ----------
def main(argv=None):
    args = parse_args(argv)
   
    scores_root = os.path.join(args.output_dir, "explainability_scores")
    qc_dir   = os.path.join(scores_root, "qc")
    stats_dir = os.path.join(scores_root, "stats")

    for d in [scores_root, qc_dir, stats_dir]:
        os.makedirs(d, exist_ok=True)

    print("\n" + "=" * 80)
    print("[compute_scores] DeepRBP TCGA - Step 1: explainability scores + Wilcoxon")
    print(f"[compute_scores]  • Category     : {args.select_category}")
    print(f"[compute_scores]  • NORMAL label : {args.normal_condition}")
    print(f"[compute_scores]  • TUMOR label  : {args.tumor_condition}")
    print(f"[compute_scores]  • Output root  : {args.output_dir}")
    print("=" * 80 + "\n")

    # --- 1) load config + model ---
    cfg = ConfigParser(args.config_path)
    cfg.update("save_per_sample_scores", True)

    print("\n[compute_scores] Loading trained model...")
    print(f"[compute_scores]  • Config YAML : {args.config_path}")
    print(f"[compute_scores]  • Checkpoint  : {args.model_ckpt_path}")
    model = PredictorModel.load_from_checkpoint(args.model_ckpt_path)

    # --- 2) NORMAL dataset ---
    print("\n[compute_scores] Loading NORMAL dataset...")
    cfg.update("select_condition", [args.normal_condition])
    args_normal = Namespace(
        **{**vars(args),
           "output_dir": os.path.join(scores_root, "Normal")}
    )
    ds_normal, getBM = data_preparation(cfg, args_normal, args.select_category)
    n_normal = ds_normal.features["scaled_rbp_df"].shape[0]
    n_rbps_total = ds_normal.features["scaled_rbp_df"].shape[1]
    print(f"[compute_scores]  • Normal samples : {n_normal}")
    print(f"[compute_scores]  • Total RBPs     : {n_rbps_total}")

    # --- 3) TUMOR dataset ---
    print("\n[compute_scores] Loading TUMOR dataset...")
    cfg.update("select_condition", [args.tumor_condition])
    args_tumor = Namespace(
        **{**vars(args),
           "output_dir": os.path.join(scores_root, "Tumor")}
    )
    ds_tumor, _ = data_preparation(cfg, args_tumor, args.select_category)
    n_tumor = ds_tumor.features["scaled_rbp_df"].shape[0]
    print(f"[compute_scores]  • Tumor samples  : {n_tumor}")

    # --- 4) Expression-based RBP filtering (for Section A & B) ---
    rbps_both, rbps_tumor_only, epsilon = classify_rbps_by_expression(ds_normal, ds_tumor)

    print("\n[compute_scores] Expression-based RBP filtering")
    print(f"[compute_scores]  • RBPs in both groups : {len(rbps_both)}")
    print(f"[compute_scores]  • RBPs tumor-only     : {len(rbps_tumor_only)}")
    print(f"[compute_scores]  • Expression epsilon  : {epsilon:.4f}")

    eps_path = os.path.join(qc_dir, "expression_threshold.json")
    with open(eps_path, "w") as f:
        json.dump(
            {
                "epsilon": float(epsilon),
                "rbps_both": rbps_both,
                "rbps_tumor_only": rbps_tumor_only,
            },
            f,
            indent=2,
        )
    print(f"[compute_scores] 💾 Saved expression summary to:\n  → {eps_path}")

    # --- 5) Explainability NORMAL/TUMOR (heavy step) ---
    print("\n[compute_scores] Running / loading DeepRBP explainability (NORMAL)...")
    res_normal = load_or_run_explainer(
        condition_label="NORMAL",
        cfg=cfg,
        args_cond=args_normal,
        dataset=ds_normal,
        model=model,
        getBM=getBM,
    )

    print("\n[compute_scores] Running / loading DeepRBP explainability (TUMOR)...")
    res_tumor = load_or_run_explainer(
        condition_label="TUMOR",
        cfg=cfg,
        args_cond=args_tumor,
        dataset=ds_tumor,
        model=model,
        getBM=getBM,
    )

    GxN_t = res_normal["df_per_sample"]
    GxT_t = res_tumor["df_per_sample"]

    print("\n[compute_scores] Per-sample DeepLIFT score matrices (Tx x SamplexRBP):")
    print(f"[compute_scores]  • NORMAL: shape={GxN_t.shape}")
    print(f"[compute_scores]  • TUMOR : shape={GxT_t.shape}")

    # --- 6) Restrict to RBPs expressed in both + QC distribution plot ---
    GxN_both = filter_score_matrix_by_rbps(GxN_t, rbps_both)
    GxT_both = filter_score_matrix_by_rbps(GxT_t, rbps_both)

    print("\n[compute_scores] Filtered matrices to RBPs expressed in both groups:")
    print(f"[compute_scores]  • NORMAL (both): shape={GxN_both.shape}")
    print(f"[compute_scores]  • TUMOR  (both): shape={GxT_both.shape}")

    # Quick QC: compare score distributions (Normal vs Tumor) on a common RBP set.
    # Note: this plot is meant to check that both groups share the same scale.
    print("\n[compute_scores] Plotting score distributions for RBPs expressed in both groups...")
    out_hist = os.path.join(qc_dir, "scores_distribution_both.png")
    plot_scores_distribution_two_groups(
        GxN_both,
        GxT_both,
        out_path=out_hist,
        title="DeepRBP scores - RBPs expressed in both groups",
        density = True,
        clip_quantiles = (0.5, 99.5)
    )
    print(f"[compute_scores] 💾 Saved QC histogram to:\n  → {out_hist}")   

    # --- 7) RBP×Transcript Wilcoxon (shared by Section A & B) ---
    print("\n[compute_scores] RBPxTranscript Wilcoxon layer (using cache if available)...")
    load_or_run_rbp_tx_wilcoxon(
        GxN_both,
        GxT_both,
        stats_dir=stats_dir,
        filename="RBPxTranscript_wilcoxon.csv",
        overwrite=False,
        verbose=False,
        show_examples=False,
    )
    print("\n[compute_scores] ✅ Done. Explainability + Wilcoxon ready for Section A/B scripts.\n")

if __name__ == "__main__":
    main()