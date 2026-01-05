
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/cli/run_tumor_only_rbp_targets.py

"""
CLI step 3/3 - Tumor-only RBP targets (Section B).

This script analyzes RNA-binding proteins (RBPs) that are expressed
exclusively in tumor samples and ranks their splicing regulatory impact
using DeepRBP explainability scores.

Unlike Section A, this step does not perform Normal vs Tumor statistical
testing. It focuses on tumor-exclusive regulatory programs that emerge
only in cancer.

Inputs (produced by steps 1 and 2)
----------------------------------
From <output_dir>/explainability_scores/:

  - qc/expression_threshold.json
       • epsilon (stored as metadata)
       • rbps_tumor_only (ENSG IDs)

  - Tumor/result_table.csv
       • one row per (RBP, Transcript) with DeepRBP Score
       • required columns:
           RBP_ID, RBP_name, Gene_ID, Transcript_ID, Score

From getBM CSV:
  - annotation table used to ensure consistent RBP_ID ↔ RBP_name mapping

Main steps
----------
1) Load tumor-only RBPs from expression_threshold.json.
2) Load tumor explainability scores and filter by RBP_ID ∈ rbps_tumor_only.
3) Validate RBP annotations using getBM.
4) Aggregate scores at the RBP level and compute:
     - mean_abs_score
     - n_genes
     - n_transcripts
     - ranking_score = mean_abs_score x log1p(n_genes)
   → Saves: tumor_only_rbp_ranking.csv

5) Export top-N target transcripts per RBP ranked by |Score|.
   → Saves: tumor_only_top_targets/<RBP_name>__top_transcripts.csv

6) Print diagnostic messages for selected RBPs (e.g. SRSF1, SRSF3, SLU7).

Outputs
-------
- tumor_only_rbp_ranking.csv
- tumor_only_top_targets/
- tumor_only_step3_metadata.json

Downstream visualization
------------------------
Designed to support:
  - Bubble plot (mean_abs_score vs n_genes; size = n_transcripts)
  - Bar plot of top-ranked tumor-only RBPs
  - Score distributions for selected RBPs
"""

# Plots que haría en Section B (tumor-only)
# - Bubble plot (ideal como figura resumen):
# X: mean_abs_score
# Y: n_genes
# tamaño: n_transcripts
# resaltar: SRSF1 / SRSF3 / SLU7

# - Barplot Top-20 por ranking_score
# etiqueta: RBP_name
# color especial (o borde) para los famosos
# Para cada RBP famoso:
# violin/histograma de Score (o abs(score)) en tumor
# y un Top targets table (genes o transcritos) exportado a CSV (ya lo generamos)

import os, json, argparse
import numpy as np
import pandas as pd

FAMOUS_RBPS_DEFAULT = ["SRSF1", "SRSF3", "SLU7"]
DEFAULT_OUTPUT_DIR = "/scratch/jsanchoz/DeepRBP/output/results/explain_tcga_compare/LIHC"
DEFAULT_GETBM_PATH = "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"

# then remove this:
args = parse_args([
    "--output_dir", "/scratch/jsanchoz/DeepRBP/output/results/explain_tcga_compare/LIHC",
    "--getBM_path", "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv",
    "--rbps_focus", "SRSF1,SRSF3,SLU7",
    "--min_abs_score", "0.0",
    "--min_n_genes", "0",
    "--top_n", "50"
])

def parse_args(args=None):
    ap = argparse.ArgumentParser(
        description="STEP 3/3 – Tumor-only RBP targets (Section B)"
    )
    ap.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR,
        help="Base output directory (default: TCGA LIHC output dir)")
    ap.add_argument("--getBM_path", default=DEFAULT_GETBM_PATH, help="Path to getBM annotation CSV")
    ap.add_argument("--rbps_focus", default=",".join(FAMOUS_RBPS_DEFAULT), help="Comma-separated list of RBPs to highlight")
    ap.add_argument("--min_abs_score", type=float, default=0.0, help="Minimum mean absolute score to keep an RBP")
    ap.add_argument("--min_n_genes", type=int, default=0, help="Minimum number of target genes per RBP")
    ap.add_argument("--top_n", type=int, default=50, help="Number of top targets to export per RBP")
    return ap.parse_args(args)

def main():
    args = parse_args()
    output_dir = args.output_dir

    # --------------------------------------------------
    # 1) Load expression-based tumor-only RBPs (STEP 2)
    # --------------------------------------------------
    expr_json = os.path.join(
        output_dir, "explainability_scores", "qc", "expression_threshold.json"
    )
    if not os.path.exists(expr_json):
        raise FileNotFoundError(f"Missing expression_threshold.json: {expr_json}")

    with open(expr_json) as f:
        expr_data = json.load(f)

    rbps_tumor_only = set(expr_data["rbps_tumor_only"])
    epsilon = expr_data.get("epsilon", None)

    # --------------------------------------------------
    # 2) Load tumor explainability scores
    # --------------------------------------------------
    tumor_scores_path = os.path.join(
        output_dir, "explainability_scores", "Tumor", "df_scores_GxRBP.csv"
    )

    if not os.path.exists(tumor_scores_path):
        raise FileNotFoundError(f"Missing df_scores_GxRBP.csv: {tumor_scores_path}")

    scores_wide = pd.read_csv(tumor_scores_path)

    scores_wide.head()
    scores_wide.shape
    scores_wide.columns[:5]

    scores_long = scores_wide.melt(
        id_vars="Gene_ID",
        var_name="RBP_ID",
        value_name="Score"
    ) 
    
    scores_long = scores_long[
        scores_long["RBP_ID"].isin(rbps_tumor_only)
    ].copy()
        
    # --------------------------------------------------
    # 3) Ranking per RBP
    # --------------------------------------------------
    summary = (
        scores
        .groupby(["RBP_ID", "RBP_name"], as_index=False)
        .agg(
            mean_abs_score=("Score", lambda x: float(np.mean(np.abs(x)))),
            n_genes=("Gene_ID", "nunique"),
            n_transcripts=("Transcript_ID", "nunique")
        )
    )

    summary["ranking_score"] = (
        summary["mean_abs_score"] * np.log1p(summary["n_genes"])
    )

    summary = summary[
        (summary["mean_abs_score"] >= args.min_abs_score) &
        (summary["n_genes"] >= args.min_n_genes)
    ].sort_values("ranking_score", ascending=False)

    summary_path = os.path.join(outdir, "tumor_only_rbp_ranking.csv")
    summary.to_csv(summary_path, index=False)

    # --------------------------------------------------
    # 4) Top targets per RBP + prints cariñosos
    # --------------------------------------------------
    focus = {x.strip() for x in args.rbps_focus.split(",") if x.strip()}
    top_dir = os.path.join(outdir, "tumor_only_top_targets")
    os.makedirs(top_dir, exist_ok=True)

    for (rbp_id, rbp_name), df in scores.groupby(["RBP_ID", "RBP_name"]):
        if rbp_name in focus:
            print("############################################")
            print(f"💙 Processing your beloved RBP: {rbp_name} ({rbp_id})")
            print("############################################")

        df = df.copy()
        df["abs_score"] = df["Score"].abs()

        top_tx = df.sort_values("abs_score", ascending=False).head(args.top_n)
        top_tx.to_csv(
            os.path.join(top_dir, f"{rbp_name}__top_transcripts.csv"),
            index=False
        )

    # --------------------------------------------------
    # 5) Metadata
    # --------------------------------------------------
    meta = {
        "epsilon_from_step2": epsilon,
        "n_rbps_tumor_only": len(rbps_tumor_only),
        "n_rows_scores": len(scores),
        "rbps_focus": sorted(list(focus))
    }

    with open(os.path.join(outdir, "tumor_only_step3_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"[OK] STEP 3 completed → {summary_path}")

if __name__ == "__main__":
    main()

