
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/cli/run_normal_vs_tumor_rbp_programs.py

"""
CLI step 2/3 - Normal vs Tumor RBP-level programs (Section A).

This script takes the heavy outputs from run_compute_scores_and_tests.py and
produces an RBP-level summary plus volcano plots.

Inputs (produced by step 1)
---------------------------
From <output_dir>/explainability_scores/:
  - qc/expression_threshold.json
       • epsilon
       • rbps_both
       • rbps_tumor_only
  - stats/RBPxTranscript_wilcoxon.csv
       • one row per (Transcript_ID, RBP) with log2_fc, p_adj_fdr, etc.

From getBM CSV:
  - Gene_ID, Gene_name, Transcript_ID, Transcript_name

Main steps
----------
1) Load expression-based RBP filtering (rbps_both, epsilon).
2) Load transcript-level Wilcoxon results (RBPxTranscript).
3) Load annotation table (getBM) and build:
     - RBP_ID → RBP_name
     - Transcript_ID → Transcript_name

4) Aggregate transcript-level tests into an RBP-level summary
   (n_transcripts, prop_sig, median_log2_fc, min_p_adj_fdr, ...).
   → Saves:  normal_vs_tumor_rbp_programs/stats/RBP_summary_from_transcripts.csv

                  RBP  n_transcripts  n_sig_tx_FDR_0_01  prop_sig_tx_FDR_0_01  median_log2_fc  median_abs_log2_fc  min_p_adj_fdr
# 0    ENSG00000004478        11459.0             6387.0              0.557378        1.477235            1.982095       0.000090
# 1    ENSG00000005022        11459.0             4062.0              0.354481       -0.359225            1.567116       0.000090
# 2    ENSG00000005075        11459.0             7411.0              0.646741       -1.797480            2.421785       0.000090
# 3    ENSG00000005100        11459.0             2619.0              0.228554       -0.383918            3.169252       0.000115
# 4    ENSG00000006047        11459.0             2211.0              0.192949       -0.311887            0.966543       0.000090
# ..               ...            ...                ...                   ...             ...                 ...            ...
# 519  ENSG00000276023        11459.0              573.0              0.050004        0.242330            0.884885       0.000090
# 520  ENSG00000277726        11459.0             1857.0              0.162056        0.163447            0.958532       0.000090
# 521  ENSG00000278619        11459.0             6704.0              0.585042        1.826004            4.028156       0.000090
# 522  ENSG00000278845        11459.0             5673.0              0.495069       -1.876273            2.108227       0.000090
# 523  ENSG00000280165        11459.0             2172.0              0.189545        0.468548            0.839201       0.000090

5) Prioritize RBPs for downstream visualizations using a set of thresholds:
     - transcript FDR (TX_FDR_THR)
     - RBP-level FDR (RBP_FDR_THR)
     - strict / moderate |log2FC| cutoffs
     - fraction of significant transcripts, etc.
   → Saves:  normal_vs_tumor_rbp_programs/stats/RBP_prioritized.csv

6) Generate:
     - One RBP-level volcano summarizing all RBPs.
       → normal_vs_tumor_rbp_programs/volcano/rbp_level/volcano_RBP_level.png
     - One transcript-level volcano per prioritized RBP.
       → normal_vs_tumor_rbp_programs/volcano/transcripts_per_RBP/volcano_tx_<symbol>_<id>.png
"""

import os
import json
import argparse
from argparse import Namespace

import pandas as pd

from ..core.wilcoxon_tests import summarize_rbp_from_transcript_tests
from ..core.rbp_prioritization import select_prioritized_rbps

from ..core.config_thresholds import (
    TX_FDR_THR,
    RBP_FDR_THR,
    PRIOR_PROP_SIG_THR_GLOBAL,
    PRIOR_MIN_ABS_LOG2FC_GLOBAL,
    PRIOR_PROP_SIG_THR_LOCAL,
    PRIOR_MAX_ABS_LOG2FC_LOCAL,
    PRIOR_EXTREME_PROP_SIG_MAX,
    PRIOR_EXTREME_MIN_ABS_LOG2FC,
    PRIOR_TOP_N_GLOBAL,
    PRIOR_TOP_N_LOCAL,
    PRIOR_TOP_N_EXTREME,
    LOG2FC_STRICT_THR,
    LOG2FC_MODERATE_THR,
)

from ..plots.volcano_plots import (
    plot_rbp_level_summary,
    plot_transcript_volcano_for_rbp,
)

_DEFAULTS = dict(
    output_dir="/scratch/jsanchoz/DeepRBP/output/results/explain_tcga_compare/LIHC",
    getBM_path="/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"    
)

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Section A: Normal vs Tumor RBP-level summary + volcano plots.")
    p.add_argument("--output_dir", default=_DEFAULTS["output_dir"], 
        help="Base output directory (must contain explainability_scores/ from the compute_scores step).")
    p.add_argument("--getBM_path", default=_DEFAULTS["getBM_path"], 
        help="Path to getBM CSV mapping Gene_ID / Gene_name / Transcript_ID / Transcript_name.")
    p.add_argument("--rbps_focus", default=None, 
        help=("Comma-separated list of RBP IDs or symbols to always prioritize (e.g. 'SLU7,SRSF3,SRSF1')."))
    
    # Thresholds: keep sensible defaults from config_thresholds, but allow overriding
    p.add_argument("--tx_fdr_thr", type=float, default=TX_FDR_THR, 
        help=f"Transcript-level FDR threshold for Wilcoxon significance (default: {TX_FDR_THR}).")
    p.add_argument("--rbp_fdr_thr", type=float, default=RBP_FDR_THR, 
        help=f"RBP-level FDR threshold (used for RBP volcano y-cut) (default: {RBP_FDR_THR}).")
    p.add_argument("--log2fc_strict_thr", type=float, default=LOG2FC_STRICT_THR, 
        help=f"Strict |log2FC| threshold for 'strong' effects and volcano cutoffs (default: {LOG2FC_STRICT_THR}).")
    p.add_argument("--log2fc_moderate_thr", type=float, default=LOG2FC_MODERATE_THR, 
        help=f"Moderate |log2FC| threshold for global program prioritization (default: {LOG2FC_MODERATE_THR}).")
    
    # --- Priorización de RBPs: parámetros global/local/extreme ---
    # GLOBAL PROGRAM
    p.add_argument("--prior_prop_sig_thr_global", type=float, default=PRIOR_PROP_SIG_THR_GLOBAL,
        help=("Minimum fraction of significant transcripts to call an RBP a "
            f"'global_program' (default: {PRIOR_PROP_SIG_THR_GLOBAL})."))
    p.add_argument("--prior_min_abs_log2fc_global", type=float, default=PRIOR_MIN_ABS_LOG2FC_GLOBAL,
        help=("Minimum |median_log2_fc| for 'global_program' RBPs "
            f"(default: {PRIOR_MIN_ABS_LOG2FC_GLOBAL})."))
    p.add_argument("--prior_top_n_global", type=int, default=PRIOR_TOP_N_GLOBAL,
        help=f"Max number of 'global_program' RBPs to keep (default: {PRIOR_TOP_N_GLOBAL}).")

    # LOCAL PROGRAM (locally rewired)
    p.add_argument("--prior_prop_sig_thr_local", type=float, default=PRIOR_PROP_SIG_THR_LOCAL,
        help=("Minimum fraction of significant transcripts to call an RBP a "
            f"'local_program' (default: {PRIOR_PROP_SIG_THR_LOCAL})."))
    p.add_argument("--prior_max_abs_log2fc_local", type=float, default=PRIOR_MAX_ABS_LOG2FC_LOCAL,
        help=("Maximum |median_log2_fc| for 'local_program' RBPs "
            f"(must stay near 0; default: {PRIOR_MAX_ABS_LOG2FC_LOCAL})."))
    p.add_argument(
        "--prior_prop_sig_thr_local_strong",
        type=float,
        default=0.40,   # este lo definimos sólo para separar strong/weak
        help=("Threshold of prop_sig to define 'local_program_strong'. "
              "Local RBPs with prop_sig between prior_prop_sig_thr_local and "
              "prior_prop_sig_thr_local_strong are treated as local_program_weak "
              "and not prioritized (default: 0.40)."),
    )
    p.add_argument(
        "--prior_top_n_local",
        type=int,
        default=PRIOR_TOP_N_LOCAL,
        help=f"Max number of 'local_program_strong' RBPs to keep (default: {PRIOR_TOP_N_LOCAL}).",
    )

    # FEW BUT EXTREME
    p.add_argument("--prior_extreme_prop_sig_max", type=float, default=PRIOR_EXTREME_PROP_SIG_MAX,
        help=("Maximum fraction of significant transcripts for the "
            f"'few_but_extreme' class (default: {PRIOR_EXTREME_PROP_SIG_MAX})."))
    p.add_argument("--prior_extreme_min_abs_log2fc", type=float, default=PRIOR_EXTREME_MIN_ABS_LOG2FC,
        help=("Minimum max |log2FC| among significant transcripts for "
            f"'few_but_extreme' RBPs (default: {PRIOR_EXTREME_MIN_ABS_LOG2FC})."))
    p.add_argument("--prior_top_n_extreme", type=int, default=PRIOR_TOP_N_EXTREME,
        help=f"Max number of 'few_but_extreme' RBPs to keep (default: {PRIOR_TOP_N_EXTREME}).")
    return p.parse_args(args=argv)

# ---------- main ----------
def main(argv=None):
    args = parse_args(argv)

    # Effective thresholds (CLI may override defaults)
    tx_fdr_thr = args.tx_fdr_thr
    rbp_fdr_thr = args.rbp_fdr_thr
    prior_prop_sig_thr_global = args.prior_prop_sig_thr_global
    prior_min_abs_log2fc_global = args.prior_min_abs_log2fc_global
    prior_top_n_global = args.prior_top_n_global
    prior_prop_sig_thr_local = args.prior_prop_sig_thr_local
    prior_max_abs_log2fc_local = args.prior_max_abs_log2fc_local
    prior_prop_sig_thr_local_strong = args.prior_prop_sig_thr_local_strong
    prior_top_n_local = args.prior_top_n_local
    prior_extreme_prop_sig_max = args.prior_extreme_prop_sig_max
    prior_extreme_min_abs_log2fc = args.prior_extreme_min_abs_log2fc
    log2fc_strict_thr = args.log2fc_strict_thr
    log2fc_moderate_thr = args.log2fc_moderate_thr
    prior_top_n_extreme = args.prior_top_n_extreme

    # Directories produced by the first main (compute_scores)
    scores_root = os.path.join(args.output_dir, "explainability_scores")
    scores_qc_dir = os.path.join(scores_root, "qc")
    scores_stats_dir = os.path.join(scores_root, "stats")

    # Section A–specific directories
    sectionA_dir = os.path.join(args.output_dir, "normal_vs_tumor_rbp_programs")
    stats_dir = os.path.join(sectionA_dir, "stats")
    volcano_dir = os.path.join(sectionA_dir, "volcano")
    volc_rbp_dir = os.path.join(volcano_dir, "rbp_level")
    volc_tx_dir = os.path.join(volcano_dir, "transcripts_per_RBP")
    rest_volc_dir = os.path.join(volc_tx_dir, "rest_of_rbps")

    for d in [sectionA_dir, stats_dir, volcano_dir, volc_rbp_dir, volc_tx_dir, rest_volc_dir]:
        os.makedirs(d, exist_ok=True)

    print("=" * 80)
    print("[sectionA] Thresholds in use:")
    print(f"[sectionA]  • TX_FDR_THR                      = {tx_fdr_thr}")
    print(f"[sectionA]  • RBP_FDR_THR                     = {rbp_fdr_thr}")
    print(f"[sectionA]  • LOG2FC_STRICT_THR               = {log2fc_strict_thr}")
    print(f"[sectionA]  • LOG2FC_MODERATE_THR             = {log2fc_moderate_thr}")

    print(f"[sectionA]  • PRIOR_PROP_SIG_THR_GLOBAL       = {prior_prop_sig_thr_global}")
    print(f"[sectionA]  • PRIOR_MIN_ABS_LOG2FC_GLOBAL     = {prior_min_abs_log2fc_global}")
    print(f"[sectionA]  • PRIOR_TOP_N_GLOBAL              = {prior_top_n_global}")

    print(f"[sectionA]  • PRIOR_PROP_SIG_THR_LOCAL        = {prior_prop_sig_thr_local}")
    print(f"[sectionA]  • PRIOR_MAX_ABS_LOG2FC_LOCAL      = {prior_max_abs_log2fc_local}")
    print(f"[sectionA]  • PRIOR_TOP_N_LOCAL               = {prior_top_n_local}")

    print(f"[sectionA]  • PRIOR_EXTREME_PROP_SIG_MAX      = {prior_extreme_prop_sig_max}")
    print(f"[sectionA]  • PRIOR_EXTREME_MIN_ABS_LOG2FC    = {prior_extreme_min_abs_log2fc}")
    print(f"[sectionA]  • PRIOR_TOP_N_EXTREME             = {prior_top_n_extreme}")
    print("=" * 80 + "\n")

    # --------------------- 1) Load expression-based RBP filtering ---------------------
    expr_json = os.path.join(scores_qc_dir, "expression_threshold.json")
    if not os.path.exists(expr_json):
        raise FileNotFoundError(
            f"[sectionA] expression_threshold.json not found at:\n"
            f"  → {expr_json}\n"
            f"Please run main_compute_scores_and_tests.py first."
        )

    with open(expr_json, "r") as f:
        expr_info = json.load(f)

    rbps_both = expr_info.get("rbps_both", [])
    epsilon = expr_info.get("epsilon", None)

    print("[sectionA] Loaded expression-based RBP filtering from compute_scores step:")
    print(f"[sectionA]  • #RBPs expressed in both groups : {len(rbps_both)}")
    if epsilon is not None:
        print(f"[sectionA]  • Expression epsilon              : {epsilon:.4f}")

    # --------------------- 2) Load transcript-level Wilcoxon tests --------------------
    wilcoxon_path = os.path.join(scores_stats_dir, "RBPxTranscript_wilcoxon.csv")
    if not os.path.exists(wilcoxon_path):
        raise FileNotFoundError(
            f"[sectionA] RBPxTranscript Wilcoxon results not found at:\n"
            f"  → {wilcoxon_path}\n"
            f"Please run main_compute_scores_and_tests.py first."
        )

    print("\n[sectionA] Loading RBPxTranscript Wilcoxon results...")
    df_results_rbptx_test = pd.read_csv(wilcoxon_path)
    print(f"[sectionA]  • #rows       : {df_results_rbptx_test.shape[0]}")
    print(f"[sectionA]  • #RBPs       : {df_results_rbptx_test['RBP'].nunique()}")
    print(f"[sectionA]  • #Transcripts: {df_results_rbptx_test['Transcript_ID'].nunique()}")

    # --------------------- 3) Load annotation table (getBM) ---------------------------
    print("\n[sectionA] Loading annotation table (getBM) to build name maps...")
    getBM = pd.read_csv(args.getBM_path)
    rbp_name_map = dict(zip(getBM["Gene_ID"], getBM["Gene_name"]))
    tx_name_map = dict(zip(getBM["Transcript_ID"], getBM["Transcript_name"]))
    print(f"[sectionA]  • #annotated RBPs       : {len(rbp_name_map)}")
    print(f"[sectionA]  • #annotated transcripts: {len(tx_name_map)}")

    # --------------------- 4) RBP-level summary from transcript tests -----------------
    print("\n[sectionA] Summarizing transcript-level tests per RBP...")
    df_rbp_summary = summarize_rbp_from_transcript_tests(
        df_results_rbptx_test,
        fdr_thr=tx_fdr_thr,
    )

    out_rbp = os.path.join(stats_dir, "RBP_summary_from_transcripts.csv")
    df_rbp_summary.to_csv(out_rbp, index=False)
    print(f"[sectionA] 💾 Saved RBP-level summary to:")
    print(f"[sectionA]   → {out_rbp}")
    print(f"[sectionA]  • #RBPs in summary: {df_rbp_summary.shape[0]}")

    # Log: top RBPs by fraction of significant transcripts
    if not df_rbp_summary.empty:
        top_rbps = (
            df_rbp_summary.sort_values("prop_sig_tx_FDR_0_01", ascending=False)
            .head(10)
        )
        print(
            f"\n[sectionA] Top 10 RBPs by fraction of significant transcripts "
            f"(Transcript FDR < {TX_FDR_THR}):"
        )
        print(
            top_rbps[
                ["RBP", "prop_sig_tx_FDR_0_01", "median_log2_fc", "min_p_adj_fdr"]
            ].to_string(index=False)
        )
    
    # --------------------- 5) Prioritize RBPs for transcript-level volcanos ----------
    if args.rbps_focus:
        rbps_focus = [x.strip() for x in args.rbps_focus.split(",") if x.strip()]
        print(f"\n[sectionA] rbps_focus (user-provided): {rbps_focus}")
    else:
        rbps_focus = None
        print("\n[sectionA] rbps_focus: None (no RBPs forced by the user).")

    print("\n[sectionA] Selecting prioritized RBPs for downstream volcano plots...")
    df_prioritized = select_prioritized_rbps(
        df_rbp=df_rbp_summary,
        df_tx=df_results_rbptx_test,
        famous_rbps=rbps_focus,
        id_col="RBP",
        name_col="RBP_name",
        rbp_name_map=rbp_name_map,
        fdr_thr=tx_fdr_thr,  # consistent with transcript-level Wilcoxon
        # GLOBAL
        prop_sig_thr_global=prior_prop_sig_thr_global,
        min_abs_log2fc_global=prior_min_abs_log2fc_global,
        top_n_global=prior_top_n_global,
        # LOCAL
        prop_sig_thr_local=prior_prop_sig_thr_local,
        max_abs_log2fc_local=prior_max_abs_log2fc_local,
        prop_sig_thr_local_strong=prior_prop_sig_thr_local_strong,
        top_n_local=prior_top_n_local,
        # FEW BUT EXTREME
        extreme_prop_sig_max=prior_extreme_prop_sig_max,
        extreme_min_abs_log2fc=prior_extreme_min_abs_log2fc,
        top_n_extreme=prior_top_n_extreme,
    )

    out_prior = os.path.join(stats_dir, "RBP_prioritized.csv")
    df_prioritized.to_csv(out_prior, index=False)
    print(f"[sectionA] 💾 Saved prioritized RBPs to:")
    print(f"[sectionA]   → {out_prior}")
    print(f"[sectionA]  • #prioritized RBPs: {df_prioritized.shape[0]}")

    if not df_prioritized.empty:
        print("\n[sectionA] Prioritized RBPs by class:")
        print(
            df_prioritized[
                [
                    "RBP",
                    "RBP_name",
                    "priority_class",
                    "prop_sig_tx_FDR_0_01",
                    "median_log2_fc",
                    "max_abs_log2_fc_sig",
                    "min_p_adj_fdr",
                ]
            ].to_string(index=False)
        )

    prioritized_ids = df_prioritized["RBP"].unique().tolist()

    # --------------------- 6) RBP-level summary plot ---------------------------------
    out_rbp_summary = os.path.join(volc_rbp_dir, "summary_RBP_level.png")
    print(
    f"\n[sectionA] Generating RBP-level summary for {len(rbps_both)} RBPs "
    f"(RBP-level p_thr={rbp_fdr_thr}, |log2FC| ≥ {log2fc_strict_thr})..."
    )
    plot_rbp_level_summary(
        df_rbp_summary,
        rbp_name_map=rbp_name_map,
        p_thr=rbp_fdr_thr,
        fc_thr=log2fc_strict_thr,
        prop_sig_local_min=prior_prop_sig_thr_local,
        prop_sig_local_strong=prior_prop_sig_thr_local_strong,
        out_path=out_rbp_summary,
        point_size=25,      # puedes ajustarlo si ves los puntos muy grandes/pequeños
        prioritized_ids=prioritized_ids
    )
    print(f"[sectionA] 💾 Saved RBP-level summary to:")
    print(f"[sectionA]   → {out_rbp_summary}")

    # --------------------- 7) Transcript-level volcanos per prioritized RBP ----------
    # Map RBP → median_log2_fc for including the RBP-level summary in the title
    all_rbps = sorted(df_results_rbptx_test["RBP"].unique())
    prioritized_set = set(prioritized_ids)
    rest_rbps = [r for r in all_rbps if r not in prioritized_set]
    
    median_map = dict(zip(df_rbp_summary["RBP"], df_rbp_summary["median_log2_fc"]))

    print(
        f"\n[sectionA] Generating transcript-level volcano plots for "
        f"{len(prioritized_ids)} prioritized RBPs "
        f"(Transcript p_thr={tx_fdr_thr}, |log2FC| ≥ {log2fc_strict_thr})..."
    )

    for rbp_id in prioritized_ids:
        rbp_symbol = rbp_name_map.get(rbp_id, rbp_id)
        out_png = os.path.join(volc_tx_dir, f"volcano_tx_{rbp_symbol}_{rbp_id}.png")

        plot_transcript_volcano_for_rbp(
            df_results_rbptx_test,
            rbp_id=rbp_id,
            rbp_name=rbp_symbol,
            tx_name_map=tx_name_map,
            p_col="p_adj_fdr",
            p_thr=tx_fdr_thr,
            fc_thr=log2fc_strict_thr,
            out_path=out_png,
            rbp_level_log2fc=median_map.get(rbp_id, None),
        )
        print(f"[sectionA]   • Saved transcript-level volcano for {rbp_symbol} to:")
        print(f"[sectionA]     → {out_png}")

    # --------------------- 7.2) Transcript-level volcanos per non-prioritized RBP ----------
    print(
        f"\n[sectionA] Generating transcript-level volcano plots for "
        f"{len(rest_rbps)} non-prioritized RBPs "
        f"(saved under 'rest_of_rbps')...")

    for rbp_id in rest_rbps:  
        rbp_symbol = rbp_name_map.get(rbp_id, rbp_id)
        out_png = os.path.join(rest_volc_dir, f"volcano_tx_{rbp_symbol}_{rbp_id}.png")

        plot_transcript_volcano_for_rbp(
            df_results_rbptx_test,
            rbp_id=rbp_id,
            rbp_name=rbp_symbol,
            tx_name_map=tx_name_map,         
            p_col="p_adj_fdr",
            p_thr=tx_fdr_thr,
            fc_thr=log2fc_strict_thr,
            out_path=out_png,
            rbp_level_log2fc=median_map.get(rbp_id, None),
        )
    print("\n[sectionA] ✅ Section A pipeline finished successfully.\n")


if __name__ == "__main__":
    main()
