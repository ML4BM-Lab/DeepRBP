
# src/deeprbp/explainability_module/differential_regulation/de_expression_conditions.py

import os
import sys
import argparse
import subprocess
from typing import List, Optional, Tuple
import pandas as pd

from ...data_loading.config_loader import ConfigParser

from .utils_reg import _parse_conditions, _slugify_condition
from .utils_de import (
    _filter_meta,
    _export_counts_for_r,
    _export_meta_for_r,
    _filter_counts_by_feature_spec,
)

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Differential expression (voom-limma) between two conditions for differential_regulation."
    )
    # Core inputs
    p.add_argument("--config_path", required=True, help="Path to config_diff_regulation.yaml")
    p.add_argument("--output_dir", required=True, help="Base output directory (DE outputs will be written under output_dir/de/...)")
    p.add_argument("--select_category", required=True, help="Tissue/category to subset samples (must exist in metadata)")
    p.add_argument("--conditions", required=True, help="Two conditions separated by comma, e.g. 'Primary_Tumor,Solid_Tissue_Normal'")
    # RBP feature specification (ONLY needed for RBP-level DE)
    p.add_argument("--feature_spec_path", default=None,
        help=(
            "Path to DeepRBP feature spec Excel. "
            "Used only for RBP-level DE. "
            "Must contain sheet 'RBPs' with column 'Ensembl_gene_id'."
        ),
    )
    # Which DE levels to run
    p.add_argument("--levels", default="genes,transcripts,rbps", help="Comma-separated levels to run: genes,transcripts,rbps. Default: genes,transcripts")
    # Counts file names (dataset-dependent but overridable)
    p.add_argument("--genes_counts_file", default="gn_counts.csv", help="Counts file for genes (samples x genes). Default: gn_counts.csv")
    p.add_argument("--transcripts_counts_file", default="trans_counts.csv", help="Counts file for transcripts (samples x transcripts). Default: trans_counts.csv")
    # R execution
    p.add_argument("--rscript_bin", default="Rscript", help="Rscript binary (default: Rscript)")
    p.add_argument("--r_runner_path", default=None, help="Path to run_voom-limma_diff_reg.R. If not provided, it is resolved relative to the package.")

    # Volcano thresholds
    p.add_argument("--logfc_thresh", type=float, default=1.0, help="Absolute log2FC threshold for volcano coloring")
    p.add_argument("--p_cut_type", choices=["fdr", "pvalue"], default="fdr", help="Significance column for volcano cut")
    p.add_argument("--p_cut_value", type=float, default=0.05, help="Threshold for p_cut_type")
    p.add_argument("--show_legend", action="store_true", help="Show volcano legend (default off).")
    return p.parse_args(args=argv)

def main(argv: Optional[List[str]] = None) -> None:
    # --------------------------------------------------------
    # Parse and validate conditions
    # --------------------------------------------------------
    args = parse_args(argv)
    conds = _parse_conditions(args.conditions)
    if len(conds) != 2:
        raise ValueError(
            "Differential expression requires exactly TWO conditions "
            "(e.g. Tumor,Normal)"
        )
    cond_a, cond_b = conds

    # --------------------------------------------------------
    # Parse and validate requested levels
    # --------------------------------------------------------
    levels = [x.strip() for x in args.levels.split(",") if x.strip()]
    valid_levels = {"genes", "transcripts", "rbps"}

    for lv in levels:
        if lv not in valid_levels:
            raise ValueError(f"Invalid level '{lv}'. Valid options: {sorted(valid_levels)}")

    # RBP-level DE requires feature spec
    if "rbps" in levels and not args.feature_spec_path:
        raise ValueError(
            "RBP-level DE requested but --feature_spec_path was not provided"
        )

    # --------------------------------------------------------
    # Prepare output structure
    # --------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)
    de_root = os.path.join(args.output_dir, "de")
    os.makedirs(de_root, exist_ok=True)

    print("\n" + "=" * 88)
    print("[diff_reg / DE] Step 2: differential expression (voom-limma)")
    print(f"[diff_reg / DE]  • category   : {args.select_category}")
    print(f"[diff_reg / DE]  • conditions : {cond_a} vs {cond_b}")
    print(f"[diff_reg / DE]  • levels     : {levels}")
    print("=" * 88 + "\n")

    # --------------------------------------------------------
    # Load YAML configuration (dataset conventions only)
    # --------------------------------------------------------
    cfg = ConfigParser(args.config_path)

    dataset_dir = cfg.get("test_path_files")
    getBM_path = cfg.get("getBM_path", default=None)
    sample_category_col = cfg.get("sample_category")
    disease_condition_col = cfg.get("disease_condition")

    if not dataset_dir or not os.path.isdir(dataset_dir):
        raise ValueError(f"Invalid dataset directory: {dataset_dir}")

    meta_path = os.path.join(dataset_dir, "phenotype_metadata.csv")
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Missing metadata file: {meta_path}")

    # --------------------------------------------------------
    # Resolve R runner
    # --------------------------------------------------------
    r_runner = args.r_runner_path or _resolve_default_r_runner()
    if not os.path.exists(r_runner):
        raise FileNotFoundError(f"Cannot find R runner: {r_runner}")

    # --------------------------------------------------------
    # Load and filter metadata (Python-side)
    # --------------------------------------------------------
    df_meta = pd.read_csv(meta_path, index_col=0)

    df_meta_filt = _filter_meta(
        df_meta,
        sample_category_col=sample_category_col,
        disease_condition_col=disease_condition_col,
        select_category=args.select_category,
        cond_a=cond_a,
        cond_b=cond_b,
    )

    # Export filtered metadata for R
    inputs_dir = os.path.join(de_root, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)

    meta_for_r_path = os.path.join(inputs_dir, "metadata_filtered.tsv")
    _export_meta_for_r(df_meta_filt, meta_for_r_path)
   
    # --------------------------------------------------------
    # Run DE independently for each requested level
    # --------------------------------------------------------
    for level in levels:
        # Select counts file and R feature mode
        if level == "genes":
            counts_file = args.genes_counts_file
            feature_level_for_r = "genes"

        elif level == "transcripts":
            counts_file = args.transcripts_counts_file
            feature_level_for_r = "transcripts"

        else: # rbps
            counts_file = args.genes_counts_file
            feature_level_for_r = "genes"

        counts_path = os.path.join(dataset_dir, counts_file)
        if not os.path.exists(counts_path):
            print(f"[diff_reg / DE] ⚠️ Missing counts file, skipping: {counts_path}")
            continue

        print(f"[diff_reg / DE] Loading counts for level '{level}': {counts_path}")
        df_counts = pd.read_csv(counts_path, index_col=0)        

        # Restrict to RBPs using feature spec (only for RBP level)
        if level == "rbps":
            df_counts = _filter_counts_by_feature_spec(
                df_counts, args.feature_spec_path
            )

        # ----------------------------------------------------
        # Align samples between counts and metadata
        # ----------------------------------------------------
        common = sorted(set(df_counts.index).intersection(df_meta_filt.index))
        if len(common) < 4:
            raise ValueError(
                f"Not enough overlapping samples after filtering for level '{level}'. "
                f"Need >=4 total and >=2 per condition. Overlap={len(common)}"
            )

        df_counts = df_counts.loc[common]
       
        # ----------------------------------------------------
        # Export counts matrix for R (features x samples)
        # ----------------------------------------------------
        counts_for_r_path = os.path.join(inputs_dir, f"counts_{level}_features_x_samples.tsv")
        _export_counts_for_r(df_counts, counts_for_r_path)

        # Output directory for this level and contrast
        out_level_dir = os.path.join(de_root, level, f"{_slugify_condition(cond_b)}_vs_{_slugify_condition(cond_a)}")
        os.makedirs(out_level_dir, exist_ok=True)

        # ----------------------------------------------------
        # Call R voom-limma runner
        # ----------------------------------------------------
        _run_one_level(
            level=level,
            counts_path=counts_for_r_path,
            meta_path=meta_for_r_path,
            condition_col=disease_condition_col,
            cond_a=cond_a,
            cond_b=cond_b,
            output_dir=out_level_dir,
            feature_level_for_r=feature_level_for_r,
            getBM_path=getBM_path, # annotation handled in R if needed
            rscript_bin=args.rscript_bin,
            r_runner_path=r_runner,
            logfc_thresh=args.logfc_thresh,
            p_cut_type=args.p_cut_type,
            p_cut_value=args.p_cut_value,
            show_legend=args.show_legend,
        )

    print("\n" + "=" * 88)
    print("[diff_reg / DE] ✅ Done.")
    print(f"[diff_reg / DE] Outputs written under:\n  → {de_root}")
    print("=" * 88 + "\n")

# ============================================================
# Helpers
# ============================================================
def _resolve_default_r_runner() -> str:
    # Resolve relative to this file: differential_regulation/ -> explainability_module/de_limma/
    here = os.path.dirname(os.path.abspath(__file__))
    # .../differential_regulation/
    de_limma_dir = os.path.normpath(os.path.join(here, "..", "de_limma"))
    return os.path.join(de_limma_dir, "run_voom-limma_diff_reg.R")

def _run_one_level(
    level: str,
    counts_path: str,
    meta_path: str,
    condition_col: str,
    cond_a: str,
    cond_b: str,
    output_dir: str,
    feature_level_for_r: str,
    getBM_path: str,
    rscript_bin: str,
    r_runner_path: str,
    logfc_thresh: float,
    p_cut_type: str,
    p_cut_value: float,
    show_legend: bool,
) -> None:
    """
    Build and execute the R voom-limma command for one feature level.
    """
    cmd = [
        rscript_bin,
        r_runner_path,
        "--counts_path", counts_path,
        "--metadata_path", meta_path,
        "--condition_col", condition_col,
        "--condition_a", cond_a,
        "--condition_b", cond_b,
        "--output_dir", output_dir,
        "--feature_level", feature_level_for_r,
        "--getBM_path", getBM_path if getBM_path else "",
        "--logfc_thresh", str(logfc_thresh),
        "--p_cut_type", p_cut_type,
        "--p_cut_value", str(p_cut_value),
        "--show_legend", "TRUE" if show_legend else "FALSE",
    ]

    print("\n" + "=" * 88)
    print(f"[diff_reg / DE] Running voom-limma for level: {level}")
    print("[diff_reg / DE] Command:")
    print("  " + " ".join(cmd))
    print("=" * 88)

    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()