
import pandas as pd
import numpy as np
import argparse
from scipy.stats import mannwhitneyu, combine_pvalues

from .utils_complex import split_to_list, map_gene_names_to_ids, build_expanded_corr, extract_block_contrasts
from .plots_complex import plot_correlation_boxplot, plot_overview_boxplots_cached

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "This script analyzes correlation patterns across genes and explainability scores of RBPs"
            "using CORUM complex data. It generates correlation heatmaps and performs "
            "statistical tests to evaluate whether RBPs within complexes show significantly higher "
            "correlations compared to RBPs outside these complexes."
        )
    )
    parser.add_argument('--corum_file_path', type=str,
                        default='/scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/complex_analysis/corum_results.txt',
                        help='Path to the CORUM complexes data file.')
    parser.add_argument('--getBM_file_path', type=str,
                        default='/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv',
                        help='Path to the getBM CSV file containing selected gene and RBP info.')
    parser.add_argument('--scores_file_path', type=str,
                        default='/scratch/jsanchoz/DeepRBP/output/results/stuff/explainability_deeplift_knock_t_stat/results/df_scores_GxRBP.csv',
                        help='Path to the CSV file with explainability scores for gene-RBP pairs.')
    parser.add_argument('--output_dir', type=str, default='',
                        help='Path to save output files such as heatmaps and statistics.')
    parser.add_argument('--include-non-family', dest='include_non_family', action='store_true', help='Include group F (non-family) in the expanded matrix.')
    return parser.parse_args()

def main():
    """
    - stat: float
        Mann-Whitney U test statistic comparing correlations inside vs outside the complex.
    - pval: float
        P-value from the Mann-Whitney U test (alternative='greater').
    """
    args = parse_args()

    # --- Load CORUM complexes ---
    print(f"Loading Corum results from {args.corum_file_path} ...")
    df_corum = pd.read_csv(args.corum_file_path, sep='\t')
    print(df_corum.head())

    # --- Load gene name -> gene id mapping (getBM) ---
    print(f"\nLoading gene mapping (getBM) from {args.getBM_file_path} ...")
    getBM = pd.read_csv(args.getBM_file_path)

    # --- Columns to lists ---
    print("\nParsing CORUM gene lists...")
    df_corum['subunits_gene_name'] = df_corum['subunits_gene_name'].apply(split_to_list)
    df_corum['subunits_gene_name_synonyms'] = df_corum['subunits_gene_name_synonyms'].apply(split_to_list)

    # --- Map names to Ensembl IDs ---
    print("\nMapping gene names to Gene_IDs for each complex...")
    df_corum['subunits_gene_id'] = df_corum.apply(lambda row: map_gene_names_to_ids(row, getBM), axis=1)

    # --- Explainability scores (abs) ---
    print(f"Loading scores from {args.scores_file_path} ...")
    scores_GxRBP = pd.read_csv(args.scores_file_path, index_col=0)
    scores_GxRBP = abs(scores_GxRBP)
    print("First rows of the loaded scores:")
    print(scores_GxRBP.head())
    print(scores_GxRBP.shape)
    print("Data preparation completed.")

    # --- Correlations across RBPs ---
    print("Computing RBP Pearson correlation matrix...")
    corr_scores = scores_GxRBP.corr() # DataFrame 1348x1348

    # --- Build expanded correlation matrix (optionally include non-family 'F') ---
    print(f"Building expanded correlation matrix (include_non_family={args.include_non_family})...")
    expanded_corr, groups = build_expanded_corr(
            corr_scores,
            df_corum,
            complex_id_col="complex_id",
            rbp_list_col="subunits_gene_id",
            complexes_order=(8369, 8370, 8371, 8372, 8391),
            include_non_family=args.include_non_family
    )

    # --- Containers for per-complex results (used for overview plot; no recomputation) ---
    rows = []            # long/tidy rows: Correlation, Group, Complex (A..E), Name (human-readable)
    per_pair_info = []   # one dict per complex with pval and y placement
    all_stats = []
    all_pvals = []

    # --- Iterate complexes in order; skip 'F' if present ---
    for i, key in enumerate(k for k in groups if k != 'F'):
        cname = df_corum.at[i, "complex_name"]  # assumes A..E follow rows 0..4
        print(f"\nAnalyzing complex {key} — {cname}")
        
        # Extract KEY-KEY (upper triangle) vs KEY-nonKEY blocks
        key_values, nonkey_values = extract_block_contrasts(
            expanded_corr, groups, group_key = key, deduplicate_nonkey_labels = True
        )
        print(f'[{key}-{key} vs {key}-not{key}]: {len(key_values)}, {len(nonkey_values)}\n')
        
        # Mann–Whitney U (one-sided: inside > outside)
        stat, pval = mannwhitneyu(key_values, nonkey_values, alternative='greater')
        print(f"Mann-Whitney U test statistic: {stat:.4f}, p-value: {pval:.4e}")
        if pval < 0.05:
            print("💥 Correlations within the complex are significantly higher.")
        else:
            print("🫠 No sufficient evidence that correlations within the complex are higher.")
        
        all_stats.append(stat)
        all_pvals.append(pval)

        # Cache data for per-complex boxplot (optional) and for the overview (no recomputation later)
        if args.output_dir:
            plot_correlation_boxplot(
                key_values, nonkey_values,
                complex_idx=key,
                complex_name=cname, 
                pval=pval,
                output_dir=args.output_dir
                )

        # Long/tidy rows for overview figure
        rows += [{"Correlation": v, "Group": "Complex", "Complex": key, "Name": cname} for v in key_values]
        rows += [{"Correlation": v, "Group": "Outside", "Complex": key, "Name": cname} for v in nonkey_values]

        # For bracket placement in the overview
        y99 = float(np.nanpercentile(np.r_[key_values, nonkey_values], 99))
        per_pair_info.append({"key": key, "name": cname, "pval": pval, "y99": y99})

    # --- Combined p-values across complexes ---
    if all_pvals:
        z_st, p_st = combine_pvalues(all_pvals, method="stouffer")
        chi_f, p_f = combine_pvalues(all_pvals, method="fisher")
        print(f"\nStouffer combined: Z={z_st:.4f}, p={p_st:.4e}")
        print(f"Fisher  combined: χ²={chi_f:.4f}, p={p_f:.4e}")

    # --- Overview figure using cached data (no re-run of experiments) ---
    stouffer_tuple = (z_st, p_st) if all_pvals else None
    fisher_tuple   = (chi_f, p_f) if all_pvals else None

    if rows:
        df_long = pd.DataFrame(rows)
        fname = "correlation_boxplot_overview_withF.png" if args.include_non_family else "correlation_boxplot_overview.png"
        plot_overview_boxplots_cached(
            df_long,
            per_pair_info,
            output_dir=args.output_dir,
            filename=fname,
            x_by="Name",                   
            stouffer=stouffer_tuple,
            fisher=fisher_tuple
        )

if __name__ == "__main__":
    main()
