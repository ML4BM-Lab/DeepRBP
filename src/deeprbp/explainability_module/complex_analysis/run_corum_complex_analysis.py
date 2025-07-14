
import pandas as pd
import argparse
from statsmodels.stats.meta_analysis import combine_pvalues


from .utils_complex import split_to_list, map_gene_names_to_ids, analyze_single_complex

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "This script analyzes correlation patterns across genes using CORUM complex data "
            "and explainability scores of RBPs. It generates correlation heatmaps and performs "
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
    return parser.parse_args()


def main():
    args = parse_args()

    # Load CORUM complex data
    print(f"Loading Corum results from {args.corum_file_path} ...")
    df_corum = pd.read_csv(args.corum_file_path, sep='\t')
    print(df_corum.head())

    # Load gene mapping data (getBM)
    print(f"\nLoading gene mapping (getBM) from {args.getBM_file_path} ...")
    getBM = pd.read_csv(args.getBM_file_path)

    # Convert 'subunits_gene_name' column from string to list
    print("\nConverting 'subunits_gene_name' to lists...")
    df_corum['subunits_gene_name'] = df_corum['subunits_gene_name'].apply(split_to_list)
    # Convert 'subunits_gene_name_synonyms' column from string to list
    print("Converting 'subunits_gene_name_synonyms' to lists...")
    df_corum['subunits_gene_name_synonyms'] = df_corum['subunits_gene_name_synonyms'].apply(split_to_list)

    # Map gene names to gene IDs using getBM mapping
    print("\nMapping gene names to Gene_IDs for each complex...")
    df_corum['subunits_gene_id'] = df_corum.apply(lambda row: map_gene_names_to_ids(row, getBM), axis=1)

    # Load explainability scores and take absolute values
    print("\nLoading explainability scores...")
    print(f"Loading scores from {args.scores_file_path} ...")
    scores_GxRBP = pd.read_csv(args.scores_file_path, index_col=0)
    print("\nTaking absolute value of scores for analysis...")
    scores_GxRBP = abs(scores_GxRBP)
    print("First rows of the loaded scores:")
    print(scores_GxRBP.head())
    print(scores_GxRBP.shape)
    print("Data preparation completed.")

    print("Calculate RBP pearson correlations across genes.")
    # Calculate RBP pearson correlations across genes
    corr_scores = scores_GxRBP.corr() # DataFrame 1348x1348

    # Analyze each complex in corum_data:
    all_stats = []
    all_pvals = []

    for idx in df_corum.index:
        print(f"\nAnalyzing complex index: {idx}")
        stat, pval = analyze_single_complex(df_corum, idx, corr_scores, args.output_dir)
        all_stats.append(stat)
        all_pvals.append(pval)

    # Combine p-values across all complexes (Fisher's method)
    stat_combined, pval_combined = combine_pvalues(all_pvals, method='stouffer')
    print(f"\nCombined Fisher's test statistic: {stat_combined:.4f}, combined p-value: {pval_combined:.4e}")

    if pval_combined < 0.05:
        print("🔥 Overall, correlations within complexes are significantly higher than outside.")
    else:
        print("🫠 Overall, no sufficient evidence that correlations within complexes are higher.")
        
if __name__ == "__main__":
    main()


















 
 