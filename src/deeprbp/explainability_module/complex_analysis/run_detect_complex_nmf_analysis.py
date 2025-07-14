# /scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/complex_analysis/run_detect_complex_nmf_analysis.py

# Negative Matrix Factorization (NMF) is a technique that decomposes a matrix into two smaller matrices with non-negativity 
# constraints, useful for uncovering latent patterns. In this case, it allows us to discover "components" (which we interpret as complexes) 
# from the gene vs RBP scores.

import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import NMF
import os
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "This script performs Non-negative Matrix Factorization (NMF) analysis "
            "to detect protein complexes from gene vs RBP explainability scores. "
            "It takes as input a CORUM complexes data file and a CSV file containing "
            "gene-RBP score matrices, runs NMF for varying numbers of components, "
            "and saves heatmaps of component loadings and reconstruction error plots "
            "to the specified output directory."
        )
    )
    parser.add_argument('--corum_file_path', type=str,
                        default='/scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/complex_analysis/corum_results.txt',
                        help='Path to the CORUM complexes data file.')
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
    n_complex = df_corum.shape[0]

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

    # perform NFM analysis to detectar los complejos de corum?? 
    nmf_complex_analysis(
            scores_GxRBP, df_corum, n_complex_max=30,
            output_dir=args.output_dir)

def nmf_complex_analysis(scores_GxRBP, df_corum, n_complex_max=None, output_dir=None):
    """
    Performs NMF analysis varying the number of complexes (components) starting from
    n_complex = number of rows in df_corum (known complexes).
    
    Saves heatmaps of the H matrix and a plot of the reconstruction error for each n_complex.
    
    Args:
        scores_GxRBP (pd.DataFrame): gene x RBP score matrix (non-negative).
        df_corum (pd.DataFrame): dataframe with complex info (used to determine initial n_complex).
        n_complex_max (int, optional): maximum number of components to test. If None, defaults to n_complex + 5.
        output_dir (str, optional): directory to save images and results.
    """
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    # Starting number of components
    n_complex_start = df_corum.shape[0]
    n_complex_max = n_complex_max or (n_complex_start + 5)

    data = scores_GxRBP.fillna(0).clip(lower=0).values

    reconstruction_errors = []
    tested_components = list(range(n_complex_start, n_complex_max + 1))

    for n in tested_components:
        print(f"Running NMF with n_components={n} ...\n")
        model = NMF(n_components=n, init='nndsvda', random_state=42, max_iter=500)
        W = model.fit_transform(data)
        H = model.components_
        reconstruction_errors.append(model.reconstruction_err_)
        # Save heatmap of H for each n_components
        plt.figure(figsize=(10, 6))
        sns.heatmap(H, cmap='viridis', cbar=True)
        plt.title(f"NMF Components vs RBPs (H matrix), n_components={n}")
        plt.xlabel("RBPs")
        plt.ylabel("Components (Complexes)")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"H_matrix_n{n}.png"))
        plt.close()

    # Save and show reconstruction error trend
    plt.figure(figsize=(8,5))
    plt.plot(tested_components, reconstruction_errors, marker='o')
    plt.title('NMF Reconstruction Error vs Number of Components')
    plt.xlabel('Number of Components (n_complex)')
    plt.ylabel('Reconstruction Error')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "nmf_reconstruction_error.png"))
    plt.show()

    print("\nSummary of reconstruction errors:")
    for n, err in zip(tested_components, reconstruction_errors):
        print(f"n_components={n} --> reconstruction error: {err:.4f}")

if __name__ == "__main__":
    main()