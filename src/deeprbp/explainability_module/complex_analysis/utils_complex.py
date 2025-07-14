# src/deeprbp/explainability_module/corum_complex_analysis/utils_complex.py

import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np

from .plots_complex import plot_upper_triangle_corr_matrix, plot_correlation_boxplot

def split_to_list(cell):
    """
    Split a semicolon-separated string into a list of trimmed strings.
    Returns an empty list if input is NaN.
    """
    if pd.isna(cell):
        return []
    return [item.strip() for item in cell.split(';')]

def map_gene_names_to_ids(row, getBM):
    """
    Map each gene name in 'subunits_gene_name' to its corresponding Gene_ID from getBM.
    If no direct match is found, try matching using each synonym (splitting by comma).
    Returns a list of Gene_IDs or None if no match is found.
    """
    gene_ids = []
    print(f"Processing complex: {row.get('complex_name', 'Unknown')}")
    for gene_name, synonym_str in zip(row['subunits_gene_name'], row['subunits_gene_name_synonyms']):
        print(f"  Trying gene name: {gene_name}")
        # Buscar match directo
        match = getBM[getBM['Gene_name'] == gene_name]
        if not match.empty:
            gene_id = match['Gene_ID'].values[0]
            print(f"    Found Gene_ID: {gene_id} for gene name: {gene_name}\n")
            gene_ids.append(gene_id)
        else:
            print(f"    Gene name '{gene_name}' not found, trying synonyms: {synonym_str}")
            # Probar cada sinónimo separado por coma
            synonyms = [s.strip() for s in synonym_str.split(',')]
            found = False
            for syn in synonyms:
                match_syn = getBM[getBM['Gene_name'] == syn]
                if not match_syn.empty:
                    gene_id = match_syn['Gene_ID'].values[0]
                    print(f"    Found Gene_ID: {gene_id} for synonym: {syn}\n")
                    gene_ids.append(gene_id)
                    found = True
                    break
            if not found:
                print(f"    No Gene_ID found for gene name '{gene_name}' or any synonym\n")
                gene_ids.append(None)
    return gene_ids

def analyze_single_complex(corum_data, complex_idx, corr_scores, output_dir=None):
    """
    Analyze correlation patterns within a single protein complex compared to RBPs outside the complex.

    Parameters:
    - corum_data: pd.DataFrame
        DataFrame containing CORUM complexes info with a 'subunits_gene_id' column listing RBPs.
    - complex_idx: int or label
        Index identifying the complex to analyze in corum_data.
    - corr_scores: pd.DataFrame
        Correlation matrix (genes x RBPs) of explainability scores.
    - output_dir: str or None
            Directory path where to save generated plots (heatmap, boxplot). If None, plots are not saved.

    Returns:
    - stat: float
        Mann-Whitney U test statistic comparing correlations inside vs outside the complex.
    - pval: float
        P-value from the Mann-Whitney U test (alternative='greater').
    """
    # Get the RBPs contained in the selected complex
    complex_rbps = corum_data.loc[complex_idx, 'subunits_gene_id']
    complex_rbps = [rbp for rbp in complex_rbps if rbp in corr_scores.columns]
    print(f"Number of RBPs in complex: {len(complex_rbps)}")
    # RBPs outside complex
    all_rbps = set(corr_scores.columns)
    outside_rbps = list(all_rbps - set(complex_rbps))
    # Reorder rows and columns in correlation matrix: complex RBPs first (top-left)
    corr_scores_ordered = corr_scores.loc[complex_rbps + outside_rbps, complex_rbps + outside_rbps]
    # Pheatmap of the reordered correlation matrix
    if output_dir is not None:
        plot_upper_triangle_corr_matrix(corr_scores_ordered, complex_idx, output_dir)
    # Submatrices: complex vs complex and outside vs outside
    complex_corrs = corr_scores_ordered.loc[complex_rbps, complex_rbps]
    outside_corrs = corr_scores_ordered.loc[outside_rbps, outside_rbps]
    # Flatten upper triangles (excluding diagonal)
    flat_complex = get_upper_triangle_flattened(complex_corrs)
    flat_outside = get_upper_triangle_flattened(outside_corrs)
    # Mann-Whitney U test
    stat, pval = mannwhitneyu(flat_complex, flat_outside, alternative='greater')
    print(f"Mann-Whitney U test statistic: {stat:.4f}, p-value: {pval:.4e}")
    if pval < 0.05:
        print("💥 Correlations within the complex are significantly higher.")
    else:
        print("🫠 No sufficient evidence that correlations within the complex are higher.")
    # Boxplot
    if output_dir is not None:
        plot_correlation_boxplot(flat_complex, flat_outside, complex_idx, output_dir)
    return stat, pval

def get_upper_triangle_flattened(corr_df):
    mask = np.triu(np.ones(corr_df.shape), k=1).astype(bool)
    return corr_df.where(mask).stack().values  # flatten sin NaNs

# ver intersecciones entre complejos
def check_complex_intersections(df):
    """
    Check for gene overlaps between complexes in the dataframe.
    
    For each pair of complexes, prints whether they share genes and which ones.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing a column 'subunits_gene_id' with lists of gene IDs per complex.
    """
    for i in range(len(df)):
        set_i = set(df.loc[i, 'subunits_gene_id'])
        for j in range(i+1, len(df)):
            set_j = set(df.loc[j, 'subunits_gene_id'])
            intersec = set_i.intersection(set_j)
            if intersec:
                print(f"Complex {i} and Complex {j} share {len(intersec)} gene(s): {intersec}")
            else:
                print(f"Complex {i} and Complex {j} share no genes.")

