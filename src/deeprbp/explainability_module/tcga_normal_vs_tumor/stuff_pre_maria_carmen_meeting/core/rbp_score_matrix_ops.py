
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/rbp_score_matrix_filter.py

import pandas as pd
from typing import List

def filter_score_matrix_by_rbps(
    score_df: pd.DataFrame,
    rbp_list: List[str]) -> pd.DataFrame:
    """
    Filter a transcript-by-(sample,RBP) score matrix and retain only the columns
    corresponding to a given set of RBPs.

    The input DataFrame is expected to have a MultiIndex on columns, with level 0
    representing samples and level 1 representing RBPs (e.g. ENSG identifiers).

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing per-transcript score values. Columns must be a 
        MultiIndex of the form (sample_id, rbp_id).
    rbp_list : List[str]
        List of RBP identifiers to retain (e.g. ENSG IDs).

    Returns
    -------
    pd.DataFrame
        A filtered version of `score_df` containing only the columns associated 
        with RBPs in `rbp_list`.
    """
    # Safety: keep only RBPs actually present in the score matrix
    rbps_available = set(score_df.columns.get_level_values("RBP"))
    rbps_to_keep = sorted(rbps_available.intersection(rbp_list))

    if not rbps_to_keep:
        raise ValueError("None of the requested RBPs are present in the score matrix.")

    # Build a filtered matrix using pandas IndexSlice
    idx = pd.IndexSlice
    filtered = score_df.loc[:, idx[:, rbps_to_keep]].copy()

    return filtered
