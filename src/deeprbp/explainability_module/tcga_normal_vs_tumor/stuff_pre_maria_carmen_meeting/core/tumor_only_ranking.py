
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/tumor_only_ranking.py

import pandas as pd

def rank_tumor_per_rbp(
    df: pd.DataFrame,
    top_n: int = 100,
    mode: str = "pos") -> pd.DataFrame:
    """
    Rank transcripts per RBP based on the 'Score' column.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain: ['RBP_ID', 'RBP_name',
                    'Transcript_ID', 'Transcript_name',
                    'Gene_ID', 'Gene_name', 'Score'].
    top_n : int
        Number of top cases to keep per RBP.
    mode : {'pos', 'neg', 'abs'}
        'pos' -> highest positive scores
        'neg' -> most negative scores
        'abs' -> largest absolute scores

    Returns
    -------
    pd.DataFrame
        Ranked table with 'rank_within_rbp' column.
    """
    df = df.copy()
    if mode == "abs":
        df["Score_for_rank"] = df["Score"].abs()
        ascending = False
    elif mode == "neg":
        df["Score_for_rank"] = df["Score"]
        ascending = True      # más negativos primero
    elif mode == "pos":
        df["Score_for_rank"] = df["Score"]
        ascending = False     # más positivos primero
    else:
        raise ValueError("mode must be one of: 'pos', 'neg', 'abs'")
    df["rank_within_rbp"] = (
        df.groupby("RBP_ID")["Score_for_rank"]
        .rank(method="first", ascending=ascending)
        .astype(int)
    )
    df = df[df["rank_within_rbp"] <= top_n]
    return df.sort_values(["RBP_ID", "rank_within_rbp"])
