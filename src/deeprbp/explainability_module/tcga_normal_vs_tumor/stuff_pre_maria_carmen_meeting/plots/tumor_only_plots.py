
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/tumor_only_plots.py

# REPASA ESTE CODE!

import os
from typing import Optional, Sequence, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------------------------------------
# Volcano-like per tumor-only RBP
# -------------------------------------------------------------

# Input: df_tx_rank_abs (o similar) filtrado a tumor-only
# x = Score
# y = abs(Score)
# punto más grande / color distinto para los top-k más extremos
# anotación de Transcript_name para top_k

def plot_tumor_only_volcano_for_rbp(
    df_tx_rank: pd.DataFrame,
    rbp_id: str,
    rbp_name: Optional[str] = None,
    score_col: str = "Score",
    transcript_name_col: str = "Transcript_name",
    top_k_annot: int = 10,
    out_path: Optional[str] = None,
    show: bool = False,
) -> None:
    """
    Volcano-like plot for a single tumor-only RBP.

    x-axis: raw DeepRBP score (can be positive or negative)
    y-axis: |score| (magnitude)
    - Points with higher |score| are more "extreme" targets.
    - Top-k most extreme targets are highlighted and annotated.

    Parameters
    ----------
    df_tx_rank : pd.DataFrame
        Table with at least columns:
        - 'RBP_ID'
        - score_col (default 'Score')
        - transcript_name_col (default 'Transcript_name')
    rbp_id : str
        RBP ID (Ensembl) to plot.
    rbp_name : str, optional
        Human-readable RBP symbol to show in title.
    score_col : str
        Column name for the score.
    transcript_name_col : str
        Column with transcript display names.
    top_k_annot : int
        Number of top extreme transcripts to annotate.
    out_path : str, optional
        If provided, saves the figure to this path.
    show : bool
        If True, calls plt.show(); otherwise closes the figure.
    """
    sub = df_tx_rank[df_tx_rank["RBP_ID"] == rbp_id].copy()
    if sub.empty:
        print(f"[tumor-only volcano] No rows found for RBP {rbp_id}. Skipping.")
        return

    rbp_label = rbp_name or rbp_id

    scores = sub[score_col].astype(float).values
    mags = np.abs(scores)

    sub["abs_score"] = mags

    # sort by |score| desc to pick top_k
    sub_sorted = sub.sort_values("abs_score", ascending=False)
    top = sub_sorted.head(top_k_annot)

    plt.figure(figsize=(6, 5))
    # all points
    plt.scatter(scores, mags, alpha=0.4, s=10)

    # highlight top_k
    plt.scatter(
        top[score_col].values,
        top["abs_score"].values,
        alpha=0.9,
        s=30,
        edgecolor="k",
    )

    # annotate top_k
    for _, row in top.iterrows():
        tx_name = str(row.get(transcript_name_col, row.get("Transcript_ID", "")))
        plt.text(
            row[score_col],
            row["abs_score"],
            tx_name,
            fontsize=6,
            ha="left",
            va="bottom",
        )

    plt.axvline(0.0, linestyle="--", linewidth=0.8)
    plt.xlabel("DeepRBP score (tumor-only)")
    plt.ylabel("|DeepRBP score|")
    plt.title(f"Tumor-only targets for {rbp_label}")

    plt.tight_layout()
    if out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"[tumor-only volcano] Saved to: {out_path}")
    plt.close()

# -------------------------------------------------------------
# Heatmap RBP × top target transcripts (tumor-only)
# -------------------------------------------------------------

# Partimos de df_tumor_only (la tabla larga).
# Para cada RBP: cogemos top_n por |Score| (o aplicamos threshold).
# Construimos una matriz transcript × RBP con los scores (o |score|).
# Dibujamos un heatmap sencillo con imshow.

def build_tumor_only_heatmap_matrix(
    df_tumor_only: pd.DataFrame,
    score_col: str = "Score",
    top_n_per_rbp: int = 30,
    use_abs: bool = True,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """
    Build a matrix (transcript x RBP) of scores for tumor-only RBPs.

    For each RBP, select top_n_per_rbp transcripts by |score| and
    pivot to a transcript x RBP matrix.

    Parameters
    ----------
    df_tumor_only : pd.DataFrame
        Tumor-only table with columns:
        - 'RBP_ID', 'Transcript_ID', 'Transcript_name', score_col
    score_col : str
        Column with the DeepRBP score.
    top_n_per_rbp : int
        Number of top transcripts per RBP to retain.
    use_abs : bool
        If True, keep the sign in the matrix but select transcripts by |score|.
        The pivot matrix will contain the raw score values.

    Returns
    -------
    mat : pd.DataFrame
        index: Transcript_ID (or name)
        columns: RBP_ID
        values: score_col
    transcripts_order : list of str
        Order of transcript IDs (rows).
    rbps_order : list of str
        Order of RBP IDs (columns).
    """
    df = df_tumor_only.copy()
    df["abs_score"] = df[score_col].astype(float).abs()

    # pick top N per RBP by |score|
    sub_list = []
    for rbp_id, grp in df.groupby("RBP_ID"):
        sub = grp.sort_values("abs_score", ascending=False).head(top_n_per_rbp)
        sub_list.append(sub)

    if not sub_list:
        return pd.DataFrame(), [], []

    top_df = pd.concat(sub_list, axis=0)

    # pivot to matrix
    mat = top_df.pivot_table(
        index="Transcript_ID",
        columns="RBP_ID",
        values=score_col,
        aggfunc="mean",
    )

    mat = mat.fillna(0.0)

    transcripts_order = mat.index.tolist()
    rbps_order = mat.columns.tolist()
    return mat, transcripts_order, rbps_order


def plot_tumor_only_heatmap(
    mat: pd.DataFrame,
    transcripts_order: Optional[Sequence[str]] = None,
    rbps_order: Optional[Sequence[str]] = None,
    out_path: Optional[str] = None,
    show: bool = False,
    title: str = "Tumor-only RBP × transcript scores",
) -> None:
    """
    Plot a heatmap for the given transcript × RBP matrix.

    Parameters
    ----------
    mat : pd.DataFrame
        transcript × RBP matrix, as returned by build_tumor_only_heatmap_matrix.
    transcripts_order : Sequence[str], optional
        Row ordering; if None, uses mat.index.
    rbps_order : Sequence[str], optional
        Column ordering; if None, uses mat.columns.
    out_path : str, optional
        Where to save the figure.
    show : bool
        If True, plt.show().
    title : str
        Figure title.
    """
    if mat.empty:
        print("[tumor-only heatmap] Empty matrix. Nothing to plot.")
        return

    if transcripts_order is not None:
        mat = mat.loc[list(transcripts_order)]
    if rbps_order is not None:
        mat = mat[rbps_order]

    data = mat.values

    plt.figure(figsize=(max(6, 0.25 * mat.shape[1]), max(6, 0.25 * mat.shape[0])))
    im = plt.imshow(data, aspect="auto", interpolation="nearest")

    plt.colorbar(im, label="DeepRBP score (tumor-only)")

    plt.xlabel("RBP (tumor-only)")
    plt.ylabel("Transcript")

    # Ticks (subsample if too many)
    if mat.shape[1] <= 40:
        plt.xticks(
            ticks=np.arange(mat.shape[1]),
            labels=mat.columns.tolist(),
            rotation=90,
            fontsize=6,
        )
    else:
        plt.xticks([])

    if mat.shape[0] <= 40:
        plt.yticks(
            ticks=np.arange(mat.shape[0]),
            labels=mat.index.tolist(),
            fontsize=6,
        )
    else:
        plt.yticks([])

    plt.title(title)
    plt.tight_layout()

    if out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"[tumor-only heatmap] Saved to: {out_path}")

    plt.close()

# -------------------------------------------------------------
# Barplots per RBP (top 10 transcripts)
# -------------------------------------------------------------

# Para cada RBP tumor-only, cogemos top 10 transcritos por |Score|
# barh, eje Y = Transcript_name, eje X = |Score|

def plot_tumor_only_barplot_for_rbp(
    df_tumor_only: pd.DataFrame,
    rbp_id: str,
    rbp_name: Optional[str] = None,
    score_col: str = "Score",
    transcript_name_col: str = "Transcript_name",
    top_n: int = 10,
    out_path: Optional[str] = None,
    show: bool = False,
) -> None:
    """
    Horizontal barplot of top-N transcripts (by |score|) for a tumor-only RBP.

    y-axis: transcript names
    x-axis: |score|

    Parameters
    ----------
    df_tumor_only : pd.DataFrame
        Must contain RBP_ID, score_col, transcript_name_col.
    rbp_id : str
        RBP ID to plot.
    rbp_name : str, optional
        Symbol to display in title.
    score_col : str
        Column with DeepRBP scores.
    transcript_name_col : str
        Column with transcript names.
    top_n : int
        Number of transcripts to show.
    out_path : str, optional
        Where to save the figure.
    show : bool
        If True, show figure.
    """
    sub = df_tumor_only[df_tumor_only["RBP_ID"] == rbp_id].copy()
    if sub.empty:
        print(f"[tumor-only barplot] No rows for RBP {rbp_id}. Skipping.")
        return

    rbp_label = rbp_name or rbp_id

    sub["abs_score"] = sub[score_col].astype(float).abs()
    sub_top = sub.sort_values("abs_score", ascending=False).head(top_n)

    # order so biggest at top
    sub_top = sub_top.iloc[::-1]

    plt.figure(figsize=(6, 4))
    plt.barh(
        sub_top[transcript_name_col].astype(str).values,
        sub_top["abs_score"].values,
    )
    plt.xlabel("|DeepRBP score| (tumor-only)")
    plt.ylabel("Transcript")
    plt.title(f"Top {top_n} targets of {rbp_label} (tumor-only)")

    plt.tight_layout()
    if out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"[tumor-only barplot] Saved to: {out_path}")

    plt.close()

# -------------------------------------------------------------
# Biotype distribution among tumor-only top targets
# -------------------------------------------------------------
# - Tomar top_k por |Score| para todos los RBPs tumor-only.
# - Mirar la distribución de Transcript_biotype en ese subconjunto.
# - Hacer un barplot de conteos o proporciones por biotype.

def plot_tumor_only_biotype_distribution(
    df_tumor_only: pd.DataFrame,
    score_col: str = "Score",
    biotype_col: str = "Transcript_biotype",
    top_n_per_rbp: int = 50,
    out_path: Optional[str] = None,
    show: bool = False,
    normalize: bool = True,
) -> None:
    """
    Plot the distribution of transcript biotypes among top tumor-only targets.

    For each tumor-only RBP, we select top_n_per_rbp transcripts by |score|,
    then aggregate their biotypes and plot global counts or proportions.

    Parameters
    ----------
    df_tumor_only : pd.DataFrame
        Table for tumor-only RBPs, with columns:
        - 'RBP_ID', score_col, biotype_col
    score_col : str
        Column with DeepRBP scores.
    biotype_col : str
        Column with transcript biotype (e.g., 'protein_coding', 'retained_intron').
    top_n_per_rbp : int
        Number of top transcripts per RBP to consider.
    out_path : str, optional
        Where to save the plot.
    show : bool
        If True, plt.show().
    normalize : bool
        If True, plot fractions (proportions). If False, raw counts.
    """
    df = df_tumor_only.copy()
    df["abs_score"] = df[score_col].astype(float).abs()

    sub_list = []
    for rbp_id, grp in df.groupby("RBP_ID"):
        sub = grp.sort_values("abs_score", ascending=False).head(top_n_per_rbp)
        sub_list.append(sub)

    if not sub_list:
        print("[tumor-only biotypes] No data to compute distribution.")
        return

    top_df = pd.concat(sub_list, axis=0)

    counts = top_df[biotype_col].value_counts().sort_values(ascending=False)
    if normalize:
        values = (counts / counts.sum()).values
        ylabel = "Fraction of top targets"
    else:
        values = counts.values
        ylabel = "Number of top targets"

    plt.figure(figsize=(6, 4))
    plt.bar(counts.index.astype(str).values, values)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel(ylabel)
    plt.title("Transcript biotypes among tumor-only RBP targets")

    plt.tight_layout()
    if out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"[tumor-only biotypes] Saved to: {out_path}")

    plt.close()
