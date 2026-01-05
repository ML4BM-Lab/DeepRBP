
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/rbp_tx_wilcoxon.py

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def compute_rbp_transcript_wilcoxon(
    GxN_t: pd.DataFrame,
    GxT_t: pd.DataFrame,
    verbose: bool = True,
    show_examples: bool = True,
    example_rows: int = 5) -> pd.DataFrame:
    """
    For each Transcript x RBP pair, collects all per-sample scores in
    NORMAL vs TUMOR and runs a Mann-Whitney U test.

    In addition to p-values, it also reports group medians and a
    signed log2 fold-change based on median score magnitudes.

    Also reports:
    - n_normal, n_tumor
    - median_normal
    - median_tumor
    - median_diff (tumor - normal)
    - log2_fc (signed, based on median magnitudes)

    Parameters
    ----------
    GxN_t : pd.DataFrame
        Scores in normal samples. Columns = MultiIndex (sample_id, RBP)
    GxT_t : pd.DataFrame
        Scores in tumor samples. Columns = MultiIndex (sample_id, RBP)
    verbose : bool
        If True, prints progress and examples of matrices.
    show_examples : bool
        If True, print N_rbp and T_rbp (head only) for each RBP.
    example_rows : int
        Number of transcript rows to print for preview.

    Expected format
    ---------------
    GxN_t, GxT_t:
        - index: Transcript_ID
        - columns: MultiIndex (sample_id, RBP_ID)
          e.g. ('S0','ENSG00000188976'), ('S1','ENSG00000188976'), ...

    Returns
    -------
    pd.DataFrame
        Columns:
        Transcript_ID | RBP | n_normal | n_tumor |
        median_normal | median_tumor | median_diff | log2_fc |
        p_value | p_adj_fdr
    """
    if not isinstance(GxN_t.columns, pd.MultiIndex):
        raise ValueError("GxN_t columns must be MultiIndex.")
    if not isinstance(GxT_t.columns, pd.MultiIndex):
        raise ValueError("GxT_t columns must be MultiIndex.")

    if not GxN_t.index.equals(GxT_t.index):
        raise ValueError("Transcript indices differ between NORMAL and TUMOR matrices.")

    rbps = sorted(
        set(GxN_t.columns.get_level_values("RBP"))
        & set(GxT_t.columns.get_level_values("RBP"))
    )

    transcripts = GxN_t.index
    results = []

    for i_rbp, rbp in enumerate(rbps):
        if verbose and i_rbp % 20 == 0:
            print(f"[wilcoxon] Processing RBP {rbp} ({i_rbp+1}/{len(rbps)})")

        # extract per-sample column submatrix for this rbp
        N_rbp = GxN_t.xs(rbp, axis=1, level="RBP")
        T_rbp = GxT_t.xs(rbp, axis=1, level="RBP")

        # Show example rows
        if verbose and show_examples:
            print(f"[wilcoxon] Example NORMAL matrix for {rbp}:")
            print(N_rbp.head(example_rows))
            print(f"[wilcoxon] Example TUMOR matrix for {rbp}:")
            print(T_rbp.head(example_rows))
        
        # Compare transcript-by-transcript
        for tx in transcripts:
            x_raw = N_rbp.loc[tx].values.astype(float)
            y_raw = T_rbp.loc[tx].values.astype(float)

            x = x_raw[~np.isnan(x_raw)]
            y = y_raw[~np.isnan(y_raw)]

            n_norm = int(x.size)
            n_tum  = int(y.size)

            median_normal = float(np.median(x)) if n_norm > 0 else np.nan
            median_tumor, median_diff, log2_fc = compute_log2fc_from_groups(x, y)

            pval = wilcoxon_two_groups(x, y)

            results.append(
                {
                    "Transcript_ID": tx,
                    "RBP": rbp,
                    "n_normal": n_norm,
                    "n_tumor": n_tum,
                    "median_normal": median_normal,
                    "median_tumor": median_tumor,
                    "median_diff": median_diff,
                    "log2_fc": log2_fc,
                    "p_value": pval,
                }
            )

    df = pd.DataFrame(results)
    df["p_adj_fdr"] = multipletests(df["p_value"].fillna(1.0), method="fdr_bh")[1] # FDR correction
    return df
    
def compute_log2fc_from_groups(x, y, eps: float = 1e-8):
    """
    Compute log2 fold-change summary between two groups of scores.

    Assumes x (Normal) and y (Tumor) are 1D NumPy arrays **without NaNs**.

    Parameters
    ----------
    x : np.ndarray
        Scores for NORMAL samples (NaN-free).
    y : np.ndarray
        Scores for TUMOR samples (NaN-free).
    eps : float
        Small constant to stabilize the ratio when medians are ~0.

    Returns
    -------
    median_tumor : float
        Median score in the tumor group.
    median_diff : float
        median_tumor - median_normal.
    log2_fc : float
        Signed log2 fold-change of median magnitudes.
        Positive → tumor has larger magnitude than normal.
        Negative → tumor has smaller magnitude than normal.
        NaN if medians are not finite.
    """
    if x.size == 0 or y.size == 0:
        return np.nan, np.nan, np.nan

    median_normal = float(np.median(x))
    median_tumor = float(np.median(y))

    if not (np.isfinite(median_normal) and np.isfinite(median_tumor)):
        return median_tumor, np.nan, np.nan

    median_diff = median_tumor - median_normal

    # log2 FC based on magnitude, sign from median_diff
    fc = (np.abs(median_tumor) + eps) / (np.abs(median_normal) + eps)
    log2_fc = float(np.log2(fc) * np.sign(median_diff))

    return median_tumor, median_diff, log2_fc

def wilcoxon_two_groups(x, y):
    """
    Mann–Whitney U test (two-sided) on two NaN-free arrays.
    Parameters
    ----------
    x : array-like
        Values for group 1 (Normal)
    y : array-like
        Values for group 2 (Tumor)

    Returns
    -------
    float
        p-value of the rank-sum test
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.size < 2 or y.size < 2:
        return np.nan
    try:
        # SciPy >= 1.7: soporta `method=`
        stat, p = mannwhitneyu(
            x, y,
            alternative="two-sided",
            method="asymptotic"
        )
        return float(p)
    except TypeError:
        # Fallback para SciPy más viejas (sin `method`)
        stat, p = mannwhitneyu(
            x, y,
            alternative="two-sided"
        )
        return float(p)
    except Exception as e:
        print("[wilcoxon] error:", e)
        return np.nan


def summarize_rbp_from_transcript_tests(
    df_rbptx: pd.DataFrame,
    fdr_thr: float = 0.01
) -> pd.DataFrame:
    """
    Aggregate transcript-level tests into an RBP-level summary.

    For each RBP, reports:
      - n_transcripts
      - n_sig_tx_FDR_0_01
      - prop_sig_tx_FDR_0_01
      - median_log2_fc
      - median_abs_log2_fc
      - min_p_adj_fdr
    """
    if df_rbptx.empty:
        return pd.DataFrame(
            columns=[
                "RBP",
                "n_transcripts",
                f"n_sig_tx_FDR_{str(fdr_thr).replace('.', '_')}",
                f"prop_sig_tx_FDR_{str(fdr_thr).replace('.', '_')}",
                "median_log2_fc",
                "median_abs_log2_fc",
                "min_p_adj_fdr",
            ]
        )

    def _agg(group: pd.DataFrame) -> pd.Series:
        n_tx = group.shape[0]
        sig_mask = group["p_adj_fdr"] < fdr_thr
        n_sig = int(sig_mask.sum())
        prop_sig = float(sig_mask.mean())
        median_log2_fc = float(group["log2_fc"].median())
        median_abs_log2_fc = float(group["log2_fc"].abs().median())
        min_fdr = float(group["p_adj_fdr"].min())
        return pd.Series(
            dict(
                n_transcripts=n_tx,
                n_sig_tx_FDR_0_01=n_sig,
                prop_sig_tx_FDR_0_01=prop_sig,
                median_log2_fc=median_log2_fc,
                median_abs_log2_fc=median_abs_log2_fc,
                min_p_adj_fdr=min_fdr,
            )
        )

    df_rbp = df_rbptx.groupby("RBP").apply(_agg).reset_index()
    return df_rbp