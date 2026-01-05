# src/deeprbp/explainability_module/postar_validation/plot_utils.py

import os, re
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from matplotlib.offsetbox import AnchoredText

def _sanitize(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')

def _fd_bins(x: np.ndarray, lo: float, hi: float, max_bins: int = 100) -> np.ndarray:
    """Freedman–Diaconis rule, acotado a [lo, hi] y con límite superior de bins."""
    x = x[(~np.isnan(x)) & (x >= lo) & (x <= hi)]
    n = x.size
    if n == 0:
        return np.linspace(lo, hi, 11)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    if iqr <= 0:
        k = int(np.sqrt(n))
    else:
        h = 2 * iqr * n ** (-1 / 3)
        k = int(np.clip(np.ceil((hi - lo) / max(h, 1e-12)), 5, max_bins))
    return np.linspace(lo, hi, k + 1)

def plot_distributions_and_roc_with_thresholds(
    df_current_rbp,
    rbp_id,
    optimal_threshold,
    fpr,
    tpr,
    optimal_idx,
    auc_score,
    path_save,
    rbp_display_name: Optional[str] = None,
    colors: Optional[dict] = None,
    score_label: str = "Explainability score",
    q_hi: float = 0.98,
    bw_adjust: float = 1.3,
    show_density_yticks: bool = False,
    stats_label: Optional[str] = None,   # e.g. "p_adj: 1.1e-12 ***"
):
    """
    KDE + ROC on white background. Minimal legends (classes only). Threshold line with a small label
    placed at the same x but near the bottom of the KDE panel. Metrics box (AUC, TPR, FPR) anchored
    to the bottom-right of the ROC.
    """
    required = {"Score", "Postar_Score"}
    if not required.issubset(df_current_rbp.columns):
        raise ValueError(f"DataFrame must contain columns: {sorted(required)}")

    # Pretty name
    if rbp_display_name is None:
        if "RBP_name" in df_current_rbp.columns:
            names = df_current_rbp["RBP_name"].dropna().astype(str)
            rbp_display_name = names.mode().iat[0] if not names.empty else rbp_id
        else:
            rbp_display_name = rbp_id

    # Colors
    if colors is None:
        colors = {1: "#7fc97f", 0: "#beaed4"}

    # Data
    scores = np.asarray(df_current_rbp["Score"], dtype=float)
    m1 = (df_current_rbp["Postar_Score"] == 1).to_numpy(bool)
    m0 = (df_current_rbp["Postar_Score"] == 0).to_numpy(bool)
    n1, n0 = int(m1.sum()), int(m0.sum())

    # X range
    lo = 0.0
    hi_all = float(np.nanmax(scores) if scores.size else 1.0)
    hi_main = float(np.nanquantile(scores, q_hi)) if np.isfinite(hi_all) else hi_all
    if np.isfinite(optimal_threshold):
        hi_main = max(hi_main, float(optimal_threshold) * 1.05)
    hi_main = min(hi_main, max(hi_all, 1e-9))

    # Medians
    med1 = float(np.nanmedian(scores[m1])) if n1 > 0 else np.nan
    med0 = float(np.nanmedian(scores[m0])) if n0 > 0 else np.nan

    # Style / sizes
    SUPTITLE_FZ = 12.5
    TITLE_FZ    = 10.5
    LABEL_FZ    = 9.5
    TICK_FZ     = 8.5
    LEGEND_FZ   = 8.5
    KDE_ALPHA   = 0.35
    TH_LW       = 1.4
    MED_LW      = 0.9
    ROC_LW      = 1.6
    RAND_LW     = 1.0
    DOT_SIZE    = 30

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.2), constrained_layout=True)
    fig.patch.set_facecolor("white")
    for ax in axes:
        ax.set_facecolor("white")
        ax.grid(False)

    # ----- KDE panel -----
    ax = axes[0]
    labels = {1: f"Binding (n={n1})", 0: f"Not-Binding (n={n0})"}

    for cls in (1, 0):
        mask = m1 if cls == 1 else m0
        x = np.clip(scores[mask], lo, hi_main)
        if x.size == 0:
            continue
        sns.kdeplot(
            x=x,
            ax=ax,
            fill=True,
            alpha=KDE_ALPHA,
            clip=(lo, hi_main),
            cut=0,
            bw_adjust=bw_adjust,
            color=colors[cls],
            label=labels[cls],
            lw=1.0,
        )

    # Threshold line + small grey label at bottom
    if np.isfinite(optimal_threshold):
        ax.axvline(optimal_threshold, color="red", linestyle="--", linewidth=TH_LW)
        y_bottom = ax.get_ylim()[0]
        ax.annotate(
            f"Th={optimal_threshold:.2f}",
            xy=(optimal_threshold, y_bottom),
            xytext=(3, 6),  # offset upward/right from the bottom edge
            textcoords="offset points",
            ha="left", va="bottom",
            fontsize=9, color="dimgray",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.6),
            clip_on=False,  # so it doesn't get cut off if it sits on the edge
        )

    # Class medians
    if np.isfinite(med1):
        ax.axvline(med1, color=colors[1], linestyle="-", linewidth=MED_LW, alpha=0.9)
    if np.isfinite(med0):
        ax.axvline(med0, color=colors[0], linestyle="-", linewidth=MED_LW, alpha=0.9)

    ax.set_xlim(lo, hi_main)
    ax.set_xlabel(score_label, fontsize=LABEL_FZ)
    ax.set_ylabel("Density (KDE)", fontsize=LABEL_FZ)
    ax.set_title("Score distribution by POSTAR label", fontsize=TITLE_FZ)
    ax.tick_params(axis="both", labelsize=TICK_FZ)
    if not show_density_yticks:
        ax.set_yticks([])

    # Legend with only classes
    leg = ax.legend(title="Postar", frameon=True)
    leg.get_title().set_fontsize(LEGEND_FZ)
    for txt in leg.get_texts():
        txt.set_fontsize(LEGEND_FZ)

    # Stats textbox (top-left)
    if stats_label:
        ax.text(
            0.02, 0.95, stats_label,
            transform=ax.transAxes, ha="left", va="top",
            fontsize=LEGEND_FZ,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none")
        )

    # ----- ROC panel -----
    ax = axes[1]
    ax.plot(fpr, tpr, lw=ROC_LW)
    if 0 <= optimal_idx < len(fpr):
        ax.scatter(float(fpr[optimal_idx]), float(tpr[optimal_idx]),
                   s=DOT_SIZE, color="red", zorder=3)
        tpr_opt = float(tpr[optimal_idx])
        fpr_opt = float(fpr[optimal_idx])
    else:
        tpr_opt = np.nan
        fpr_opt = np.nan

    ax.plot([0, 1], [0, 1], "--", color="gray", lw=RAND_LW)
    ax.set_xlabel("False Positive Rate", fontsize=LABEL_FZ)
    ax.set_ylabel("True Positive Rate", fontsize=LABEL_FZ)
    ax.set_title("ROC curve", fontsize=TITLE_FZ)
    ax.tick_params(axis="both", labelsize=TICK_FZ)

    # Anchored metrics box at bottom-right
    metrics_text = f"AUC = {auc_score:.2f}\nTPR = {tpr_opt:.2f}\nFPR = {fpr_opt:.2f}"
    at = AnchoredText(metrics_text, loc="lower right",
                      prop=dict(size=LEGEND_FZ), frameon=True, borderpad=0.4)
    at.patch.set_alpha(0.85)
    at.patch.set_edgecolor("none")
    ax.add_artist(at)

    # Supertitle & save
    fig.suptitle(f"RBP: {rbp_display_name} ({rbp_id})", fontsize=SUPTITLE_FZ)
    os.makedirs(path_save, exist_ok=True)
    out_path = os.path.join(path_save, f"figure_{_sanitize(rbp_display_name)}.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
