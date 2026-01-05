
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import transforms
import re
from matplotlib import transforms
from typing import Optional, Tuple, List, Dict, Any
from matplotlib.ticker import FixedLocator
from matplotlib import ticker as mticker
from matplotlib.patches import Patch

from .utils_complex import _infer_tumor_from_output_dir

def _p_to_stars(p):
    return "ns" if p >= 0.05 else ("*" if p >= 1e-2 else ("**" if p >= 1e-3 else ("***" if p >= 1e-4 else "****")))

def _annotate_sig(ax, x1, x2, y, text, h=0.015, lw=1.6):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=lw, c="black", clip_on=False)
    ax.text((x1+x2)/2, y+h*1.35, text, ha="center", va="bottom", fontsize=12)

def plot_correlation_boxplot(flat_complex, flat_outside, complex_idx=None, output_dir=None, pval=None, complex_name=None, p_yfrac=0.9):
    """
    Boxplot (no grid) comparing 'Complex' vs 'Outside' and annotating the p-value.

    Parameters
    ----------
    flat_complex : array-like
        Flattened correlation values within the complex (upper triangle, no diagonal).
    flat_outside : array-like
        Flattened correlation values for the complex vs non-complex block.
    complex_idx : str | int | None
        Label used to identify the complex in the filename/title (e.g., 'A', 'B', ...).
    output_dir : str | None
        If provided, the figure is saved there; otherwise it is shown.
    pval : float | None
        If provided, an annotation with significance stars and the p-value is drawn.
    complex_name : str | None
        Human-readable complex name; if given, it is appended to the plot title.
    """
    
    data = pd.DataFrame({
        'Correlation': list(flat_complex) + list(flat_outside),
        'Group': ['Complex'] * len(flat_complex) + ['Outside'] * len(flat_outside)
    })

    # Labels with sample sizes
    labels = [f"Complex (n={len(flat_complex)})", f"Outside (n={len(flat_outside)})"]
    palette = ["#0072B2", "#E69F00"]  # azul, naranja

    # Minimal look, no grid
    sns.set_theme(style="white", rc={"axes.grid": False})
    fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=300)

    ax = sns.boxplot(
        x="Group", y="Correlation", hue="Group", data=data, ax=ax,
        width=0.55, showfliers=False, whis=(5, 95), palette=palette, dodge=False,
        medianprops={"linewidth": 2.2, "color": "black"},
        boxprops={"linewidth": 1.4}, whiskerprops={"linewidth": 1.4}, capprops={"linewidth": 1.4},
    )
    # Quitamos la leyenda automática que aparece al usar hue
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    ax.xaxis.set_major_locator(FixedLocator([0, 1]))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels, fontsize=11)

    # Title
    ttl = complex_name if complex_name else (f"Complex {complex_idx}" if complex_idx is not None else "")
    if ttl:
        ax.set_title(ttl, fontsize=16, pad=6, weight="bold")

    ax.set_xlabel("")
    ax.set_ylabel("Correlation", fontsize=12)
    sns.despine(ax=ax) # Remove top/right spines

    # --- Tight Y limits using robust percentiles + modest headroom for bracket ---
    combined = np.r_[flat_complex, flat_outside].astype(float)
    q1, q99  = np.nanpercentile(combined, [1, 99])
    span     = max(q99 - q1, 1e-3)
    ax.set_ylim(q1 - 0.02*span, q99 + 0.06*span)

    # P-value annotation with significance stars
    if pval is not None:
        stars = _p_to_stars(pval)
        txt   = f"{stars}  (p = {pval:.2e})"
        
        # Bracket at y = p_yfrac (e.g., 0.83), independent of data scale
        trans_line = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        x1, x2 = 0, 1
        h = 0.015  # bracket height in axes units
        ax.plot([x1, x1, x2, x2],
                [p_yfrac - h, p_yfrac, p_yfrac, p_yfrac - h],
                transform=trans_line, color="black", lw=1.4, clip_on=False)

        # Text centered between boxes, also in axes coords
        ax.text(0.5, p_yfrac + 0.006, txt,
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=11)

    plt.tight_layout()
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        fname = f"correlation_boxplot_complex_{complex_idx if complex_idx is not None else 'all'}.png"
        fig.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches="tight")
        plt.close(fig)

def plot_overview_boxplots(
    df_long: pd.DataFrame,
    per_pair_info: List[Dict[str, Any]],
    output_dir: Optional[str] = None,
    filename: str = "correlation_boxplot_overview.png",
    x_by: str = "Complex",
    # Robust by default; set both to None to plot the full range
    y_lower_q: Optional[float] = 5.0,
    y_upper_q: Optional[float] = 99.5,
    inner_sep: Optional[float] = None,       # None => auto from box_width
    box_width: float = 0.22,
    colors: Tuple[str, str] = ("#4C72B0", "#DD8452"),
    violin: bool = False,
    violin_alpha: float = 0.18,
    star_pad_frac: float = 0.012,            # vertical pad above the max (fraction of span)
    bracket_height_frac: float = 0.008,      # bracket height (fraction of span)
):
    import os, re
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib import ticker as mticker

    # ---- Order and data ----
    order = [d["key"] for d in per_pair_info] if x_by == "Complex" else [d["name"] for d in per_pair_info]
    x_col = "Complex" if x_by == "Complex" else "Name"
    n_cat = len(order)

    complex_vals, outside_vals = [], []
    for cat in order:
        vc = df_long[(df_long[x_col] == cat) & (df_long["Group"] == "Complex")]["Correlation"].dropna().to_numpy()
        vo = df_long[(df_long[x_col] == cat) & (df_long["Group"] == "Outside")]["Correlation"].dropna().to_numpy()
        complex_vals.append(vc)
        outside_vals.append(vo)

    if inner_sep is None:
        inner_sep = box_width * 0.60
    centers = np.arange(n_cat, dtype=float)
    pos_complex = centers - inner_sep
    pos_outside = centers + inner_sep

    # ---- Y limits (percentiles lado a lado; full-range sólo donde sea None) ----
    y_all = df_long["Correlation"].to_numpy(float)

    def _q_or_extreme(arr, q, extreme_fn):
        return float(extreme_fn(arr)) if q is None else float(np.nanpercentile(arr, q))

    y_lo = _q_or_extreme(y_all, y_lower_q, np.nanmin)   # usa percentil si y_lower_q != None
    y_hi = _q_or_extreme(y_all, y_upper_q, np.nanmax)   # usa percentil si y_upper_q != None

    # Evita rango degenerado
    if not np.isfinite(y_lo): y_lo = -1.0
    if not np.isfinite(y_hi): y_hi = 1.0
    if y_hi <= y_lo:
        y_hi = y_lo + 1e-3

    span = y_hi - y_lo

    # ---- Figure ----
    fig_w = max(6.0, 0.9 + 0.9 * n_cat)
    with plt.rc_context({
        "font.size": 9, "axes.titlesize": 18, "axes.labelsize": 11,
        "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10
    }):
        fig, ax = plt.subplots(figsize=(fig_w, 4.8), dpi=300)

        col_c, col_o = colors

        # Optional violins (drawn first)
        if violin:
            vwidth = box_width * 1.6
            vp1 = ax.violinplot(complex_vals, positions=pos_complex, widths=vwidth, showextrema=False)
            for b in vp1['bodies']:
                b.set_facecolor(col_c); b.set_edgecolor('none'); b.set_alpha(violin_alpha)
            vp2 = ax.violinplot(outside_vals, positions=pos_outside, widths=vwidth, showextrema=False)
            for b in vp2['bodies']:
                b.set_facecolor(col_o); b.set_edgecolor('none'); b.set_alpha(violin_alpha)

        # Boxplots on top
        bp1 = ax.boxplot(
            complex_vals, positions=pos_complex, widths=box_width,
            whis=(5, 95), showfliers=False, patch_artist=True, manage_ticks=False
        )
        bp2 = ax.boxplot(
            outside_vals, positions=pos_outside, widths=box_width,
            whis=(5, 95), showfliers=False, patch_artist=True, manage_ticks=False
        )

        def _style(bp, face):
            for box in bp["boxes"]:
                box.set_facecolor(face); box.set_alpha(0.92)
                box.set_edgecolor("black"); box.set_linewidth(1.0)
            for med in bp["medians"]:
                med.set_color("black"); med.set_linewidth(1.8)
            for w in bp["whiskers"]:
                w.set_linewidth(1.0)
            for cap in bp["caps"]:
                cap.set_linewidth(1.0)
        _style(bp1, col_c); _style(bp2, col_o)

        # X ticks
        ax.set_xticks(centers)
        def _short(lbl: str) -> str:
            return re.sub(r'^\s*Spliceosome,\s*', '', lbl).replace(" complex", "")
        ax.set_xticklabels([_short(c) for c in order])

        # Legend (top center, compact)
        patches = [
            Patch(facecolor=col_c, edgecolor="black", label="Complex"),
            Patch(facecolor=col_o, edgecolor="black", label="Non-complex"),
        ]
        leg = ax.legend(handles=patches, loc="lower center",
                        bbox_to_anchor=(0.5, 1.06), ncol=2, frameon=True)
        leg.get_frame().set_alpha(0.96)
        leg.get_frame().set_linewidth(0.8)

        # Title: TCGA short (if we can infer it)
        if output_dir:
            tcga, _ = _infer_tumor_from_output_dir(output_dir)
            fig.suptitle(tcga, weight="bold", y=0.98)
        else:
            fig.suptitle("TCGA", weight="bold", y=0.98)

        ax.set_xlabel("Spliceosome Complexes")
        ax.set_ylabel(r"Pearson correlation ($r$)")

        # Final Y-lims with small padding to ensure stars/brackets fit
        ax.set_ylim(y_lo, y_hi + 0.08 * span)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune="both"))
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))

        # Tight horizontal margins
        ax.margins(x=0.02)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

        # Stars above the highest observed point per complex
        pad = star_pad_frac * span
        bh  = bracket_height_frac * span
        for i, info in enumerate(per_pair_info):
            p = float(info["pval"])
            max_val = float(np.nanmax(np.r_[complex_vals[i], outside_vals[i]]))
            y_base = max_val + pad
            x0, x1 = pos_complex[i], pos_outside[i]
            ax.plot([x0, x0, x1, x1], [y_base, y_base + bh, y_base + bh, y_base],
                    color="black", lw=1.0, clip_on=False)
            ax.text((x0 + x1) / 2, y_base + bh * 0.95, _p_to_stars(p),
                    ha="center", va="bottom", fontsize=12, fontweight="bold")

        plt.tight_layout()
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            fig.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()