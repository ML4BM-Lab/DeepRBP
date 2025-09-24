
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import transforms
import re
from typing import Optional, Tuple, List, Dict, Any
from matplotlib.ticker import FixedLocator

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

    # # Boxplot without outliers (whis 5–95) and emphasized medians
    # sns.boxplot(
    #     x="Group", y="Correlation", data=data, ax=ax, width=0.55,
    #     showfliers=False, whis=(5, 95), palette=palette, saturation=1,
    #     medianprops={"linewidth": 2.2, "color": "black"},
    #     boxprops={"linewidth": 1.4}, whiskerprops={"linewidth": 1.4}, capprops={"linewidth": 1.4},
    # )
    # ax.set_xticklabels(labels, fontsize=11)

    # ✅ FIX: añade hue="Group" (y quitamos cualquier legend kw)
    ax = sns.boxplot(
        x="Group", y="Correlation", hue="Group", data=data, ax=ax,
        width=0.55, showfliers=False, whis=(5, 95), palette=palette, dodge=False,
        medianprops={"linewidth": 2.2, "color": "black"},
        boxprops={"linewidth": 1.4}, whiskerprops={"linewidth": 1.4}, capprops={"linewidth": 1.4},
    )
    # Quitamos la leyenda automática que aparece al usar hue
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # ✅ FIX: fija los ticks antes de poner etiquetas
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

def plot_overview_boxplots_cached(
    df_long: pd.DataFrame,
    per_pair_info: List[Dict[str, Any]],
    output_dir: Optional[str] = None,
    filename: str = "correlation_boxplot_overview.png",
    x_by: str = "Complex",
    stouffer: Optional[Tuple[float, float]] = None,   # (z_st, p_st)
    fisher: Optional[Tuple[float, float]] = None,     # (chi_f, p_f)
    p_yfrac: float = 0.88,                            # altura (fracción eje Y) para las barras/p
    w_per_cat: float = 1.00,
    box_width: float = 0.30,
    xtick_rotation: int = 0
) -> None:
    order = [d["key"] for d in per_pair_info] if x_by == "Complex" else [d["name"] for d in per_pair_info]
    n_cat = max(1, len(order))
    fig_w = max(7.0, 1.2 + w_per_cat * n_cat)

    # ↓↓↓ Tipografías más pequeñas en todo el plot
    with plt.rc_context({
        "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 10,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9
    }):
        sns.set_theme(style="white", rc={"axes.grid": False})
        fig, ax = plt.subplots(figsize=(fig_w, 5.0), dpi=300)

        palette = ["#0072B2", "#E69F00"]  # Okabe–Ito
        x_col = "Complex" if x_by == "Complex" else "Name"
        sns.boxplot(
            data=df_long, x=x_col, y="Correlation", hue="Group", order=order,
            width=box_width, dodge=True, showfliers=False, whis=(5, 95),
            palette=palette,
            medianprops={"linewidth": 1.6, "color": "black"},
            boxprops={"linewidth": 1.1}, whiskerprops={"linewidth": 1.1}, capprops={"linewidth": 1.1},
            ax=ax
        )

        # Leyenda arriba, más compacta
        h, l = ax.get_legend_handles_labels()
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        leg = ax.legend(h, l, title="", loc="lower center",
                        bbox_to_anchor=(0.5, 1.10), ncol=2, frameon=True, prop={"size": 9})
        leg.get_frame().set_alpha(0.96)
        leg.get_frame().set_facecolor("white")
        leg.get_frame().set_edgecolor("#999999")
        leg.get_frame().set_linewidth(0.8)

        # Título (TCGA si se puede inferir)
        if output_dir:
            tcga, long_name = _infer_tumor_from_output_dir(output_dir)
            fig.suptitle(f"{tcga} — {long_name.replace('_',' ')}", fontsize=14, weight="bold", y=0.98)
        else:
            fig.suptitle("TCGA overview", fontsize=14, weight="bold", y=0.98)

        # Subtítulo
        subtitle = []
        if stouffer is not None:
            z_st, p_st = stouffer
            subtitle.append(f"Stouffer: Z = {z_st:.2f}, p = {p_st:.2e}")
        if fisher is not None:
            chi_f, p_f = fisher
            subtitle.append(f"Fisher: χ² = {chi_f:.2f}, p = {p_f:.2e}")
        if subtitle:
            ax.set_title("   |   ".join(subtitle), fontsize=9.5, pad=6)

        ax.set_xlabel("Complexes", labelpad=6, fontsize=10)
        ax.set_ylabel(r"Pearson correlation ($r$)", fontsize=10)
        sns.despine(ax=ax)

        # X labels sin “Spliceosome, ” ni “ complex”
        def _short(lbl: str) -> str:
            return re.sub(r'^\s*Spliceosome,\s*', '', lbl).replace(" complex", "")
        
        # ax.set_xticklabels([_short(t.get_text()) for t in ax.get_xticklabels()],
        #                    rotation=xtick_rotation, ha="center")

        # ✅ FIX: usa las posiciones esperadas (0..len(order)-1)
        tick_locs = np.arange(len(order))
        ax.xaxis.set_major_locator(FixedLocator(tick_locs))
        ax.set_xticks(tick_locs)
        ax.set_xticklabels([_short(lbl) for lbl in order], rotation=xtick_rotation, ha="center")

        # Límites Y robustos
        y = df_long["Correlation"].to_numpy(float)
        q1, q99 = np.nanpercentile(y, [1, 99])
        span    = max(q99 - q1, 1e-3)
        ax.set_ylim(q1 - 0.02 * span, q99 + 0.06 * span)

        # Barras + ⭐ y p en dos líneas
        sep = (box_width / 2.0) * 1.05
        left_off, right_off = -sep, +sep
        h_axes = 0.014
        trans_line = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        tick_lookup = {lab: i for i, lab in enumerate(order)}

        for d in per_pair_info:
            cat = d["key"] if x_by == "Complex" else d["name"]
            if cat not in tick_lookup:
                continue
            i = tick_lookup[cat]
            p = float(d["pval"])
            x0, x1 = i + left_off, i + right_off

            # bracket
            ax.plot([x0, x0, x1, x1],
                    [p_yfrac - h_axes, p_yfrac, p_yfrac, p_yfrac - h_axes],
                    transform=trans_line, color="black", lw=1.0, clip_on=False)
            # ⭐ arriba
            ax.text((x0 + x1) / 2, p_yfrac + h_axes * 1.15, _p_to_stars(p),
                    transform=trans_line, ha="center", va="bottom", fontsize=10)
            # (p=...) debajo de la estrella (pero aún encima del bracket)
            ax.text((x0 + x1) / 2, p_yfrac + h_axes * 0.10, f"(p={p:.1e})",
                    transform=trans_line, ha="center", va="bottom", fontsize=7.5)

        plt.tight_layout()
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            fig.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()


# def plot_upper_triangle_corr_matrix(corr_matrix, complex_idx=None, output_dir=None, subset_size=100): # muy mejorable
#     """
#     Plot a simple heatmap showing only the upper triangle (including diagonal) of a subset
#     of a large correlation matrix, without labels or clustering.

#     Parameters:
#     - corr_matrix: pd.DataFrame, square correlation matrix.
#     - complex_idx: str or int or None (default None)
#         Identifier of the protein complex to show in the plot title.
#     - output_dir: str or None, path to save the figure (if None, no save)
#     - subset_size: int, number of genes to keep for plotting (default 100)
#     """
#     # Subset matrix
#     corr_small = corr_matrix.iloc[:subset_size, :subset_size].copy()
#     # Create mask for lower triangle
#     mask = np.tril(np.ones_like(corr_small, dtype=bool), k=-1)
#     plt.figure(figsize=(8,8))
#     sns.set_theme(style="white")
#     ax = sns.heatmap(
#         corr_small,
#         mask=mask,
#         cmap="RdYlBu_r",
#         square=True,
#         cbar_kws={"label": "Correlation"},
#         xticklabels=False,
#         yticklabels=False,
#         linewidths=0,
#         vmin=0, vmax=1
#     )
#     plt.title("Upper Triangle Correlation Heatmap (subset)", fontsize=14, fontweight='bold')
#     if output_dir:
#         os.makedirs(output_dir, exist_ok=True)
#         path = os.path.join(output_dir, f"upper_triangle_heatmap_{complex_idx}.png")
#         plt.savefig(path, dpi=300, bbox_inches='tight')
#         plt.close()