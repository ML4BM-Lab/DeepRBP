# src/deeprbp/util/plots.py

import os 
import math
from math import ceil
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import FuncFormatter, MaxNLocator
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List, Union, Iterable, Optional, Tuple

def scatter_real_vs_pred(
    category: str,
    metrics: Dict[str, float],
    pred: Union[List[float], np.ndarray],
    labels: Union[List[float], np.ndarray],
    output_dir: str
) -> None:
    """
    Creates a scatter plot comparing predicted vs. real values for a regression model.
    Includes regression line and performance metrics for clear visualization.

    Parameters:
    -----------
    category : str
        The category of the samples (e.g., cancer type or dataset name).
    metrics : Dict[str, float]
        Dictionary of evaluation metrics: Spearman correlation, Pearson correlation, MSE, R².
    pred : Union[List[float], np.ndarray]
        List or array of predicted values.
    labels : Union[List[float], np.ndarray]
        List or array of true/real values.
    output_dir : str
        Base path to save the plot as a PNG image.

    Returns:
    --------
    None
        The function saves a PNG image at the specified path.
    """
    sns.set_style("white")
    plt.figure(figsize=(8, 6))
    plt.xlabel('Predicted Values', fontsize=16, fontweight='bold')
    plt.ylabel('Real Values', fontsize=16, fontweight='bold')
    plt.title(f'Real vs Predicted: {category})', fontsize=18, fontweight='bold')
    sns.regplot(
        x=pred, y=labels,
        scatter_kws={'alpha': 0.3, 'color': 'blue'},
        line_kws={'color': 'red', 'lw': 2}
    )
    legend_text = (
        f"Spearman Corr: {metrics['spearman_corr']:.4f}\n"
        f"Pearson Corr: {metrics['pearson_corr']:.4f}\n"
        f"MSE: {metrics['mse']:.4f}\n"
        f"R²: {metrics['r2']:.4f}\n"
        f"Spearman Corr per Gene: {metrics['mean_corr_per_gene']:.4f} \n"
        f"Spearman Corr per Gene (max trans): {metrics['mean_corr_max_trans_per_gene']:.4f} \n"
    )
    plt.text(
        0.05, 0.95, legend_text, transform=plt.gca().transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle="round", edgecolor="black", facecolor="white")
    )
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{category}.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved at: {output_path}")
 
def plot_small_multiples_real_vs_pred_grid(
    panels: Iterable[dict],
    set_name: str,
    output_dir: str,
    order_codes: list[str],
    grid: tuple[int, int] = (6, 6),
    use_hexbin: bool = True,
    gridsize: int = 30,
    show_title: bool = False,        # default OFF
    show_colorbar: bool = False,     # default OFF
    figsize: tuple[float, float] = (11, 11),
    n_major_ticks: int = 6,          # nº aprox. de ticks
    xlabel: Optional[str] = "Predicted log2(TPM + 1)",
    ylabel: Optional[str] = "Observed log2(TPM + 1)",
    axis_range: Optional[tuple[float, float]] = None,  # (lo, hi) común
    # --- NUEVO: control de densidad/colores ---
    cmap: str = "inferno",
    density_scale: str = "log",      # "log" o "linear"
    share_density_norm: bool = True, # misma normalización en todos los paneles
):
    """
    Small-multiples Observed vs. Predicted log2(TPM+1) por tejido.

    • Límites comunes y línea identidad (y = x).
    • Misma rejilla de ticks (sin decimales) en X e Y.
    • Hexbin opcional con escala log y normalización compartida (color más
      intenso donde hay mayor densidad de puntos).
    • Etiquetas globales de ejes; títulos por panel con el código TCGA.
    """
    panels = list(panels)
    if not panels:
        return
    # ---- Orden por código TCGA (estable; desconocidos al final)
    pos = {code: i for i, code in enumerate(order_codes)}
    panels = sorted(panels, key=lambda d: (pos.get(d.get("short", ""), 10**9)))
    # ---- Límites comunes
    if axis_range is None:
        all_pred = np.concatenate([p["pred"] for p in panels])
        all_true = np.concatenate([p["true"] for p in panels])
        lo = float(np.nanmin([all_pred.min(), all_true.min()]))
        hi = float(np.nanmax([all_pred.max(), all_true.max()]))
    else:
        lo, hi = axis_range
    # ---- Ticks compartidos (enteros)
    locator = MaxNLocator(nbins=n_major_ticks, steps=[1, 2, 2.5, 5, 10])
    tick_values = locator.tick_values(lo, hi)
    tick_values = tick_values[(tick_values >= lo - 1e-9) & (tick_values <= hi + 1e-9)]
    int_formatter = FuncFormatter(lambda x, pos: f"{x:.0f}")
    # ---- Normalización de densidad común (para hexbin)
    norm = None
    if use_hexbin and share_density_norm:
        # Aproximamos la densidad con un hist2d para fijar vmax (coherente entre paneles)
        max_count = 1
        for p in panels:
            H, _, _ = np.histogram2d(
                p["pred"], p["true"],
                bins=gridsize,
                range=[[lo, hi], [lo, hi]],
            )
            m = int(H.max())
            if m > max_count:
                max_count = m
        norm = LogNorm(vmin=1, vmax=max_count) if density_scale == "log" else Normalize(vmin=0, vmax=max_count)
    # ---- Figura
    nrows, ncols = grid
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    axes = axes.ravel()
    out_dir = os.path.join(output_dir, "scat_plot_real_vs_pred_grid", set_name)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"grid_{set_name}.png")
    last_hb = None
    for i, ax in enumerate(axes):
        if i < len(panels):
            p = panels[i]
            x, y = p["pred"], p["true"]
            if use_hexbin:
                last_hb = ax.hexbin(
                    x, y,
                    gridsize=gridsize,
                    extent=(lo, hi, lo, hi),  # asegura bins alineados en todos los ejes
                    mincnt=1,
                    linewidths=0,
                    cmap=cmap,
                    norm=norm,                # <- misma escala de color
                )
                # rasteriza la colección para PDFs/SVG más ligeros
                last_hb.set_rasterized(True)
            else:
                ax.scatter(x, y, s=2, alpha=0.08, edgecolors="none")
            # identidad + límites
            ax.plot([lo, hi], [lo, hi], color="0.7", lw=0.8)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_aspect("equal", adjustable="box")
            # misma rejilla de ticks
            ax.set_xticks(tick_values)
            ax.set_yticks(tick_values)
            ax.xaxis.set_major_formatter(int_formatter)
            ax.yaxis.set_major_formatter(int_formatter)
            # título corto
            ax.set_title(p.get("short", p.get("category", "")), fontsize=9, pad=2)
            # oculta etiquetas interiores
            if (i // ncols) < (nrows - 1):
                ax.set_xticklabels([])
            else:
                ax.tick_params(axis="x", labelsize=6, length=2)
            if (i % ncols) != 0:
                ax.set_yticklabels([])
            else:
                ax.tick_params(axis="y", labelsize=6, length=2)
            # spines ligeros
            for spine in ax.spines.values():
                spine.set_linewidth(0.6)
                spine.set_alpha(0.8)
        else:
            ax.axis("off")
    # Etiquetas globales
    if xlabel:
        fig.supxlabel(xlabel, fontsize=11)
    if ylabel:
        fig.supylabel(ylabel, fontsize=11)
    # Colorbar opcional (comparte la misma norm)
    if show_colorbar and use_hexbin and last_hb is not None:
        cbar = fig.colorbar(last_hb, ax=axes.tolist(), shrink=0.65, pad=0.01)
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label("Point density", fontsize=8)
    if show_title:
        fig.suptitle(f"Observed vs Predicted — {set_name.upper()} (all tissues)", fontsize=12, y=0.996)
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[grid] Saved: {out_path}")


def plot_all_metrics_history(metrics_df: pd.DataFrame, output_dir: str) -> None:
    metrics_pairs = {
        'Correlation Pearson History': ('train_corr_pearson', 'validation_corr_pearson'),
        'Correlation Spearman History': ('train_corr_spearman', 'validation_corr_spearman'),
        'Loss History': ('train_loss', 'validation_loss'),
        'R2 Score History': ('train_r2', 'validation_r2')
    }
    for title, (train_metric, val_metric) in metrics_pairs.items():
        plot_metric_history(
            train_history=metrics_df[train_metric].tolist(),
            val_history=metrics_df[val_metric].tolist(),
            title=title,
            output_dir=output_dir,
            plot_name=title.lower().replace(' ', '_')
        )

def plot_metric_history(train_history: List[float], val_history: List[float], 
                        title: str, 
                        output_dir: str, 
                        plot_name: str) -> None:
    """
    Plots the training and validation metrics history to visualize the performance of the model over time.

    Parameters:
    - train_history (List[float]): List of training metric values for each epoch.
    - val_history (List[float]): List of validation metric values for each epoch.
    - title (str): Title of the plot.
    - output_dir (str): Path to save the plot as a .png file.
    - plot_name (str): Name of the plot file.

    Returns:
    - None: The function will save the plot.
    """
    os.makedirs(output_dir, exist_ok=True)
    plt.style.use('seaborn-v0_8-muted')
    plt.figure(figsize=(8, 5))
    # Plot both training and validation metrics
    plt.plot(train_history, label='Training', color='royalblue', linestyle='-', linewidth=2.5)
    plt.plot(val_history, label='Validation', color='darkorange', linestyle='--', linewidth=2.5)
    plt.title(title, fontsize=18, fontweight='bold', pad=10)
    plt.xlabel('Epoch', fontsize=14, labelpad=8)
    plt.ylabel('Performance', fontsize=14, labelpad=8)  
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12, loc='upper left', bbox_to_anchor=(1.04, 1), frameon=True, shadow=True, fancybox=True)
    plt.tight_layout(pad=2)
    # Save the plot with the specified name
    plt.savefig(os.path.join(output_dir, f'{plot_name}.png'), dpi=300)
    plt.close()
    print(f"Plot saved to {output_dir}/{plot_name}.png")

def plot_transcript_to_gene_ratio_distributions(
    ratios_pred: np.ndarray,
    ratios_label: np.ndarray,
    output_path: str,
):
    """
    Plots histograms of predicted and labeled transcript-to-gene expression ratios.
    Saves the plot as a high-quality PNG file.
    """
    # Set figure size and style
    plt.figure(figsize=(14, 7))
    plt.style.use('ggplot')  # Elegant grid style
    # Define color palette
    color_pred = "#1f77b4"  # Blue
    color_label = "#ff7f0e"  # Orange
    # Predicted ratio histogram
    plt.subplot(1, 2, 1)
    plt.hist(ratios_pred, bins=40, color=color_pred, edgecolor='black', linewidth=1.2)
    mean_pred = np.mean(ratios_pred)
    std_pred = np.std(ratios_pred)
    plt.title(f"Predicted Gene Ratio (Mean > 5 TPM)", fontsize=12, fontweight='bold')
    plt.xlabel("Transcript-to-Gene Ratio", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.axvline(mean_pred, color=color_pred, linestyle='dashed', linewidth=1.5)
    plt.text(mean_pred + 0.1, plt.ylim()[1] * 0.9, f'Mean: {mean_pred:.3f}\nSTD: {std_pred:.3f}', color=color_pred, fontsize=9)
    # Labeled ratio histogram
    plt.subplot(1, 2, 2)
    plt.hist(ratios_label, bins=40, color=color_label, edgecolor='black', linewidth=1.2)
    mean_label = np.mean(ratios_label)
    std_label = np.std(ratios_label)
    plt.title(f"Labeled Gene Ratio (Mean > 5 TPM)", fontsize=12, fontweight='bold')
    plt.xlabel("Transcript-to-Gene Ratio", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.axvline(mean_label, color=color_label, linestyle='dashed', linewidth=1.5)
    plt.text(mean_label + 0.01, plt.ylim()[1] * 0.9, f'Mean: {mean_label:.3f}\nSTD: {std_label:.3f}', color=color_label, fontsize=9)
    x_min = 0.8
    x_max = 1.2
    plt.xlim((x_min, x_max))
    plt.xticks(np.arange(x_min, x_max + 0.01, 0.05))
    # Adjust layout
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_path}")
 
############################################################
# # Load training log
# metrics_df = pd.read_csv("/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/csv_logs/deep_rbp_predictor/version_0/metrics.csv")

# # Where to save the grid figure
# out_dir = "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor/history"
# os.makedirs(out_dir, exist_ok=True)

# plot_selected_history_metrics_grid(
#     df=metrics_df,
#     metric_keys=["loss", "spearman", "r2", "pearson"],
#     out_path=os.path.join(out_dir, "training_dynamics_grid.png"),
#     best_epoch=124,
#     ncols=2,
#     figsize=(10, 7),
#     suptitle="Training and validation curves (checkpoint: epoch 124)",
#     use_auto_ylim=True,
#     auto_ylim_pad=0.02,
#     mse_log=True,             # loss en log
#     annotate_best=True,
# )

############################################################


def plot_selected_history_metrics_grid(
    df: pd.DataFrame,
    metric_keys: List[str],
    out_path: str,
    best_epoch: Optional[int] = None,
    ncols: int = 2,
    figsize: Tuple[float, float] = (10, 7),
    suptitle: Optional[str] = None,
    use_auto_ylim: bool = True,
    auto_ylim_pad: float = 0.02,
    mse_log: bool = True,
    annotate_best: bool = True,
) -> None:
    METRIC_SPECS = {
        "loss":     {"train": "train_loss",              "val": "validation_loss",              "ylabel": "MSE"},
        "pearson":  {"train": "train_corr_pearson",      "val": "validation_corr_pearson",      "ylabel": "Pearson r"},
        "spearman": {"train": "train_corr_spearman",     "val": "validation_corr_spearman",     "ylabel": "Spearman \u03C1"},
        "spearman_per_gene": {
            "train": "train_corr_spearman_per_gene",
            "val":   "validation_corr_spearman_per_gene",
            "ylabel": "Mean per-gene Spearman \u03C1",
        },
        "r2":       {"train": "train_r2",                "val": "validation_r2",                "ylabel": "R\u00B2"},
    }
    TRAIN_COLOR = "#303030"
    VAL_COLOR   = "#D98C00"
    TRAIN_STYLE = "-"
    VAL_STYLE   = "--"
    def _fmt_value(v: float, key: str) -> str:
        if np.isnan(v):
            return "NA"
        return f"{v:.4f}" if key == "loss" else f"{v:.3f}"
    nplots = len(metric_keys)
    nrows = ceil(nplots / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True)
    axes = np.atleast_1d(axes).ravel()
    legend_handles, legend_labels = [], ["Training", "Validation"]
    x = np.arange(len(df), dtype=int)
    for i, key in enumerate(metric_keys):
        ax = axes[i]
        spec = METRIC_SPECS[key]
        tr = np.asarray(df[spec["train"]].values, dtype=float)
        va = np.asarray(df[spec["val"]].values,   dtype=float)
        ln_tr, = ax.plot(x, tr, color=TRAIN_COLOR, linestyle=TRAIN_STYLE, linewidth=1.8, alpha=0.95, zorder=2)
        ln_va, = ax.plot(x, va, color=VAL_COLOR,   linestyle=VAL_STYLE,   linewidth=2.0, alpha=0.95, zorder=3)
        if not legend_handles:
            legend_handles = [ln_tr, ln_va]
        # --- BEST EPOCH: marker + vline + VALUE LABEL (placement tweak for loss)
        if best_epoch is not None and 0 <= best_epoch < len(va) and not np.isnan(va[best_epoch]):
            ax.axvline(best_epoch, color=VAL_COLOR, linestyle="--", linewidth=1.0, alpha=0.35, zorder=1)
            if annotate_best:
                ax.scatter(best_epoch, va[best_epoch], s=24, color=VAL_COLOR,
                        edgecolors="white", linewidths=0.8, zorder=4)
                # Place labels: below for all metrics, but a bit higher-above for loss
                if key == "loss":
                    dy_points = 12    # was 6 → “a little higher”
                    va_text = "bottom"
                else:
                    dy_points = -18   # keep below for non-loss metrics
                    va_text = "top"
                ax.annotate(
                    f"{_fmt_value(va[best_epoch], key)} @ e{best_epoch}",
                    xy=(best_epoch, va[best_epoch]),
                    xytext=(6, dy_points), textcoords="offset points",
                    ha="left", va=va_text,
                    fontsize=8, color=VAL_COLOR,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=VAL_COLOR, lw=0.5, alpha=0.9),
                    zorder=5
                )
        ax.set_ylabel(spec["ylabel"], labelpad=6)
        if key == "loss" and mse_log:
            ax.set_yscale("log")
        if use_auto_ylim:
            y_all = np.concatenate([tr[~np.isnan(tr)], va[~np.isnan(va)]])
            if y_all.size > 0 and np.isfinite(y_all.min()) and np.isfinite(y_all.max()):
                ymin, ymax = y_all.min(), y_all.max()
                pad = auto_ylim_pad * (ymax - ymin if ymax > ymin else 1.0)
                ax.set_ylim(ymin - pad, ymax + pad)
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.2)
    for j in range(nplots, nrows * ncols):
        axes[j].axis("off")
    fig.supxlabel("Epoch", fontsize=11)
    if suptitle:
        fig.suptitle(suptitle, fontsize=12, y=0.985)
        legend_y = 0.952
    else:
        legend_y = 0.985
    fig.legend(legend_handles, legend_labels, loc="upper center",
               bbox_to_anchor=(0.5, legend_y), ncol=2, frameon=False, handlelength=2.8)
    plt.tight_layout(rect=(0.03, 0.03, 0.97, 0.90))
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
