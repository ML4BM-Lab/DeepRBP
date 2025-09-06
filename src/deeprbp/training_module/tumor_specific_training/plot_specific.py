
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm

from ..tcga_codes import TCGA_CODE

ABBREVIATION_DICT = {**TCGA_CODE, "all": "all"}

RANGES = {
    "pearson_corr": dict(cmap="viridis", vmin=0.80, vmax=0.98, norm=None),
    "spearman_corr": dict(cmap="viridis", vmin=0.70, vmax=0.89, norm=None),
    "r2": dict(cmap="viridis", vmin=0.60, vmax=0.96, norm=None),
    "mean_corr_per_gene": dict(cmap="viridis", vmin=0.65, vmax=0.85, norm=None),
    "mse": dict(cmap="magma_r", vmin=0.06, vmax=1.30, norm="log"),
}

_METRIC_LABELS = {
    "pearson_corr": (
        r"Pearson $r$ across tumor types",
        r"Pearson $r$",
    ),
    "spearman_corr": (
        r"Spearman $\rho$ across tumor types",
        r"Spearman $\rho$",
    ),
    "r2": (
        r"Coefficient of determination $R^2$ across tumor types",
        r"$R^2$",
    ),
    "mean_corr_per_gene": (
        r"Mean per-gene Spearman ($\overline{\rho}_{\mathrm{gene}}$) across tumor types",
        r"$\overline{\rho}_{\mathrm{gene}}$",
    ),
    "mse": (
        "Mean squared error (MSE) across tumor types",
        "MSE",
    ),
}

def plot_metric_confusion_matrix(
    df,
    metric_name,
    output_path="confusion_matrix.png",
    title=None
):
    metric_name = metric_name.lower()
    # --- Title & colorbar label (paper-style)
    if title is None:
        if metric_name in _METRIC_LABELS:
            title, cbar_label = _METRIC_LABELS[metric_name]
        else:
            base = metric_name.replace('_', ' ')
            title = (base.upper() if metric_name == "mse" else base.capitalize()) + " across tumor types"
            cbar_label = base.upper() if metric_name == "mse" else base.capitalize()
    else:
        cbar_label = _METRIC_LABELS.get(metric_name, (None, metric_name.replace('_', ' ').capitalize()))[1]
    # --- Map long names → TCGA short codes and enforce plotting order
    df_plot = df.copy()
    df_plot.index = [ABBREVIATION_DICT.get(idx, idx) for idx in df_plot.index]
    df_plot.columns = [ABBREVIATION_DICT.get(col, col) for col in df_plot.columns]
    df_plot.index.name = 'Train Tumor Type'
    df_plot.columns.name = 'Eval Tumor Type'
    desired_order = list(ABBREVIATION_DICT.values())
    df_plot = df_plot.reindex(index=[c for c in desired_order if c in df_plot.index],
                              columns=[c for c in desired_order if c in df_plot.columns])
    # --- Figure sizing
    fig_width = min(max(4, len(df_plot.columns) * 0.6), 20)
    fig_height = min(max(3, len(df_plot) * 0.45), 15)
    plt.figure(figsize=(fig_width, fig_height))
    base_fontsize = 10
    scaling_factor = min(fig_width / 10, fig_height / 6)
    font_size = max(6, base_fontsize * scaling_factor)
    # --- Apply RANGES spec (cmap, bounds, optional log scale)
    cfg = RANGES.get(metric_name, dict(cmap="viridis", vmin=None, vmax=None, norm=None))
    use_log = (cfg.get("norm") == "log")
    norm = LogNorm(vmin=cfg.get("vmin"), vmax=cfg.get("vmax")) if use_log else None
    ax = sns.heatmap(
        df_plot,
        annot=True,
        fmt=".2f",
        cmap=cfg.get("cmap", "viridis"),
        vmin=None if norm else cfg.get("vmin"),
        vmax=None if norm else cfg.get("vmax"),
        norm=norm,
        cbar_kws={"label": cbar_label},
        annot_kws={"fontsize": font_size * 0.6},
        linewidths=0.5,
        linecolor='lightgrey',
    )
    ax.set_title(title, fontsize=font_size + 2, pad=20)
    plt.xticks(rotation=90, fontsize=font_size)
    plt.yticks(rotation=0, fontsize=font_size)
    plt.xlabel(df_plot.columns.name, fontsize=font_size, labelpad=15)
    plt.ylabel(df_plot.index.name, fontsize=font_size, labelpad=15)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=font_size)
    cbar.set_label(cbar_label, fontsize=font_size, labelpad=6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
