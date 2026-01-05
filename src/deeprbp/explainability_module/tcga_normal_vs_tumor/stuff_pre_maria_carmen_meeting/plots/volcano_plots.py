
# src/deeprbp/explainability_module/tcga_normal_vs_tumor/volcano_plots.py

import os
from typing import Optional, Iterable, Tuple, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from adjustText import adjust_text
    _HAVE_ADJUSTTEXT = True
except ImportError:
    _HAVE_ADJUSTTEXT = False

def plot_rbp_level_summary(
    df_rbp: pd.DataFrame,
    rbp_name_map: dict = None,
    log2fc_col: str = "median_log2_fc",
    p_col: str = "min_p_adj_fdr",
    prop_sig_col: str = "prop_sig_tx_FDR_0_01",
    p_thr: float = 0.01,
    fc_thr: float = 1.0,
    prop_sig_local_min: float = 0.30,
    prop_sig_local_strong: float = 0.40,
    title: str = "RBP-level DeepRBP reprogramming",
    out_path: str = None,
    point_size: float = 40.0,
    prioritized_ids: Optional[Iterable[str]] = None,
    max_labels: int = 45):
    """
    RBP-level reprogramming map summarizing transcript-level DeepRBP dysregulation.

    Each point represents one RBP (~524 RBPs). For each RBP, the plot
    summarizes its transcript-level changes using three aggregates:
      - median_log2_fc        → global direction of change (tumor vs normal),
      - prop_sig_tx_FDR_0_01  → fraction of targets with significant shifts,
      - min_p_adj_fdr         → minimum transcript FDR (used to call RBPs).

    Clasificación visual:
      - Background (gris):
          RBPs sin evidencia clara de reprogramación, más
          los 'local_program_weak' (prop_sig_local_min ≤ prop_sig <
          prop_sig_local_strong).
      - Locally rewired (morado):
          'local_program_strong': min_p_adj_fdr < p_thr,
          |median_log2_fc| < fc_thr y prop_sig ≥ prop_sig_local_strong.
      - Global down / up (azul / rojo):
          min_p_adj_fdr < p_thr y |median_log2_fc| ≥ fc_thr,
          con el signo determinando direction (tumor < normal vs tumor > normal).

    Biological questions addressed
    ------------------------------
    - Which RBPs undergo global gain / loss of regulatory influence in tumors?
    - Which RBPs show local rewiring of subsets of targets without a net shift?

    Plot definition
    ---------------
    x-axis: median_log2_fc  (tumor vs normal)
    y-axis: prop_sig_tx_FDR_0_01  (fraction of significant targets per RBP)

    RBP categories (based on median_log2_fc and min_p_adj_fdr)
    ----------------------------------------------------------
      - Global up:
            min_p_adj_fdr < p_thr  and  median_log2_fc ≥ +fc_thr
      - Global down:
            min_p_adj_fdr < p_thr  and  median_log2_fc ≤ -fc_thr
      - Locally rewired:
            min_p_adj_fdr < p_thr  and  |median_log2_fc| < fc_thr
      - Background:
            min_p_adj_fdr ≥ p_thr

    Global up/down RBPs show coherent and large-magnitude shifts in DeepRBP
    scores across many transcripts. Locally rewired RBPs also pass the FDR
    threshold but have median_log2_fc ≈ 0, indicating mixed (up and down)
    transcript-level effects and pathway-specific rewiring.

    Parameters
    ----------
    df_rbp : pd.DataFrame
        RBP-level summary table produced after aggregating transcript-level
        Wilcoxon results.
    rbp_name_map : dict, optional
        Map from RBP ID → RBP symbol. If None, the RBP ID is used.
    log2fc_col : str
        Column name for the median transcript-level log2 fold-change.
    p_col : str
        Column name for the minimum adjusted p-value across transcripts.
    prop_sig_col : str
        Column with the proportion of significant transcripts (FDR < 0.01).
    p_thr : float
        Adjusted p-value threshold at RBP level.
    fc_thr : float
        Threshold in log2 units to classify large-magnitude global effects.
    title : str
        Title for the plot.
    out_path : str
        If not None, save the figure to this path.
     point_size : float
        Marker size for all points.
    """
    if df_rbp.empty:
        print("[summary-RBP] Empty df_rbp; nothing to plot.")
        return

    df = df_rbp.copy()

    # IDs → symbols
    if rbp_name_map is not None:
        df["RBP_name"] = df["RBP"].map(rbp_name_map).fillna(df["RBP"])
    else:
        df["RBP_name"] = df["RBP"]

    x = df[log2fc_col].astype(float).values
    prop = df[prop_sig_col].astype(float).values
    pvals = df[p_col].astype(float).values

    sig = pvals < p_thr

    # Global up/down (según mediana)
    global_up = sig & (x >= fc_thr)
    global_down = sig & (x <= -fc_thr)

    # Local en sentido amplio
    local_any = sig & (np.abs(x) < fc_thr) & (prop >= prop_sig_local_min)

    # Local strong vs weak
    local_strong = local_any & (prop >= prop_sig_local_strong)
    local_weak = local_any & (prop < prop_sig_local_strong)

    # Background: todo lo que no es global ni local strong/weak
    background = ~(global_up | global_down | local_any)

    # colores
    col_bg = "lightgrey"
    col_local = "#984ea3"   # morado
    col_down = "#377eb8"    # azul
    col_up = "#e41a1c"      # rojo

    fig, ax = plt.subplots(figsize=(6.4, 6.2))

    # background + local_weak → gris
    bg_mask = background | local_weak
    ax.scatter(
        x[bg_mask],
        prop[bg_mask],
        s=point_size,
        color=col_bg,
        alpha=0.35,
        edgecolor="none",
        label="Background",
    )

    # local strong
    ax.scatter(
        x[local_strong],
        prop[local_strong],
        s=point_size,
        color=col_local,
        alpha=0.8,
        edgecolor="none",
        label="Locally rewired",
    )

    # global down
    ax.scatter(
        x[global_down],
        prop[global_down],
        s=point_size,
        color=col_down,
        alpha=0.9,
        edgecolor="none",
        label="Loss (tumor < normal)",
    )

    # global up
    ax.scatter(
        x[global_up],
        prop[global_up],
        s=point_size,
        color=col_up,
        alpha=0.9,
        edgecolor="none",
        label="Gain (tumor > normal)",
    )

    # líneas verticales en ±fc_thr
    ax.axvline(-fc_thr, linestyle="--", color="black", linewidth=0.8)
    ax.axvline(+fc_thr, linestyle="--", color="black", linewidth=0.8)

    # límites simétricos en X
    if np.any(np.isfinite(x)):
        max_abs = np.nanmax(np.abs(x[np.isfinite(x)]))
        ax.set_xlim(-max_abs * 1.05, max_abs * 1.05)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Median log2FC of DeepRBP score\n(tumor vs normal)")
    ax.set_ylabel("Fraction of significant targets")
    ax.set_title(title)

    # resumen de recuentos (solo global + local strong)
    n_gup = int(global_up.sum())
    n_gdown = int(global_down.sum())
    n_local_strong = int(local_strong.sum())
    txt = f"Shifted: {n_gdown + n_gup} (↓{n_gdown}, ↑{n_gup}) \n |  Local strong: {n_local_strong}"
    
    # anotación de RBPs priorizados (si se pasa la lista)
    if prioritized_ids is not None and len(prioritized_ids) > 0 and max_labels > 0:
        mask_prior = df["RBP"].isin(prioritized_ids).values
        df_prior = df[mask_prior].copy()

        if not df_prior.empty:
            # criterio de orden: primero por fracción de dianas, luego por |log2FC|
            df_prior["abs_log2FC"] = df_prior[log2fc_col].abs()
            df_prior = df_prior.sort_values(
                [prop_sig_col, "abs_log2FC"],
                ascending=[False, False],
            ).head(max_labels)

            texts = []
            for _, row in df_prior.iterrows():
                xx = float(row[log2fc_col])
                yy = float(row[prop_sig_col])
                label = row["RBP_name"]
                t = ax.text(
                    xx,
                    yy,
                    label,
                    fontsize=8,
                    ha="right" if xx < 0 else "left",
                    va="bottom",
                )
                texts.append(t)

            # recolocar etiquetas para reducir solapamientos (si adjustText está instalado)
            if _HAVE_ADJUSTTEXT and texts:
                adjust_text(
                    texts,
                    ax=ax,
                    only_move={"points": "y", "texts": "y"},
                    arrowprops=dict(arrowstyle="-", lw=0.5),
                )

    # leyenda centrada a la derecha
    ax.legend(
        frameon=True,           # ← activa la caja
        fancybox=True,          # esquinas redondeadas (opcional)
        framealpha=0.9,         # ligera transparencia
        edgecolor="black",      # borde negro fino
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        borderaxespad=0.0,
        fontsize=8,
    )

    ax.text(
        1.02,
        0.10,              # algo por debajo del centro; ajusta si hace falta
        txt,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )

    plt.tight_layout()

    if out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"[summary-RBP] Saved RBP-level summary to: {out_path}")

    plt.close()

def plot_transcript_volcano_for_rbp(  
    df_tx: pd.DataFrame,
    rbp_id: str,
    rbp_name: Optional[str] = None,
    tx_name_map: Optional[Dict[str, str]] = None,
    p_col: str = "p_adj_fdr",
    fc_col: str = "log2_fc",
    p_thr: float = 0.01,
    fc_thr: float = 1.0,  # |log2FC| >= 1 → ≥ 2-fold change
    out_path: Optional[str] = None,
    max_labels: int = 20,
    rbp_level_log2fc: Optional[float] = None
):
    """
    Volcano plot for all transcripts of a single RBP (global driver or locally rewired).

    One point = one transcript (Transcript_ID x RBP).

    df_tx must contain, for each row:
      - 'RBP' (ID),
      - 'Transcript_ID',
      - fc_col (log2FC of DeepRBP score, tumor vs normal),
      - p_col (FDR-adjusted p-value).

    Parameters
    ----------
    df_rbptx : pd.DataFrame
        Output of compute_rbp_transcript_wilcoxon
        (one row per Transcript_ID x RBP).
    rbp_id : str
        RBP Ensembl ID to plot (e.g. 'ENSG00000004478').
    rbp_name : str, optional
        Human-readable symbol for title. If None, uses rbp_id.
    log2fc_col : str
        Name of the column with log2 fold-change.
    p_col : str
        Name of the column with (adjusted) p-values.
    p_thr : float
        Significance threshold for (adjusted) p-values.
    fc_thr : float
        Threshold in log2 units for calling “strong” change.
    out_path : str
        If not None, save the figure to this path.
    """
    sub = df_tx[df_tx["RBP"] == rbp_id].copy()
    if sub.empty:
        print(f"[volcano-tx] No rows for RBP {rbp_id}. Skipping.")
        return

    if rbp_name is None:
        rbp_name = rbp_id

    sub["log10p"] = -np.log10(np.clip(sub[p_col].astype(float), 1e-300, 1.0))
    sub["abs_fc"] = sub[fc_col].astype(float).abs()

    sig_mask = (sub[p_col] < p_thr) & (sub["abs_fc"] >= fc_thr)

    # colores
    color_ns   = "lightgrey"
    color_up   = "#e41a1c"   # rojo
    color_down = "#377eb8"   # azul
    thr_color  = "dimgray"

    fig, ax = plt.subplots(figsize=(6.2, 6.0))

    # no significativos (más pequeños y transparentes)
    ax.scatter(
        sub.loc[~sig_mask, fc_col],
        sub.loc[~sig_mask, "log10p"],
        s=5,
        alpha=0.15,
        color=color_ns,
    )

    # significativos up
    up = sub[sig_mask & (sub[fc_col] > 0)]
    ax.scatter(
        up[fc_col],
        up["log10p"],
        s=20,
        alpha=0.9,
        color=color_up,
        label="Tumor > Normal",
    )

    # significativos down
    down = sub[sig_mask & (sub[fc_col] < 0)]
    ax.scatter(
        down[fc_col],
        down["log10p"],
        s=20,
        alpha=0.9,
        color=color_down,
        label="Tumor < Normal",
    )

    # líneas de corte
    ax.axhline(-np.log10(p_thr), linestyle="--", linewidth=1, color="black")
    ax.axvline(+fc_thr, linestyle="--", linewidth=1, color="black")
    ax.axvline(-fc_thr, linestyle="--", linewidth=1, color="black")

    # límites simétricos en X (opcional pero mejora la lectura)
    if np.any(np.isfinite(sub[fc_col].values)):
        max_abs = float(np.nanmax(np.abs(sub[fc_col].values)))
        ax.set_xlim(-max_abs * 1.05, max_abs * 1.05)

    # etiquetas (balanceadas entre up y down)
    if tx_name_map is not None and sig_mask.any() and max_labels > 0:
        n_each = max_labels // 2 if max_labels >= 2 else max_labels

        sig_up = up.sort_values("abs_fc", ascending=False).head(n_each)
        sig_down = down.sort_values("abs_fc", ascending=False).head(n_each)
        sig_points = pd.concat([sig_up, sig_down], axis=0)

        texts = []
        for _, row in sig_points.iterrows():
            tx_id = row["Transcript_ID"]
            label = tx_name_map.get(tx_id, tx_id)
            t = ax.text(
                row[fc_col],
                row["log10p"],
                label,
                fontsize=8,
                ha="center",
                va="bottom",
            )
            texts.append(t)

        if _HAVE_ADJUSTTEXT and texts:
            adjust_text(
                texts,
                ax=ax,
                only_move={"points": "y", "texts": "y"},
                arrowprops=dict(arrowstyle="-", lw=0.5),
            )

    # información N total / N significativos
    n_all = len(sub)
    n_sig = int(sig_mask.sum())
    info_text = f"N = {n_all} (sig = {n_sig})"
    # ax.text(
    #     0.02,
    #     0.98,
    #     info_text,
    #     transform=ax.transAxes,
    #     ha="left",
    #     va="top",
    #     fontsize=9,
    # )

    ax.text(
    0.5,
    -0.28,            # más abajo que la leyenda
    info_text,
    transform=ax.transAxes,
    ha="center",
    va="top",
    fontsize=9,
    )

    # título conectando con resumen a nivel RBP (si se pasa)
    base_title = rbp_name
    if rbp_level_log2fc is not None and np.isfinite(rbp_level_log2fc):
        base_title += f" (RBP-level median log2FC = {rbp_level_log2fc:+.2f})"
    ax.set_title(base_title)

    # ejes
    ax.set_xlabel("Transcript log2FC of DeepRBP score\n(tumor vs normal)")
    ax.set_ylabel(r"$-\log_{10}(\mathrm{FDR})$")

    # estilizar ejes
    ax.tick_params(labelsize=9)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # leyenda centrada a la derecha
    ax.legend(
        frameon=True,
        fancybox=True,
        framealpha=0.9,
        edgecolor="black",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),   # debajo del xlabel, centrado
        ncol=2,                         # queda mejor en una fila horizontal
        fontsize=8,
    )
    # ax.legend(
    #     frameon=True,           # ← activa la caja
    #     fancybox=True,          # esquinas redondeadas (opcional)
    #     framealpha=0.9,         # ligera transparencia
    #     edgecolor="black",      # borde negro fino
    #     bbox_to_anchor=(1.02, 0.5),
    #     loc="center left",
    #     borderaxespad=0.0,
    #     fontsize=8,
    # )

    plt.tight_layout()

    if out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"[volcano-tx] Saved volcano for {rbp_name} to:\n  → {out_path}")

    plt.close()