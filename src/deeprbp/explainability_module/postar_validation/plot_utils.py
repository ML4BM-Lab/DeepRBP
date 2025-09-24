# src/deeprbp/explainability_module/postar_validation/plot_utils.py

import os, re
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

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
    q_hi: float = 0.98,          # cuantíl para fijar el límite X superior del panel principal
    bw_adjust: float = 1.3, #0.7        # suavizado KDE (↑ más suave, ↓ más picudo)
    show_density_yticks: bool = False
):
    """
    KDE + ROC con fondo blanco y sin rejilla. Muestra threshold y medianas por clase.

    - Recorta eje X por cuantíl `q_hi` para que outliers no aplasten la vista.
    - Clip a [0, hi] y cut=0 para impedir “colas” negativas en la KDE.
    - `bw_adjust` controla el suavizado de la KDE.
    """
    required = {"Score", "Postar_Score"}
    if not required.issubset(df_current_rbp.columns):
        raise ValueError(f"DataFrame must contain columns: {sorted(required)}")

    # Nombre legible: usa RBP_name si existe; si no, el ID
    if rbp_display_name is None:
        if "RBP_name" in df_current_rbp.columns:
            names = df_current_rbp["RBP_name"].dropna().astype(str)
            rbp_display_name = names.mode().iat[0] if not names.empty else rbp_id
        else:
            rbp_display_name = rbp_id

    # Colores (tus originales)
    if colors is None:
        colors = {1: "#7fc97f", 0: "#beaed4"}

    # Datos base
    scores = np.asarray(df_current_rbp["Score"], dtype=float)
    m1 = (df_current_rbp["Postar_Score"] == 1).to_numpy(bool)
    m0 = (df_current_rbp["Postar_Score"] == 0).to_numpy(bool)
    n1, n0 = int(m1.sum()), int(m0.sum())

    # Rango X: no-negativo, recortado por cuantíl; asegura que entre el threshold
    lo = 0.0
    hi_all = float(np.nanmax(scores) if scores.size else 1.0)
    hi_main = float(np.nanquantile(scores, q_hi)) if np.isfinite(hi_all) else hi_all
    if np.isfinite(optimal_threshold):
        hi_main = max(hi_main, float(optimal_threshold) * 1.05)
    hi_main = min(hi_main, max(hi_all, 1e-9))

    # Medianas de cada clase
    med1 = float(np.nanmedian(scores[m1])) if n1 > 0 else np.nan
    med0 = float(np.nanmedian(scores[m0])) if n0 > 0 else np.nan

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
    
    # ======= FIGURA =======
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.2), constrained_layout=True)
    fig.patch.set_facecolor("white")
    for ax in axes:
        ax.set_facecolor("white")
        ax.grid(False)  # asegurar sin rejilla

    # ----- Panel Izquierdo: KDE por clase -----
    ax = axes[0]
    labels = {1: f"Binding (n={n1})", 0: f"Not binding (n={n0})"}

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

    # Threshold + medianas
    if np.isfinite(optimal_threshold):
        ax.axvline(optimal_threshold, color="red", linestyle="--", linewidth=TH_LW, label=f"Th={optimal_threshold:.2f}")
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

    # Leyenda compacta
    handles, texts = ax.get_legend_handles_labels()
    has_th = any(lbl.startswith("Th=") for lbl in texts)
    if not has_th and np.isfinite(optimal_threshold):
        from matplotlib.patches import Patch
        handles.append(Patch(facecolor="none", edgecolor="red", linestyle="--"))
        texts.append(f"Th={optimal_threshold:.2f}")
    leg = ax.legend(handles, texts, title="Postar", frameon=True)
    leg.get_title().set_fontsize(LEGEND_FZ)
    for txt in leg.get_texts():
        txt.set_fontsize(LEGEND_FZ)

    # ----- Panel Derecho: ROC -----
    ax = axes[1]
    ax.plot(fpr, tpr, lw=ROC_LW, label=f"AUC = {auc_score:.2f}")
    if 0 <= optimal_idx < len(fpr):
        ax.scatter(float(fpr[optimal_idx]), float(tpr[optimal_idx]),
                   s=DOT_SIZE, color="red", zorder=3, label=f"Th={optimal_threshold:.2f}")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=RAND_LW, label="Random")
    ax.set_xlabel("False Positive Rate", fontsize=LABEL_FZ)
    ax.set_ylabel("True Positive Rate", fontsize=LABEL_FZ)
    ax.set_title("ROC curve", fontsize=TITLE_FZ)
    ax.tick_params(axis="both", labelsize=TICK_FZ)
    leg2 = ax.legend(frameon=True)
    leg2.get_title().set_fontsize(LEGEND_FZ)
    for txt in leg2.get_texts():
        txt.set_fontsize(LEGEND_FZ)

    # Supertítulo
    fig.suptitle(f"RBP: {rbp_display_name} ({rbp_id})", fontsize=SUPTITLE_FZ)

    # Guardado
    os.makedirs(path_save, exist_ok=True)
    out_path = os.path.join(path_save, f"figure_{_sanitize(rbp_display_name)}.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



# def plot_distributions_and_roc_with_thresholds(
#     df_current_rbp,
#     rbp_id,
#     optimal_threshold,
#     fpr,
#     tpr,
#     optimal_idx,
#     auc_score,
#     path_save,
#     left_panel: str = "hist",          # "hist" | "kde" | "both"
#     rbp_display_name: Optional[str] = None,  
#     colors: Optional[dict] = None,
#     score_label: str = "Explainability score",
#     x_max: Optional[float] = None,            # tope duro del eje X (panel izq.)
#     x_max_percentile: float = 99.0,           # si x_max=None, usar percentil global
#     show_inset: bool = True                   # inset con el rango completo
#     ):
#     """
#     Plots the distributions of scores and the ROC curve with the optimal threshold.

#     Args:
#         df_current_rbp (DataFrame): DataFrame containing the absolute scores and POSTAR labels.
#         rbp_id (str): The ID of the RNA Binding Protein (RBP).
#         optimal_threshold (float): The optimal threshold for classification.
#         fpr (array-like): False positive rates for the ROC curve.
#         tpr (array-like): True positive rates for the ROC curve.
#         optimal_idx (int): Index of the optimal threshold in the fpr and tpr arrays.
#         auc_score (float): Calculated auc score between a particular RBP Postar_Score and RBP explainability score.
#         path_save (str): Path where the figure will be saved.
#         left_panel: "hist" (recomendado), "kde" o "both"
#         rbp_display_name: si None, intenta usar df['RBP_name'] y cae a rbp_id
#         colors: dict como {1: "#...", 0: "#..."}; si None usa buenos default
#     """
#     required = {'Score', 'Postar_Score'}
#     if not required.issubset(df_current_rbp.columns):
#         raise ValueError(f"DataFrame must contain columns: {sorted(required)}")
#     # --- nombre a mostrar: intenta usar RBP_name del DF; si no, cae al id ---
#     if rbp_display_name is None:
#         if 'RBP_name' in df_current_rbp.columns:
#             names = df_current_rbp['RBP_name'].dropna().astype(str)
#             rbp_display_name = names.mode().iat[0] if not names.empty else rbp_id
#         else:
#             rbp_display_name = rbp_id

#     # NOTA: usamos los scores tal cual llegan (ya en valor absoluto en tu pipeline)
#     scores = df_current_rbp["Score"].to_numpy(dtype=float)
#     # Conteos por clase para la leyenda y medianas por clase
#     m1 = (df_current_rbp["Postar_Score"] == 1).to_numpy(bool)
#     m0 = (df_current_rbp["Postar_Score"] == 0).to_numpy(bool)
#     n1, n0 = int(m1.sum()), int(m0.sum())
    
#     # rango completo y rango de visualización
#     lo_full = 0.0
#     hi_full = float(np.nanmax(scores) if scores.size else 1.0)
#     if x_max is not None:
#         hi_disp = float(x_max)
#     else:
#         hi_disp = float(np.nanpercentile(scores, x_max_percentile)) if np.isfinite(hi_full) else 1.0
#         hi_disp = max(1e-8, min(hi_disp, hi_full))  # evita hi=0 y no exceder el máx real

#     bins_disp = _fd_bins(scores, lo_full, hi_disp)
#     bins_full = _fd_bins(scores, lo_full, hi_full)

#     # conteos recortados (para avisar)
#     clipped_1 = int((scores[m1] > hi_disp).sum())
#     clipped_0 = int((scores[m0] > hi_disp).sum())

#     # medianas por clase (sobre todo el rango)
#     med1 = float(np.nanmedian(scores[m1])) if n1 > 0 else np.nan
#     med0 = float(np.nanmedian(scores[m0])) if n0 > 0 else np.nan

#     # estilo
#     sns.set_style("whitegrid")
#     fig, axes = plt.subplots(1, 2, figsize=(12, 4))

#     # -------- Panel izquierdo (distribución) ----------
#     ax = axes[0]
#     color_binding   = "#7fc97f"  # Class-1
#     color_notbind   = "#beaed4"  # Class-0
#     labels = {1: f"Binding (n={n1})", 0: f"Not binding (n={n0})"}

#     def _plot_hist(ax_, mask, color, label, bins):
#         x = np.clip(scores[mask], lo_full, hi_disp)  # recorta a [0, hi_disp]
#         if x.size == 0:
#             return
#         sns.histplot(x=x, bins=bins, stat="density",
#                      element="step", fill=True, alpha=0.35,
#                      edgecolor="k", linewidth=0.2,
#                      color=color, ax=ax_, label=label)

#     def _plot_kde(ax_, mask, color, label):
#         x = np.clip(scores[mask], lo_full, hi_disp)
#         if x.size == 0:
#             return
#         sns.kdeplot(x=x, ax=ax_,
#                     fill=(left_panel == "kde"), alpha=0.5 if left_panel == "kde" else 1.0,
#                     clip=(lo_full, hi_disp), cut=0, warn_singular=False,
#                     color=color, label=label)

#     if left_panel in ("hist", "both"):
#         _plot_hist(ax, m1, color_binding,   labels[1] if left_panel != "both" else None, bins_disp)
#         _plot_hist(ax, m0, color_notbind,   labels[0] if left_panel != "both" else None, bins_disp)
#     if left_panel in ("kde", "both"):
#         _plot_kde(ax,  m1, color_binding,   labels[1] if left_panel == "kde" else None)
#         _plot_kde(ax,  m0, color_notbind,   labels[0] if left_panel == "kde" else None)

#     # umbral y medianas
#     ax.axvline(optimal_threshold, color="red", linestyle="--", linewidth=1.5)
#     if np.isfinite(med1):
#         ax.axvline(med1, color=color_binding, linestyle="-", linewidth=1, alpha=0.9)
#     if np.isfinite(med0):
#         ax.axvline(med0, color=color_notbind, linestyle="-", linewidth=1, alpha=0.9)

#     ax.set_xlim(lo_full, hi_disp)
#     ax.set_xlabel(score_label)
#     ax.set_ylabel("Density")
#     ax.set_title("Score distribution by POSTAR label")

#     # leyenda compacta (evita duplicados)
#     if left_panel == "both":
#         handles = [
#             Patch(facecolor=color_binding, edgecolor="k", alpha=0.35, label=labels[1]),
#             Patch(facecolor=color_notbind, edgecolor="k", alpha=0.35, label=labels[0]),
#         ]
#         ax.legend(handles=handles, title="Postar", frameon=True)
#     else:
#         ax.legend(title="Postar", frameon=True)

#     # anotación de recorte si aplica
#     if clipped_1 or clipped_0:
#         msg = f"clipped > {hi_disp:.2g}: 1→{clipped_1}, 0→{clipped_0}"
#         ax.text(0.99, 0.98, msg, transform=ax.transAxes,
#                 ha="right", va="top", fontsize=9,
#                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, lw=0.0))

#     # inset con rango completo (opcional)
#     if show_inset and hi_full > hi_disp:
#         axins = inset_axes(ax, width="40%", height="40%", loc="upper right", borderpad=1.0)
#         _plot_hist(axins, m1, color_binding, None, bins_full)
#         _plot_hist(axins, m0, color_notbind, None, bins_full)
#         axins.axvline(optimal_threshold, color="red", linestyle="--", linewidth=1.0)
#         axins.set_xlim(lo_full, hi_full)
#         axins.set_ylim(0, None)
#         axins.set_xticks([])
#         axins.set_yticks([])
#         axins.set_title("full range", fontsize=9, pad=2)

#     # med1 = float(np.nanmedian(scores[m1])) if n1 > 0 else np.nan
#     # med0 = float(np.nanmedian(scores[m0])) if n0 > 0 else np.nan
    
#     # # Asegura rango no-negativo en el eje X (evita colas negativas visuales)
#     # clip_lo, clip_hi = 0.0, float(np.nanmax(scores) if scores.size else 1.0)
#     # bins = _fd_bins(scores, clip_lo, clip_hi)

#     # # ===== Panel izquierdo: histograma binned (evita colas “negativas” del KDE) =====
#     # # Reemplaza NaN por 0.0 solo para poder graficar sin romper el hist/KDE (no deberían haberlos igualmente).
#     # # Forzamos el rango del eje X a 0,𝑚𝑎𝑥(𝑠𝑐𝑜𝑟𝑒) para reflejar que trabajamos con valores absolutos.
#     # # Con _fd_bins damos una binned view robusta.
#     # # Separación por clases POSTAR (0/1)
    
#     # sns.set_style("whitegrid")
#     # fig, axes = plt.subplots(1, 2, figsize=(12, 4))
#     # # --- Distribución por clases ---
#     # ax = axes[0]
#     # color_group1 = "#7fc97f"  # Binding
#     # color_group0 = "#beaed4"  # Not binding
#     # labels = {1: f"Binding (n={n1})", 0: f"Not binding (n={n0})"}

#     # if left_panel in ("hist", "both"):
#     #     # Hist por clase (densidad), barras rellenas con contorno suave
#     #     for cls, color in ((1, color_group1), (0, color_group0)):
#     #         mask = m1 if cls == 1 else m0
#     #         x = np.clip(scores[mask], clip_lo, clip_hi)
#     #         if x.size == 0:
#     #             continue
#     #         sns.histplot(
#     #             x=x,
#     #             bins=bins,
#     #             stat="density",
#     #             element="step",
#     #             fill=True,
#     #             alpha=0.35,
#     #             edgecolor="k",
#     #             linewidth=0.2,
#     #             color=color,
#     #             ax=ax,
#     #             label=labels[cls] if left_panel != "both" else None,
#     #         )

#     # if left_panel in ("kde", "both"):
#     #     # KDE acotada (cut=0, clip=(0, hi)) para evitar colas negativas
#     #     for cls, color in ((1, color_group1), (0, color_group0)):
#     #         mask = m1 if cls == 1 else m0
#     #         x = np.clip(scores[mask], clip_lo, clip_hi)
#     #         if x.size == 0:
#     #             continue
#     #         sns.kdeplot(
#     #             x=x,
#     #             ax=ax,
#     #             fill=(left_panel == "kde"),
#     #             alpha=0.5 if left_panel == "kde" else 1.0,
#     #             clip=(clip_lo, clip_hi),
#     #             cut=0,
#     #             warn_singular=False,
#     #             color=color,
#     #             label=labels[cls] if left_panel == "kde" else None,
#     #         )

#     # # Umbral óptimo
#     # ax.axvline(optimal_threshold, color="red", linestyle="--", linewidth=1.5)

#     # # Medianas de cada clase
#     # if not np.isnan(med1):
#     #     ax.axvline(med1, color=color_group1, linestyle="-", linewidth=1, alpha=0.9)
#     # if not np.isnan(med0):
#     #     ax.axvline(med0, color=color_group0, linestyle="-", linewidth=1, alpha=0.9)

#     # ax.set_xlim(clip_lo, clip_hi)
#     # ax.set_xlabel(score_label)
#     # ax.set_ylabel('Density')
#     # ax.set_title('Score distribution by POSTAR label')
#     # #ax.legend(title='Postar')

#     # # Leyenda compacta
#     # if left_panel == "both":
#     #     # construimos leyenda manual para evitar duplicados
#     #     handles = [
#     #         Patch(facecolor=color_group1, edgecolor="k", alpha=0.35, label=labels[1]),
#     #         Patch(facecolor=color_group0, edgecolor="k", alpha=0.35, label=labels[0]),
#     #     ]
#     #     ax.legend(handles=handles, title="Postar", frameon=True)
#     # else:
#     #     ax.legend(title="Postar", frameon=True)

#     # --- ROC ---
#     ax = axes[1]
#     ax.plot(fpr, tpr, lw=2, label=f'AUC = {auc_score:.2f}')
#     ax.scatter(float(fpr[optimal_idx]), float(tpr[optimal_idx]),
#                s=30, color='red', zorder=3, label=f'Th={optimal_threshold:.2f}')
#     ax.plot([0, 1], [0, 1], '--', color='gray', lw=1, label='Random')
#     ax.set_xlabel('False Positive Rate')
#     ax.set_ylabel('True Positive Rate')
#     ax.set_title('ROC curve')
#     ax.legend(frameon=True)
#     # Supertítulo profesional con gene name y el ID entre paréntesis
#     fig.suptitle(f'RBP: {rbp_display_name} ({rbp_id})', fontsize=13)
#     fig.tight_layout()
#     os.makedirs(path_save, exist_ok=True)
#     out_path = os.path.join(path_save, f'figure_{_sanitize(rbp_display_name)}.png')
#     fig.savefig(out_path, dpi=300, bbox_inches='tight', transparent=False)
#     plt.close(fig)


#def plot_distributions_and_roc_with_thresholds(df_current_rbp, rbp_id, optimal_threshold, fpr, tpr, optimal_idx, auc_score, path_save):
#     """
#     Plots the distributions of scores and the ROC curve with the optimal threshold.

#     Args:
#         df_current_rbp (DataFrame): DataFrame containing the absolute scores and POSTAR labels.
#         rbp_id (str): The ID of the RNA Binding Protein (RBP).
#         optimal_threshold (float): The optimal threshold for classification.
#         fpr (array-like): False positive rates for the ROC curve.
#         tpr (array-like): True positive rates for the ROC curve.
#         optimal_idx (int): Index of the optimal threshold in the fpr and tpr arrays.
#         auc_score (float): Calculated auc score between a particular RBP Postar_Score and RBP explainability score.
#         path_save (str): Path where the figure will be saved.
#     """
#     # Validate input DataFrame
#     if not {'Score', 'Postar_Score'}.issubset(df_current_rbp.columns):
#         raise ValueError("DataFrame must contain 'Score' and 'Postar_Score' columns.")
#     plt.figure(figsize=(12, 4))
#     sns.set(style="whitegrid")
#     color_group1 = '#7fc97f'   
#     color_group0 = '#beaed4'  
#     # Plot distribution of 0s and 1s
#     plt.subplot(1, 2, 1)
#     sns.kdeplot(data=df_current_rbp, x='Score', hue='Postar_Score', fill=True, 
#                 palette={1: color_group1, 0: color_group0}, common_norm=False)
#     plt.axvline(x=optimal_threshold, color='red', linestyle='--', 
#                 label=f'Threshold = {optimal_threshold:.2f}')
#     plt.title(f'Distribution of 0s and 1s for RBP: {rbp_id}')
#     plt.xlabel('Scores')
#     plt.ylabel('Density')
#     plt.legend(title='Postar', labels=['Class-1', 'Class-0'])
#     plt.grid(False)
#     # Plot ROC curve
#     plt.subplot(1, 2, 2)
#     plt.plot(fpr, tpr, label=f'AUC = {auc_score:.2f}')
#     plt.scatter(fpr[optimal_idx], tpr[optimal_idx], marker='o', color='red', 
#                 label=f'Threshold = {optimal_threshold:.2f}')
#     plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random')
#     plt.title(f'ROC Curve for RBP: {rbp_id}')
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     plt.legend()
#     plt.grid(False)
#     plt.tight_layout()
#     # Save the figure
#     plt.savefig(f'{path_save}/figure_{rbp_id}.png', transparent=True)
#     plt.close()
