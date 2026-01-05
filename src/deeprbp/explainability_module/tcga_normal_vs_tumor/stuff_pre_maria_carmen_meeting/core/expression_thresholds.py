# src/deeprbp/explainability_module/tcga_normal_vs_tumor/rbp_expression_tcga.py

import pandas as pd
from typing import Tuple, List
import numpy as np
from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks

def get_rbp_expr_df(dataset) -> pd.DataFrame:
    """
    Convierte dataset.features['scaled_rbp_df'] (B x R) en un DataFrame
    con columnas = rbp_names y filas = muestras S.
    """
    arr = dataset.features['scaled_rbp_df'].numpy()   # (B, R)
    cols = dataset.rbp_names                          # lista de ENSG...
    return pd.DataFrame(arr, columns=cols)

def classify_rbps_by_expression(
    ds_normal,
    ds_tumor,
    epsilon: float = None
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], List[str]]:
    """
    Devuelve:
    - rbps_both: RBPs expresados en ambos grupos
    - rbps_tumor_only: RBPs expresados solo en tumor

    'Expresado' se define como media > epsilon.
    """
    df_rbp_norm  = get_rbp_expr_df(ds_normal)
    df_rbp_tumor = get_rbp_expr_df(ds_tumor)
    mean_norm  = df_rbp_norm.mean(axis=0)
    mean_tumor = df_rbp_tumor.mean(axis=0)
    if epsilon is None:
        df_all = pd.concat([df_rbp_norm, df_rbp_tumor], axis=0)
        epsilon = estimate_epsilon_from_scaled(df_all)
        print(f"[classify_rbps] epsilon óptimo (data-driven): {epsilon:.4f}")
    rbps_both = mean_norm[(mean_norm > epsilon) & (mean_tumor > epsilon)].index.tolist()
    rbps_tumor_only = mean_tumor[(mean_tumor > epsilon) & (mean_norm <= epsilon)].index.tolist()
    return rbps_both, rbps_tumor_only, epsilon

def estimate_epsilon_from_scaled(df_scaled, bandwidth=0.02):
    """
    df_scaled: DataFrame con valores escalados+clippeados [0,1].
    Devuelve epsilon basado en el mínimo entre los dos modos principales > 0.
    """
    vals = df_scaled.values.flatten()
    xs = np.linspace(0, 1, 500)[:, None]
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    kde.fit(vals[:, None])
    logdens = kde.score_samples(xs)
    peaks, _ = find_peaks(logdens)
    # ordenar picos por posición y escoger el de x≈0 y el siguiente
    # (primer pico cercano a 0, segundo > 0.1, por ejemplo)
    peak_xs = xs[peaks].ravel()
    order = np.argsort(peak_xs)
    peaks = peaks[order]
    # pico cerca de 0
    i0 = peaks[0]
    # primer pico con x > 0.1
    candidates = [p for p, x in zip(peaks[1:], peak_xs[1:]) if x > 0.1]
    if not candidates:
        return 0.05
    i1 = candidates[0]
    valley_idx = np.argmin(logdens[i0:i1]) + i0
    return float(xs[valley_idx])