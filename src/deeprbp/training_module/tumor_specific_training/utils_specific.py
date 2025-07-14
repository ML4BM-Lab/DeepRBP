
import os
import pandas as pd

from .main_specific_vs_general import LIST_TUMOR_TYPES

def collect_cross_tumor_metrics(base_output_dir, metric_name, all_output_dir=None):
    """
    Builds a cross-tumor metric matrix (as a DataFrame) using the specified metric.
    Each row represents the tumor type used for training (ttype),
    and each column represents the tumor type used for evaluation (category).
    Adds an extra row 'all' for the model trained with all tumor types.

    Parameters:
        base_output_dir (str): Base directory containing subdirectories per tumor type.
        metric_name (str): Name of the metric to extract (e.g., 'spearman_corr').
        all_output_dir (str, optional): Directory containing 'test_tumor_category_results.csv' 
                                        for the model trained with all tumor types.
    Returns:
        pd.DataFrame: Cross-tumor metric matrix.
    """
    df_matrix = pd.DataFrame(index=LIST_TUMOR_TYPES, columns=LIST_TUMOR_TYPES, dtype=float)
    # Tumor-specific training rows
    for ttype in LIST_TUMOR_TYPES:
        file_path = os.path.join(base_output_dir, ttype, "test_tumor_category_results.csv")
        if not os.path.exists(file_path):
            print(f"⚠️  Warning: File not found for {ttype}: {file_path}")
            continue
        try:
            df = pd.read_csv(file_path)
            metric_series = df.set_index('category')[metric_name]
            for category, value in metric_series.items():
                if category in LIST_TUMOR_TYPES:
                    df_matrix.loc[ttype, category] = value
        except Exception as e:
            print(f"❌ Error reading {file_path}: {e}")
    # Model trained with all tumor types row
    if all_output_dir is not None:
        all_path = os.path.join(all_output_dir, "test_tumor_category_results.csv")
        if os.path.exists(all_path):
            try:
                df_all = pd.read_csv(all_path)
                metric_series_all = df_all.set_index('category')[metric_name]
                df_matrix.loc['all'] = None  # Inicializa la nueva fila
                for category, value in metric_series_all.items():
                    if category in LIST_TUMOR_TYPES:
                        df_matrix.loc['all', category] = value
            except Exception as e:
                print(f"❌ Error reading {all_path}: {e}")
        else:
            print(f"⚠️  Warning: File not found in all_output_dir: {all_path}")
    df_matrix.index.name = 'Tumor type (train)'
    df_matrix.columns.name = 'Tumor type (eval)'
    return df_matrix


