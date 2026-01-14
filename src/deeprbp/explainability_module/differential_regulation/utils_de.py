
# src/deeprbp/explainability_module/differential_regulation/utils_de.py

import pandas as pd

def _filter_counts_by_feature_spec(
    counts: pd.DataFrame,
    feature_spec_path: str,
    sheet_name: str = "RBPs"
) -> pd.DataFrame:
    """
    Filter counts matrix (samples x genes) to RBPs using DeepRBP feature spec.

    The Excel file is expected to contain a sheet named 'RBPs' with
    a column 'Ensembl_gene_id'.
    """
    spec = pd.read_excel(feature_spec_path, sheet_name=sheet_name)

    if "Ensembl_gene_id" not in spec.columns:
        raise ValueError(
            f"Feature spec sheet '{sheet_name}' must contain column 'Ensembl_gene_id'. "
            f"Found: {list(spec.columns)}"
        )

    rbp_genes = (
        spec["Ensembl_gene_id"]
        .astype(str)
        .unique()
        .tolist()
    )

    common = sorted(set(rbp_genes).intersection(counts.columns))

    if not common:
        raise ValueError(
            "No overlapping genes between RBP feature spec and counts matrix"
        )

    return counts.loc[:, common]

def _filter_meta(
    meta: pd.DataFrame,
    sample_category_col: str,
    disease_condition_col: str,
    select_category: str,
    cond_a: str,
    cond_b: str) -> pd.DataFrame:
    
    if sample_category_col not in meta.columns:
        raise ValueError(f"Metadata missing sample_category column: {sample_category_col}. Found: {list(meta.columns)}")
    if disease_condition_col not in meta.columns:
        raise ValueError(f"Metadata missing disease_condition column: {disease_condition_col}. Found: {list(meta.columns)}")

    m = meta.copy()
    m = m[m[sample_category_col] == select_category]

    if m.empty:
        raise ValueError(
            f"No samples found for select_category='{select_category}' using column '{sample_category_col}'."
        )

    m = m[m[disease_condition_col].isin([cond_a, cond_b])]
    if m.empty:
        raise ValueError(
            f"No samples left after filtering conditions '{cond_a}'/'{cond_b}' in column '{disease_condition_col}'."
        )
    return m


def _export_counts_for_r(counts_samples_x_features: pd.DataFrame, out_path: str) -> None:
    """
    R runner expects: features x samples, with first column as Feature_ID.
    We'll write TSV for speed and to avoid comma/quote edge cases.
    """
    # transpose to features x samples
    mat = counts_samples_x_features.T
    mat.insert(0, "Feature_ID", mat.index)
    mat.to_csv(out_path, sep="\t", index=False)


def _export_meta_for_r(meta: pd.DataFrame, out_path: str) -> None:
    """
    R runner expects Sample_ID as first column.
    """
    m = meta.copy()
    m.insert(0, "Sample_ID", m.index)
    m.to_csv(out_path, sep="\t", index=False)