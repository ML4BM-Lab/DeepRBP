# /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/utils_datauser.py

import re
from pathlib import Path
import pandas as pd

# Kallisto utilities (already used in RealKD)
from ..explainability_module.real_knockdowns.preprocess_data.utils_preprocess import (
    load_expression_from_abundance,
)

def load_kallisto_expression( 
    kallisto_output_dir: str,
    sample_ids: list,
    load_counts: bool = True,
) -> dict:
    """
    Load gene- and transcript-level TPM matrices from kallisto output
    and return a dict compatible with clean_expression_data().
    """
    transcript_tpm_dfs = []
    gene_tpm_dfs = []
    transcript_counts_dfs = []
    gene_counts_dfs = []

    missing_samples = []
    kallisto_output_dir = Path(kallisto_output_dir)

    for sample_id in sample_ids:
        abundance_path = kallisto_output_dir / sample_id / "abundance.tsv"

        if not abundance_path.is_file():
            missing_samples.append(sample_id)
            continue

        if load_counts:
            df_trans, df_genes, df_trans_counts, df_genes_counts = load_expression_from_abundance(
                kallisto_output_dir,
                sample_id,
                load_counts,
            )
        else:
            df_trans, df_genes = load_expression_from_abundance(
                kallisto_output_dir,
                sample_id,
                load_counts,
            )

        df_trans.set_index("sample", inplace=True)
        df_genes.set_index("sample", inplace=True)

        transcript_tpm_dfs.append(df_trans)
        gene_tpm_dfs.append(df_genes)

        if load_counts:
            df_trans_counts.set_index("sample", inplace=True)
            df_genes_counts.set_index("sample", inplace=True)

            transcript_counts_dfs.append(df_trans_counts)
            gene_counts_dfs.append(df_genes_counts)

    if len(transcript_tpm_dfs) == 0:
        raise RuntimeError(
            f"No valid kallisto samples found. Missing examples: {missing_samples[:5]}"
        )

    if missing_samples:
        print(
            f"[load_kallisto_expression] ⚠️ Missing abundance.tsv for "
            f"{len(missing_samples)} samples (e.g. {missing_samples[:3]})"
        )

    # Concatenate samples
    df_trans_tpm = pd.concat(transcript_tpm_dfs, axis=1).reset_index()
    df_genes_tpm = pd.concat(gene_tpm_dfs, axis=1).reset_index()

    data = {
        "df_genes": df_genes_tpm,
        "df_trans": df_trans_tpm,
    }

    if load_counts:
        df_trans_counts = pd.concat(transcript_counts_dfs, axis=1).reset_index()
        df_gene_counts  = pd.concat(gene_counts_dfs, axis=1).reset_index()
        data["df_counts"] = df_gene_counts
        data["df_trans_counts"] = df_trans_counts
    return data

def _sanitize(value) -> str:
    """Make filesystem-safe strings."""
    if pd.isna(value):
        return "NA"
    s = str(value).strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^\w\.-]+", "_", s)
    return s if s else "NA"

def build_group_name(group_cols, group_values) -> str:
    """
    Build a readable and filesystem-safe group name
    e.g. tissue_type-Primary__tumor_stage-II
    """
    if not isinstance(group_values, tuple):
        group_values = (group_values,)

    chunks = []
    for col, val in zip(group_cols, group_values):
        chunks.append(f"{col}-{_sanitize(val)}")

    return "__".join(chunks)