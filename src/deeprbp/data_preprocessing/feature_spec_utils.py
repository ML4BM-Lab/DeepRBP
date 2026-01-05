
import pandas as pd
from typing import List, Tuple
import warnings

def load_feature_spec(feature_spec_path: str) -> Tuple[List[str], List[str], List[str]]:
    """
    Load ordered RBPs, Genes and Transcripts from a feature-spec Excel file.

    Required sheets:
      - RBPs            → Ensembl_gene_id
      - Genes           → Ensembl_gene_id
      - Transcripts     → Ensembl_transcript_id

    The order in each sheet is preserved and treated as canonical.
    """
    spec = pd.read_excel(feature_spec_path, sheet_name=None)

    required_sheets = ["RBPs", "Genes", "Transcripts"]
    for sheet in required_sheets:
        if sheet not in spec:
            raise ValueError(
                f"[feature_spec] Missing required sheet '{sheet}' in {feature_spec_path}"
            )
    
    def _require_column(df: pd.DataFrame, col: str, sheet: str) -> List[str]:
        if col not in df.columns:
            raise ValueError(
                f"[feature_spec] Sheet '{sheet}' must contain column '{col}'"
            )
        return df[col].astype(str).tolist()

    list_rbps = _require_column(spec["RBPs"], "Ensembl_gene_id", "RBPs")
    list_genes = _require_column(spec["Genes"], "Ensembl_gene_id", "Genes")
    list_trans = _require_column(spec["Transcripts"], "Ensembl_transcript_id", "Transcripts")
    return list_rbps, list_genes, list_trans


def enforce_feature_order(
    df: pd.DataFrame,
    expected_ids: List[str],
    name: str,
    fill_missing: bool = False
) -> pd.DataFrame:
    """
    Enforce exact feature order on a DataFrame indexed by feature IDs.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame indexed by feature IDs (genes / transcripts).
    expected_ids : List[str]
        Canonical ordered list of expected feature IDs.
    name : str
        Name used for error messages.
    fill_missing : bool
        If True, missing features are added as zero rows.

    Returns
    -------
    pd.DataFrame
        Reordered DataFrame with exact feature alignment.
    """
    df = df.copy()
    df.index = df.index.astype(str)

    missing = [x for x in expected_ids if x not in df.index]

    if missing:
        msg = f"[{name}] Missing {len(missing)} features"
        if fill_missing:
            warnings.warn(msg + " → filling with zeros")
            zeros = pd.DataFrame(0, index=missing, columns=df.columns)
            df = pd.concat([df, zeros], axis=0)
        else:
            raise ValueError(msg)
    return df.loc[expected_ids]
