# src/deeprbp/explainability_module/corum_complex_analysis/utils_complex.py

import pandas as pd
import numpy as np
from typing import Iterable, Dict, List, Sequence, Tuple
from collections import Counter
from pathlib import PurePath
from string import ascii_uppercase

from ...training_module.tcga_codes import get_tcga_code, TCGA_CODE

def split_to_list(cell):
    """
    Split a semicolon-separated string into a list of trimmed strings.
    Returns an empty list if input is NaN.
    """
    if pd.isna(cell):
        return []
    return [item.strip() for item in cell.split(';')]

def map_gene_names_to_ids(row, getBM):
    """
    Map each gene name in 'subunits_gene_name' to its corresponding Gene_ID from getBM.
    If no direct match is found, try matching using each synonym (splitting by comma).
    Returns a list of Gene_IDs or None if no match is found.
    """
    gene_ids = []
    print(f"Processing complex: {row.get('complex_name', 'Unknown')}")
    for gene_name, synonym_str in zip(row['subunits_gene_name'], row['subunits_gene_name_synonyms']):
        print(f"  Trying gene name: {gene_name}")
        # Buscar match directo
        match = getBM[getBM['Gene_name'] == gene_name]
        if not match.empty:
            gene_id = match['Gene_ID'].values[0]
            print(f"    Found Gene_ID: {gene_id} for gene name: {gene_name}\n")
            gene_ids.append(gene_id)
        else:
            print(f"    Gene name '{gene_name}' not found, trying synonyms: {synonym_str}")
            # Probar cada sinónimo separado por coma
            synonyms = [s.strip() for s in synonym_str.split(',')]
            found = False
            for syn in synonyms:
                match_syn = getBM[getBM['Gene_name'] == syn]
                if not match_syn.empty:
                    gene_id = match_syn['Gene_ID'].values[0]
                    print(f"    Found Gene_ID: {gene_id} for synonym: {syn}\n")
                    gene_ids.append(gene_id)
                    found = True
                    break
            if not found:
                print(f"    No Gene_ID found for gene name '{gene_name}' or any synonym\n")
                gene_ids.append(None)
    return gene_ids

def build_expanded_corr(
    corr_scores: pd.DataFrame,
    df_corum: pd.DataFrame,
    complex_id_col: str = "complex_id",
    rbp_list_col: str = "subunits_gene_id",
    complexes_order: Sequence[int] = (8369, 8370, 8371, 8372, 8391),  # A..E
    include_non_family: bool = True,  # F
    sort_within_groups: bool = True) -> pd.DataFrame:
    """
    Returns the expanded correlation matrix reordered by block rows/cols:
      Row/col A: [A-A | A-B | A-C | A-D | A-E | (A-F if include_non_family)]
      ...
      Row/col F (if include_non_family): [F-A | F-B | F-C | F-D | F-E | F-F]

    NOTES:
    - Genes present in multiple complexes are NOT removed; they may appear in multiple blocks.
    - Only genes present in 'corr_scores' are kept.
    """
    if set(corr_scores.index) != set(corr_scores.columns):
        raise ValueError("The matrix must be square with identical row/column IDs.")
    universe = set(corr_scores.index)
    # Map A..E from complex_id in the given order
    group_names = ["A", "B", "C", "D", "E"]
    groups: Dict[str, List[str]] = {}
    for name, cid in zip(group_names, complexes_order):
        row = df_corum.loc[df_corum[complex_id_col] == cid, rbp_list_col]
        if row.empty:
            genes = []
        else:
            # Keep group membership, but filter to genes present in the matrix
            genes = [g for g in row.iloc[0] if g in universe]
        if sort_within_groups:
            genes = sorted(genes)
        groups[name] = genes
    # F = genes not assigned to any of the A..E complexes
    if include_non_family:
        in_any = set(g for lst in groups.values() for g in lst)
        F_genes = [g for g in corr_scores.index if g not in in_any]
        if sort_within_groups:
            F_genes = sorted(F_genes)
        groups["F"] = F_genes
    # Final ordering by blocks (A..E and optionally F)
    ordered_labels: List[str] = []
    for name in ["A", "B", "C", "D", "E"] + (["F"] if include_non_family else []):
        ordered_labels.extend(groups.get(name, []))
    # Reindex rows and columns with that order (labels may repeat)
    expanded_corr = corr_scores.loc[ordered_labels, ordered_labels]
    return expanded_corr, groups

def extract_block_contrasts(
    expanded_corr: pd.DataFrame,
    groups: Dict[str, List[str]],
    group_key: str = "A",
    deduplicate_nonkey_labels: bool = True) -> Tuple[List[float], List[float]]:
    """
    Returns (key_values, nonkey_values) for the selected group.

    IMPORTANT:
    - Rows (KEY group) are ALWAYS deduplicated by label, keeping the FIRST
      occurrence only. This behavior is intentional and not configurable. 

    Definitions:
    - key_values: upper-triangle (k=1, no diagonal) of the 'KEY vs KEY' block
                  using ROWS = unique KEY RBPs and COLS = unique KEY RBPs.

    - nonkey_values: ALL values from the 'KEY vs NON-KEY' block using:
        * ROWS  = unique KEY RBPs (first occurrence kept)
        * COLS  = gene_ids NOT in KEY.
                  If deduplicate_nonkey_labels=True (default), NON-KEY labels
                  are also deduplicated by keeping their FIRST occurrence only.

    Notes:
    - The matrix must be square with identical row/column IDs.
    - Works even if labels repeat; filtering is by label (gene_id).
    """
    if set(expanded_corr.index) != set(expanded_corr.columns):
        raise ValueError("The matrix must be square with identical row/column IDs.")
    if group_key not in groups:
        raise ValueError(f"'{group_key}' is not in 'groups'.")
    key_labels = set(groups[group_key])
    # --- ROWS (KEY group) ---
    # Always deduplicate KEY rows by label: keep first occurrence only.
    row_labels = expanded_corr.index.to_list()
    row_is_key = pd.Series(row_labels).isin(key_labels).to_numpy()
    seen_rows = set()
    row_pos_key: List[int] = []
    for i, (lab, ok) in enumerate(zip(row_labels, row_is_key)):
        if ok and lab not in seen_rows:
            row_pos_key.append(i)
            seen_rows.add(lab)
    row_pos_key = np.array(row_pos_key, dtype=int)
    # --- COLUMNS (NON-KEY group) ---
    col_labels = expanded_corr.columns.to_list()
    col_is_nonkey = ~pd.Series(col_labels).isin(key_labels).to_numpy()
    if deduplicate_nonkey_labels:
        # Deduplicate NON-KEY columns by label: keep first occurrence only.
        seen_cols = set()
        col_pos_nonkey: List[int] = []
        for j, (lab, ok) in enumerate(zip(col_labels, col_is_nonkey)):
            if ok and lab not in seen_cols:
                col_pos_nonkey.append(j)
                seen_cols.add(lab)
        col_pos_nonkey = np.array(col_pos_nonkey, dtype=int)
    else:
        col_pos_nonkey = np.where(col_is_nonkey)[0]
    # --- KEY vs KEY → upper triangle (k=1) ---
    if row_pos_key.size >= 2:
        kk_block = expanded_corr.to_numpy()[np.ix_(row_pos_key, row_pos_key)]
        iu = np.triu_indices(kk_block.shape[0], k=1)
        key_values = kk_block[iu].astype(float).tolist()
    else:
        key_values = []
    # --- KEY vs NON-KEY → flattened ---
    if row_pos_key.size >= 1 and col_pos_nonkey.size >= 1:
        kn_block = expanded_corr.to_numpy()[np.ix_(row_pos_key, col_pos_nonkey)]
        nonkey_values = kn_block.ravel().astype(float).tolist()
    else:
        nonkey_values = []
    return key_values, nonkey_values

def _infer_tumor_from_output_dir(output_dir: str) -> tuple[str, str]:
    """
    Devuelve (tcga_code, long_name) intentando encontrar un componente de la ruta
    que exista en TCGA_CODE. Si no lo encuentra, usa el padre del directorio.
    """
    parts = list(PurePath(output_dir).parts)
    long_name = None
    for p in reversed(parts):
        if p in TCGA_CODE:        # coincide con las claves tipo 'Acute_Myeloid_Leukemia'
            long_name = p
            break
    if long_name is None:
        long_name = PurePath(output_dir).parent.name  # p.ej., 'Acute_Myeloid_Leukemia'
    tcga = get_tcga_code(long_name)                   # p.ej., 'LAML'
    return tcga, long_name

def build_design_matrix_with_F(emb_df: pd.DataFrame, df_corum: pd.DataFrame,
                               col='subunits_gene_id', family_labels=None) -> pd.DataFrame:
    """
    Build a design matrix (RBP × Families) with an extra 'F' (others) column.

    emb_df: DataFrame (X×R) with columns = ENSG IDs of RBPs
    df_corum[col]: column with list-like gene IDs per family (list/tuple/ndarray)
    family_labels: optional list of names for families; defaults to A,B,C,...
    Returns: DataFrame (R × (F+1)) with columns = [families..., 'F'] and {0,1}
    """
    # 1) Family labels
    fam_ids = list(df_corum.index)
    F = len(fam_ids)
    if family_labels is None:
        base = list(ascii_uppercase)
        family_labels = [base[i] if i < len(base) else f"Family_{i}" for i in range(F)]
    assert len(family_labels) == F, "family_labels must match number of families (rows in df_corum)."
    # 2) Normalize to lists (inline, no extra helper)
    series_lists = []
    for x in df_corum[col].tolist():
        if x is None:
            series_lists.append([])
        elif isinstance(x, list):
            series_lists.append(x)
        elif isinstance(x, tuple):
            series_lists.append(list(x))
        elif isinstance(x, np.ndarray):
            series_lists.append(x.tolist())
        else:
            # treat any other scalar/string as single-element membership
            series_lists.append([x])
    # 3) Map ENSG -> set of family indices
    gene2fams = {}
    for fam_idx, genes in zip(fam_ids, series_lists):
        for g in genes:
            g = str(g).strip()
            if g:
                gene2fams.setdefault(g, set()).add(fam_idx)
    # 4) Base matrix without 'F'
    rbps = emb_df.columns.astype(str)
    design = pd.DataFrame(0, index=rbps, columns=family_labels, dtype=np.int8)
    famid_to_label = {fid: family_labels[i] for i, fid in enumerate(fam_ids)}
    for g in rbps:
        for fid in gene2fams.get(g, ()):
            design.at[g, famid_to_label[fid]] = 1
    # 5) Add 'F' (others): 1 if RBP not assigned to any family
    design['F'] = (design.sum(axis=1) == 0).astype(np.int8)
    return design
