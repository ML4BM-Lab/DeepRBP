
# src/deeprbp/data_preprocessing/preprocess_data_user.py

import argparse
import timeit
import pandas as pd

from .feature_spec_utils import load_feature_spec
from .preprocess_data import (
    clean_expression_data,
    generate_model_input_matrices,
    transform_expression_data,
    transpose_dataframes,
    save_processed_data,
)

from .utils_datauser import (
    build_group_name,
    load_kallisto_expression,
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare DeepRBP input matrices from kallisto output (user datasets)."
    )
    parser.add_argument("--kallisto_output", required=True)
    parser.add_argument("--metadata_csv", required=True)
    parser.add_argument("--feature_spec", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--group_col", default=None, help="Optional column used to group samples. If not provided, all samples are processed together.")
    parser.add_argument("--split_by_tumor_stage", action="store_true")
    parser.add_argument("--split_by_cell_line", action="store_true") 
    return parser.parse_args()

# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()
    start_time = timeit.default_timer()

    # 1) Load metadata
    metadata = pd.read_csv(args.metadata_csv)
    metadata["Run"] = metadata["Run"].astype(str)

    # 2) Define grouping columns
    if args.group_col is None:
        group_cols = None
    else:
        group_cols = [args.group_col]
        if args.split_by_tumor_stage:
            group_cols.append("tumor_stage")
        if args.split_by_cell_line:
            group_cols.append("cell_line")

    print("[preprocess_data_user]")
    print(f"  kallisto_output : {args.kallisto_output}")
    print(f"  metadata_csv    : {args.metadata_csv}")
    print(f"  feature_spec    : {args.feature_spec}")
    print(f"  output_dir      : {args.output_dir}")
    print(f"  group_cols      : {group_cols}")

    # 3) Load feature spec
    list_rbps, list_genes, list_transcripts = load_feature_spec(args.feature_spec)
    print("[feature_spec]")
    print(f"  RBPs       : {len(list_rbps)}")
    print(f"  Genes      : {len(list_genes)}")
    print(f"  Transcripts: {len(list_transcripts)}")

    # 4) Iterate over groups
    if group_cols is None:
        # Single dataset
        process_one_group(
            sample_ids=metadata["Run"].tolist(),
            group_name="ALL",
            kallisto_output_dir=args.kallisto_output,
            output_dir=args.output_dir,
            list_rbps=list_rbps,
            list_genes=list_genes,
            list_transcripts=list_transcripts,
            phenotype_df=metadata,
        )
    else:
        # Group-wise processing
        for group_values, meta_group in metadata.groupby(group_cols, dropna=False):
            group_name = build_group_name(group_cols, group_values)
            sample_ids = meta_group["Run"].tolist()

            process_one_group(
                sample_ids=sample_ids,
                group_name=group_name,
                kallisto_output_dir=args.kallisto_output,
                output_dir=args.output_dir,
                list_rbps=list_rbps,
                list_genes=list_genes,
                list_transcripts=list_transcripts,
                phenotype_df=meta_group,
            )

    # Timing
    elapsed = timeit.default_timer() - start_time
    print(f"\n✅ Finished in {elapsed/60:.1f} min")

def process_one_group(
    sample_ids: list,
    group_name: str,
    kallisto_output_dir: str,
    output_dir: str,
    list_rbps: list,
    list_genes: list,
    list_transcripts: list,
    phenotype_df: pd.DataFrame,
) -> None:
    """
    End-to-end preprocessing for one user-defined group.

    Preprocess a user-defined group of RNA-seq samples quantified with kallisto
    into DeepRBP-compatible input matrices.

    This function:
      1) loads transcript- and gene-level TPMs from kallisto `abundance.tsv` files,
      2) cleans and standardizes expression matrices,
      3) builds DeepRBP input matrices using the feature specification,
      4) applies the same transformations used during training,
      5) saves the resulting matrices and phenotype metadata to disk.

    Parameters
    ----------
    sample_ids : list
        Sample identifiers to include (run accessions).
    group_name : str
        Name of the group being processed (used as output subfolder and study label).
    kallisto_output_dir : str
        Directory containing kallisto output folders (one per sample).
    output_dir : str
        Base directory where processed files will be written.
    list_rbps : list
        Ordered list of RBP Ensembl gene IDs.
    list_genes : list
        Ordered list of gene Ensembl IDs.
    list_transcripts : list
        Ordered list of transcript Ensembl IDs.
    phenotype_df : pandas.DataFrame
        Sample metadata for the current group (must contain a `Run` column).

    Returns
    -------
    None
        Processed matrices are written to disk.
    """
    print(f"\n[process_one_group] ▶ {group_name}")
    print(f"  Samples requested: {len(sample_ids)}")

    # --------------------------------------------------
    # 1) Load kallisto expression
    # --------------------------------------------------
    data = load_kallisto_expression(
        kallisto_output_dir=kallisto_output_dir,
        sample_ids=sample_ids,
        load_counts=True
    )

    # --------------------------------------------------
    # 2) Clean expression data
    # --------------------------------------------------
    data = clean_expression_data(data)

    # --------------------------------------------------
    # 3) Prepare phenotype table
    # --------------------------------------------------
    df_phenotype = phenotype_df.copy()
    df_phenotype["Run"] = df_phenotype["Run"].astype(str)
    df_phenotype = df_phenotype.set_index("Run")
    df_phenotype["study"] = group_name
    
    # --------------------------------------------------
    # 4) Generate model input matrices
    # --------------------------------------------------
    data_study, df_pheno_study = generate_model_input_matrices(
        data=data,
        df_phenotype=df_phenotype,
        study_name=group_name,
        list_rbps=list_rbps,
        list_genes=list_genes,
        list_transcripts=list_transcripts,
    )

    if not data_study:
        print(f"[process_one_group] ⚠️ No valid samples after intersection, skipping.")
        return

    # --------------------------------------------------
    # 5) Transform expression data
    # --------------------------------------------------
    # Kallisto TPMs are NOT log2p
    data_transformed = transform_expression_data(
        data_study,
        from_log2p=False,
        counts_are_log2p=False,
    )

    # --------------------------------------------------
    # 6) Transpose (samples as rows)
    # --------------------------------------------------
    data_transformed = transpose_dataframes(data_transformed)

    # --------------------------------------------------
    # 7) Save outputs
    # --------------------------------------------------
    save_processed_data(
        data=data_transformed,
        output_dir=output_dir,
        df_phenotype_study=df_pheno_study,
        study_name=group_name,
    )

    print(f"[process_one_group] ✅ Finished {group_name}")

if __name__ == "__main__":
    main()