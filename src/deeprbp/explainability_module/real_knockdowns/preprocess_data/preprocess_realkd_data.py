# src/deeprbp/explainability_module/real_knockdowns/preprocess_data/preprocess_realkd_data.py

# 22-dic esto hay que volver a probarlo porque es nueva version del code y utilizalo para la info nueva de Maria!
import os 
import pandas as pd
import argparse
import timeit

from .utils_preprocess import identify_samples_from_tissue, load_expression_from_abundance
from ....data_preprocessing.feature_spec_utils import load_feature_spec
from ....data_preprocessing.preprocess_data import (
    clean_expression_data,
    generate_model_input_matrices,
    transform_expression_data,
    transpose_dataframes,
    save_processed_data,
)

def process_data_chunk(
    id_samples: list,
    study_name: str,
    output_dir: str,
    path_dataset: str,
    list_rbps: list,
    list_genes: list,
    list_transcripts: list) -> None:
    """
    Preprocess a group of RNA-seq samples corresponding to one experimental
    condition (e.g. knockdown or control).

    This function:
      1) loads Kallisto TPMs for each sample,
      2) cleans expression matrices,
      3) builds DeepRBP-compatible inputs using the feature specification,
      4) applies the same transformations as in training,
      5) saves the resulting matrices to disk.

    Parameters
    ----------
    id_samples : list
        Sample identifiers to include.
    study_name : str
        Label for the experimental condition (used as output subfolder).
    output_dir : str
        Base directory where processed data will be written.
    path_dataset : str
        Root dataset directory (must contain `kallisto_output/`).
    list_rbps : list
        Ordered list of RBP Ensembl gene IDs.
    list_genes : list
        Ordered list of gene Ensembl IDs.
    list_transcripts : list
        Ordered list of transcript Ensembl IDs.

    Returns
    -------
    None
        Processed matrices are written to disk.
    """
   
    # --------------------------------------------------
    # 1) Load expression from kallisto
    # --------------------------------------------------
    transcript_tpm_dfs = [] # List of DataFrames with Transcript_ID as index and TPMs per sample
    gene_tpm_dfs = [] # List of DataFrames with Gene_ID as index and TPMs per sample
    
    for sample_id in id_samples:    
        df_trans, df_genes = load_expression_from_abundance(
            os.path.join(path_dataset, 'kallisto_output'),
            sample_id
        )
        df_trans.set_index('sample', inplace=True)
        df_genes.set_index('sample', inplace=True)

        transcript_tpm_dfs.append(df_trans)
        gene_tpm_dfs.append(df_genes)
        
    df_trans_tpm = pd.concat(transcript_tpm_dfs, axis=1).reset_index()
    df_genes_tpm = pd.concat(gene_tpm_dfs, axis=1).reset_index()

    data = {
        'df_genes': df_genes_tpm,
        'df_trans': df_trans_tpm,
    }
        
    # --------------------------------------------------
    # 2) Clean expression data
    # --------------------------------------------------
    data = clean_expression_data(data)

    # --------------------------------------------------
    # 3) Build model inputs (SAME PATTERN AS TRAINING)
    # --------------------------------------------------
    # Fake phenotype (all samples belong to this "study")
    df_phenotype = pd.DataFrame(
        {'study': study_name},
        index=id_samples
    )
    df_phenotype.index = df_phenotype.index.astype(str)

    data_study, _ = generate_model_input_matrices(
        data=data,
        df_phenotype=df_phenotype,
        study_name=study_name,
        list_rbps=list_rbps,
        list_genes=list_genes,
        list_transcripts=list_transcripts,
    )

    if not data_study:
        return

    # --------------------------------------------------
    # 4) Transform + transpose
    # --------------------------------------------------
    # Kallisto TPMs are not log2p
    data_transformed = transform_expression_data(data_study, from_log2p=False)
    data_transformed = transpose_dataframes(data_transformed)
    
    # --------------------------------------------------
    # 5) Save
    # --------------------------------------------------
    save_processed_data(data = data_transformed, output_dir = output_dir, study_name = study_name)
                                
def process_multiple_conditions(
    df_info_samples: pd.DataFrame,
    conditions: dict,
    output_dir: str,
    path_dataset: str,
    list_rbps: list,
    list_genes: list,
    list_transcripts: list,
    ) -> None:
    """
    Process multiple experimental conditions (e.g. KD vs control)
    in a single run.

    Parameters
    ----------
    df_info_samples : pd.DataFrame
        Sample metadata table.
    conditions : dict
        Mapping {study_name -> condition_name}.
    output_dir : str
        Base directory for processed outputs.
    path_dataset : str
        Root dataset directory.
    list_rbps, list_genes, list_transcripts : list
        Feature specification lists.

    Returns
    -------
    None
    """
    print("[process_multiple_conditions] Starting batch processing...\n")
    for study_name, condition_name in conditions.items():
        print(
            f"[process_multiple_conditions] "
            f"Processing {study_name} (condition: {condition_name})"
        )
        id_samples = identify_samples_from_tissue(df_info_samples, condition_name)
        
        process_data_chunk(
            id_samples=id_samples,
            study_name=study_name,
            output_dir=output_dir,
            path_dataset=path_dataset,
            list_rbps=list_rbps,
            list_genes=list_genes,
            list_transcripts=list_transcripts,
        )
    print('\n[process_multiple_conditions] ✅ All conditions processed.\n')

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess KD/control RNA-seq data for model input generation.')
    parser.add_argument('--path_dataset', type=str, required=True,
                        help='Path to the root directory of the dataset (contains kallisto_output and info_samples.txt).')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory to save the processed data.')
    parser.add_argument('--feature_spec', type=str, required=True, help='Excel file defining ordered RBPs, Genes and Transcripts')
    parser.add_argument('--condition_control', type=str, required=True,
                        help='Condition name for control samples (default: Rescued_tdp43).')
    parser.add_argument('--condition_knockdown', type=str, required=True,
                        help='Condition name for knockdown samples (default: tdp43_ko).')
    return parser.parse_args()

def main():
    """
    Entry point for the RealKD preprocessing pipeline.
    """
    args = parse_args()
    start_time = timeit.default_timer()

    # --------------------------------------------------
    # Load metadata
    # --------------------------------------------------
    df_info_samples = pd.read_csv(
        os.path.join(args.path_dataset, 'info_samples.txt'), 
        sep='\t',
    ) # before delimiter
    
    # --------------------------------------------------
    # Load feature specification
    # --------------------------------------------------
    list_rbps, list_genes, list_transcripts = load_feature_spec(args.feature_spec)
    print('[feature_spec]')
    print(f'  RBPs       : {len(list_rbps)}')
    print(f'  Genes      : {len(list_genes)}')
    print(f'  Transcripts: {len(list_transcripts)}')
    
    # --------------------------------------------------
    # Define experimental conditions
    # --------------------------------------------------
    conditions = {
        'knockdown': args.condition_knockdown,
        'control': args.condition_control
    }

    # --------------------------------------------------
    # Run preprocessing
    # --------------------------------------------------
    process_multiple_conditions(
        df_info_samples=df_info_samples,
        conditions=conditions,
        output_dir=args.output_dir,
        path_dataset=args.path_dataset,
        list_rbps=list_rbps,
        list_genes=list_genes,
        list_transcripts=list_transcripts,
    )

    end_time = timeit.default_timer()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\n✅ Total execution time: {int(hours)}h {int(minutes)}m {int(seconds)}s")

if __name__ == '__main__':
    main()