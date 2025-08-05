# src/deeprbp/explainability_module/real_knockdowns/preprocess_realkd_data.py

import os 
import pandas as pd
import argparse
import timeit

from .utils_kds import identify_samples_from_tissue, load_expression_from_abundance
from ...data_preprocessing.preprocess_data import (clean_expression_data, filter_transcripts, 
                                                   select_genes_for_modeling, select_rna_binding_proteins,
                                                   transform_expression_data, transpose_dataframes,
                                                   save_processed_data)

def process_data_chunk(
    id_samples: list,
    study_name: str,
    output_dir: str,
    path_dataset: str,
    getBM: pd.DataFrame,
    selected_genes_dir: str,
    splicing_genes_file: str,
    cancer_genes_file: str,
    gene_census_file: str,
    rbp_genes_file: str,
    gene_selection: bool = True) -> None:
    """
    Full preprocessing pipeline for a group of RNA-seq samples. Loads expression data, filters and selects genes/transcripts/RBPs,
    transforms and saves model input matrices.

    Parameters:
    - id_samples (list): List of sample IDs to include in this data chunk.
    - study_name (str): Name of the study (used to name output folder).
    - output_dir (str): Directory to save processed data.
    - path_dataset (str): Path to the root directory of the dataset (contains Kallisto output, info_samples with sample conditions, etc.).
    - getBM (pd.DataFrame): DataFrame with gene and transcript annotations from BioMart.
    - selected_genes_dir (str): Directory containing the list of RNA-binding proteins (RBPs) and other selected genes for modeling.
    - gene_selection (bool): Specifies which gene set to use. Boolean flag to indicate whether gene selection should be performed. 
                            If True, genes related to cancer and alternative splicing are used.
    - gene_transcript_mapping_file (str): Filename for the output that relates transcript IDs and names to gene IDs, names, and their 
                            biotypes for protein coding genes.
    - splicing_genes_file (str): Filename for the Excel file containing genes whose alternative splicing has been characterized to contribute to cancer.
    - cancer_genes_file (str): Filename for the Excel file containing 900 genes predicted to be potential cancer drivers based on mutations and/or 
                            copy number alterations.
    - gene_census_file (str): Filename for the tab-separated values (TSV) file containing the Cancer Gene Census data.
    - rbp_genes_file (str): Filename for the Excel file containing a list of RNA-binding proteins (RBPs).
    Returns:
    - None. Processed expression matrices are saved to disk.
    """
    # Step 1. Merge samples
    transcript_tpm_dfs = [] # List of DataFrames with Transcript_ID as index and TPMs per sample
    gene_tpm_dfs = [] # List of DataFrames with Gene_ID as index and TPMs per sample
    for sample_id in id_samples:    
        df_sample_trans, df_sample_genes = load_expression_from_abundance(os.path.join(path_dataset, 'kallisto_output'), sample_id)
        # Set index for proper alignment during concatenation
        df_sample_trans.set_index('sample', inplace=True)
        df_sample_genes.set_index('sample', inplace=True)
        transcript_tpm_dfs.append(df_sample_trans)
        gene_tpm_dfs.append(df_sample_genes)
    # Concatenate all samples into final expression matrices
    df_trans_tpm = pd.concat(transcript_tpm_dfs, axis=1).reset_index()
    df_genes_tpm = pd.concat(gene_tpm_dfs, axis=1).reset_index()
    data = {
        'df_genes': df_genes_tpm,
        'df_trans': df_trans_tpm
        #'df_counts': df_counts
    }
    # Define and validate the column names used for gene matching.
    # This dictionary contains the default names of the columns used to match genes between
    # the input DataFrames and the getBM DataFrame.
    parm_match_columns = {
        "df_gene_name_col" : "HGNC symbol",    # Column in df_genes_names containing gene names to match.
        "getBM_gene_name_col" : 'Gene_name',   # Column in getBM containing gene names for matching.
        "gene_id_col" : "Gene_ID",             # Column in getBM containing gene IDs corresponding to the gene names.
        "synonyms_col" : "Synonyms",           # Column in df_genes_names containing synonyms for gene names.
    }
    # Step 2: Clean and filter expression data
    data = clean_expression_data(data)  
    getBM_copy = getBM.copy()
    data, getBM_copy = filter_transcripts(data, getBM_copy)  
    # Step 3: Select genes for modeling
    _, list_transcripts, list_genes = select_genes_for_modeling(   
            getBM_copy, 
            selected_genes_dir, 
            splicing_genes_file, 
            cancer_genes_file, 
            gene_census_file,
            parm_match_columns,
            gene_selection
        )
    # Step 4: Filter a list of RNA-binding proteins (RBPs) for modeling   
    list_rbps = select_rna_binding_proteins(selected_genes_dir, rbp_genes_file, getBM, parm_match_columns)
    # Step 5: Obtain the RBP expression, Gene expression and Transcript expression datasets
    output_data = {
        'df_rbp_gene': data['df_genes'].loc[list_rbps, :].copy(),      # get the RBP gene expression dataset
        'df_trans': data['df_trans'].loc[list_transcripts, :].copy(),  # get the transcript expression dataset 
        'df_genes': data['df_genes'].loc[list_genes, :].copy()         # get the dataset for selected genes
    }
    data.update(output_data)
    print(f'[generate_model_input_matrices] RBP matrix shape: {data["df_rbp_gene"].shape}')
    print(f'[generate_model_input_matrices] Transcript matrix shape: {data["df_trans"].shape}')
    print(f'[generate_model_input_matrices] Genes matrix shape: {data["df_genes"].shape}')
    #print(f'[generate_model_input_matrices] Genes matrix count shape: {df_counts.shape}')
    print('\n')
    # Step 6: Transform expression data
    data_transformed = transform_expression_data(data, from_log2p=False)
    # Step 7: Traspose dataframes
    data_transformed = transpose_dataframes(data_transformed)
    # Step 8: Save processed data
    save_processed_data(data = data_transformed, output_dir = output_dir, study_name = study_name)
                                
def process_multiple_conditions(
    df_info_samples: pd.DataFrame,
    conditions: dict,
    output_dir: str,
    path_dataset: str,
    getBM: pd.DataFrame,
    selected_genes_dir: str,
    splicing_genes_file: str,
    cancer_genes_file: str,
    gene_census_file: str,
    rbp_genes_file: str,
    gene_selection: bool = True) -> None:
    """
    Wrapper to process multiple experimental conditions in one go.

    Parameters:
    - conditions (dict): Dictionary mapping study_name -> condition_name in phenotype table.
                         Example: {'knockdown': 'tdp43_ko', 'control': 'Rescued_tdp43'}
    """
    print('[process_multiple_conditions] Starting batch processing of conditions...\n')
    for study_name, condition_name in conditions.items():
        print(f'[process_multiple_conditions] Processing study: {study_name} (condition: {condition_name})')
        id_samples = identify_samples_from_tissue(df_info_samples, condition_name)
        process_data_chunk(
            id_samples=id_samples,
            study_name=study_name,
            output_dir=output_dir,
            path_dataset=path_dataset,
            getBM=getBM,
            selected_genes_dir=selected_genes_dir,
            splicing_genes_file=splicing_genes_file,
            cancer_genes_file=cancer_genes_file,
            gene_census_file=gene_census_file,
            rbp_genes_file=rbp_genes_file,
            gene_selection=gene_selection
        )
    print('\n[process_multiple_conditions] ✅ All conditions processed.\n')

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess KD/control RNA-seq data for model input generation.')
    parser.add_argument('--path_dataset', type=str, required=True,
                        help='Path to the root directory of the dataset (contains kallisto_output and info_samples.txt).')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory to save the processed data.')
    parser.add_argument('--selected_genes_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps', 
                        help='Directory containing the selected genes and RNA-binding proteins (RBPs).')
    parser.add_argument('--gene_transcript_mapping_file', type=str, default='getBM.csv',
                        help='Filename for transcript-to-gene mapping data.')
    parser.add_argument('--splicing_genes_file', type=str, default='Table_S5_Cancer_splicing_gene_eyras.xlsx', 
                        help='Filename for alternative splicing genes.')
    parser.add_argument('--cancer_genes_file', type=str, default='Table_S6_Cancer_gene_eyras.xlsx', 
                        help='Filename for potential cancer driver genes.')
    parser.add_argument('--gene_census_file', type=str, default='Table_Cancer_Gene_Census.tsv', 
                        help='Filename for Cancer Gene Census data.')
    parser.add_argument('--rbp_genes_file', type=str, default='Table_S2_list_RBPs_eyras.xlsx', 
                        help='Filename for RNA-binding proteins (RBPs) list.')
    parser.add_argument('--gene_selection', type=bool, default=True, 
                        help='Boolean flag to indicate whether gene selection should be performed. Otherwise all protein-coding genes with more than one isoform will be used')
    parser.add_argument('--condition_control', type=str, required=True,
                        help='Condition name for control samples (default: Rescued_tdp43).')
    parser.add_argument('--condition_knockdown', type=str, required=True,
                        help='Condition name for knockdown samples (default: tdp43_ko).')
    return parser.parse_args()

def main():
    args = parse_args()
    start_time = timeit.default_timer()

    # Load sample metadata
    df_info_samples = pd.read_csv(os.path.join(args.path_dataset, 'info_samples.txt'), delimiter='\t')
    # Load gene-transcript mapping (getBM)
    getBM_path = os.path.join(args.selected_genes_dir, args.gene_transcript_mapping_file)
    getBM = pd.read_csv(getBM_path)

    # Define experimental conditions
    conditions = {
        'knockdown': args.condition_knockdown,
        'control': args.condition_control
    }

    # Process all conditions in batch
    process_multiple_conditions(
        df_info_samples=df_info_samples,
        conditions=conditions,
        output_dir=args.output_dir,
        path_dataset=args.path_dataset,
        getBM=getBM,
        selected_genes_dir=args.selected_genes_dir,
        splicing_genes_file=args.splicing_genes_file,
        cancer_genes_file=args.cancer_genes_file,
        gene_census_file=args.gene_census_file,
        rbp_genes_file=args.rbp_genes_file,
        gene_selection=args.gene_selection
    )

    # Print elapsed time
    end_time = timeit.default_timer()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\n✅ Total execution time: {int(hours)}h {int(minutes)}m {int(seconds)}s")

if __name__ == '__main__':
    main()



