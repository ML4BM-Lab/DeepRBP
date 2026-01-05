# HAY QUE ACTUALIZAR EN CONSECUENCIA TB el preprocess_realkd_data.py!!!!

# src/deeprbp/data_preprocessing/preprocess_data.py

import os
import copy
import pandas as pd
import numpy as np
from typing import Dict, Union, List
from tqdm import tqdm
import warnings
import argparse
import timeit

from .feature_spec_utils import (
    load_feature_spec,
    enforce_feature_order
)

def preprocessing( 
    raw_data_dir: str, 
    transcript_expression_file: str, 
    gene_expression_file: str,
    gene_counts_file: str,
    transcript_counts_file: str,
    phenotype_data_file: str,
    feature_spec_file: str,
    chunk_size: int,
    output_dir: str,
    ): 
    """
    Preprocess raw RNA-seq matrices into DeepRBP-ready inputs using a feature specification file to
    enforce a fixed feature set and ordering.

    The pipeline loads gene/transcript expression (and expected counts) in sample-wise chunks for
    memory efficiency, cleans identifiers and phenotype metadata, and writes per-study outputs
    (e.g., TCGA and GTEX) with standardized filenames used across the repository.

    Parameters
    ----------
    raw_data_dir : str
        Directory containing the raw input files.
    transcript_expression_file : str
        Transcript-level expression matrix (typically gzipped TSV; features in the `sample` column).
    gene_expression_file : str
        Gene-level expression matrix (typically gzipped TSV; features in the `sample` column).
    gene_counts_file : str
        Gene-level expected counts matrix (typically gzipped TSV).
    transcript_counts_file : str
        Transcript-level expected counts matrix (typically gzipped TSV).
    phenotype_data_file : str
        Sample metadata table (TSV). Must include `sample` (ID), `detailed_category`, and `study`.
    feature_spec_file : str
        Excel feature manifest defining RBPs, genes, and transcripts (and their exact order).
    chunk_size : int
        Number of sample columns processed per iteration.
    output_dir : str
        Output directory where processed files are written (typically into `TCGA/` and `GTEX/` subfolders).

    Returns
    -------
    None
        Results are written to disk under `output_dir`.
    """

    # Load patient IDs and phenotype data
    df_phenotype = pd.read_csv(f"{raw_data_dir}/{phenotype_data_file}", sep='\t', encoding='ISO-8859-1')
    df_phenotype = clean_and_format_phenotype_data(df_phenotype)
    patient_ids = df_phenotype.index.tolist()

    # Load all columns from df_counts once to avoid repeated reads
    available_cols_counts = pd.read_csv(f"{raw_data_dir}/{gene_counts_file}", compression='gzip', sep='\t', nrows=0).columns.tolist()
    available_cols_set = set(available_cols_counts)  # Use a set for faster membership testing

    # Load feature specification  
    list_rbps_spec, list_genes_spec, list_trans_spec = load_feature_spec(feature_spec_file)
    print('[feature_spec]')
    print(f'  RBPs       : {len(list_rbps_spec)}')
    print(f'  Genes      : {len(list_genes_spec)}')
    print(f'  Transcripts: {len(list_trans_spec)}')

    # Initialize list for processed patients
    processed_patients = {"Patient_ID": [], "Chunk": []}
   
    # Process data in chunks
    col_start = 0
    chunk_idx = 0
    total_cols = len(patient_ids)

    with tqdm(total=total_cols, desc="Prepare model inputs") as pbar:
        while col_start < total_cols:
            
            col_end = min(col_start + chunk_size, total_cols)
            selected_cols = patient_ids[col_start:col_end]
          
            if 'sample' not in selected_cols:
                selected_cols.insert(0, 'sample')

            df_genes = pd.read_csv(f"{raw_data_dir}/{gene_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)
            df_trans = pd.read_csv(f"{raw_data_dir}/{transcript_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)
            
            # Filter selected_cols to include only those that exist in df_counts
            valid_selected_cols = [col for col in selected_cols if col in available_cols_set]

            # Ensure 'sample' is at the beginning of valid_selected_cols
            if 'sample' not in valid_selected_cols:
                valid_selected_cols.insert(0, 'sample')

            # Read df_counts using valid selected columns
            df_counts = pd.read_csv(f"{raw_data_dir}/{gene_counts_file}", compression='gzip', sep='\t', usecols=valid_selected_cols)
            df_trans_counts = pd.read_csv(f"{raw_data_dir}/{transcript_counts_file}", compression='gzip', sep='\t', usecols=valid_selected_cols)

            process_data_chunk(
                df_genes, df_trans, df_counts, df_trans_counts, df_phenotype,
                list_rbps_spec, list_genes_spec, list_trans_spec, output_dir)
            
            # Store processed patients
            processed_patients["Patient_ID"].extend(selected_cols[1:])  
            processed_patients["Chunk"].extend([chunk_idx] * (len(selected_cols) - 1))  # Remove 'sample' 
            chunk_idx += 1

            # Update the progress bar
            pbar.update(len(selected_cols) - 1)  # Use the actual number of processed patients
            col_start = col_end

    # Convert processed patients to a DataFrame
    df_processed_patients = pd.DataFrame(processed_patients)

    # Save the processed patients to an Excel file in the output directory
    output_excel_path = f"{output_dir}/processed_patients.xlsx"
    df_processed_patients.to_excel(output_excel_path, index=False)
    print(f"Processed patients saved to {output_excel_path}")

# ----------------------------
# Phenotype + expression cleaning
# ----------------------------
def clean_and_format_phenotype_data(dictionary_data: pd.DataFrame) -> pd.DataFrame:
    """
    Cleans and standardizes phenotype data by performing the following transformations:
    
    - Sets the 'sample' column as the index.
    - Replaces hyphens ('-') with spaces in the dataset.
    - Replaces parentheses ('(' and ')') with underscores ('_').
    - Replaces spaces with underscores in the data entries.
    - Collapses multiple consecutive underscores ('__', '___', etc.) into a single underscore ('_').
    - Strips leading and trailing spaces from column names, and replaces spaces with underscores.
    - Replaces tissue/disease names for 'GTEX' study in 'primary_disease_or_tissue' column based on patterns.

    Args:
    - dictionary_data (pd.DataFrame): A DataFrame containing phenotype data.

    Returns:
    - pd.DataFrame: A cleaned and standardized version of the input phenotype data.
    """
    dictionary_data = dictionary_data[dictionary_data['detailed_category'].notna()]
    dictionary_data = dictionary_data.set_index("sample")

    obj_cols = dictionary_data.select_dtypes(include=["object"]).columns
    dictionary_data[obj_cols] = dictionary_data[obj_cols].apply(
        lambda x: x.str.replace("-", " ", regex=False)
    )
    dictionary_data[obj_cols] = dictionary_data[obj_cols].apply(
        lambda x: x.str.replace(r"[\(\)]", "_", regex=True)
    )
    dictionary_data[obj_cols] = dictionary_data[obj_cols].apply(
        lambda x: x.str.replace(" ", "_", regex=False)
    )
    dictionary_data[obj_cols] = dictionary_data[obj_cols].apply(
        lambda x: x.str.replace("_{2,}", "_", regex=True)
    )

    dictionary_data.columns = [
        col.strip().replace(" ", "_") if not col.startswith("_")
        else col.strip().replace(" ", "_")[1:]
        for col in dictionary_data.columns
    ]

    # Collapse some GTEX categories
    replacements = {
        '^Skin_.+': 'Skin',
        '^Brain_.+': 'Brain',
        '^Adipose_.+': 'Adipose',
        '^Artery_.+': 'Artery',
        '^Colon_.+': 'Colon',
        '^Cervix_.+': 'Cervix',
        '^Esophagus_.+': 'Esophagus'
    }
    
    if "study" in dictionary_data.columns and "primary_disease_or_tissue" in dictionary_data.columns:
        mask = dictionary_data["study"] == "GTEX"
        dictionary_data.loc[mask, "primary_disease_or_tissue"] = dictionary_data.loc[
            mask, "primary_disease_or_tissue"
        ].replace(replacements, regex=True)
    return dictionary_data

def clean_expression_data(data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Cleans the gene and transcript expression data by removing genome version and aggregating loci.
    
    Args:
    - data (dict): Dictionary containing 'df_genes', 'df_trans', 'df_counts' DataFrames.
    
    Returns:
    - dict: Dictionary with cleaned gene, transcript, and count DataFrames.
    """
    print('Cleaning expression data...')

    # Internal function to clean and aggregate data
    def clean_df(df):
        df['sample'] = df['sample'].str.rsplit('.').str[0].str.replace('R', '0')
        return df.iloc[:, 1:].groupby(df['sample']).sum()
    
    # Apply the cleaning function to each DataFrame in the dictionary
    for key in data.keys():
        data[key] = clean_df(data[key])
    print('[clean_expression_data] Cleaning done\n')
    return data

def generate_model_input_matrices(
    data: dict, 
    df_phenotype: pd.DataFrame, 
    study_name: str, 
    list_rbps: list, 
    list_transcripts: list, 
    list_genes: list
    ) -> tuple[dict, pd.DataFrame]:
    """
    Generate model input matrices for a given study (TCGA / GTEX).

    Feature spaces are treated independently:
      - RBPs        → gene-level regulatory inputs
      - Genes       → gene-level outputs / context
      - Transcripts → transcript-level targets
      - Counts      → raw counts (genes + transcripts)

    All feature selection and sample intersection is handled here.
 
    This function filters and organizes expression data into separate matrices for RBP, transcript, and gene expression based on the provided study name. 
    It ensures that only samples common to both the phenotype and expression data are included, enabling accurate modeling for downstream analysis.
    """
    if "study" not in df_phenotype.columns:
        raise ValueError(
            "[phenotype] Expected column 'study' in phenotype_metadata "
            "to split TCGA vs GTEX."
        )
    print(f'[generate_model_input_matrices] Generating input matrices for study: {study_name}...')
    df_phenotype_study = df_phenotype[df_phenotype.study == study_name].copy()

    if df_phenotype_study.empty:
        print(f"[generate_model_input_matrices] No samples for {study_name}")
        return {}, df_phenotype_study

    # --- sample intersection (deterministic order)---
    samples_expr = sorted(
        set(df_phenotype_study.index)
        & set(data['df_genes'].columns)
        & set(data['df_trans'].columns)
    )

    samples_counts = sorted(
        set(df_phenotype_study.index)
        & set(data['df_counts'].columns)
        & set(data['df_trans_counts'].columns) 
    )

    print(f"  • Expression samples: {len(samples_expr)}")
    print(f"  • Count samples     : {len(samples_counts)}")

    # --- Slice samples first (never features yet) ---
    df_genes_raw = data["df_genes"].loc[:, samples_expr]
    df_trans_raw = data["df_trans"].loc[:, samples_expr]
    
    df_gene_counts = data["df_counts"].loc[:, samples_counts]
    df_trans_counts = data["df_trans_counts"].loc[:, samples_counts]

    # --- Build feature spaces (INDEPENDENT) ---
    # RBPs (regulatory input space)
    df_rbp_gene = enforce_feature_order(
        df_genes_raw,
        list_rbps,
        name="RBPs",
        fill_missing=False,
    )

    # Genes (gene-level space)
    df_genes = enforce_feature_order(
        df_genes_raw,
        list_genes,
        name="Genes",
        fill_missing=False,
    )

    # Transcripts (transcript-level space)
    df_trans = enforce_feature_order(
        df_trans_raw,
        list_transcripts,
        name="Transcripts",
        fill_missing=False,
    )

    # Final output
    output_data = {
        "df_rbp_gene": df_rbp_gene,       # (n_RBPs × samples)
        "df_genes": df_genes,             # (n_genes × samples)
        "df_trans": df_trans,             # (n_transcripts × samples)
        "df_counts": df_gene_counts,      # (genes × samples)
        "df_trans_counts": df_trans_counts,  # (transcripts × samples)
    }

    print(f"  • RBPs        : {df_rbp_gene.shape}")
    print(f"  • Genes       : {df_genes.shape}")
    print(f"  • Transcripts : {df_trans.shape}")
    print(f"  • Gene counts : {df_gene_counts.shape}")
    print(f"  • Trans counts: {df_trans_counts.shape}")
    print()
    return output_data, df_phenotype_study

def transform_expression_data(data: dict, 
                              from_log2p: bool = True,
                              epsilon: float = 0.001,
                              epsilon_counts: float = 1.0,
                              counts_are_log2p: bool = True) -> dict:
    """
    Transforms expression data (TPM or log2p(TPM)) for RBPs, transcripts, genes, and optionally counts.

    Parameters:
    - data (dict): Dictionary containing expression DataFrames. Expected keys include:
        - 'df_rbp_gene', 'df_trans', 'df_genes' (required)
        - 'df_counts' (optional)
    - from_log2p (bool): Whether to first transform from log2p(TPM + epsilon) to TPM. Default is True.
    - epsilon (float): Value that was added before log2. Used in inverse log2p. Default is 0.001.
    - epsilon_counts (float): Value used in log2p(COUNT + epsilon_counts) for counts inversion.
    - counts_are_log2p (bool): Whether df_counts is in log2p. If False, no transformation is applied.

    Returns:
    - dict: Updated dictionary with transformed expression matrices.
    """
    print('[transform_expression_data] Starting transformation...')
    data = copy.deepcopy(data)
    # Define which keys to transform normally (TPM/log2p)
    expression_keys = [k for k in data.keys() if k != 'df_counts']
    # Step 1: If from_log2p, convert back to TPM
    if from_log2p:
        print('[transform_expression_data] Inverting log2p(TPM + ε) to TPM...')
        print(f'  • TPM from_log2p: {from_log2p}  (epsilon = {epsilon})')
        for key in expression_keys:
            print(f'  ↪️ Inverting {key}')
            data[key] = np.power(2, data[key]) - epsilon
    # Step 2: Clip negative values (just in case)
    for key in expression_keys:
        data[key] = data[key].clip(lower=0)
    # Step 3: Apply log2p(x + 1) to df_trans y df_rbp_gene
    for key in ['df_trans', 'df_rbp_gene']:
        if key in data:
            print(f'  🔁 Applying log2p(x + 1) to {key}')
            data[key] = np.log2(data[key] + 1)
    # Step 4: Process df_counts if present
    if 'df_counts' in data:
        if counts_are_log2p:
            print(f'[transform_expression_data] Inverting log2p(COUNTS + {epsilon_counts})...')
            data['df_counts'] = np.power(2, data['df_counts']) - epsilon_counts
            data['df_counts'] = data['df_counts'].clip(lower=0)
        else:
            print('[transform_expression_data] df_counts is already in raw count scale (no inverse log2p).')
        print('[transform_expression_data] Rounding and converting df_counts to integers...')
        data['df_counts'] = data['df_counts'].round().astype(int)
    print('[transform_expression_data] ✅ Transformation complete.\n')
    return data

def transpose_dataframes(data: dict) -> dict:
    """
    Transposes the expression DataFrames to have patients (samples) as index and genes/transcripts as columns.

    Parameters:
    - data (dict): Dictionary containing expression DataFrames.

    Returns:
    - dict: Dictionary with transposed DataFrames.
    """
    print('[transpose_dataframes] Transposing expression data...')
    data = data.copy()
    for key in ['df_rbp_gene', 'df_trans', 'df_genes', 'df_counts']:
        if key in data:
            print(f'  🔁 Transposing {key}')
            data[key] = data[key].T
    print('[transpose_dataframes] ✅ Transposition complete.\n')
    return data

def save_processed_data(
    data: dict,
    output_dir: str,
    df_phenotype_study: pd.DataFrame = None,
    study_name: str = None
) -> None:
    """
    Saves processed expression matrices and optional phenotype metadata to CSV files.

    Parameters:
    - data (dict): Dictionary containing expression matrices (e.g. df_rbp_gene, df_trans, df_genes, df_counts).
    - output_dir (str): Base directory to save files.
    - df_phenotype_study (pd.DataFrame, optional): Phenotype metadata to save (default: None).
    - study_name (str, optional): Name of the study to create a subdirectory (default: None).
    """
    print('[save_processed_data] Saving processed data...')
    # Determine output path
    path = os.path.join(output_dir, study_name) if study_name else output_dir
    os.makedirs(path, exist_ok=True)
    # Mapping of keys to filenames
    save_map = {
        'df_rbp_gene': 'RBPs_log2p_tpm.csv',
        'df_trans': 'trans_log2p_tpm.csv',
        'df_genes': 'gn_tpm.csv',
        'df_counts': 'gn_counts.csv',
        'df_trans_counts': 'trans_counts.csv'    
    }
    # Save each DataFrame if it exists
    for key, filename in save_map.items():
        if key in data:
            file_path = os.path.join(path, filename)
            print(f'[save_processed_data] Saving {key} to {file_path}...')
            data[key].to_csv(file_path, mode='a', header=not os.path.exists(file_path))
        else:
            print(f'[save_processed_data] Skipping {key} (not found in data).')
    # Save phenotype metadata if provided
    if df_phenotype_study is not None and 'df_rbp_gene' in data:
        path_phenotype = os.path.join(path, 'phenotype_metadata.csv')
        print(f'[save_processed_data] Saving phenotype metadata to {path_phenotype}...')
        df_phenotype_study.loc[data['df_rbp_gene'].index].to_csv(
            path_phenotype, mode='a', header=not os.path.exists(path_phenotype)
        )
    elif df_phenotype_study is None:
        print('[save_processed_data] Skipping phenotype metadata (not provided).')
    print('[save_processed_data] ✅ Data saving completed.\n')

def process_data_chunk(  
    df_genes: pd.DataFrame,
    df_trans: pd.DataFrame,
    df_counts: pd.DataFrame,
    df_trans_counts: pd.DataFrame,
    df_phenotype: pd.DataFrame,
    list_rbps_spec: list,
    list_genes_spec: list,
    list_trans_spec: list,
    output_dir: str) -> None:
    """
    Processes data chunks including gene and transcript expression, phenotype data, and gene selection to generate the model inputs

    Args:
    - df_genes (DataFrame): DataFrame containing gene expression data.
    - df_trans (DataFrame): DataFrame containing transcript expression data.
    - df_counts (DataFrame): DataFrame containing gene counts data.
    - df_phenotype (DataFrame): DataFrame containing phenotype data, such as detailed category, 
      primary disease, primary site, sample type, gender, and study information.
   
    Returns:
    - None
    """
    data = {
        'df_genes': df_genes,
        'df_trans': df_trans,
        'df_counts': df_counts,
        'df_trans_counts': df_trans_counts
    }

    # Step 1: clean expression data ONLY
    data = clean_expression_data(data)
    
    # Step 2: generate per-study inputs
    for study_name in ['TCGA', 'GTEX']:
        print(f"Processing study: {study_name}")
        data_study, df_phenotype_study = generate_model_input_matrices(
            data=data, 
            df_phenotype=df_phenotype, 
            study_name=study_name, 
            list_rbps=list_rbps_spec, 
            list_genes=list_genes_spec,
            list_transcripts=list_trans_spec
        )

        if not data_study:
            continue

        data_transformed = transform_expression_data(data_study)
        data_transformed = transpose_dataframes(data_transformed)

        save_processed_data(
            data = data_transformed, 
            df_phenotype_study = df_phenotype_study, 
            output_dir = output_dir, 
            study_name = study_name)
        print('\n')

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess gene and transcript expression data for model input generation.')
    parser.add_argument('--raw_data_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/raw', 
                        help='Directory containing the raw data files.')
    parser.add_argument('--transcript_expression_file', type=str, default='TcgaTargetGtex_rsem_isoform_tpm.gz', 
                        help='Filename for transcript-level expression data.')
    parser.add_argument('--gene_expression_file', type=str, default='TcgaTargetGtex_rsem_gene_tpm.gz', 
                        help='Filename for gene-level expression data.')
    parser.add_argument('--gene_counts_file', type=str, default='TcgaTargetGTEX_gene_expected_count.gz', 
                        help='Filename for gene-level counts data.')
    parser.add_argument('--transcript_counts_file', type=str, default='TcgaTargetGtex_isoform_expected_count.gz', 
                        help='Filename for gene-level counts data.')
    parser.add_argument('--phenotype_data_file', type=str, default='TcgaTargetGTEX_phenotype.txt', 
                        help='Filename for phenotype information.')
    parser.add_argument('--feature_spec_file', type=str, required=True, help='Excel file defining ordered RBPs, Genes and Transcripts')                   
    parser.add_argument('--chunk_size', type=int, default=1000, 
                        help='Number of rows to process per chunk for memory efficiency.')
    parser.add_argument('--output_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/processed', 
                        help='Directory to save the processed data.')
    return parser.parse_args()

def main():
    args = parse_args()
    start_time = timeit.default_timer()
    
    preprocessing(
        raw_data_dir=args.raw_data_dir, 
        transcript_expression_file=args.transcript_expression_file, 
        gene_expression_file=args.gene_expression_file, 
        gene_counts_file=args.gene_counts_file, 
        transcript_counts_file=args.transcript_counts_file,
        phenotype_data_file=args.phenotype_data_file,
        feature_spec_file=args.feature_spec_file,
        chunk_size=args.chunk_size, 
        output_dir=args.output_dir, 
    )
    
    end_time = timeit.default_timer()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Total execution time: {int(hours)}h {int(minutes)}m {int(seconds)}s")
    
if __name__ == '__main__':
    main()