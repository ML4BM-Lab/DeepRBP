# src/deeprbp/explainability_module/real_knockdowns/preprocess_data/utils_preprocess.py

import pandas as pd

def identify_samples_from_tissue(df_info_samples, condition):  
    """
    Identify sample IDs from a given condition using the sample_title column.

    Parameters:
    - df_info_samples (DataFrame): DataFrame with 'run_accession' and 'sample_title' columns.
    - condition (str): The condition to filter on (must match a value in 'sample_title').

    Returns:
    - id_samples (list): List of sample IDs (run_accession) matching the condition.
    """
    print(f"[identify_samples_from_tissue] Selected condition: {condition}")
    if condition in df_info_samples['sample_title'].unique():
        id_samples = df_info_samples[df_info_samples['sample_title'] == condition]['run_accession'].tolist()
    else:
        print(f"[identify_samples_from_tissue] Unknown condition: {condition}")
        id_samples = []
    print(f"[identify_samples_from_tissue] Found samples: {id_samples}\n")
    return id_samples

def clean_abundance_dataframe(df_abundance):
    """
    Cleans and parses a Kallisto abundance DataFrame by:
    - Splitting the target_id field into transcript/gene metadata
    - Removing version numbers from Transcript_ID and Gene_ID
    - Reordering columns for clarity

    Parameters:
    - df_abundance (pd.DataFrame): Raw DataFrame read from abundance.tsv

    Returns:
    - df_clean (pd.DataFrame): Cleaned and structured DataFrame
    """
    print('[clean_abundance_dataframe] Parsing target_id into annotated fields...')
    df_abundance[['Transcript_ID', 'Gene_ID', 'OTTHUMG', 'OTTHUMT', 'Gene_Symbol', 'Gene_Name', 'Transcript_Length', 'Biotype']] = (
        df_abundance['target_id'].str.rstrip('|').str.split('|', expand=True)
    )
    print('[clean_abundance_dataframe] Dropping unused columns and reordering...')
    df_abundance.drop(columns=['target_id'], inplace=True)
    # Reorder columns logically
    df_abundance = df_abundance[
        ['Transcript_ID', 'Gene_ID', 'OTTHUMT', 'OTTHUMG', 'Gene_Symbol', 'Gene_Name', 'Biotype',
         'length', 'eff_length', 'est_counts', 'Transcript_Length', 'tpm']
    ]
    print('[clean_abundance_dataframe] Removing transcript/gene version suffixes...')
    df_abundance['Transcript_ID'] = df_abundance['Transcript_ID'].str.split('.').str[0]
    df_abundance['Gene_ID'] = df_abundance['Gene_ID'].str.split('.').str[0]
    return df_abundance

def load_expression_from_abundance(path, sample_id, load_counts: bool = False): 
    """
    Load transcript- and gene-level expression vectors (TPM) from a sample's Kallisto abundance file.

    Parameters:
    - path (str): Path to the dataset directory.
    - sample_id (str): ID of the sample (i.e. run_accession folder name).

    Returns:
    - df_trans_tpm (pd.DataFrame): DataFrame with transcript-level TPM values.
    - df_genes_tpm (pd.DataFrame): DataFrame with gene-level TPM values.
    """
    # Read the abundance.tsv file
    df_abundance = pd.read_csv(f'{path}/{sample_id}/abundance.tsv', sep='\t')
    print(f'[{load_expression_from_abundance.__name__}] ✔️ Loaded abundance.tsv for sample {sample_id}')
    
    # Clean and parse the raw abundance data to extract transcript/gene-level metadata and structure the DataFrame
    df_abundance = clean_abundance_dataframe(df_abundance)
    print(f'[{load_expression_from_abundance.__name__}] ✔️ Cleaned and structured abundance data for sample {sample_id}')
    
    # Create transcript-level TPM dataframe
    df_trans_tpm = pd.DataFrame(
         df_abundance['tpm'].values,
         index=df_abundance['Transcript_ID'].values,
         columns=[sample_id]
    ).reset_index().rename(columns={'index': 'sample'})
    print(f'[{load_expression_from_abundance.__name__}] ✔️ Created transcript-level TPM DataFrame for sample {sample_id}')
    
    # Create gene-level TPM dataframe
    tpm_sum_by_gene = df_abundance.groupby('Gene_ID')['tpm'].sum().reset_index()
    df_genes_tpm = pd.DataFrame(
         tpm_sum_by_gene['tpm'].values,
         index=tpm_sum_by_gene['Gene_ID'].values,
         columns=[sample_id]
    ).reset_index().rename(columns={'index': 'sample'})
    print(f'[{load_expression_from_abundance.__name__}] ✔️ Created gene-level TPM DataFrame for sample {sample_id}\n')
    
    if load_counts:
        # transcript-level counts
        df_trans_counts = pd.DataFrame(
            df_abundance["est_counts"].values,
            index=df_abundance["Transcript_ID"].values,
            columns=[sample_id]
        ).reset_index().rename(columns={"index": "sample"})

        # gene-level counts
        counts_sum_by_gene = (
            df_abundance.groupby("Gene_ID")["est_counts"]
            .sum()
            .reset_index()
        )

        df_genes_counts = pd.DataFrame(
            counts_sum_by_gene["est_counts"].values,
            index=counts_sum_by_gene["Gene_ID"].values,
            columns=[sample_id]
        ).reset_index().rename(columns={"index": "sample"})

    if load_counts:
        return df_trans_tpm, df_genes_tpm, df_trans_counts, df_genes_counts
    else:
        return df_trans_tpm, df_genes_tpm
