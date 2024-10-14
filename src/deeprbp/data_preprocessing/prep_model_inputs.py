# /src/deeprbp/data_preprocessing/prep_model_inputs.py

import os
import pandas as pd
import numpy as np
from typing import Union, List
from tqdm import tqdm
import warnings
import argparse
import timeit
from verification_output import save_processed_patients_to_excel, compare_patients_with_excel

def prepare_inputs( 
    raw_data_dir: str, 
    selected_genes_dir: str,
    output_dir: str,
    transcript_expression_file: str, 
    gene_expression_file: str,
    phenotype_data_file: str,
    chunk_size: int,
    gene_selection: bool, # Boolean to determine whether to filter cancer genes
    gene_transcript_mapping_file: str,
    splicing_genes_file: str,
    cancer_genes_file: str,
    gene_census_file: str,
    rbp_genes_file: str, 
    **kwargs): # el getBM reducido puede no guardarse en un fichero sino que al cargar las clases el absoluto sea su input y en self.getBM se guarde el subset con los trans que estemos usando
    
    """
    Loads and preprocesses gene and transcript expression data in chunks for memory efficiency. 
    Additionally, filters the gene expression matrix to create a subset containing RNA-binding proteins (RBPs), 
    which serves as the primary input for the model. The final output includes three matrices: 
    one for RBPs, one for all genes, and another for transcripts for each study.

    Args:
    - raw_data_dir (str): Directory containing the raw data files.
    - selected_genes_dir (str): Directory containing the list of RNA-binding proteins (RBPs) and other selected genes for modeling.
    - transcript_expression_file (str): Filename for transcript-level expression data ('trans x n_patients' matrix) in log2(tpm+0.001)).
    - gene_expression_file (str): Filename for gene-level expression data ('genes x n_patients' matrix) in log2(tpm+0.001).
    - phenotype_data_file (str): Filename for phenotype information ('genes x n_patients' matrix).
    - chunk_size (int): The number of rows to process per chunk (for memory efficiency).
    - gene_selection (bool): Specifies which gene set to use. Boolean flag to indicate whether gene selection should be performed.
        If true genes related with cancer and alternative splicing are used.
    - gene_transcript_mapping_file (str): Filename for the output that relates transcript IDs and names to gene IDs, names, and their biotypes for protein coding genes.
    - splicing_genes_file (str): Filename for the Excel file containing genes whose alternative splicing has been characterized to contribute to cancer.
    - cancer_genes_file (str): Filename for the Excel file containing 900 genes predicted to be potential cancer drivers based on mutations and/or copy number alterations.
    - gene_census_file (str): Filename for the tab-separated values (TSV) file containing the Cancer Gene Census data.
    - rbp_genes_file (str): Filename for the Excel file containing a list of RNA-binding proteins (RBPs).
    - output_dir (str): Directory to save the processed files.
    
    Returns:
    None
    """

    # Load patient IDs and phenotype data
    patient_ids = pd.read_csv(f"{raw_data_dir}/{gene_expression_file}", compression='gzip', sep='\t', nrows=1).columns.tolist()
    df_phenotype = pd.read_csv(f"{raw_data_dir}/{phenotype_data_file}", sep='\t', encoding='ISO-8859-1')
    df_phenotype = clean_and_format_phenotype_data(df_phenotype)
 
    # Load the getBM data
    getBM = pd.read_csv(f"{selected_genes_dir}/{gene_transcript_mapping_file}")

    # Track whether phenotype has been saved for each study
    #phenotype_saved = {"TCGA": False, "GTEX": False} # esto creo que no va a hacer falta.

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
            if col_start != 0:
                selected_cols.insert(0, 'sample')

            df_gene = pd.read_csv(f"{raw_data_dir}/{gene_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)
            df_trans = pd.read_csv(f"{raw_data_dir}/{transcript_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)

            process_data_chunk(df_gene, df_trans, df_phenotype, getBM, selected_genes_dir, gene_selection, 
                               rbp_genes_file, output_dir, #phenotype_saved,
                               splicing_genes_file if gene_selection else None,
                               cancer_genes_file if gene_selection else None,
                               gene_census_file if gene_selection else None,    
                               **kwargs)
            
            col_start = col_end
            pbar.update(len(selected_cols) - 1)

            # Store processed patients
            processed_patients["Patient_ID"].extend(selected_cols[1:])  # Excluimos 'sample' (esto lo puedo comentar una vez que el código se haya verificado otra vez)
            processed_patients["Chunk"].extend([chunk_idx] * (len(selected_cols) - 1))  # Excluimos 'sample' (esto lo puedo comentar una vez que el código se haya verificado otra vez)
            chunk_idx += 1
             
    # Save the processed patient information to Excel
    save_processed_patients_to_excel(processed_patients, output_dir) # (esto lo puedo comentar una vez que el código se haya verificado otra vez)
    # Compare processed patients with patient_ids to see if all of them were processed
    compare_patients_with_excel(output_dir, set(patient_ids)) # (esto lo puedo comentar una vez que el código se haya verificado otra vez)
    # Align all matrices according to a common set of patient IDs
    #sort_indices(output_dir)

def clean_and_format_phenotype_data(dictionary_data: pd.DataFrame):
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
    dictionary_data = dictionary_data.set_index("sample")
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace("-", " "), axis=1)
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace(r'[\(\)]', '_', regex=True))
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace(" ", "_"),  axis=1)
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace("_{2,}", "_", regex=True),  axis=1)
    dictionary_data.columns = [col.strip().replace(" ", "_") if not col.startswith("_") else col.strip().replace(" ", "_")[1:] for col in dictionary_data.columns]
    replacements = {
    '^Skin_.+': 'Skin',
    '^Brain_.+': 'Brain',
    '^Adipose_.+': 'Adipose',
    '^Artery_.+': 'Artery',
    '^Colon_.+': 'Colon',
    '^Cervix_.+': 'Cervix',
    '^Esophagus_.+': 'Esophagus'
    }
    mask = dictionary_data['study'] == 'GTEX'
    dictionary_data.loc[mask, 'primary_disease_or_tissue'] = dictionary_data.loc[mask, 'primary_disease_or_tissue'].replace(replacements, regex=True)
    return dictionary_data

def clean_expression_data(data: dict) -> dict:
    """
    Cleans the gene and transcript expression data by removing genome version and aggregating loci.
    Args:
    - data (dict): Dictionary containing 'df_gene' and 'df_trans' DataFrames.
    Returns:
    - dict: Dictionary with cleaned gene and transcript DataFrames.
    """
    print('Cleaning expression data...')
    data['df_gene']['sample'] = data['df_gene']['sample'].str.rsplit('.').str[0].str.replace('R', '0')
    data['df_trans']['sample'] = data['df_trans']['sample'].str.rsplit('.').str[0].str.replace('R', '0')
    data['df_gene'] = data['df_gene'].iloc[:, 1:].groupby(data['df_gene']['sample']).sum()
    data['df_trans'] = data['df_trans'].iloc[:, 1:].groupby(data['df_trans']['sample']).sum()
    print('[clean_expression_data] Cleaning done\n')
    return data

def filter_transcripts(data: dict, getBM: pd.DataFrame) -> dict:
    """
    Filters out transcripts of genes that have only one isoform.
    Args:
    - data (dict): Dictionary containing 'df_trans_log2p_tpm'.
    - getBM (DataFrame): DataFrame containing gene-to-transcript mappings.
    Returns:
    - tuple:
        - data (dict): Updated dictionary with the filtered 'df_trans' DataFrame, excluding transcripts 
                       from genes that have only one isoform.
        - getBM_filtered (DataFrame): DataFrame containing only genes with multiple isoforms and 
                                      their corresponding transcripts.
    """
    print('Filtering transcripts with single isoforms...')
    num_iso_per_gn = getBM['Gene_ID'].value_counts()
    list_gn_one_iso = num_iso_per_gn[num_iso_per_gn == 1].index.tolist()
    set_gn_mult_iso = set(getBM['Gene_ID'].unique()) - set(list_gn_one_iso)
    getBM_filtered = getBM[getBM['Gene_ID'].isin(set_gn_mult_iso)]
    getBM_filtered = getBM_filtered[getBM_filtered.Transcript_ID.isin(data['df_trans'].index)]
    data['df_trans'] = data['df_trans'].loc[getBM_filtered.Transcript_ID]
    print('[filter_transcripts] Filtering done\n')
    return data, getBM_filtered

def find_gene_ids_from_synonyms(dfs: list[pd.DataFrame], synonyms_col: str, getBM: pd.DataFrame, gene_id_col: str) -> list:
    """
    Identifies gene IDs from one or more DataFrames containing a 'Synonyms' column.

    Args:
    - dfs (list[pd.DataFrame]): List of DataFrames, each containing a 'Synonyms' column.
    - synonyms_col (str): Column name in the DataFrames that contains gene synonyms as comma-separated strings.
    - getBM (pd.DataFrame): DataFrame containing 'Gene_name' and 'Gene_ID' columns for gene lookup.
    
    Returns:
    - List of unique Gene_IDs found in the 'Synonyms' columns across all DataFrames.
    """
    gene_ids = set()
    for df in dfs:
        synonyms_series = df[synonyms_col].dropna().str.split(',').explode()
        ensg_synonyms = synonyms_series[synonyms_series.str.startswith('ENSG')].str.split('.').str[0].unique()
        matched_gene_ids = getBM[getBM[gene_id_col].isin(ensg_synonyms)][gene_id_col].unique()
        gene_ids.update(matched_gene_ids)  # Use a set to avoid duplicates
    return list(gene_ids)

def match_dataframe_with_getBM(df_genes_names: Union[pd.DataFrame, List[pd.DataFrame]], 
                               getBM: pd.DataFrame, 
                               parm_match_colums: dict) -> pd.DataFrame:
    """
    Matches the genes from a given DataFrame (or a list of DataFrames) with the corresponding Gene IDs in the getBM DataFrame, 
    including those found by matching synonyms (if applicable).

    Args:
    - df_genes_names (DataFrame or List[DataFrame]): DataFrame or list of DataFrames containing gene names to be matched.
    - getBM (DataFrame): DataFrame containing the mapping of gene names and gene IDs.
    - parm_match_colums (dict): Dictionary containing column name mappings for the various datasets:
        - "df_gene_name_col": Column name in df_genes_names that contains the gene names to be matched (default: 'HGNC symbol').
        - "getBM_gene_name_col": Column name in getBM that contains the gene names to be matched (default: 'Gene_name').
        - "gene_id_col": Column name in getBM that contains the gene IDs (default: 'Gene_ID').
        - "synonyms_col": Column name in df_genes_names that contains synonyms for the genes (default: 'Synonyms').

    Returns:
    - DataFrame: Filtered getBM DataFrame containing only the matched gene IDs.
    """
    if isinstance(df_genes_names, pd.DataFrame):
        df_genes_names = [df_genes_names]
    df_gene_name_col = parm_match_colums['df_gene_name_col']
    getBM_gene_name_col = parm_match_colums['getBM_gene_name_col']
    gene_id_col = parm_match_colums['gene_id_col']
    synonyms_col = parm_match_colums['synonyms_col']  
    gene_names = pd.concat([df[df_gene_name_col] for df in df_genes_names]).unique().tolist()
    matched_gene_ids = getBM[getBM[getBM_gene_name_col].isin(gene_names)][gene_id_col].unique().tolist()
    synonym_ids = find_gene_ids_from_synonyms(dfs=df_genes_names, 
                                                synonyms_col=synonyms_col, 
                                                getBM=getBM,
                                                gene_id_col=gene_id_col)   
    matched_gene_ids.extend(synonym_ids)
    matched_gene_ids = list(set(matched_gene_ids))
    getBM_filtered = getBM[getBM[gene_id_col].isin(matched_gene_ids)].reset_index(drop=True)
    return getBM_filtered

def select_genes_for_modeling(getBM: pd.DataFrame, 
                              selected_genes_dir: str,
                              splicing_genes_file: str,
                              cancer_genes_file: str,
                              gene_census_file: str,
                              parm_match_colums: dict) -> tuple[pd.DataFrame, list, list]:
    """
    Selects genes and their transcripts for modeling, either based on cancer-related genes or all protein-coding genes.
    Args:
    - getBM (DataFrame): DataFrame containing gene-to-transcript mappings.
    - selected_genes_dir (str): Directory path where the raw data is stored.
    - splicing_genes_file (str): Filename for the splicing genes dataset.
    - cancer_genes_file (str): Filename for the cancer genes dataset.
    - gene_census_file (str): Filename for the gene census dataset.
    - parm_match_colums (dict): Dictionary containing column name mappings for the various datasets. 
        Refer to the `match_dataframe_with_getBM` function for details on the expected keys and their meanings.
    Returns:
    - tuple:
        - getBM_filtered (DataFrame): Filtered DataFrame based on the selected genes.
        - list_transcripts (list): List of Transcript_IDs corresponding to the selected genes.
        - list_genes_mapped_to_trans (list): List of selected Gene_IDs mapped to each of Transcripts_IDs
    """
    print('[select_genes_for_modeling] Selecting the genes for modeling ...')
    genes_s5_eyras = pd.read_excel(f'{selected_genes_dir}/{splicing_genes_file}', skiprows=2)
    genes_s6_eyras = pd.read_excel(f'{selected_genes_dir}/{cancer_genes_file}', skiprows=2)
    genes_cosmic = pd.read_csv(f'{selected_genes_dir}/{gene_census_file}', sep="\t")
    selected_genes_dfs = [genes_s5_eyras, genes_s6_eyras, genes_cosmic]
    getBM_filtered = match_dataframe_with_getBM(df_genes_names=selected_genes_dfs, 
                                                getBM=getBM, 
                                                parm_match_colums=parm_match_colums)
    list_genes = getBM_filtered['Gene_ID'].unique().tolist()
    list_transcripts = getBM_filtered['Transcript_ID'].tolist()
    list_genes_mapped_to_trans = [
            getBM_filtered.loc[getBM_filtered['Transcript_ID'] == trans_id, 'Gene_ID'].values[0]
            for trans_id in list_transcripts
            ]
    print('[select_genes_for_modeling] Number of Genes:', len(list_genes))
    print('[select_genes_for_modeling] Number of Transcripts:', len(list_transcripts))
    print('\n')
    return getBM_filtered, list_transcripts, list_genes_mapped_to_trans

def generate_model_input_matrices(data: dict, 
                                  df_phenotype: pd.DataFrame, 
                                  study_name: str, 
                                  list_rbps: list, 
                                  list_transcripts: list, 
                                  list_genes_mapped_to_trans: list) -> tuple[dict, pd.DataFrame]:
    """
    Generates model input matrices from expression data for the specified study.
    Args:
    - data (dict): A dictionary containing gene and transcript expression data. 
                   It should have the following keys:
                   - 'df_gene': DataFrame with gene expression values, where columns represent patients.
                   - 'df_trans': DataFrame with transcript expression values, with patients as columns.
    - df_phenotype (DataFrame): DataFrame containing phenotype data, including patient identifiers and study information.
    - study_name (str): The name of the study ('TCGA' or 'GTEX') used to filter relevant samples.
    - list_rbps (list): List of gene IDs representing RNA-binding proteins (RBPs) to be included in the analysis.
    - list_transcripts (list): List of transcript IDs for which expression data is required.
    - list_genes_mapped_to_trans (list): List of gene IDs that correspond to the transcripts of interest.
    Returns:
    - tuple: A tuple containing:
        - output_data (dict): A dictionary with the following key-value pairs:
            - 'df_rbp_gene' (DataFrame): A matrix of RBP expression values with patients as columns and RBPs as rows.
            - 'df_trans' (DataFrame): A matrix of transcript expression values with patients as columns and transcripts as rows.
            - 'df_genes_mapped_to_trans' (DataFrame): A matrix of gene expression values with patients as columns and genes as rows, where gene indices are set to correspond to the mapped transcripts.
        - df_phenotype_study (DataFrame): A filtered DataFrame containing phenotype information for the patients relevant to the specified study.

    This function filters and organizes expression data into separate matrices for RBP, transcript, and gene expression based on the provided study name. It ensures that only samples common to both the phenotype and expression data are included, enabling accurate modeling for downstream analysis.
    """
    print(f'[generate_model_input_matrices] Generating input matrices for study: {study_name}...')
    df_phenotype_study = df_phenotype[df_phenotype.study == study_name].copy()
    patients_df_gene = data['df_gene'].columns
    patients_df_trans = data['df_trans'].columns  
    # Get samples for the specific study
    samples = list(set(df_phenotype_study.index) & 
                   set(patients_df_gene) & 
                   set(patients_df_trans))
    print(f'[generate_model_input_matrices] Number of samples found: {len(samples)}')
    df_gene = data['df_gene'].loc[:, samples]
    # Obtain the RBP expression matrix
    df_rbp_gene = df_gene.loc[list_rbps, :]
    # Get the dataset with transcript expression
    df_trans = data['df_trans'].loc[list_transcripts, samples]
    # Get the dataset for genes mapped to transcripts
    df_genes_mapped_to_trans = df_gene.loc[list_genes_mapped_to_trans, :]
    # Set index names for genes mapped to transcripts
    df_genes_mapped_to_trans.index = list_transcripts
    output_data = {
        'df_rbp_gene': df_rbp_gene,
        'df_trans': df_trans,
        'df_genes_mapped_to_trans': df_genes_mapped_to_trans
    }
    print(f'[generate_model_input_matrices] RBP matrix shape: {df_rbp_gene.shape}')
    print(f'[generate_model_input_matrices] Transcript matrix shape: {df_trans.shape}')
    print(f'[generate_model_input_matrices] Genes mapped to transcripts matrix shape: {df_genes_mapped_to_trans.shape}')
    print('\n')
    return output_data, df_phenotype_study

def transform_expression_data(data: dict) -> dict:
    """
    Transforms gene expression to TPM and RBP and transcript expression to log2(tpm+1).
    Args:
    - data (dict): Dictionary containing cleaned 'df_gene' and 'df_trans' DataFrames.
    Returns:
    - dict: Dictionary with transformed gene and transcript DataFrames.
    """
    print('Transforming expression data...')
    data['df_rbp_gene'] = np.power(2, data['df_rbp_gene']) - 0.001
    data['df_genes_mapped_to_trans'] = np.power(2, data['df_genes_mapped_to_trans']) - 0.001
    data['df_trans'] = np.power(2, data['df_trans']) - 0.001
    data['df_rbp_gene'] = data['df_rbp_gene'].clip(lower=0)
    data['df_genes_mapped_to_trans'] = data['df_genes_mapped_to_trans'].clip(lower=0)
    data['df_trans'] = data['df_trans'].clip(lower=0)
    data['df_rbp_gene'] = np.log2(data['df_rbp_gene'] + 1)
    data['df_trans'] = np.log2(data['df_trans'] + 1)
    print('[transform_expression_data] Transformation done\n')
    return data

def transpose_dataframes(data: dict) -> dict:
    """
    Transposes the expression DataFrames to have patients as index and genes (or transcript IDs) as columns.

    Args:
    - data (dict): Dictionary containing DataFrames to transpose:
        - 'df_rbp_gene': RBP expression data.
        - 'df_trans': Transcript expression data.
        - 'df_genes_mapped_to_trans': Gene expression data mapped to transcripts.
    Returns:
    - dict: Dictionary with transposed DataFrames for RBP, transcripts, and mapped genes.
    """
    print('Transposing expression data...')
    data['df_rbp_gene'] = data['df_rbp_gene'].T
    data['df_trans'] = data['df_trans'].T
    data['df_genes_mapped_to_trans'] = data['df_genes_mapped_to_trans'].T
    print('[transpose_dataframes] Transposition done\n')
    return data

def save_processed_data(data: dict, df_phenotype_study: pd.DataFrame, output_dir: str, study_name: str) -> None:
    """
    Saves the processed data chunks and phenotype metadata to CSV files in the specified output directory.

    Args:
    - data (dict): Dictionary containing processed data matrices.
    - df_phenotype_study (DataFrame): DataFrame containing phenotype data for the specific study.
    - output_dir (str): Directory where output files will be saved.
    - study_name (str): Name of the study (e.g., 'TCGA', 'GTEX').
    """
    print(f'[save_processed_data] Saving processed data for study: {study_name}...')
    path = os.path.join(output_dir, study_name)
    os.makedirs(path, exist_ok=True)
    path_rbp = os.path.join(path, 'RBPs_log2p_tpm.csv')
    path_trans = os.path.join(path, 'trans_log2p_tpm.csv')
    path_gn = os.path.join(path, 'gn_expr_each_iso_tpm.csv')
    path_phenotype = os.path.join(path, 'phenotype_metadata.csv')
    print(f'[save_processed_data] Saving RBP expression matrix to {path_rbp}...')
    data['df_rbp_gene'].to_csv(path_rbp, mode='a', header=not os.path.exists(path_rbp))
    print(f'[save_processed_data] Saving transcript expression matrix to {path_trans}...')
    data['df_trans'].to_csv(path_trans, mode='a', header=not os.path.exists(path_trans))
    print(f'[save_processed_data] Saving gene expression matrix mapped to transcripts to {path_gn}...')
    data['df_genes_mapped_to_trans'].to_csv(path_gn, mode='a', header=not os.path.exists(path_gn))
    print(f'[save_processed_data] Saving phenotype metadata to {path_phenotype}...')
    df_phenotype_study.loc[data['df_rbp_gene'].index].to_csv(path_phenotype, mode='a', header=not os.path.exists(path_phenotype))
    print(f'[save_processed_data] Saving Phenotype metadata to {path_phenotype}.')
    print('[save_processed_data] Data saving completed.\n')

    # if not phenotype_saved.get(study_name, False):
    #     print(f'[save_processed_data] Saving phenotype metadata to {path_phenotype}...')
    #     df_phenotype_study.to_csv(path_phenotype, mode='w', header=True)
    #     phenotype_saved[study_name] = True

def process_data_chunk(df_gene: pd.DataFrame,  
                        df_trans: pd.DataFrame,
                        df_phenotype: pd.DataFrame,
                        getBM: pd.DataFrame,
                        selected_genes_dir: str,
                        gene_selection: bool,
                        rbp_genes_file: str,
                        output_dir: str,
                        #phenotype_saved: dict,
                        splicing_genes_file: str = None,
                        cancer_genes_file: str = None,
                        gene_census_file: str = None, 
                        **kwargs) -> None:
    """
    Processes data chunks including gene and transcript expression, phenotype data, and gene selection to generate the model inputs

    Args:
    - df_gene (DataFrame): DataFrame containing gene expression data.
    - df_trans (DataFrame): DataFrame containing transcript expression data.
    - df_phenotype (DataFrame): DataFrame containing phenotype data, such as detailed category, 
      primary disease, primary site, sample type, gender, and study information.
    - getBM (DataFrame): DataFrame containing gene and transcript information.
   
    Returns:
    - None
    """
    getBM_copy = getBM.copy()
    
    data = {
        'df_gene': df_gene,
        'df_trans': df_trans
    }

    # Step 0: Define and validate the column names used for gene matching.
    # This dictionary contains the default names of the columns used to match genes between
    # the input DataFrames and the getBM DataFrame.
    parm_match_colums = {
        "df_gene_name_col" : "HGNC symbol",    # Column in df_genes_names containing gene names to match.
        "getBM_gene_name_col" : 'Gene_name',   # Column in getBM containing gene names for matching.
        "gene_id_col" : "Gene_ID",             # Column in getBM containing gene IDs corresponding to the gene names.
        "synonyms_col" : "Synonyms",           # Column in df_genes_names containing synonyms for gene names.
    }

    for key, value in kwargs.items():
        if key in parm_match_colums:
            parm_match_colums[key] = value
        else:
            warnings.warn("Parameter " + key + " is not defined.")

    # For a more detailed description of the parameter names and their usage, 
    # please refer to the match_dataframe_with_getBM function.

    # Step 1: Clean and filter expression data
    data = clean_expression_data(data)
    data, getBM_copy = filter_transcripts(data, getBM_copy)

    # Step 2: Select genes for modeling (if applicable)
    if gene_selection:
        _, list_transcripts, list_genes_mapped_to_trans = select_genes_for_modeling(
            getBM_copy, 
            selected_genes_dir, 
            splicing_genes_file, 
            cancer_genes_file, 
            gene_census_file,
            parm_match_colums
        )

    # Step 3: Filter a list of RNA-binding proteins (RBPs) for modeling
    df_rbps_names = pd.read_excel(f'{selected_genes_dir}/{rbp_genes_file}', skiprows=2)
    getBM_filtered = match_dataframe_with_getBM(df_genes_names=df_rbps_names, 
                                                getBM=getBM, 
                                                parm_match_colums=parm_match_colums)
    list_rbps = getBM_filtered['Gene_ID'].unique().tolist()

    # Step 4: Generate input matrices for the model for each study (TCGA and GTEX) and save
    for study_name in ['TCGA', 'GTEX']:
        print(f"Processing study: {study_name}")
        data_study, df_phenotype_study = generate_model_input_matrices(
            data, df_phenotype, study_name=study_name, list_rbps=list_rbps, 
            list_transcripts=list_transcripts, list_genes_mapped_to_trans=list_genes_mapped_to_trans)

        data_transformed = transform_expression_data(data_study)
        data_transformed = transpose_dataframes(data_transformed)
        
        # Save the processed data and phenotype metadata
        save_processed_data(data_transformed, df_phenotype_study, output_dir, study_name) #, phenotype_saved)
        print('\n')

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess gene and transcript expression data for model input generation.')

    parser.add_argument('--raw_data_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/raw', 
                        help='Directory containing the raw data files.')
    parser.add_argument('--selected_genes_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps', 
                        help='Directory containing the selected genes and RNA-binding proteins (RBPs).')
    parser.add_argument('--output_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/processed', 
                        help='Directory to save the processed data.')
    parser.add_argument('--transcript_expression_file', type=str, default='TcgaTargetGtex_rsem_isoform_tpm.gz', 
                        help='Filename for transcript-level expression data.')
    parser.add_argument('--gene_expression_file', type=str, default='TcgaTargetGtex_rsem_gene_tpm.gz', 
                        help='Filename for gene-level expression data.')
    parser.add_argument('--phenotype_data_file', type=str, default='TcgaTargetGTEX_phenotype.txt', 
                        help='Filename for phenotype information.')
    parser.add_argument('--chunk_size', type=int, default=1000, 
                        help='Number of rows to process per chunk for memory efficiency.')
    parser.add_argument('--gene_selection', type=bool, default=True, 
                        help='Boolean flag to indicate whether gene selection should be performed.')
    parser.add_argument('--gene_transcript_mapping_file', type=str, default='getBM.csv', 
                        help='Filename for gene-transcript mapping data.')
    parser.add_argument('--splicing_genes_file', type=str, default='Table_S5_Cancer_splicing_gene_eyras.xlsx', 
                        help='Filename for alternative splicing genes.')
    parser.add_argument('--cancer_genes_file', type=str, default='Table_S6_Cancer_gene_eyras.xlsx', 
                        help='Filename for potential cancer driver genes.')
    parser.add_argument('--gene_census_file', type=str, default='Table_Cancer_Gene_Census.tsv', 
                        help='Filename for Cancer Gene Census data.')
    parser.add_argument('--rbp_genes_file', type=str, default='Table_S2_list_RBPs_eyras.xlsx', 
                        help='Filename for RNA-binding proteins (RBPs) list.')
    return parser.parse_args()

def main():
    args = parse_args()
    start_time = timeit.default_timer()
    
    prepare_inputs(
                raw_data_dir=args.raw_data_dir, 
                selected_genes_dir=args.selected_genes_dir, 
                output_dir=args.output_dir, 
                transcript_expression_file=args.transcript_expression_file, 
                gene_expression_file=args.gene_expression_file, 
                phenotype_data_file=args.phenotype_data_file,
                chunk_size=args.chunk_size, 
                gene_selection=args.gene_selection, 
                gene_transcript_mapping_file=args.gene_transcript_mapping_file, 
                splicing_genes_file=args.splicing_genes_file, 
                cancer_genes_file=args.cancer_genes_file, 
                gene_census_file=args.gene_census_file, 
                rbp_genes_file=args.rbp_genes_file
            )
    
    end_time = timeit.default_timer()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Total execution time: {int(hours)}h {int(minutes)}m {int(seconds)}s")
    
if __name__ == '__main__':
    main()

# cosas que aun podría meter: 
# - hacer la explicacion del README

# 3. Error Handling
# Add error handling using try-except blocks to handle issues like missing files, incorrect formats, or other data issues gracefully. For example, when reading files:
# python
# Copiar código
# try:
#     df_gene = pd.read_csv(f"{raw_data_dir}/{gene_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)
# except FileNotFoundError as e:
#     print(f"Error: File not found {e.filename}")
#     return

# 9. Unit Testing
# Add unit tests to validate individual components of the preprocessing pipeline. This will help ensure that changes to the code do not break functionality. For example, test the clean_expression_data and filter_transcripts functions separately.

# old: si funciona bien lo puedo quitar ya:
 # n_chunks = (n_patients + chunk_size - 1) // chunk_size

    # for chunk_idx, col_start in tqdm(enumerate(range(0, n_patients, chunk_size + 1)), total=n_chunks, desc='Prepare model inputs'):
    #     print(f'[prepare_inputs] Chunk Number : {chunk_idx}')
    #     col_end = min(col_start + chunk_size, n_patients)
    #     selected_cols = patient_ids[col_start:col_end]
    #     #print(selected_cols)
    #     print('\n')
    #     selected_cols.insert(0, 'sample')  # Include the sample information
    #     #print(selected_cols)
    #     print('\n')
    #     #df_gene = pd.read_csv(f"{raw_data_dir}/{gene_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)
    #     #df_trans = pd.read_csv(f"{raw_data_dir}/{transcript_expression_file}", compression='gzip', sep='\t', usecols=selected_cols)
    #     # process_data_chunk(df_gene, df_trans, df_phenotype, getBM, selected_genes_dir, gene_selection, 
    #     #                     rbp_genes_file, output_dir, phenotype_saved,
    #     #                     splicing_genes_file if gene_selection else None,
    #     #                     cancer_genes_file if gene_selection else None,
    #     #                     gene_census_file if gene_selection else None,    
    #     #                     **kwargs)
    #     # Store processed patients
    #     processed_patients["Patient_ID"].extend(selected_cols[1:])  # Exclude 'sample'
    #     processed_patients["Chunk"].extend([chunk_idx] * (len(selected_cols) - 1))  # Exclude 'sample'
        
    # # Save the processed patient information to Excel
    # save_processed_patients_to_excel(processed_patients, output_dir)
    # # Compare processed patients with patient_ids to see if all of them were processed
    # compare_patients_with_excel(output_dir, set(patient_ids))