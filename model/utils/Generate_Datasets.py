import os
import shutil
import pandas as pd
import numpy as np
from glob import glob
from collections import namedtuple
from sklearn.model_selection import train_test_split
from tqdm import tqdm

########################################## Generate Datasets for TCGA and GTeX data ################################
def check_create_new_directory(path):
    """ 
    Function that checks if a directory is created and if not create the new directory
    Parameters
    ----------
    path: string with the path
    """
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print(f'[Utils] The new directory {path} is created!')

def replace_comma(filename):
    with open(filename, 'r', encoding='ISO-8859-1') as file:
        lines = file.readlines()
    with open(filename, 'w') as file:
        for line in lines:
            line = line.replace(',', '_')
            file.write(line)   

def process_phenotype_data(dictionary_data):
    """
    Processes phenotype data stored in a dictionary format.
    """
    dictionary_data = dictionary_data.set_index("sample")
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace("-", " "), axis=1)
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace(r'[\(\)]', '_', regex=True))
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace(" ", "_"),  axis=1)
    dictionary_data = dictionary_data.apply(lambda x: x.str.replace("_{2,}", "_", regex=True),  axis=1)
    dictionary_data.columns = [col.strip().replace(" ", "_") if not col.startswith("_") else col.strip().replace(" ", "_")[1:] for col in dictionary_data.columns]
    return dictionary_data

def process_data_chunk(df_gene, df_trans, dictionary_data, getBM, path_data, select_genes):  
    """
    Processes data chunks including gene and transcript expression, phenotype data, and gene selection.

    Args:
    - df_gene (DataFrame): DataFrame containing gene expression data.
    - df_trans (DataFrame): DataFrame containing transcript expression data.
    - dictionary_data (DataFrame): DataFrame containing phenotype data.
    - getBM (DataFrame): DataFrame containing gene and transcript information.
    - path_data (str): Path to data directory.
    - select_genes (bool): Indicator whether to select specific genes.

    Returns:
    - None
    """
    print('Processing the chunk...')
    ##### 1) Remove the genoma version from the variable 'sample'
    df_gene['sample'] = df_gene['sample'].str.rsplit('.').str[0]
    df_trans['sample'] = df_trans['sample'].str.rsplit('.').str[0]
    print('[process_data_chunk] Remove the genoma version from the variable ... -> DONE')

    ##### 2) Give to different loci the same gen_id and trans_id 
    df_gene['sample'] = df_gene['sample'].str.replace('R', '0')
    df_trans['sample'] = df_trans['sample'].str.replace('R', '0')
    print('[process_data_chunk] Give to different loci the same gen_id and trans_id ... -> DONE')

    ##### 3) Aggregate different version information and loci.
    df_gene = df_gene.iloc[:, 1:].groupby(df_gene['sample']).sum()      
    df_trans = df_trans.iloc[:, 1:].groupby(df_trans['sample']).sum()
    print('[process_data_chunk] Aggregate different version information and loci ... -> DONE')

    ##### 4)  Convert gene expression to TPMs and transcript expression to log2(tpm+1)
    df_gene_tpm = np.power(2, df_gene) - 0.001
    df_trans_tpm =  np.power(2, df_trans) - 0.001 # new
    # Replace low expression samples with negative values with zeros
    df_gene_tpm = df_gene_tpm.clip(lower=0)
    df_trans_tpm = df_trans_tpm.clip(lower=0)

    df_trans_log2p_tpm = np.log2(df_trans_tpm + 1)
    print('[process_data_chunk] Convert gene expression to TPMs and transcript expression to log2(tpm+1) ... -> DONE')

    ##### 5) Phenotype subset data for TCGA and GTEX sources.
    df_tcga_phtype = dictionary_data[dictionary_data.study == 'TCGA'].copy()
    #df_tcga_phtype = df_tcga_phtype[df_tcga_phtype.sample_type == 'Primary_Tumor'] # ahora mismo entrenamos con todo tcga, bien tumoral, bien normal
    df_gtex_phtype = dictionary_data[dictionary_data.study == 'GTEX'].copy()

    ##### 5.1) Combine some category types in GTEX
    # Define the pattern-replacement dictionary
    replacements = {'^Skin_.+': 'Skin', '^Brain_.+': 'Brain', '^Adipose_.+': 'Adipose', '^Artery_.+': 'Artery', '^Colon_.+': 'Colon', 
                    '^Cervix_.+': 'Cervix', '^Esophagus_.+': 'Esophagus'}
    df_gtex_phtype['primary_disease_or_tissue'] = df_gtex_phtype['primary_disease_or_tissue'].replace(replacements, regex=True)
    print('[process_data_chunk] Combine some category types in GTEX ... -> DONE')

    ##### 6) List of RBPs:
    df_RBPs = pd.read_excel(f'{path_data}/selected_genes_rbps/Table_S2_list_RBPs_eyras.xlsx', skiprows=2) # List of splicing genes
    ids_RBPs = getBM[getBM['Gene_name'].isin(df_RBPs['HGNC symbol'].tolist())]['Gene_ID'].unique().tolist()
    print('[process_data_chunk] List of RBPs created ... -> DONE')
    
    ##### 7) Drop the genes with just one transcript 
    getBM = pd.read_csv(f"{path_data}/extra/getBM_total.csv")
    
    # Obtain the number of isoforms for each gene and create a list of genes with only one isoform.
    num_iso_per_gn = getBM['Gene_ID'].value_counts()
    list_gn_uniqueiso = num_iso_per_gn[num_iso_per_gn == 1].index.tolist()
    set_gn_more_one_iso = set(getBM['Gene_ID'].unique()) - set(list_gn_uniqueiso)

    # Filter the getBM DataFrame to include only genes with > than 1 isoform
    getBM = getBM.loc[getBM['Gene_ID'].isin(set_gn_more_one_iso)]

    # 7.1) Filter the Transcript data matching the getBM information:
    getBM = getBM[getBM.Transcript_ID.isin(df_trans_log2p_tpm.index)]
    df_trans_log2p_tpm = df_trans_log2p_tpm.loc[list(getBM.Transcript_ID),:]
    print('[process_data_chunk] Drop the unique transcripts ... -> DONE')

    ##### 8) Select Genes:
    print('[process_data_chunk] Selecting the genes for modelling ... -> DONE')

    if select_genes is not None: # select_genes == cancer_genes
        print('[process_data_chunk] We are using just cancer related genes!')
        genes_s5_eyras = pd.read_excel(f'{path_data}/selected_genes_rbps/Table_S5_Cancer_splicing_gene_eyras.xlsx', skiprows=2)
        genes_s6_eyras = pd.read_excel(f'{path_data}/selected_genes_rbps/Table_S6_Cancer_gene_eyras.xlsx', skiprows=2)
        genes_cosmic = pd.read_csv(f'{path_data}/selected_genes_rbps/Table_Cancer_Gene_Census.tsv', sep="\t")

        selected_genes = pd.concat([genes_s5_eyras['HGNC symbol'], genes_s6_eyras['HGNC symbol'], genes_cosmic['Gene Symbol']])
        selected_genes = selected_genes.unique()
        ids_selected_genes = getBM[getBM['Gene_name'].isin(selected_genes)]['Gene_ID'].unique().tolist()
        
        getBM = getBM.iloc[[a in ids_selected_genes for a in getBM.Gene_ID], :] 
        getBM = getBM.reset_index(drop=True)

        ##### 9) Save reduced getBM set
        getBM.to_csv(f'{path_data}/extra/getBM_reduced.csv')  
        
    else: # else selected_genes = '', we are using all the protein coding genes
        print('We are using all protein-coding genes!')
        
        ##### 9) Save reduced getBM set
        getBM.to_csv(f'{path_data}/extra/getBM_reduced.csv') 
        
    list_transcripts = list(getBM.Transcript_ID)
    list_genes = list(getBM.Gene_ID)
    print('[process_data_chunk] Number of Genes:', len(getBM.Gene_ID.unique()))
    print('[process_data_chunk] Number of Transcripts:', len(list_transcripts))
    print('[process_data_chunk] Save reduced getBM set ... -> DONE')
    
    ##### 10) Filter categories from_phenotype and process datasets
    print('[process_data_chunk] Filter categories from_phenotype and process datasets...')
    divide_from_category_phenotype(df_tpm = df_gene_tpm, df_trans_log2p_tpm = df_trans_log2p_tpm, d_phtype = df_tcga_phtype,
                                    patients_df_tpm = df_gene_tpm.columns, patients_df_trans_log2p_tpm = df_trans_log2p_tpm.columns,
                                    ids_RBPs = ids_RBPs, selected_genes = list_genes, selected_transcripts = list_transcripts,
                                    path_data=path_data) # TCGA (tumor & normal together)
    divide_from_category_phenotype(df_tpm = df_gene_tpm, df_trans_log2p_tpm = df_trans_log2p_tpm, d_phtype = df_gtex_phtype,
                                    patients_df_tpm = df_gene_tpm.columns, patients_df_trans_log2p_tpm = df_trans_log2p_tpm.columns,
                                    ids_RBPs = ids_RBPs, selected_genes = list_genes, selected_transcripts = list_transcripts,
                                    path_data=path_data)  # GTEX (normal samples)
    print('[process_data_chunk] Filter categories from_phenotype and process datasets ... -> DONE')

def divide_from_category_phenotype(df_tpm, df_trans_log2p_tpm, d_phtype, patients_df_tpm, patients_df_trans_log2p_tpm, ids_RBPs, 
                                   selected_genes, selected_transcripts, path_data):
    """
    Divides data into categories based on phenotype information and saves the datasets.

    Args:
    - df_tpm (DataFrame): DataFrame containing gene expression data.
    - df_trans_log2p_tpm (DataFrame): DataFrame containing transcript expression data.
    - d_phtype (DataFrame): DataFrame containing phenotype data.
    - patients_df_tpm (list): List of patients for gene expression data.
    - patients_df_trans_log2p_tpm (list): List of patients for transcript expression data.
    - ids_RBPs (list): List of RNA-binding proteins IDs.
    - selected_genes (list): List of selected genes.
    - selected_transcripts (list): List of selected transcripts.
    - path_data (str): Path to data directory.

    Returns:
    - None
    """
    study_name = d_phtype['study'].unique().tolist()[0]
    categories = d_phtype['primary_disease_or_tissue'].dropna().unique().tolist()
    for i in tqdm(categories, desc='Processing categories'):
        print('[divide_from_category_phenotype] category_name:', i)
        samples = list(set(d_phtype[d_phtype.detailed_category == i].index) & set(patients_df_tpm) & set(patients_df_trans_log2p_tpm))
        if not samples:
            continue
        else:
            #### Obtain the sets:
            ### 1) df expression of all genes for type i
            df_tpm_gn_i = df_tpm.loc[:, samples]
            print(f'df_tpm_gn_i: {df_tpm_gn_i}')

            ### 2) Obtaining the datafame of the RBP expression from the general gene with
            # a log2+1 transformation of the gene expression
            df_sf_tpm_i = df_tpm_gn_i.loc[ids_RBPs,:]
            df_sf_log2p_tpm_i = np.log2(1+df_sf_tpm_i)
            print(f'df_sf_log2p_tpm_i: {df_sf_log2p_tpm_i}')

            ### 3) Get the dataset with the expression of the transcripts
            df_trans_log2p_tpm_i = df_trans_log2p_tpm.loc[selected_transcripts, samples]
            print(f'df_trans_log2p_tpm_i: {df_trans_log2p_tpm_i}')

            ### 4) Get the dataset with the gene expression for each transcript
            df_gn_each_iso_tpm_i = df_tpm_gn_i.loc[selected_genes,:]
            print(f'df_gn_each_iso_tpm_i: {df_gn_each_iso_tpm_i}')

            ### 5) Get the transpose of the sets so the information of the patients is in the index and the expression in 
            #the columns
            df_sf_log2p_tpm_i = df_sf_log2p_tpm_i.T
            df_trans_log2p_tpm_i = df_trans_log2p_tpm_i.T
            df_gn_each_iso_tpm_i = df_gn_each_iso_tpm_i.T

            ### 6) We put to the gene_expr_each_iso dataframe the colnames of the transcripts they are related to
            df_gn_each_iso_tpm_i.columns = selected_transcripts

            ### 7) Change the colnames of rbp set by gene_name
            getBM = pd.read_csv(f"{path_data}/extra/getBM_total.csv")
            gene_id_name_map = {row['Gene_ID']: row['Gene_name'] for index, row in getBM[getBM.Gene_ID.isin(ids_RBPs)].drop_duplicates(subset='Gene_name').iterrows()}
            df_sf_log2p_tpm_i = df_sf_log2p_tpm_i.loc[:,gene_id_name_map.keys()] #the are 4 rbp_ids that have the same rbp_name
            df_sf_log2p_tpm_i = df_sf_log2p_tpm_i.rename(columns=gene_id_name_map)
            print(f'[divide_from_category_phenotype] Obtain the RBPs, Transcript and Gene set for {i} ... -> DONE')

            ### 8) Save the whole sets in files:
            path = f'{path_data}/processed/splitted_datasets/{study_name}/{i}'
            check_create_new_directory(path)
            path_sf = f'{path}/{i}_{study_name}_RBPs_log2p_tpm.csv'
            path_trans =  f'{path}/{i}_{study_name}_trans_log2p_tpm.csv'
            path_gn = f'{path}/{i}_{study_name}_gn_expr_each_iso_tpm.csv'
            
            df_sf_log2p_tpm_i.to_csv(path_sf, mode='a', header=not os.path.exists(path_sf)) #sfs
            df_trans_log2p_tpm_i.to_csv(path_trans, mode='a', header=not os.path.exists(path_trans)) #trans
            df_gn_each_iso_tpm_i.to_csv(path_gn, mode='a', header=not os.path.exists(path_gn)) #genes for each iso
            print(f'[divide_from_category_phenotype] Saving the {study_name} datasets for {i} in {path}!')

def split_training_test_sets(df_sfs, df_trans, df_gns_each_trans, test_size=0.2):
    """
    Splits the datasets into training and testing sets.

    Args:
    - df_sfs (DataFrame): DataFrame containing splicing factor expression data.
    - df_trans (DataFrame): DataFrame containing transcript expression data.
    - df_gns_each_trans (DataFrame): DataFrame containing gene expression for each transcript.
    - test_size (float): Proportion of the dataset to include in the test split.

    Returns:
    - namedtuple: Namedtuple containing the training and testing sets.
    """
    df_train, df_test = train_test_split(df_sfs, test_size=test_size, random_state=0)  
    # labels (we need the same patients so we use the same index selection)
    train_labels = df_trans.loc[df_train.index]
    test_labels = df_trans.loc[df_test.index]
    # gene_expr for each isoform
    train_gn = df_gns_each_trans.loc[df_train.index]
    test_gn = df_gns_each_trans.loc[df_test.index]
    DataSplit = namedtuple('DataSplit', ['df_train', 'train_labels', 'train_gn',
                                         'df_test', 'test_labels', 'test_gn'])
    data_split = DataSplit(df_train, train_labels, train_gn, df_test, test_labels, test_gn)          
    return data_split   

def from_directory_train_test_split(path_data, categories_source, study_name):
    """
    Splits the datasets into training and testing sets based on directories.

    Args:
    - path_data (str): Path to the data directory.
    - categories_source (list): List of categories to process.
    - study_name (str): Name of the study.

    Returns:
    None
    """

    for i in categories_source: 
        df_sf_log2p_tpm_i = pd.read_csv(f'{path_data}/processed/splitted_datasets/{study_name}/{i}/{i}_{study_name}_RBPs_log2p_tpm.csv', index_col=0)
        df_trans_log2p_tpm_i = pd.read_csv(f'{path_data}/processed/splitted_datasets/{study_name}/{i}/{i}_{study_name}_trans_log2p_tpm.csv', index_col=0)
        df_gn_each_iso_tpm_i = pd.read_csv(f'{path_data}/processed/splitted_datasets/{study_name}/{i}/{i}_{study_name}_gn_expr_each_iso_tpm.csv', index_col=0)

        if study_name == 'TCGA': # new
            ### Divide in training and test
            splitted_datasets = split_training_test_sets(df_sf_log2p_tpm_i, df_trans_log2p_tpm_i, df_gn_each_iso_tpm_i)
            df_train = splitted_datasets.df_train
            df_test = splitted_datasets.df_test
            train_labels = splitted_datasets.train_labels
            test_labels = splitted_datasets.test_labels
            train_gn = splitted_datasets.train_gn
            test_gn = splitted_datasets.test_gn 

            ### Create Training paths and folder
            path_train = f'{path_data}/processed/splitted_datasets/{study_name}/training/{i}'
            check_create_new_directory(path_train)

            path_train_sf = f'{path_train}/{i}_training_{study_name}_RBPs_log2p_tpm.csv'
            path_train_trans = f'{path_train}/{i}_training_{study_name}_trans_log2p_tpm.csv'
            path_train_gn = f'{path_train}/{i}_training_{study_name}_gn_expr_each_iso_tpm.csv'

            # Save Training samples
            df_train.to_csv(path_train_sf, mode='a', header=not os.path.exists(path_train_sf)) #sfs
            train_labels.to_csv(path_train_trans, mode='a', header=not os.path.exists(path_train_trans)) #trans
            train_gn.to_csv(path_train_gn, mode='a', header=not os.path.exists(path_train_gn)) #genes for each iso

        else: 
            ### if GTEX ALL the samples will be considered as 'test'
            df_test = df_sf_log2p_tpm_i
            test_labels = df_trans_log2p_tpm_i
            test_gn = df_gn_each_iso_tpm_i

        ### Create Test paths and folder
        path_test = f'{path_data}/processed/splitted_datasets/{study_name}/test/{i}'
        check_create_new_directory(path_test)

        path_test_sf = f'{path_test}/{i}_test_{study_name}_RBPs_log2p_tpm.csv'
        path_test_trans = f'{path_test}/{i}_test_{study_name}_trans_log2p_tpm.csv'
        path_test_gn = f'{path_test}/{i}_test_{study_name}_gn_expr_each_iso_tpm.csv'

        # Test
        df_test.to_csv(path_test_sf, mode='a', header=not os.path.exists(path_test_sf)) #sfs
        test_labels.to_csv(path_test_trans, mode='a', header=not os.path.exists(path_test_trans)) #trans
        test_gn.to_csv(path_test_gn, mode='a', header=not os.path.exists(path_test_gn)) #genes for each iso
        
        # Drop the original folder with the entire set
        shutil.rmtree(f'{path_data}/processed/splitted_datasets/{study_name}/{i}')

def create_datasets(path_data, chunksize, select_genes):
    """
    Reads and processes gene data by chunks and divides the generated data into training and test sets.

    Args:
    - path_data (str): Path to the data directory.
    - chunksize (int): Size of the data chunks to process.
    - select_genes (bool): Whether to select specific genes for modeling.

    Returns:
    None
    """

#FileNotFoundError: [Errno 2] No such file or directory: #'/Users/joseba/Downloads/ML4BMLab2/DeepRBP/data/input_create_model/../raw/TcgaTargetGtex_rsem_gene_tpm.gz'

#FileNotFoundError: [Errno 2] No such file or directory: './raw/TcgaTargetGtex_rsem_gene_tpm.gz'
#data/input_create_model/raw/TcgaTargetGTEX_phenotype.txt
#data/input_create_model/processed/create_data.sh
#data/input_create_model/raw
    
    columns_patients = pd.read_csv(f"{path_data}/raw/TcgaTargetGtex_rsem_gene_tpm.gz", compression='gzip', sep='\t', nrows=1).columns.tolist()
    n_patients = len(columns_patients)-1
    #### 1) Read the general Data 
    dictionary_data = pd.read_csv(f"{path_data}/raw/TcgaTargetGTEX_phenotype.txt", sep='\t', encoding='ISO-8859-1')
    dictionary_data = process_phenotype_data(dictionary_data)
    print('[create_datasets] Read the dictionary_data ... -> DONE')
    # 1.1) Load the total getBM (information that relates Gene-Transcript information for genes that are protein-coding - the transcripts 
    # of these genes can have also other functions)
    getBM = pd.read_csv(f"{path_data}/extra/getBM_total.csv")
    print('[create_datasets] Load the getBM data ... -> DONE')

    #### 2) Read the GENE DATA by chunks and process it.
    print('[create_datasets] Read the GENE DATA by chunks and process it')
    # Calculate the number of chunks
    n_chunks = (n_patients + chunksize - 1) // chunksize

    for chunk_idx, col_start in tqdm(enumerate(range(1, n_patients, chunksize+1)), total=n_chunks, desc='Processing chunks'): # for chunk:
        if col_start + chunksize > n_patients:
            col_end = col_start+n_patients-col_start
            print('[create_datasets] col_start:', col_start)
            print('[create_datasets] col_end:', col_end)
        else:
            col_end = col_start + chunksize
            print('[create_datasets] col_start:', col_start)
            print('[create_datasets] col_end:', col_end)
        
        selected_cols = columns_patients[col_start:col_end]
        print('[create_datasets] len(selected_cols):', len(selected_cols))
        selected_cols.insert(0, 'sample') # To selected cols insert the sample information

        # Read the selected cols:
        df_gene = pd.read_csv(f"{path_data}/raw/TcgaTargetGtex_rsem_gene_tpm.gz", compression='gzip', sep='\t', usecols=selected_cols)
        df_trans = pd.read_csv(f"{path_data}/raw/TcgaTargetGtex_rsem_isoform_tpm.gz", compression='gzip', sep='\t', usecols=selected_cols)
        print('[create_datasets] df_gene & df_trans chunks read ... -> DONE')
        process_data_chunk(df_gene, df_trans, dictionary_data, getBM, path_data, select_genes) #create_datasets
        
    print('[create_datasets] Divide the generated data in the folders in training and test...')
    #### 3) Divide the generated data in the folders in training and test:
    categories_tcga = [os.path.basename(folder) for folder in glob(f'{path_data}/processed/splitted_datasets/TCGA/*')]
    categories_gtex = [os.path.basename(folder) for folder in glob(f'{path_data}/processed/splitted_datasets/GTEX/*')]

    print(categories_tcga)
    from_directory_train_test_split(path_data, categories_source = categories_tcga, study_name='TCGA')
    from_directory_train_test_split(path_data, categories_source = categories_gtex, study_name='GTEX')
    print('[create_datasets] Divide the generated data in the folders in training and test ... -> DONE')

########################################## Generate Datasets for in-silico validation of DeepRBP using real kds ################################

def create_vec_exp_matrix_from_abundance(path_data, id_sample):
    """
    Reads abundance.tsv file, expands detailed information from target_id column, removes version names from Transcript_ID and Gene_ID,
    saves the final abundance file, and creates final vectors of transcript expression and gene expression for the sample.

    Args:
    - path_data (str): Path to the data directory.
    - id_sample (str): ID of the sample.

    Returns:
    - df_trans_tpm (DataFrame): DataFrame containing transcript TPM values.
    - df_trans_counts (DataFrame): DataFrame containing transcript counts.
    - df_genes_tpm (DataFrame): DataFrame containing gene TPM values.
    """
    ### Read the abundance.tsv file
    df = pd.read_csv(f'{path_data}/{id_sample}/abundance.tsv', sep='\t')
    print(f'[create_vec_exp_matrix_from_abundance] Read the abundance.tsv file of sample: {id_sample} ... -> DONE')

    ### Expand detailed information from target_id column and reorder columns
    df[['Transcript_ID', 'Gene_ID', 'OTTHUMG', 'OTTHUMT', 'Gene_Symbol', 'Gene_Name', 'Length', 'Biotype']] = df['target_id'].str.rstrip('|').str.split('|', expand=True)
    df.drop(columns=['target_id', 'Length'], inplace=True)
    df = df[['Transcript_ID', 'Gene_ID', 'OTTHUMT', 'OTTHUMG', 'Gene_Symbol', 'Gene_Name', 'Biotype', 'length', 'eff_length', 'est_counts', 'tpm']]
    print('[create_vec_exp_matrix_from_abundance] Expand detailed information from target_id column and reorder columns ... -> DONE')
    
    ### Remove version name from Transcript_ID and Gene_ID
    df['Transcript_ID'] = df['Transcript_ID'].str.split('.').str[0]
    df['Gene_ID'] = df['Gene_ID'].str.split('.').str[0]
    print(f'Are there unique version of the Transcripts? : {(len(df.Transcript_ID.unique()) == df.shape[0])}')
    print('[create_vec_exp_matrix_from_abundance] Remove version name from Transcript_ID and Gene_ID ... -> DONE')
    
    ### Save the final abundance file (with counts)
    df.to_csv(f'{path_data}/{id_sample}/processed_abundance.csv')

    ### Create the final vectors of transcript expression and gene expression for this sample
    # trans (tpm & counts)
    df_trans_tpm = pd.DataFrame(df['tpm'].values, index=df['Transcript_ID'].values, columns=[f'{id_sample}']).reset_index()
    df_trans_counts = pd.DataFrame(df['est_counts'].values, index=df['Transcript_ID'].values, columns=[f'{id_sample}']).reset_index()
    
    # genes
    sum_tpm_df = df.groupby('Gene_ID')['tpm'].sum().reset_index() # Create the grouping by 'Gene_ID' and calculate the sum of the 'tpm' column
    df_genes_tpm = pd.DataFrame(sum_tpm_df['tpm'].tolist(), index=sum_tpm_df['Gene_ID'].tolist(), columns=[f'{id_sample}']).reset_index()
    
    df_trans_tpm.rename(columns={'index': 'sample'}, inplace=True)
    df_trans_counts.rename(columns={'index': 'sample'}, inplace=True)
    df_genes_tpm.rename(columns={'index': 'sample'}, inplace=True)
    print('[create_vec_exp_matrix_from_abundance] Create the final vectors of transcript expression and gene expression for this sample ... -> DONE')
    print('\n')
    return df_trans_tpm, df_trans_counts, df_genes_tpm

def process_and_save_expression_data(path_data, id_samples, condition, path_model): # This code now maybe doesnt work
    """
    Processes expression data for samples, merges them, performs data cleaning and aggregation,
    log2 transformation of transcript data, obtains RBP expression dataframe,
    selects genes depending on the model used, transposes DataFrames, and saves the created files.

    Args:
    - path_data (str): Path to the data directory.
    - id_samples (list): List of sample IDs.
    - condition (str): Condition of the samples.
    - path_model (str): Path to the model directory.
    """
    # Initialize cumulative DataFrames
    df_trans_tpm_cumulative = pd.DataFrame()
    df_trans_counts_cumulative = pd.DataFrame()
    df_genes_tpm_cumulative = pd.DataFrame()

    # 1. Merge each sample
    print(f'[process_samples_and_merge] 1. Merging the samples ... ')
    for sample in id_samples:
        print(f'[process_samples_and_merge] We are processing the sample: {sample} from tissue: {condition}')
        # Obtain data for the current sample
        df_trans_tpm, df_trans_counts, df_genes_tpm = create_vec_exp_matrix_from_abundance(path_data, sample)
        # Merge in the cumulative DataFrames on the 'sample' column
        if df_trans_tpm_cumulative.empty:
            # If the cumulative DataFrame is empty, simply copy the data
            df_trans_tpm_cumulative = df_trans_tpm.copy()
            df_trans_counts_cumulative = df_trans_counts.copy()
            df_genes_tpm_cumulative = df_genes_tpm.copy()
        else:
            # Merge in the cumulative DataFrames on the 'sample' column
            df_trans_tpm_cumulative = pd.merge(df_trans_tpm_cumulative, df_trans_tpm, on='sample', how='outer')
            df_trans_counts_cumulative = pd.merge(df_trans_counts_cumulative, df_trans_counts, on='sample', how='outer')
            df_genes_tpm_cumulative = pd.merge(df_genes_tpm_cumulative, df_genes_tpm, on='sample', how='outer')
    print(f'[process_samples_and_merge] 1. Merging the samples ... -> DONE')
    
    # 2. Data Cleaning and aggregation of different loci
    # Give to different loci the same gen_id and trans_id 
    df_trans_tpm_cumulative['sample'] = df_trans_tpm_cumulative['sample'].str.replace('R', '0')
    df_trans_counts_cumulative['sample'] = df_trans_counts_cumulative['sample'].str.replace('R', '0')
    df_genes_tpm_cumulative['sample'] = df_genes_tpm_cumulative['sample'].str.replace('R', '0')
    df_trans_tpm_cumulative = df_trans_tpm_cumulative.iloc[:, 1:].groupby(df_trans_tpm_cumulative['sample']).sum()
    df_trans_counts_cumulative = df_trans_counts_cumulative.iloc[:, 1:].groupby(df_trans_counts_cumulative['sample']).sum()
    df_genes_tpm_cumulative = df_genes_tpm_cumulative.iloc[:, 1:].groupby(df_genes_tpm_cumulative['sample']).sum()      
    print('[process_samples_and_merge] 2. Data Cleaning and aggregation of different loci ... -> DONE')
    
    # 3. Log2 transformation of Transcript data
    df_trans_log2p_tpm_cumulative = np.log2(df_trans_tpm_cumulative + 1)
    print('[process_samples_and_merge] 3. Log2 transformation of Transcript data ... -> DONE')
    
    # 4. Obtaining RBP expression dataframe in log2(tpm+1)
    getBM = pd.read_csv(f"{path_data}/extra/getBM_total.csv")

    df_RBPs = pd.read_excel(f'{path_model}/selected_genes_rbps/Table_S2_list_RBPs_eyras.xlsx', skiprows=2) # List of splicing genes
    ids_RBPs = getBM[getBM['Gene_name'].isin(df_RBPs['HGNC symbol'].tolist())]['Gene_ID'].unique().tolist()
    df_rbp_tpm_cumulative = df_genes_tpm_cumulative.loc[ids_RBPs,:]
    df_rbp_log2p_tpm_cumulative = np.log2(1+df_rbp_tpm_cumulative)
    gene_id_name_map = {row['Gene_ID']: row['Gene_name'] for index, row in getBM[getBM.Gene_ID.isin(ids_RBPs)].drop_duplicates(subset='Gene_name').iterrows()} # Change the colnames of rbp set by gene_name
    df_rbp_log2p_tpm_cumulative = df_rbp_log2p_tpm_cumulative.loc[gene_id_name_map.keys(), :] #the are 4 rbp_ids that have the same rbp_name
    df_rbp_log2p_tpm_cumulative = df_rbp_log2p_tpm_cumulative.rename(index=gene_id_name_map)
    print('[process_samples_and_merge] 4. Obtaining RBP expression dataframe in log2(tpm+1) ... -> DONE')
    
    # 5. Gene selection depending on model used (cancer genes or protein coding genes) and get the dataset with the gene expression for each transcript
    getBM_mini = pd.read_csv(f'{path_model}/extra/getBM_reduced.csv', index_col=0)
    list_transcripts = list(getBM_mini.Transcript_ID)
    list_genes = list(getBM_mini.Gene_ID)
    df_trans_log2p_tpm_cumulative = df_trans_log2p_tpm_cumulative.loc[list_transcripts, :] # Get the dataset with the expression of the transcripts
    df_gn_each_iso_tpm_cumulative = df_genes_tpm_cumulative.loc[list_genes,:] # Get the dataset with the gene expression for each transcript
    print('[process_samples_and_merge] 5. Gene selection depending on model used and get the dataset with the gene expression for each transcript ... -> DONE')
    
    # 6. Transpose DataFrames
    df_rbp_log2p_tpm_cumulative = df_rbp_log2p_tpm_cumulative.T
    df_trans_log2p_tpm_cumulative = df_trans_log2p_tpm_cumulative.T
    df_gn_each_iso_tpm_cumulative = df_gn_each_iso_tpm_cumulative.T
    df_gn_each_iso_tpm_cumulative.columns = list_transcripts # We put to the gene_expr_each_iso dataframe the colnames of the transcripts they are related to
    print('[process_samples_and_merge] 6. Transpose DataFrames ... -> DONE')
    
    # 7. Save created files in folder
    path = f'{path_data}/datasets/{condition}'
    check_create_new_directory(path)
    path_rbp = f'{path}/{condition}_RBPs_log2p_tpm.csv'
    path_trans =  f'{path}/{condition}_trans_log2p_tpm.csv'
    path_trans_counts = f'{path}/{condition}_trans_est_counts.csv'
    path_gn = f'{path}/{condition}_gn_expr_each_iso_tpm.csv'
    df_rbp_log2p_tpm_cumulative.to_csv(path_rbp) #rbps
    df_trans_log2p_tpm_cumulative.to_csv(path_trans) #trans
    df_trans_counts_cumulative.to_csv(path_trans_counts)
    df_gn_each_iso_tpm_cumulative.to_csv(path_gn) #genes for each iso
    print('[process_samples_and_merge] 7. Save created files in folder ... -> DONE \n')
    print(f'[process_samples_and_merge] {condition} files were saved in {path}')

