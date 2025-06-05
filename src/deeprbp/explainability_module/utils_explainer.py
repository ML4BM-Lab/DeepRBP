# src/deeprbp/explainability_module/utils_explainer.py

import pandas as pd
from collections import namedtuple
from typing import List, Optional

from .deeplift_handler import DeepLiftHandler
from .pseudoknock_handler import PseudoKnockHandler

def initialize_explainer_handler(config, dataset, model):
    """Initialize the ExplainerHandler with the loaded model and data.
    
    This function creates an instance of the appropriate explainer handler 
    based on the configuration provided. Supported explanation methods include 
    'DeepLIFT' and 'Pseudoknockdown'. If an unknown method is specified, 
    a ValueError is raised.

    Args:
        config (ConfigParser): An instance of the ConfigParser class containing settings for the explainer.
        dataset (DeepRBPExpressionDataset): An instance of the DeepRBPExpressionDataset class that contains the necessary data 
                        for the explanation process.
        model (PredictorModel): An instance of the PredictorModel class, which is a trained PyTorch 
                                Lightning model that will be used for calculating 
        explainability scores.

    Returns:
        ExplainerHandler: An instance of the appropriate explainer handler 
                        (DeepLiftHandler or PseudoKnockHandler) initialized 
                        with the provided dataset and model.

    Raises:
        ValueError: If the provided explanation method in the configuration is unknown 
                    or unsupported. 
    """
    print("🔍 Initializing Explainer handler...")
    explain_method = config.get('explanation_method') 
    if explain_method == "DeepLIFT":
        explainer_handler = DeepLiftHandler(config, dataset, model)
    elif explain_method == "Pseudoknockdown":   
        explainer_handler = PseudoKnockHandler(config, dataset, model)
    else:
        print(f"❌ Unknown explanation method: {explain_method}", level=1)
        raise ValueError(f"Unknown explanation method: {explain_method}")
    print(" ✅ Explainer handler initialized successfully.")
    return explainer_handler  

def filter_scores_for_low_expressed_transcripts(scores, dataset):
    """
    Filters RBP x transcript explainability scores for transcripts that are never expressed.

    This function identifies transcripts that have zero total expression 
    and sets their corresponding scores to zero. It is useful for 
    ensuring that only relevant transcripts are included in the analysis, 
    as transcripts that are never expressed do not contribute meaningful 
    information to the model's explanations.

    Args:
        scores (pd.DataFrame): A DataFrame containing explainability scores at transcript level.
        dataset (DeepRBPExpressionDataset): An instance of the DeepRBPExpressionDataset class,
                                             which contains expression data, specifically 
                                             the 'isoform_df' DataFrame that holds the 
                                             expression levels of each transcript.

    Returns:
        pd.DataFrame: An updated DataFrame of explainability scores with scores for 
                       low-expressed transcripts (those with total expression 
                       equal to zero) set to zero.
    """
    print("🔍 Filtering scores for low-expressed transcripts...")
    # Convert log2-transcripts per million (log2p(tpm)) expression data to transcripts per million (tpm)
    trans_expr = 2 ** dataset.features['isoform_df'] - 1
    trans_expr_df = pd.DataFrame(trans_expr.numpy(), columns=dataset.trans_names)
    # Identify transcripts that have a total expression of 0 (never expressed)
    transcripts_never_expressed = (trans_expr_df.sum(axis=0) == 0)  # Transcripts that have a total expression of 0
    # Set scores to 0 for the transcripts that are never expressed
    scores.loc[transcripts_never_expressed, :] = 0
    # Count the number of transcripts that are never expressed
    num_never_expressed = transcripts_never_expressed.sum()
    print(f"✅ Filter results: {num_never_expressed} transcripts never expressed (total expression = 0).")
    return scores

def filter_scores_for_low_expressed_genes(scores, dataset, threshold=1):
    """
    Set low-expressed genes (mean expression < threshold) to 0 in the TxRBP scores DataFrame.

    Args:
        scores (pd.DataFrame): A DataFrame containing explainability scores at transcript level.
        dataset (DeepRBPExpressionDataset): An instance of the DeepRBPExpressionDataset class,
                                             which contains expression data, specifically 
                                             the 'gene_df' DataFrame that holds the 
                                             expression levels of each gene.
        threshold (float): Expression threshold (default=1 TPM).

    Returns:
        pd.DataFrame: Updated DeepLIFT scores with low-expressed genes set to 0.
    """
    print("🔍 Filtering scores for low-expressed genes...")
    gene_expr_df = pd.DataFrame(dataset.features['gene_df'].numpy(), columns=dataset.trans_names)  # DataFrame with gene expression values in TPM. 
    # colnames are the transcript ids are in this case we are working with the expanded matrix
    low_expr_genes = gene_expr_df.mean() < threshold
    scores.loc[low_expr_genes, :] = 0
    print(f"✅ Low-expressed genes (mean expression < {threshold} TPM) have been excluded from the TxRBP scores.")
    return scores

def collapse_transcript_scores_to_genes(scores, getBM, gene_collapse_method, dataset):
    """
    Collapse explainability scores from transcript level to gene level by aggregating the scores for each gene-RBP pair. 
    Specifically, this method transforms the given DataFrame of scores into a long format, merges it with gene information, 
    and then aggregates to find the maximum absolute score for each gene-RBP pair. 
    Finally, it creates a pivot table to summarize the results.

    Args:
        scores (pd.DataFrame): DataFrame containing e scores at the transcript level, where rows correspond to transcripts and columns 
                                correspond to RBPs. 
        getBM (pd.DataFrame): A DataFrame containing mappings of Transcript_IDs to Gene_IDs.
        gene_collapse_method (str): The method to use for collapsing scores to the gene level. 
                                     Supported methods include:
                                     - 'max_absolute_value': Selects the maximum absolute value score for 
                                       each gene-RBP pair.
        dataset (DeepRBPExpressionDataset): An instance of the DeepRBPExpressionDataset class, which 
                                             contains information about RBPs and transcripts, including 
                                             the names of RBPs to be used in the output.
    Returns:
    GeneScores: A namedtuple containing:
        - result_table (pd.DataFrame): A DataFrame with RBP and transcript details.
        - df_scores_GxRBP (pd.DataFrame): DataFrame with DeepLIFT scores (GxRBP).
    """ 
    print("🔍 Collapsing transcripts scores to genes...")
    GeneScores = namedtuple('GeneScores', ['result_table', 'df_scores_GxRBP']) 
    # Transform the wide format DataFrame into a long format
    deeplift_scores_long = scores.stack().reset_index()
    deeplift_scores_long.columns = ['Transcript_ID', 'RBP_ID', 'Score']  
    # Get RBP names from their IDs
    deeplift_scores_long['RBP_name'] = get_gene_info(deeplift_scores_long['RBP_ID'], getBM, return_type='names')
    # Merge with gene information to get Gene_IDs and additional metadata
    deeplift_scores_long = deeplift_scores_long.merge(getBM, on='Transcript_ID', how='left')
    # Determine collapse method based on configuration
    print(f'Determining gene collapse method {gene_collapse_method} based on configuration.')
    if gene_collapse_method == 'max_absolute_value':
        deeplift_scores_long['Score_abs'] = deeplift_scores_long['Score'].abs() 
        # Find the index of the maximum score for each Gene-RBP pair
        max_indices = deeplift_scores_long.loc[deeplift_scores_long.groupby(['Gene_ID', 'RBP_ID'])['Score_abs'].idxmax()]  
        result_table = max_indices[['RBP_ID', 'RBP_name', 'Gene_ID', 'Gene_name', 
                                        'Transcript_ID', 'Transcript_name', 
                                        'Transcript_biotype', 'Score']].reset_index(drop=True)
        # Count the number of transcripts per Gene_ID
        num_transcripts_per_gene = getBM['Gene_ID'].value_counts().reset_index()
        num_transcripts_per_gene.columns = ['Gene_ID', 'Num_trans_per_gene']
        # Merge the count of transcripts with result_table
        result_table = result_table.merge(num_transcripts_per_gene, on='Gene_ID', how='left')
        # Create a pivot table to summarize scores by Gene_ID and RBP_ID
        df_scores_GxRBP = result_table.pivot_table(
            index='Gene_ID', 
            columns='RBP_ID', 
            values='Score', 
            aggfunc='first'
        )
    print("✅ Collapsing transcripts scores to genes...")
    # Return the results as a namedtuple
    return GeneScores(result_table=result_table, df_scores_GxRBP=df_scores_GxRBP[dataset.rbp_names])

def get_gene_info(gene_ids_or_names: List[str], getBM: pd.DataFrame, return_type: str = 'names') -> List[Optional[str]]:
    """
    Retrieves Gene_names from Gene_IDs or Gene_IDs from Gene_names using the getBM DataFrame.

    Parameters:
    gene_ids_or_names (List[str]): A list of Gene_IDs or Gene_names.
    getBM (pd.DataFrame): A DataFrame that contains the relationship between Gene_ID and Gene_name.
    return_type (str): Indicates what to return:
                       - 'names': return Gene_names for given Gene_IDs
                       - 'ids': return Gene_IDs for given Gene_names

    Returns:
    List[Optional[str]]: A list of Gene_names or Gene_IDs corresponding to the provided input.
                         If a Gene_ID or Gene_name does not have a corresponding entry, None will be returned.
    """
    getBM_subset = getBM[['Gene_ID', 'Gene_name']].drop_duplicates()
    if return_type == 'names':
        # Convert Gene_IDs to Gene_names
        getBM_subset.set_index('Gene_ID', inplace=True)
        return getBM_subset.loc[gene_ids_or_names]['Gene_name'].values.tolist()
    elif return_type == 'ids':
        # Convert Gene_names to Gene_IDs
        getBM_subset.set_index('Gene_name', inplace=True)
        #gene_ids = getBM_subset.loc[gene_ids_or_names]
        # Prepare to retrieve Gene_IDs and print duplicates
        gene_ids = []
        for gene_name in gene_ids_or_names:
            if gene_name in getBM_subset.index:
                occurrences = getBM_subset.loc[gene_name]
                if len(occurrences) > 1:
                    print(f"Duplicate entries found for Gene_name '{gene_name}': {occurrences['Gene_ID'].values.tolist()}")
                # Use the first occurrence as the final result
                if isinstance(occurrences, pd.DataFrame):
                    first_gene_id = occurrences['Gene_ID'].iloc[0]  # Safe access for DataFrame
                else:
                    first_gene_id = occurrences['Gene_ID']  # Direct access for Series
                print(f"Using '{first_gene_id}' as Gene_ID for Gene_name '{gene_name}'")
                gene_ids.append(first_gene_id)
            else:
                gene_ids.append(None)
        return gene_ids
    else:
        raise ValueError("Invalid return_type. Use 'names' or 'ids'.")