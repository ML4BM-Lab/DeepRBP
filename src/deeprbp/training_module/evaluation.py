# src/deeprbp/training_module/evaluation.py

import os
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from typing import List, Dict, Any
from tqdm import tqdm 

from ..util.utils import CustomTensorDataset, adjust_batch_size, filter_data_by_sample_ids
from ..util.plots import scatter_real_vs_pred, plot_transcript_to_gene_ratio_distributions

def calculate_spearmanr(predictions, true_values):
    """
    Calculate Spearman correlation coefficient value between predicted and true values.

    Parameters:
    predictions (ndarray): Array of predicted values (flattened).
    true_values (ndarray): Array of true values (flattened).

    Returns:
    float: Spearman correlation coefficient.
    """
    spearman_corr = spearmanr(predictions, true_values)[0]
    return spearman_corr

def calculate_metrics(predictions, true_values):
    """
    Calculates general metrics like Spearman Correlation, MSE, and Pearson Correlation between flattened model 
    log2(tpm+1) predictions vs real values
   
    Args:
        predictions: Predicted values (in log2(tpm+1))
        true_values: True labels (in log2(tpm+1))
    Returns:
        A dictionary with the calculated metrics
    """
    spearman_corr = calculate_spearmanr(predictions, true_values)
    pearson_corr = pearsonr(predictions, true_values)[0]
    mse = mean_squared_error(predictions, true_values)
    r2 = r2_score(predictions, true_values)
    return {
        'spearman_corr': spearman_corr,
        'pearson_corr': pearson_corr,
        'mse': mse,
        'r2': r2
    }

def calculate_spearman_corr_per_gene(dataset, pred, label, getBM):
    """
    Calculate Spearman correlation between flattened model log2(tpm+1) predictions and real values for 
    the transcripts within each gene. This function calculates the correlation based on the ranking of 
    transcripts specific to each gene, rather than across all transcripts globally.

    Parameters:
    dataset (CustomTensorDataset): A dataset object containing gene and transcript information.
        - attribute gene_names (list): A list of gene IDs.
        - attribute trans_names (list): A list of transcript IDs.
        
    pred (ndarray): Array of predicted values (log2(tpm+1)).
    label (ndarray): Array of true values (log2(tpm+1)).
    getBM (DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.

    Returns:
    dict: A dictionary containing:
        - mean_corr (float): The mean Spearman correlation across all genes.
        - gene_corrs_dict (dict): A dictionary where keys are gene IDs and values are their 
          corresponding Spearman correlation coefficients.
        - mean_corr_max_trans (float): The mean Spearman correlation for the most expressed 
          transcripts of each gene.
        - gene_corrs_dict_max_trans (dict): A dictionary where keys are gene-transcript pairs 
          and values are their Spearman correlation coefficients for the most expressed transcripts.
    """
    gene_corrs_dict = {}
    gene_corrs_dict_max_trans = {}
    for gene_id in dataset.gene_names:
        related_transcripts = getBM[getBM['Gene_ID'] == gene_id]['Transcript_ID'].values
        if len(related_transcripts) > 0:
            # Find the indices of the related transcripts in the dataset and filter the predictions and true values for the related transcripts
            index = [np.where(np.array(dataset.trans_names) == transcript)[0][0] 
                     for transcript in related_transcripts if transcript in dataset.trans_names]
            filtered_preds = pred[:, index]
            filtered_true = label[:, index]
            # Step 1: Calculate Spearman's correlation for the current gene
            spearman_corr = calculate_spearmanr(filtered_preds.flatten(), filtered_true.flatten())
            gene_corrs_dict[gene_id] = spearman_corr
            # Step 2: Calculate the Spearman correlation for the most expressed transcript
            mean_expressions = np.mean(filtered_true, axis=0) # Calculate the true mean expression for each transcript
            max_trans_index = np.argmax(mean_expressions) # Find the index of the transcript with the highest mean expression
            spearman_corr_max_trans = calculate_spearmanr(filtered_preds[:, max_trans_index], filtered_true[:, max_trans_index])
            gene_corrs_dict_max_trans[f'{gene_id}-{related_transcripts[max_trans_index]}'] = spearman_corr_max_trans
    # Step 3: Calculate the mean correlation across all genes
    mean_corr = np.mean(list(gene_corrs_dict.values()))   
    # Step 4: Calculate the mean correlation for the most expressed transcripts
    mean_corr_max_trans = np.mean(list(gene_corrs_dict_max_trans.values()))  
    print("Mean Spearman Correlation across all genes:", mean_corr)
    print("Mean Spearman Correlation for the most expressed transcripts:", mean_corr_max_trans)
    return {
        'mean_corr_per_gene': mean_corr,
        'gene_corrs_dict': gene_corrs_dict,
        'mean_corr_max_trans_per_gene': mean_corr_max_trans,
        'gene_corrs_dict_max_trans': gene_corrs_dict_max_trans,
    }   

# def save_metrics_summary(
#     set_names_list: List[str], 
#     metrics_list: Dict[str, List[float]], 
#     output_path: str = None) -> pd.DataFrame:
#     """
#     Generates and saves a DataFrame from metrics for different keys.

#     Args:
#         set_names_list (List[str]): List of names for the key index.
#         metrics_list (Dict[str, List[float]]): Dictionary containing lists of metrics for each set.
#         output_path (str, optional): Path to save the DataFrame as a CSV file. Default is None.

#     Returns:
#         pd.DataFrame: DataFrame with metrics as columns and sets as rows.
#     """
#     #all_metrics = {set_name: metrics for set_name, metrics in zip(set_names_list, metrics_list)}
#     #df = pd.DataFrame(all_metrics).T
#     df = pd.DataFrame(metrics_list, index=set_names_list)
#     df.index.name = "Set"
#     if output_path:
#         df.to_csv(output_path, index=True)

def map_transcripts_and_filter_gene_matrix(
    gene_df: pd.DataFrame, getBM: pd.DataFrame, min_mean_expr: float = 5) -> pd.DataFrame:
    """
    Processes the gene expression dataframe by renaming columns with gene IDs,
    removing duplicates, and filtering genes based on mean expression.
    """
    gene_df.columns = getBM["Gene_ID"]
    gene_df = gene_df.loc[:, ~gene_df.columns.duplicated()]
    mean_expr = gene_df.mean()
    return gene_df[mean_expr[mean_expr > min_mean_expr].index]

def aggregate_transcripts_to_genes(
    transcript_df: pd.DataFrame, getBM: pd.DataFrame, filtered_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates transcript-level expression into gene-level by summing transcript values per gene.
    Ensures consistency with the filtered gene list.
    """
    transcript_df = np.power(2, transcript_df) - 1  # Convert log-transformed values to TPM
    transcript_df.columns = getBM["Gene_ID"]
    #aggregated_df = transcript_df.groupby(transcript_df.columns, axis=1).sum()
    aggregated_df = transcript_df.T.groupby(transcript_df.columns).sum().T
    return aggregated_df.loc[:, filtered_genes.columns]

def calculate_ratios(
    aggregated_df: pd.DataFrame, real_gene_df: pd.DataFrame) -> np.ndarray:
    """
    Calculates the ratio between the aggregated transcript-level expression and the real gene-level expression.
    Filters NaN and infinite values.
    """
    ratio = aggregated_df / real_gene_df
    ratio_values = ratio.values.flatten()
    return ratio_values[~np.isnan(ratio_values) & np.isfinite(ratio_values)]

def analyze_transcript_to_gene_ratios(
    getBM: pd.DataFrame,
    pred_df: pd.DataFrame,
    labels_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    category: str,
    output_dir: str,
    source_name: str):
    """
    Calculates and plots histograms of the predicted and labeled transcript-to-gene expression ratios.
    Saves the histogram in the specified output directory.
    """
    # 1. Process gene expression data
    getBM = getBM.set_index("Transcript_ID").loc[labels_df.columns.tolist()].reset_index()
    genes_filtered = map_transcripts_and_filter_gene_matrix(genes_df, getBM)
    # 2. Aggregate transcript expressions to gene level
    pred_agg = aggregate_transcripts_to_genes(pred_df, getBM, genes_filtered)
    label_agg = aggregate_transcripts_to_genes(labels_df, getBM, genes_filtered)
    # 3. Calculate ratios
    pred_ratios = calculate_ratios(pred_agg, genes_filtered)
    label_ratios = calculate_ratios(label_agg, genes_filtered)
    # 4. Plot histograms
    output_path = os.path.join(output_dir, f"histogram_ratio_{category}-{source_name}.png")
    plot_transcript_to_gene_ratio_distributions(pred_ratios, label_ratios, source_name, output_path)

def calculate_metrics_per_category( # ESTA FUNCION HAY QUE ACTUALIZAR!!
    test_data: Dict[str, pd.DataFrame], 
    trainer: Any,
    output_dir: str,
    set_name: str,
    sample_category: str,
    batch_size: int,
    source_name: str,
    plot_results: bool = False,
    getBM: str = None) -> List[Dict[str, float]]:
    """
    Calculates metrics for each category in the any dataset, including:
    - Performance metrics for each category  
    - Visualization of predicted vs actual values with scatter plots.
    - Distribution of transcript-to-gene ratios via histograms.

    Args:
        test_data (Dict[str, pd.DataFrame]): dataset to test including metadata, scaled features, and true labels.
        trainer (Any): Trainer object with `generate_predictions`.
        output_dir (str): Directory path to save results and plots.
        set_name (str): Name of the dataset being processed.
        sample_category (str): Column name for the sample categories in the metadata.
        batch_size (int): Batch size for the DataLoader.
        source_name (str): Name of the data source, used for saving results.
        plot_results (bool): Whether to generate plots for predictions and ratios.
        getBM (str): `getBM` dataframe for mapping transcripts to genes.

    Returns:
        List[Dict[str, float]]: List of metrics dictionaries, one for each category.
    """
    categories = test_data['metadata_df'].detailed_category.unique().tolist()
    metadata_df = test_data['metadata_df']
    metrics_list = []
    set_names_list = []
    # Use tqdm to create a progress bar for processing categories
    for category in tqdm(categories, desc="Processing categories", unit="category"):
        # Filter samples for the current category
        category_samples = metadata_df.loc[metadata_df[sample_category] == category].index
        # Filter test data for selected samples
        test_subset = filter_data_by_sample_ids(data=test_data, selected_sample_ids=category_samples)
        # Create DataLoader
        test_loader = DataLoader(
            CustomTensorDataset(test_subset), 
            batch_size=adjust_batch_size(test_subset['scaled_rbp_expr_log2p_tpm_df'], batch_size * 2)
        )
        # Generate predictions and calculate metrics
        predictions, true_values, predictions_matrix = trainer.generate_predictions(test_loader)
        metrics = calculate_metrics(predictions, true_values)
        metrics_list.append(metrics)
        set_names_list.append(category)
        # Prepare data for visualization
        labels_df = test_subset['trans_expr_log2p_tpm_df']
        pred_df = pd.DataFrame(
                    predictions_matrix, 
                    columns=list(test_subset['trans_expr_log2p_tpm_df'].columns), 
                    index=category_samples)
        genes_df = test_subset['gn_expr_each_iso_tpm_df'] #gn_expr_tpm_each_iso_df
        ###
        # Optional: Plot results
        if plot_results:
            scatter_real_vs_pred(
                category=category,
                source_name=source_name,
                metrics=metrics,
                pred=predictions,
                labels=true_values,
                output_dir=os.path.join(output_dir, 'scat_plot_real_vs_pred_value', source_name, set_name, category),
            )
            if getBM is not None:  # Ensure getBM is not None
                analyze_transcript_to_gene_ratios(
                    getBM=getBM,
                    pred_df=pred_df,
                    labels_df=labels_df,
                    genes_df=genes_df,
                    category=category,
                    output_dir=os.path.join(output_dir, 'pred_label_expression_ratio_histogram', source_name, set_name, category),
                    source_name=source_name,
                )
    return metrics_list, set_names_list
