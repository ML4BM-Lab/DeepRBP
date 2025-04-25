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

def calculate_spearman_corr_per_gene(gene_names, trans_names, pred, label, getBM):
    """
    Calculate Spearman correlation between flattened model log2(tpm+1) predictions and real values for 
    the transcripts within each gene. This function calculates the correlation based on the ranking of 
    transcripts specific to each gene, rather than across all transcripts globally.

    Parameters:
    gene_names (list): A list of gene IDs used in our model.
    trans_names (list): A list of transcript IDs used in our model.
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
    nan_count_genes = 0 
    nan_count_max_trans = 0 
    for gene_id in gene_names:
        related_transcripts = getBM[getBM['Gene_ID'] == gene_id]['Transcript_ID'].values
        if len(related_transcripts) > 0:
            # Find the indices of the related transcripts in the dataset and filter the predictions and true values for the related transcripts
            index = [np.where(np.array(trans_names) == transcript)[0][0] 
                     for transcript in related_transcripts if transcript in trans_names]
            filtered_preds = pred[:, index]
            filtered_true = label[:, index]
            # Step 1: Calculate Spearman's correlation for the current gene
            spearman_corr = calculate_spearmanr(filtered_preds.flatten(), filtered_true.flatten())
            gene_corrs_dict[gene_id] = spearman_corr
            # Check if the correlation is NaN
            if np.isnan(spearman_corr):
                nan_count_genes += 1
            # Step 2: Calculate the Spearman correlation for the most expressed transcript
            mean_expressions = np.mean(filtered_true, axis=0) # Calculate the true mean expression for each transcript
            max_trans_index = np.argmax(mean_expressions) # Find the index of the transcript with the highest mean expression
            spearman_corr_max_trans = calculate_spearmanr(filtered_preds[:, max_trans_index], filtered_true[:, max_trans_index])
            gene_corrs_dict_max_trans[f'{gene_id}-{related_transcripts[max_trans_index]}'] = spearman_corr_max_trans
            if np.isnan(spearman_corr_max_trans):
                nan_count_max_trans += 1  # Increment NaN count for max trans correlations
    # Step 3: Calculate the mean correlation across all genes, excluding NaN values
    valid_corrs = [corr for corr in gene_corrs_dict.values() if not np.isnan(corr)]
    mean_corr = np.mean(valid_corrs) if valid_corrs else np.nan  # Avoid NaN if there are no valid correlations
    # Step 4: Calculate the mean correlation for the most expressed transcripts
    valid_corrs_max_trans = [corr for corr in gene_corrs_dict_max_trans.values() if not np.isnan(corr)]
    mean_corr_max_trans = np.mean(valid_corrs_max_trans) if valid_corrs_max_trans else np.nan  # Avoid NaN if there are no valid correlations
    # Log the results
    if nan_count_genes > 0:
        print(f"Warning: Excluded {nan_count_genes} gene correlations due to constant inputs (NaN values).")
    if nan_count_max_trans > 0:
        print(f"Warning: Excluded {nan_count_max_trans} maximum transcript correlations due to constant inputs (NaN values).")
    return {
        'mean_corr_per_gene': mean_corr,
        'gene_corrs_dict': gene_corrs_dict,
        'mean_corr_max_trans_per_gene': mean_corr_max_trans,
        'gene_corrs_dict_max_trans': gene_corrs_dict_max_trans,
    }   

def filter_low_expressed_genes(
    gene_df: pd.DataFrame, 
    min_mean_expr: int = 5
) -> pd.DataFrame:
    """
    Filters out genes with mean expression below a specified threshold
    and returns the filtered dataframe.

    Args:
        gene_df (pd.DataFrame): DataFrame containing gene expression data.
        min_mean_expr (float): Minimum mean expression threshold for filtering genes.

    Returns:
        pd.DataFrame: Filtered gene expression DataFrame.
    """
    # Calculate mean expression for each gene
    mean_expr = gene_df.mean()
    # Filter genes based on mean expression
    filtered_genes = gene_df[mean_expr[mean_expr > min_mean_expr].index]
    # Print the filtering details
    print(f"Filtering out low-expressed genes with a minimum mean expression of: {min_mean_expr}")
    print(f"Number of genes remaining after filtering: {filtered_genes.shape[1]}")
    return filtered_genes

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
    genes_filtered = filter_low_expressed_genes(genes_df)
    # 2. Aggregate transcript expressions to gene level
    pred_agg = aggregate_transcripts_to_genes(pred_df, getBM, genes_filtered)
    label_agg = aggregate_transcripts_to_genes(labels_df, getBM, genes_filtered)
    # 3. Calculate ratios
    pred_ratios = calculate_ratios(pred_agg, genes_filtered)
    label_ratios = calculate_ratios(label_agg, genes_filtered)
    # 4. Plot histograms
    output_path = os.path.join(output_dir, f"histogram_ratio_{category}-{source_name}.png")
    plot_transcript_to_gene_ratio_distributions(pred_ratios, label_ratios, source_name, output_path)

def calculate_metrics_per_category(  
    test_data: Dict[str, pd.DataFrame], 
    trainer: Any,
    output_dir: str,
    set_name: str,
    config: object,
    source_name: str,
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
        config
        source_name (str): Name of the data source, used for saving results.
        getBM (str): `getBM` dataframe for mapping transcripts to genes.

    Returns:
        List[Dict[str, float]]: List of metrics dictionaries, one for each category.
    """
    sample_category = config.get('sample_category')
    categories = test_data['metadata_df'][sample_category].unique().tolist()
    metadata_df = test_data['metadata_df']
    # Initialize a list to store results
    results_list = []
    # Use tqdm to create a progress bar for processing categories
    for category in tqdm(categories, desc="Processing categories", unit="category"):
        print(f"Currently processing category: {category}")
        # Filter samples for the current category
        category_samples = metadata_df.loc[metadata_df[sample_category] == category].index
        # Create a copy of the test data to avoid modifying the original data
        test_data_copy = {key: df.copy() for key, df in test_data.items()}
        # Filter test data for selected samples
        test_subset = filter_data_by_sample_ids(data=test_data_copy, selected_sample_ids=category_samples)
        # Create DataLoader
        test_subdataset = CustomTensorDataset(test_subset, getBM)
        test_loader = DataLoader(
            test_subdataset,
            batch_size=adjust_batch_size(test_subdataset, config.get('val_batch_size'))
        )
        # Generate predictions
        preds, labels = trainer.generate_predictions(test_loader) 
        # Calculate metrics
        metrics_general = calculate_metrics(preds.flatten(), labels.flatten())  
        # Calculate correlations per gene
        corr_per_gene = calculate_spearman_corr_per_gene(test_subdataset.gene_names, 
                                                         test_subdataset.trans_names, 
                                                         preds, labels, getBM)  
        # Collect results for the current category
        metrics = {
            'category': category,
            'spearman_corr': metrics_general['spearman_corr'],
            'pearson_corr': metrics_general['pearson_corr'],
            'mse': metrics_general['mse'],
            'r2': metrics_general['r2'],
            'mean_corr_per_gene': corr_per_gene['mean_corr_per_gene'],
            'mean_corr_max_trans_per_gene': corr_per_gene['mean_corr_max_trans_per_gene']
        }
        # Append the metrics dictionary to the results list
        results_list.append(metrics)
        # Optional: Plot results
        if config.get('plot_results'):
            scatter_real_vs_pred(
                category=category,
                source_name=source_name,
                metrics=dict(list(metrics.items())[1:]),
                pred=preds.flatten(),
                labels=labels.flatten(),
                output_dir=os.path.join(output_dir, 'scat_plot_real_vs_pred_value', source_name, set_name, category),
            )
            if getBM is not None:  # Ensure getBM is not None
                analyze_transcript_to_gene_ratios(
                    getBM=getBM,
                    pred_df=pd.DataFrame(preds, columns=list(test_subset['isoform_df'].columns), index=category_samples),
                    labels_df=test_subset['isoform_df'],
                    genes_df=test_subset['gene_df'],  
                    category=category,
                    output_dir=os.path.join(output_dir, 'pred_label_expression_ratio_histogram', source_name, set_name, category),
                    source_name=source_name
                )
    # Save results
    results_df = pd.DataFrame(results_list)
    if config.get('save_results'):
        print(f"📂 Saving metrics per category summary for the {set_name} dataset...")
        results_df.to_csv(os.path.join(output_dir, f'{set_name}_tumor_category_results.csv'), index=False)
        print(f"✅ Results for the {set_name} dataset have been successfully saved to files.")

    

  