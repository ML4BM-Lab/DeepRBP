# src/deeprbp/training_module/evaluation.py

import os
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from typing import List, Dict
import torch
import lightning as L

from ..util.utils import filter_data_by_sample_ids, log_section_separator, print_if_main
from .preds2visualization import scatter_real_vs_pred, plot_transcript_to_gene_ratio_distributions

from ..data_preparation.prepare_data import DeepRBPExpressionDataset

def calculate_general_metrics(true_values, predictions):  
    """
    Calculates general metrics like Spearman Correlation, MSE, and Pearson Correlation between flattened model 
    log2(tpm+1) predictions vs real values
   
    Args:
        predictions: Predicted values (in log2(tpm+1))
        true_values: True labels (in log2(tpm+1))
    Returns:
        A dictionary with the calculated metrics
    """
    spearman_corr = spearmanr(true_values, predictions)[0]
    pearson_corr = pearsonr(true_values, predictions)[0]
    mse = mean_squared_error(true_values, predictions)
    r2 = r2_score(true_values, predictions)
    return {
        'spearman_corr': spearman_corr,
        'pearson_corr': pearson_corr,
        'mse': mse,
        'r2': r2}

def spearmanr_per_gene(gene_names, getBM, trans_names, outputs, labels): # new
    """
    Calculate Spearman correlation between flattened model log2(tpm+1) predictions and real values for 
    the transcripts within each gene. This function calculates the correlation based on the ranking of 
    transcripts specific to each gene, rather than across all transcripts globally.

    Parameters:
    gene_names (list): A list of gene IDs used in our model.
    getBM (DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.
    trans_names (list): A list of transcript IDs used in our model.
    outputs (ndarray): Array of predicted values (log2(tpm+1)).
    labels (ndarray): Array of true values (log2(tpm+1)).
    
    Returns:
    dict: A dictionary containing:
        - mean_corr (float): The mean Spearman correlation across all genes.
        - gene_corrs_dict (dict): A dictionary where keys are gene IDs and values are their 
          corresponding Spearman correlation coefficients.
        - mean_corr_max (float): The mean Spearman correlation for the most expressed 
          transcripts of each gene.
        - gene_corrs_dict_max (dict): A dictionary where keys are gene-transcript pairs 
          and values are their Spearman correlation coefficients for the most expressed transcripts.
    """
    gene_corrs_dict = {}
    gene_corrs_dict_max = {}
    nan_count_genes = 0  # Counter for NaN correlations
    nan_count_max_trans = 0  # Counter for NaN max transcript correlations
    # Loop through each gene to calculate Spearman correlations
    for gene_id in gene_names:
        related_transcripts = getBM[getBM['Gene_ID'] == gene_id]['Transcript_ID'].values
        # Find and filter the indices of the related transcripts in predictions and actual values tensor
        index = [np.where(np.array(trans_names) == transcript)[0][0] for transcript in related_transcripts if transcript in trans_names]
        filtered_preds = outputs[:, index]
        filtered_true = labels[:, index]
        # Step 1: Calculate Spearman's correlation for the current gene
        corr = spearmanr(filtered_preds.flatten(), filtered_true.flatten())[0]
        ###
        # Step 2: Calculate the Spearman correlation for the most expressed transcript
        mean_expressions = np.mean(filtered_true, axis=0) # Calculate the true mean expression for each transcript
        max_trans_index = np.argmax(mean_expressions) # Find the index of the transcript with the highest mean expression
        corr_max = spearmanr(filtered_preds[:, max_trans_index], filtered_true[:, max_trans_index])[0]
        ###
        # Step 3: Add the calculations of the current gene
        gene_corrs_dict[gene_id] = corr
        gene_corrs_dict_max[f'{gene_id}-{related_transcripts[max_trans_index]}'] = corr_max
        # Increment NaN counters if correlations are NaN
        if np.isnan(corr):
            nan_count_genes += 1
        if np.isnan(corr_max):
            nan_count_max_trans += 1
        ###
    # Step 4: Calculate the mean correlation across all genes
    corrs = [corr for corr in gene_corrs_dict.values() if not np.isnan(corr)]
    mean_corr = np.mean(corrs)
    # Step 5: Calculate the mean correlation for the most expressed transcripts
    corrs_max = [corr for corr in gene_corrs_dict_max.values() if not np.isnan(corr)]
    mean_corr_max = np.mean(corrs_max) 
    # Print the results about NaN exclusions
    if nan_count_genes > 0:
        print_if_main(f"\n[spearmanr_per_gene] Excluded {nan_count_genes} gene correlations due to NaN values.")
    if nan_count_max_trans > 0:
        print_if_main(f"[spearmanr_per_gene] Excluded {nan_count_max_trans} maximum transcript correlations due to NaN values.\n")
    return {
        'mean_corr': mean_corr,
        'gene_corrs_dict': gene_corrs_dict,
        'mean_corr_max': mean_corr_max,
        'gene_corrs_dict_max': gene_corrs_dict_max
    }

def evaluate_and_visualize_metrics_by_category(  
    test_data: Dict[str, pd.DataFrame], 
    trainer: L.Trainer,
    model: L.LightningModule,
    dm: L.LightningDataModule,
    output_dir: str,
    set_name: str,
    plot_results: bool) -> List[Dict[str, float]]:
    """
    Evaluates metrics for each category in a dataset and optionally plots the results.
    - Performance metrics for each category  
    - Visualization of predicted vs actual values with scatter plots.
    - Distribution of transcript-to-gene ratios (pred/real value) via histograms.

    Args:
        test_data (Dict[str, pd.DataFrame]): A dictionary containing the dataset to test, 
                                              including metadata, scaled features, and true labels.
        trainer (L.Trainer): PyTorch Lightning Trainer object to handle the prediction process.
        model (L.LightningModule): PyTorch Lightning model object to be used for making predictions.
        dm (L.LightningDataModule): An instance of DeepRBPDataModule used for managing data loading and preparation, 
                                providing a consistent interface for accessing datasets.
        output_dir (str): Directory path to save results and plots.
        set_name (str): Name of the dataset being processed (e.g., 'validation', 'test').
        plot_results (bool): A flag indicating whether to generate and save visualizations of the results. 
                             If True, scatter plots and histograms will be created and saved in the output directory.
    Returns:
        List[Dict[str, float]]: List of metrics dictionaries, one for each category.
    """
    getBM = dm.getBM.copy()
    metadata_df = test_data['metadata_df'].copy()
    categories = metadata_df[dm.sample_category].unique().tolist()
    # Initialize a list to store results
    results_list = []
    total_categories = len(categories)
    for index, category in enumerate(categories):
        print_if_main('\n')
        log_section_separator(f"Processing Category: {category} ({index + 1}/{total_categories})")
        print_if_main('\n')
        # Filter samples for the current category
        category_samples = metadata_df.loc[metadata_df[dm.sample_category] == category].index
        # Create a copy of the test data to avoid modifying the original data
        test_data_copy = {key: df.copy() for key, df in test_data.items()}
        # Filter test data for selected samples
        test_subset = filter_data_by_sample_ids(data=test_data_copy, selected_sample_ids=category_samples)
        # Create DataLoader
        test_subdataset = DeepRBPExpressionDataset(test_subset, getBM )
        test_loader = DataLoader(
            test_subdataset,
            batch_size=len(test_subdataset) #adjust_batch_size(test_subdataset, batch_size)
        )
        # Generate predictions
        print_if_main(f"[evaluate_and_visualize_metrics_by_category] 🔮 Generating predictions for category '{category}'...")
        results = trainer.predict(model, dataloaders=dm.predict_dataloader(mode='predict', custom_loader=test_loader))
        true_values = []
        predictions = []
        for preds, labels in results:
            predictions.append(preds.detach().cpu().numpy())
            true_values.append(labels.detach().cpu().numpy())
        # Concatenate predictions and true values
        predictions = np.concatenate(predictions)
        true_values = np.concatenate(true_values)
        # Calculate metrics
        metrics_general = calculate_general_metrics(true_values.flatten(), predictions.flatten())  
        metrics_per_gene = spearmanr_per_gene(test_subdataset.gene_names,
                                              getBM ,
                                              test_subdataset.trans_names, 
                                              predictions, true_values)
        # Collect results for the current category
        metrics = {
            'category': category,
            'spearman_corr': metrics_general['spearman_corr'],
            'pearson_corr': metrics_general['pearson_corr'],
            'mse': metrics_general['mse'],
            'r2': metrics_general['r2'],
            'mean_corr_per_gene': metrics_per_gene['mean_corr'],
            'mean_corr_max_trans_per_gene': metrics_per_gene['mean_corr_max']
        }
        print_if_main(f"[evaluate_and_visualize_metrics_by_category] 📈 Metrics for category '{category}': {metrics}\n")
        # Append the metrics dictionary to the results list
        results_list.append(metrics)
        # Optional: Plot results
        # Do plot only for rank 0 (GPU-0 or CPU))
        if trainer.global_rank == 0 or not torch.cuda.is_available():
            if plot_results:
                print_if_main(f"[evaluate_and_visualize_metrics_by_category] 📊 Plotting log2(tpm+1) predictions vs real values scatter plot for category '{category}'.")
                scatter_real_vs_pred(
                    category=category,
                    metrics=dict(list(metrics.items())[1:]),
                    pred=predictions.flatten(),
                    labels=true_values.flatten(),
                    output_dir=os.path.join(output_dir, 'scat_plot_real_vs_pred_value', set_name, category))
                print_if_main(f"\n[evaluate_and_visualize_metrics_by_category] 🔍 Analyzing transcript-to-gene ratios for category '{category}'...")
                analyze_transcript_to_gene_ratios(
                    getBM=getBM ,
                    pred_df=pd.DataFrame(predictions, columns=list(test_subset['isoform_df'].columns), index=category_samples),
                    labels_df=test_subset['isoform_df'],
                    genes_df=test_subset['gene_df'],  
                    category=category,
                    output_dir=os.path.join(output_dir, 'pred_label_expression_ratio_histogram', set_name, category))
    # Save results
    results_df = pd.DataFrame(results_list)
    results_df.to_csv(os.path.join(output_dir, f'{set_name}_tumor_category_results.csv'), index=False)
    print_if_main(f"[evaluate_and_visualize_metrics_by_category] ✅ Results for the {set_name} dataset have been successfully saved to files.")

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
    print(f"\nFiltering out low-expressed genes with a minimum mean expression of: {min_mean_expr}")
    print(f"Number of genes remaining after filtering: {filtered_genes.shape[1]}\n")
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
    output_dir: str):
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
    output_path = os.path.join(output_dir, f"histogram_ratio_{category}.png")
    plot_transcript_to_gene_ratio_distributions(pred_ratios, label_ratios, output_path)
