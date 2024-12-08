# evaluation_utils.py
import os
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from typing import Tuple, List, Dict, Any, Union

from utils import CustomTensorDataset, adjust_batch_size, filter_data_by_sample_ids
from plots import scatter_real_vs_pred, plot_transcript_to_gene_ratio_distributions

def calculate_metrics(predictions, true_values):
    """
    Calculates general metrics like Spearman Correlation, MSE, and Pearson Correlation.
    Args:
        predictions: Predicted values
        true_values: True labels
    Returns:
        A dictionary with the calculated metrics
    """
    spearman_corr = spearmanr(predictions, true_values)[0]
    pearson_corr = pearsonr(predictions, true_values)[0]
    mse = mean_squared_error(predictions, true_values)
    r2 = r2_score(predictions, true_values)
    return {
        'spearman_corr': spearman_corr,
        'pearson_corr': pearson_corr,
        'mse': mse,
        'r2': r2
    }

def save_metrics_summary(
    set_names_list: List[str], 
    metrics_list: List[Dict[str, float]], 
    output_path: str = None) -> pd.DataFrame:
    """
    Generates and saves a DataFrame from metrics for different keys.

    Args:
        set_names_list (List[str]): List of names for the key index.
        metrics_list (List[Dict[str, float]]): List of dictionaries containing metrics for each set.
        output_path (str, optional): Path to save the DataFrame as a CSV file. Default is None.

    Returns:
        pd.DataFrame: DataFrame with metrics as columns and sets as rows.
    """
    all_metrics = {set_name: metrics for set_name, metrics in zip(set_names_list, metrics_list)}
    df = pd.DataFrame(all_metrics).T
    df.index.name = "Set"
    if output_path:
        df.to_csv(output_path, index=True)

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
    config: Dict[str, Union[str, int]],
):
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
    output_path = os.path.join(
        output_dir, f"histogram_ratio_{category}-{config['source_name']}.png"
    )
    plot_transcript_to_gene_ratio_distributions(pred_ratios, label_ratios, config["source_name"], output_path)


def calculate_metrics_per_category(
    test_data: Dict[str, pd.DataFrame], 
    config: Dict[str, Any], 
    trainer: Any,
    output_dir: str,
    set_name: str) -> List[Dict[str, float]]:
    """
    Calculates metrics for each category in the test dataset, including:
    - Performance metrics for each category (e.g., accuracy, precision, recall).
    - Visualization of predicted vs actual values with scatter plots.
    - Distribution of transcript-to-gene ratios via histograms.

    Args:
        test_data (Dict[str, pd.DataFrame]): Test dataset including metadata, scaled features, and true labels.
        config (Dict[str, Any]): Configuration dictionary containing category and batch size information.
        trainer (Any): Trainer object with `generate_predictions`.
        output_dir (str): Directory path to save results and plots.

    Returns:
        List[Dict[str, float]]: List of metrics dictionaries, one for each category.
    """
    categories = test_data['metadata_df'].detailed_category.unique().tolist()
    metadata_df = test_data['metadata_df']
    metrics_list = []
    set_names_list = []

    for category in categories:
        print(f"Processing category: {category}")
        # Filter samples for the current category
        category_samples = metadata_df.loc[metadata_df[config['sample_category']] == category].index
        
        # Filter test data for selected samples
        test_subset = filter_data_by_sample_ids(data=test_data, selected_sample_ids=category_samples)
        
        # Create DataLoader
        test_loader = DataLoader(
            CustomTensorDataset(test_subset), 
            batch_size=adjust_batch_size(test_subset['scaled_rbp_expr_df'], config['training']['batch_size'] * 2)
        )
        
        # Generate predictions and calculate metrics
        predictions, true_values, predictions_matrix = trainer.generate_predictions(test_loader)
        metrics = calculate_metrics(predictions, true_values)
        print(f"Metrics for category {category}: {metrics}")
        metrics_list.append(metrics)
        set_names_list.append(category)

        labels_df = test_subset['trans_expr_df']
        pred_df = pd.DataFrame(
                    predictions_matrix, 
                    columns=list(test_subset['trans_expr_df'].columns), 
                    index=category_samples)
        genes_df = test_subset['gene_expr_df']

        if config.get('plot_results', False):   
            # Plot real vs pred values
            scatter_real_vs_pred(
                category=category,
                config=config,
                metrics=metrics,
                pred=predictions,
                labels=true_values,
                output_dir=os.path.join(output_dir, 'scat_plot_real_vs_pred_value', config["source_name"], set_name, category)
            )

            if config['data_paths'].get('getBM_path'):
                getBM = pd.read_csv(config['data_paths']['getBM_path']) 
                analyze_transcript_to_gene_ratios(
                    getBM=getBM,
                    pred_df=pred_df,
                    labels_df=labels_df,
                    genes_df=genes_df,
                    category=category,
                    output_dir=os.path.join(output_dir, 'pred_label_expression_ratio_histogram', config["source_name"], set_name, category),
                    config=config,
                )
    return metrics_list, set_names_list
