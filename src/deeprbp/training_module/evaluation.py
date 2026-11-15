# src/deeprbp/training_module/evaluation.py

import os
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from typing import Dict
import torch
import lightning as L

from .tcga_codes import TCGA_CODE
from ..util.utils import filter_data_by_sample_ids, log_section_separator, print_if_main
from .preds2visualization import scatter_real_vs_pred, plot_transcript_to_gene_ratio_distributions, plot_small_multiples_real_vs_pred_grid

from ..data_preparation.prepare_data import DeepRBPExpressionDataset

def generate_predictions(trainer, model, dataloader):
    """Generates predictions and true values using the given model and dataloader."""
    results = trainer.predict(model, dataloaders=dataloader)
    predictions, true_values = [], []
    for preds, labels in results:
        predictions.append(preds.detach().cpu().numpy())
        true_values.append(labels.detach().cpu().numpy())
    return np.concatenate(predictions), np.concatenate(true_values)

def save_prediction_outputs(
    predictions: np.ndarray,
    true_values: np.ndarray,
    test_data: Dict[str, pd.DataFrame],
    getBM: pd.DataFrame,
    output_dir: str,
    set_name: str = "test",
) -> None:
    """
    Save full-split prediction matrices and companion data needed for downstream analyses.

    Saved matrices:
    - predicted transcript abundance in log2(TPM + 1)
    - true transcript abundance in log2(TPM + 1)
    - gene expression in TPM
    - sample metadata
    - transcript-to-gene mapping used by the evaluator
    """
    isoform_df = test_data["isoform_df"]
    gene_df = test_data["gene_df"]
    metadata_df = test_data["metadata_df"]

    if predictions.shape != isoform_df.shape:
        raise ValueError(
            f"Prediction shape {predictions.shape} does not match isoform_df shape {isoform_df.shape}"
        )

    if true_values.shape != isoform_df.shape:
        raise ValueError(
            f"True-value shape {true_values.shape} does not match isoform_df shape {isoform_df.shape}"
        )

    predictions_dir = os.path.join(output_dir, "prediction_matrices", set_name)
    os.makedirs(predictions_dir, exist_ok=True)

    transcript_ids = isoform_df.columns.tolist()
    sample_ids = isoform_df.index.tolist()

    pred_df = pd.DataFrame(
        predictions,
        index=sample_ids,
        columns=transcript_ids,
    )

    true_df = pd.DataFrame(
        true_values,
        index=sample_ids,
        columns=transcript_ids,
    )

    mapping_df = (
        getBM
        .drop_duplicates(subset=["Transcript_ID"])
        .set_index("Transcript_ID")
        .reindex(transcript_ids)
        .reset_index()
    )

    pred_df.to_csv(
        os.path.join(predictions_dir, "pred_trans_log2p_tpm.csv.gz"),
        compression="gzip",
    )

    true_df.to_csv(
        os.path.join(predictions_dir, "true_trans_log2p_tpm.csv.gz"),
        compression="gzip",
    )

    gene_df.to_csv(
        os.path.join(predictions_dir, "gene_tpm.csv.gz"),
        compression="gzip",
    )

    metadata_df.to_csv(
        os.path.join(predictions_dir, "metadata.csv"),
    )

    mapping_df.to_csv(
        os.path.join(predictions_dir, "transcript_to_gene_mapping.csv"),
        index=False,
    )

    print_if_main(
        f"[save_prediction_outputs] Saved prediction matrices to: {predictions_dir}"
    )

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

def calculate_metrics(predictions, true_values, gene_names, getBM, trans_names):
    """Calculates general and per-gene metrics."""
    metrics_general = calculate_general_metrics(true_values.flatten(), predictions.flatten())
    metrics_per_gene = spearmanr_per_gene(gene_names, getBM, trans_names, predictions, true_values)
    return {
        'spearman_corr': metrics_general['spearman_corr'],
        'pearson_corr': metrics_general['pearson_corr'],
        'mse': metrics_general['mse'],
        'r2': metrics_general['r2'],
        'mean_corr_per_gene': metrics_per_gene['mean_corr'],
        'mean_corr_max_trans_per_gene': metrics_per_gene['mean_corr_max']
    }

def spearmanr_per_gene(gene_names, getBM, trans_names, outputs, labels):  
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

def evaluate_and_visualize(
    trainer: L.Trainer,
    model: L.LightningModule,
    dm: L.LightningDataModule,
    output_dir: str,
    set_name: str):
    """
    Evaluates metrics for a dataset with no categories and plots the results.
    - Performance metrics
    - Visualization of predicted vs actual values with scatter plots.
    - Distribution of transcript-to-gene ratios (pred/real value) via histograms.

    Args:
        trainer (L.Trainer): PyTorch Lightning Trainer object to handle the prediction process.
        model (L.LightningModule): PyTorch Lightning model object to be used for making predictions.
        dm (L.LightningDataModule): An instance of DeepRBPDataModule used for managing data loading and preparation, 
                                providing a consistent interface for accessing datasets.
        output_dir (str): Directory path to save results and plots.
        set_name (str): Name of the dataset being processed (e.g., 'validation', 'test').
    This function does not return anything but saves the plots to the specified output directory.
    """
    getBM = dm.getBM.copy()
    print_if_main(f"[evaluate_and_visualize] 🔮 Generating predictions for '{set_name}'...")
    predictions, true_values = generate_predictions(trainer, model, dm.predict_dataloader(mode='test'))
    metrics = calculate_metrics(predictions, true_values, dm.test_dataset.gene_names, getBM, dm.test_dataset.trans_names)
    print_if_main(f"[evaluate_and_visualize] 📈 Metrics for '{set_name}': {metrics}\n")
    print_if_main(f"[evaluate_and_visualize] 📊 Plotting log2(tpm+1) predictions vs real values scatter plot for '{set_name}'.")
    scatter_real_vs_pred(
        category=set_name,
        metrics=dict(list(metrics.items())),
        pred=predictions.flatten(),
        labels=true_values.flatten(),
        output_dir=os.path.join(output_dir, 'scat_plot_real_vs_pred_value', set_name))
    print_if_main(f"\n[evaluate_and_visualize] 🔍 Analyzing transcript-to-gene ratios for '{set_name}'...")
    analyze_transcript_to_gene_ratios(
        getBM=getBM ,
        pred_df=pd.DataFrame(predictions, columns=list(dm.test_data['isoform_df'].columns), index=dm.test_data['isoform_df'].index),
        labels_df=dm.test_data['isoform_df'],
        genes_df=dm.test_data['gene_df'],  
        category=set_name,
        output_dir=os.path.join(output_dir, 'pred_label_expression_ratio_histogram', set_name))

def evaluate_and_visualize_metrics_by_category(  
    test_data: Dict[str, pd.DataFrame], 
    trainer: L.Trainer,
    model: L.LightningModule,
    dm: L.LightningDataModule,
    output_dir: str,
    set_name: str,
    plot_results: bool):
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
    
    This function does not return anything but saves the evaluation results and plots to the specified output directory.
    """
    getBM = dm.getBM.copy()
    metadata_df = test_data['metadata_df'].copy()
    categories = metadata_df[dm.sample_category].unique().tolist()
    results_list = [] # Initialize a list to store results
    grid_panels = []

    for index, category in enumerate(categories):
        print_if_main('\n')
        log_section_separator(f"Processing Category: {category} ({index + 1}/{len(categories)})")
        print_if_main('\n')

        # Filter samples for the current category
        category_samples = metadata_df.loc[metadata_df[dm.sample_category] == category].index
        test_data_copy = {key: df.copy() for key, df in test_data.items()} # Create a copy of the test data to avoid modifying the original data
        test_subset = filter_data_by_sample_ids(data=test_data_copy, selected_sample_ids=category_samples) # Filter test data for selected samples
        test_subdataset = DeepRBPExpressionDataset(test_subset, getBM) # Create DataLoader
        test_loader = DataLoader(test_subdataset, batch_size=len(test_subdataset)) #adjust_batch_size(test_subdataset, batch_size)
        print_if_main(f"[evaluate_and_visualize_metrics_by_category] 🔮 Generating predictions for category '{category}'...")
        predictions, true_values = generate_predictions(trainer, model, dm.predict_dataloader(mode='predict', custom_loader=test_loader))

        # Calculate metrics
        metrics = calculate_metrics(predictions, true_values, test_subdataset.gene_names, getBM, test_subdataset.trans_names)
        metrics['category'] = category
        print_if_main(f"[evaluate_and_visualize_metrics_by_category] 📈 Metrics for category '{category}': {metrics}\n")
        results_list.append(metrics)

        # Fill data for grid plot
        short = TCGA_CODE.get(category, category)  # usa código si existe
        grid_panels.append({
            "category": category,
            "short": short,
            "pred": predictions.flatten(),
            "true": true_values.flatten()
        })
        # Optional: Plot results  
        if trainer.global_rank == 0 or not torch.cuda.is_available(): # Do plot only for rank 0 (GPU-0 or CPU))
            if plot_results:
                print_if_main(f"[evaluate_and_visualize_metrics_by_category] 📊 Plotting log2(tpm+1) predictions vs real values scatter plot for category '{category}'.")
                scatter_real_vs_pred(
                    category=category,
                    metrics= {k: v for k, v in metrics.items() if k != 'category'}, #dict(list(metrics.items())[1:]),
                    pred=predictions.flatten(),
                    labels=true_values.flatten(),
                    output_dir=os.path.join(output_dir, 'scat_plot_real_vs_pred_value', set_name, category))
                print_if_main(f"\n[evaluate_and_visualize_metrics_by_category] 🔍 Analyzing transcript-to-gene ratios for category '{category}'...")
                analyze_transcript_to_gene_ratios(
                    getBM=getBM,
                    pred_df=pd.DataFrame(predictions, columns=list(test_subset['isoform_df'].columns), index=category_samples),
                    labels_df=test_subset['isoform_df'],
                    genes_df=test_subset['gene_df'],  
                    category=category,
                    output_dir=os.path.join(output_dir, 'pred_label_expression_ratio_histogram', set_name, category))     

    # Save results
    results_df = pd.DataFrame(results_list)
    results_df.to_csv(os.path.join(output_dir, f'{set_name}_tumor_category_results.csv'), index=False)
    print_if_main(f"[evaluate_and_visualize_metrics_by_category] ✅ Results for the {set_name} dataset have been successfully saved to files.")

    # Generate unique figure for scatter plots
    if (trainer.global_rank == 0 or not torch.cuda.is_available()) and plot_results:
        print_if_main(f"[evaluate_and_visualize_metrics_by_category] 🧩 Building small-multiples grid for '{set_name}'.")
        plot_small_multiples_real_vs_pred_grid(
            panels=grid_panels,
            set_name=set_name,
            output_dir=output_dir,
            order_codes=list(TCGA_CODE.values()),
            axis_range=(0, 15),          # si quieres rango idéntico en todos
            cmap="plasma",              # o "magma", "plasma"
            density_scale="log",         # más contraste en zonas densas
            share_density_norm=True,     # mismo mapeo de color en todos los paneles
            show_colorbar=True          # lo puedes activar si quieres comprobar la escala
        )
    print_if_main(f"[evaluate_and_visualize_metrics_by_category] ✅ Done for split '{set_name}'.")


#### #### #### #### #### #### #### #### #### #### #### ####  remove this when plot is definitive  
# output_dir = '/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor'
# panels = load_panels(os.path.join(output_dir, "test_panels.pkl.gz"))

# for p in panels:
#     if p.get("category") in (
#        "Pheochromocytoma_&_Paraganglioma",
#          "Pheochromocytoma_and_Paraganglioma",
#      ) or p.get("short") in (
#          "Pheochromocytoma_&_Paraganglioma",
#          "Pheochromocytoma_and_Paraganglioma",
#      ):
#         p["short"] = "PCPG"

# plot_small_multiples_real_vs_pred_grid(
#      panels=panels,
#      set_name="test",
#      output_dir=output_dir,
#      order_codes=list(TCGA_CODE.values()),
#      axis_range=(0, 15),          # si quieres rango idéntico en todos
#      cmap="plasma",              # o "magma", "plasma"
#      density_scale="log",         # más contraste en zonas densas
#      share_density_norm=True,     # mismo mapeo de color en todos los paneles
#      show_colorbar=True,          # lo puedes activar si quieres comprobar la escala
#      font_scale=1.40
# )

#### #### #### #### #### #### #### #### #### #### #### #### 

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

# ============================================================
# Isoform-usage evaluation
# Reconstructed for PLOS manuscript reproducibility
# ============================================================


def calculate_isoform_jsd(
    pred_df: pd.DataFrame,
    true_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    input_scale: str = "log2p1",
    min_mean_expr: float = 5.0,
    min_sample_expr: float = 5.0,
) -> dict:
    """Compare observed and predicted isoform distributions per gene/sample.

    Inputs: transcript abundances in log2(TPM+1) (or TPM if input_scale='tpm'),
    observed gene TPM, and a Transcript_ID/Gene_ID annotation mapping.
    Only genes with >=2 modelled transcripts and mean gene TPM > min_mean_expr
    are considered. Individual gene/sample pairs require gene TPM strictly
    > min_sample_expr. Each transcript distribution is normalized separately.

    Returns a dictionary with `pairs` (sample_id, Gene_ID, js_divergence,
    js_distance), `n_eligible_genes`, `n_candidate_pairs`, `n_invalid_pairs`.
    Jensen-Shannon divergence uses log base 2, and distance = sqrt(divergence).
    Invalid/nonfinite transcript distributions and zero-sum distributions are
    excluded and counted; nonfinite gene TPM is ignored for gene means and
    does not qualify for the per-sample threshold. No figures/files are made.
    """
    from scipy.spatial.distance import jensenshannon

    if input_scale not in ("log2p1", "tpm"):
        raise ValueError("input_scale must be 'log2p1' or 'tpm'")
    if not isinstance(pred_df, pd.DataFrame) or not isinstance(true_df, pd.DataFrame):
        raise TypeError("pred_df and true_df must be pandas DataFrames")
    if not pred_df.index.equals(true_df.index) or not pred_df.columns.equals(true_df.columns):
        raise ValueError("Prediction and observation sample/transcript axes differ")
    if not pred_df.index.equals(genes_df.index):
        raise ValueError("Gene TPM sample index differs from transcript matrices")
    if (not pred_df.index.is_unique or not pred_df.columns.is_unique
            or not genes_df.columns.is_unique or not genes_df.index.is_unique):
        raise ValueError("Duplicate sample, transcript or gene identifiers")
    if not {"Transcript_ID", "Gene_ID"}.issubset(mapping.columns):
        raise ValueError("mapping must contain Transcript_ID and Gene_ID")
    if mapping["Transcript_ID"].isna().any() or mapping["Transcript_ID"].duplicated().any():
        raise ValueError("Missing/duplicate Transcript_ID in annotation mapping")
    if not pred_df.columns.isin(mapping["Transcript_ID"]).all():
        raise ValueError("Unmapped modelled transcript")
    aligned_mapping = mapping.set_index("Transcript_ID").loc[pred_df.columns]
    if aligned_mapping["Gene_ID"].isna().any():
        raise ValueError("Missing Gene_ID for a modelled transcript")
    if not aligned_mapping["Gene_ID"].isin(genes_df.columns).all():
        raise ValueError("Mapped gene is absent from gene TPM matrix")

    # The original filtering function uses a strict > mean-TPM criterion.
    clean_genes = genes_df.where(np.isfinite(genes_df))
    filtered = filter_low_expressed_genes(clean_genes, min_mean_expr)
    transcript_groups = aligned_mapping.groupby("Gene_ID", sort=False).indices
    eligible = [g for g in filtered.columns if len(transcript_groups.get(g, ())) >= 2]
    if input_scale == "log2p1":
        with np.errstate(over="ignore", invalid="ignore"):
            pred_tpm = np.maximum(np.exp2(pred_df.to_numpy(dtype=float)) - 1, 0)
            true_tpm = np.maximum(np.exp2(true_df.to_numpy(dtype=float)) - 1, 0)
    else:
        pred_tpm = np.maximum(pred_df.to_numpy(dtype=float), 0)
        true_tpm = np.maximum(true_df.to_numpy(dtype=float), 0)

    parts = []
    candidate_pairs = 0
    for gene_id in eligible:
        indices = transcript_groups[gene_id]
        sample_pos = np.flatnonzero(filtered[gene_id].to_numpy() > min_sample_expr)
        candidate_pairs += len(sample_pos)
        if not len(sample_pos):
            continue
        p = pred_tpm[np.ix_(sample_pos, indices)]
        q = true_tpm[np.ix_(sample_pos, indices)]
        # Do not treat any NaN/Inf transcript as zero or partially sum it.
        finite = np.isfinite(p).all(axis=1) & np.isfinite(q).all(axis=1)
        p_sum, q_sum = p.sum(axis=1), q.sum(axis=1)
        finite &= np.isfinite(p_sum) & np.isfinite(q_sum) & (p_sum > 0) & (q_sum > 0)
        if not finite.any():
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            p = p[finite] / p_sum[finite, None]
            q = q[finite] / q_sum[finite, None]
            # scipy.spatial.distance.jensenshannon returns the DISTANCE.
            distances = jensenshannon(p, q, base=2, axis=1)
        keep = np.isfinite(distances)
        if keep.any():
            parts.append(pd.DataFrame({
                "sample_id": pred_df.index.to_numpy()[sample_pos[finite][keep]],
                "Gene_ID": gene_id,
                "js_divergence": distances[keep] ** 2,
                "js_distance": distances[keep],
            }))

    columns = ["sample_id", "Gene_ID", "js_divergence", "js_distance"]
    pairs = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    return {
        "pairs": pairs,
        "n_eligible_genes": len(eligible),
        "n_candidate_pairs": candidate_pairs,
        "n_invalid_pairs": candidate_pairs - len(pairs),
    }


def calculate_isoform_usage_sum(
    pred_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    input_scale: str = "log2p1",
    min_mean_expr: float = 5.0,
) -> dict:
    """Predicted summed transcript TPM / observed gene TPM per gene/sample.

    Select genes with observed mean gene TPM > min_mean_expr. Retain
    gene/sample observations whose observed gene TPM is positive, including
    those <= 5 TPM. Do NOT normalize predicted isoforms beforehand. Reuse
    the existing filter, log-to-TPM aggregation, and finite-ratio functions.

    Invalid/nonfinite ratios (including zero denominators) and transcript
    groups containing NaN/Inf are excluded; valid zero ratios are retained.
    Return `pairs` (sample_id, Gene_ID, ratio) plus gene/pair counts.
    """
    if input_scale not in ("log2p1", "tpm"):
        raise ValueError("input_scale must be 'log2p1' or 'tpm'")
    if not isinstance(pred_df, pd.DataFrame) or not isinstance(genes_df, pd.DataFrame):
        raise TypeError("pred_df and genes_df must be pandas DataFrames")
    if not pred_df.index.equals(genes_df.index):
        raise ValueError("Gene TPM sample index differs from transcript predictions")
    if (not pred_df.index.is_unique or not pred_df.columns.is_unique
            or not genes_df.index.is_unique or not genes_df.columns.is_unique):
        raise ValueError("Duplicate sample, transcript or gene identifiers")
    if not {"Transcript_ID", "Gene_ID"}.issubset(mapping.columns):
        raise ValueError("mapping must contain Transcript_ID and Gene_ID")
    if mapping["Transcript_ID"].isna().any() or mapping["Transcript_ID"].duplicated().any():
        raise ValueError("Missing/duplicate Transcript_ID in annotation mapping")
    if not pred_df.columns.isin(mapping["Transcript_ID"]).all():
        raise ValueError("Unmapped modelled transcript")
    aligned_mapping = mapping.set_index("Transcript_ID").loc[pred_df.columns].reset_index()
    if aligned_mapping["Gene_ID"].isna().any():
        raise ValueError("Missing Gene_ID for a modelled transcript")
    if not aligned_mapping["Gene_ID"].isin(genes_df.columns).all():
        raise ValueError("Mapped gene is absent from gene TPM matrix")

    clean_genes = genes_df.where(np.isfinite(genes_df))
    filtered = filter_low_expressed_genes(clean_genes, min_mean_expr)
    # The historical aggregator requires every retained gene to be modelled.
    # Unmodelled gene-expression columns cannot have a predicted isoform sum.
    filtered = filtered.loc[:, filtered.columns.isin(aligned_mapping["Gene_ID"])]
    # Preserve historical aggregate_transcripts_to_genes() conversion route,
    # while clamping negative transcript TPM as in the reconstructed baseline.
    if input_scale == "log2p1":
        pred_log = pred_df.clip(lower=0)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            pred_log = np.log2(1.0 + pred_df.clip(lower=0))
    predicted_by_gene = aggregate_transcripts_to_genes(pred_log, aligned_mapping, filtered)

    # Pandas groupby.sum() silently skips NaN; avoid turning a partially
    # nonfinite transcript distribution into an apparently valid gene sum.
    invalid_transcripts = pd.DataFrame(
        ~np.isfinite(pred_log.to_numpy(dtype=float)),
        index=pred_df.index,
        columns=aligned_mapping["Gene_ID"].to_numpy(),
    )
    invalid_gene_pairs = invalid_transcripts.T.groupby(level=0).any().T
    predicted_by_gene = predicted_by_gene.mask(invalid_gene_pairs.loc[:, filtered.columns])

    positive_genes = filtered.where(filtered > 0)
    # The historical function owns the numeric ratio calculation and its
    # finite-value filtering; reconstruct identifiers from the same mask.
    values = calculate_ratios(predicted_by_gene, positive_genes)
    with np.errstate(divide="ignore", invalid="ignore"):
        valid = np.isfinite(predicted_by_gene.to_numpy() / positive_genes.to_numpy())
    sample_pos, gene_pos = np.nonzero(valid)
    if len(values) != len(sample_pos):
        raise RuntimeError("Mismatch with historical calculate_ratios() ordering")
    pairs = pd.DataFrame({
        "sample_id": pred_df.index.to_numpy()[sample_pos],
        "Gene_ID": filtered.columns.to_numpy()[gene_pos],
        "ratio": values,
    })
    candidate_pairs = len(pred_df) * len(filtered.columns)
    return {
        "pairs": pairs,
        "n_eligible_genes": len(filtered.columns),
        "n_candidate_pairs": candidate_pairs,
        "n_invalid_pairs": candidate_pairs - len(pairs),
    }


def evaluate_isoform_usage_by_category(
    pred_df: pd.DataFrame,
    true_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    category_col: str = "detailed_category",
    input_scale: str = "log2p1",
    min_mean_expr: float = 5.0,
    min_sample_expr: float = 5.0,
) -> dict:
    """Compute per-category JSD and predicted isoform usage sums.

    Matrices and metadata must have identical sample indices in identical
    order; transcript predictions/labels must have identical column order.
    Mapping is aligned explicitly by Transcript_ID, independent of its row
    order. No paths, checkpoints, figures or output files are involved.

    Returns `jsd_by_category` and `sum_by_category` summary DataFrames,
    `jsd_pairs` and `sum_pairs` detailed DataFrames, and `pooled` global
    summaries (population SD, ddof=0). The JSD summaries contain distance
    statistics for comparison with the reconstructed PLOS Table S9; the
    divergence is separately reported to avoid conflating the two metrics.
    """
    if not (pred_df.index.equals(true_df.index)
            and pred_df.index.equals(genes_df.index)
            and pred_df.index.equals(metadata_df.index)):
        raise ValueError("Sample index/order is not identical across input matrices")
    if not pred_df.columns.equals(true_df.columns):
        raise ValueError("Transcript column/order differs between predictions and labels")
    if category_col not in metadata_df.columns:
        raise ValueError(f"Missing category column: {category_col}")
    if metadata_df[category_col].isna().any():
        raise ValueError("Missing category labels")
    if not pred_df.index.is_unique or not pred_df.columns.is_unique:
        raise ValueError("Duplicate sample or transcript identifiers")

    jsd_results, sum_results, jsd_parts, sum_parts = [], [], [], []
    for category in sorted(metadata_df[category_col].unique()):
        mask = metadata_df[category_col].eq(category)
        p, q, g = pred_df.loc[mask], true_df.loc[mask], genes_df.loc[mask]
        js = calculate_isoform_jsd(
            p, q, g, mapping, input_scale=input_scale,
            min_mean_expr=min_mean_expr, min_sample_expr=min_sample_expr,
        )
        sm = calculate_isoform_usage_sum(
            p, g, mapping, input_scale=input_scale,
            min_mean_expr=min_mean_expr,
        )
        j, r = js["pairs"], sm["pairs"]
        jsd_results.append({
            "category": category,
            "n_samples": int(mask.sum()),
            "n_eligible_genes": js["n_eligible_genes"],
            "n_candidate_pairs": js["n_candidate_pairs"],
            "n_invalid_pairs": js["n_invalid_pairs"],
            "n_jsd": len(j),
            "jsd_mean": j["js_distance"].mean(),
            "jsd_sd": j["js_distance"].std(ddof=0),
            "js_divergence_mean": j["js_divergence"].mean(),
            "js_divergence_sd": j["js_divergence"].std(ddof=0),
        })
        sum_results.append({
            "category": category,
            "n_samples": int(mask.sum()),
            "n_genes": sm["n_eligible_genes"],
            "n_valid_pairs": len(r),
            "n_invalid_pairs": sm["n_invalid_pairs"],
            "ratio_mean": r["ratio"].mean(),
            "ratio_sd": r["ratio"].std(ddof=0),
        })
        jsd_parts.append(j.assign(category=category))
        sum_parts.append(r.assign(category=category))

    j_all = pd.concat(jsd_parts, ignore_index=True)
    r_all = pd.concat(sum_parts, ignore_index=True)
    pooled = pd.DataFrame([{
        "n_categories": len(jsd_results),
        "n_jsd": len(j_all),
        "jsd_mean": j_all["js_distance"].mean(),
        "jsd_sd": j_all["js_distance"].std(ddof=0),
        "js_divergence_mean": j_all["js_divergence"].mean(),
        "js_divergence_sd": j_all["js_divergence"].std(ddof=0),
        "n_valid_pairs": len(r_all),
        "ratio_mean": r_all["ratio"].mean(),
        "ratio_sd": r_all["ratio"].std(ddof=0),
    }])
    return {
        "jsd_by_category": pd.DataFrame(jsd_results),
        "sum_by_category": pd.DataFrame(sum_results),
        "jsd_pairs": j_all,
        "sum_pairs": r_all,
        "pooled": pooled,
    }
