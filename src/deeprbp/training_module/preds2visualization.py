# src/deeprbp/util/plots.py

import os 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List, Union

def scatter_real_vs_pred(
    category: str,
    metrics: Dict[str, float],
    pred: Union[List[float], np.ndarray],
    labels: Union[List[float], np.ndarray],
    output_dir: str
) -> None:
    """
    Creates a scatter plot comparing predicted vs. real values for a regression model.
    Includes regression line and performance metrics for clear visualization.

    Parameters:
    -----------
    category : str
        The category of the samples (e.g., cancer type or dataset name).
    metrics : Dict[str, float]
        Dictionary of evaluation metrics: Spearman correlation, Pearson correlation, MSE, R².
    pred : Union[List[float], np.ndarray]
        List or array of predicted values.
    labels : Union[List[float], np.ndarray]
        List or array of true/real values.
    output_dir : str
        Base path to save the plot as a PNG image.

    Returns:
    --------
    None
        The function saves a PNG image at the specified path.
    """
    plt.figure(figsize=(8, 6))
    plt.xlabel('Predicted Values', fontsize=16, fontweight='bold')
    plt.ylabel('Real Values', fontsize=16, fontweight='bold')
    plt.title(f'Real vs Predicted: {category})', fontsize=18, fontweight='bold')
    sns.regplot(
        x=pred, y=labels,
        scatter_kws={'alpha': 0.3, 'color': 'blue'},
        line_kws={'color': 'red', 'lw': 2}
    )
    legend_text = (
        f"Spearman Corr: {metrics['spearman_corr']:.2f}\n"
        f"Pearson Corr: {metrics['pearson_corr']:.2f}\n"
        f"MSE: {metrics['mse']:.2f}\n"
        f"R²: {metrics['r2']:.2f}\n"
        f"Spearman Corr per Gene: {metrics['mean_corr_per_gene']:.2f} \n"
        f"Spearman Corr per Gene (max trans): {metrics['mean_corr_max_trans_per_gene']:.2f} \n"
    )
    plt.text(
        0.05, 0.95, legend_text, transform=plt.gca().transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle="round", edgecolor="black", facecolor="white")
    )
    sns.set_style("whitegrid")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{category}.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved at: {output_path}")
 
def plot_all_metrics_history(metrics_df: pd.DataFrame, output_dir: str) -> None:
    metrics_pairs = {
        'Correlation Pearson History': ('train_corr_pearson', 'val_corr_pearson'),
        'Correlation Spearman History': ('train_corr_spearman', 'val_corr_spearman'),
        'Loss History': ('train_loss', 'val_loss'),
        'R2 Score History': ('train_r2', 'val_r2')
    }
    for title, (train_metric, val_metric) in metrics_pairs.items():
        print(metrics_df[train_metric].tolist())
        print(metrics_df[val_metric].tolist())
        plot_metric_history(
            train_history=metrics_df[train_metric].tolist(),
            val_history=metrics_df[val_metric].tolist(),
            title=title,
            output_dir=output_dir,
            plot_name=title.lower().replace(' ', '_')
        )

def plot_metric_history(train_history: List[float], val_history: List[float], 
                        title: str, 
                        output_dir: str, 
                        plot_name: str) -> None:
    """
    Plots the training and validation metrics history to visualize the performance of the model over time.

    Parameters:
    - train_history (List[float]): List of training metric values for each epoch.
    - val_history (List[float]): List of validation metric values for each epoch.
    - title (str): Title of the plot.
    - output_dir (str): Path to save the plot as a .png file.
    - plot_name (str): Name of the plot file.

    Returns:
    - None: The function will save the plot.
    """
    os.makedirs(output_dir, exist_ok=True)
    plt.style.use('seaborn-v0_8-muted')
    plt.figure(figsize=(8, 5))
    # Imprimir valores a graficar
    print(f"Plotting {title}:")
    print(f"Train values: {train_history}")
    print(f"Validation values: {val_history}")
    # Plot both training and validation metrics
    plt.plot(train_history, label='Training', color='royalblue', linestyle='-', linewidth=2.5)
    plt.plot(val_history, label='Validation', color='darkorange', linestyle='--', linewidth=2.5)
    plt.title(title, fontsize=18, fontweight='bold', pad=10)
    plt.xlabel('Epoch', fontsize=14, labelpad=8)
    print('\n\n')
    plt.ylabel('Performance', fontsize=14, labelpad=8)  # Dynamically set the ylabel
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12, loc='upper right', frameon=True, shadow=True, fancybox=True)
    plt.tight_layout(pad=2)
    # Save the plot with the specified name
    plt.savefig(os.path.join(output_dir, f'{plot_name}.png'), dpi=300)
    plt.close()
    print(f"Plot saved to {output_dir}/{plot_name}.png")

def plot_transcript_to_gene_ratio_distributions(
    ratios_pred: np.ndarray,
    ratios_label: np.ndarray,
    output_path: str,
):
    """
    Plots histograms of predicted and labeled transcript-to-gene expression ratios.
    Saves the plot as a high-quality PNG file.
    """
    # Set figure size and style
    plt.figure(figsize=(14, 7))
    plt.style.use('ggplot')  # Elegant grid style
    # Define color palette
    color_pred = "#1f77b4"  # Blue
    color_label = "#ff7f0e"  # Orange
    # Predicted ratio histogram
    plt.subplot(1, 2, 1)
    plt.hist(ratios_pred, bins=40, color=color_pred, edgecolor='black', linewidth=1.2)
    mean_pred = np.mean(ratios_pred)
    std_pred = np.std(ratios_pred)
    plt.title(f"Predicted Gene Ratio (Mean > 5 TPM)", fontsize=12, fontweight='bold')
    plt.xlabel("Transcript-to-Gene Ratio", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.axvline(mean_pred, color=color_pred, linestyle='dashed', linewidth=1.5)
    plt.text(mean_pred + 0.1, plt.ylim()[1] * 0.9, f'Mean: {mean_pred:.3f}\nSTD: {std_pred:.3f}', color=color_pred, fontsize=9)
    # Labeled ratio histogram
    plt.subplot(1, 2, 2)
    plt.hist(ratios_label, bins=40, color=color_label, edgecolor='black', linewidth=1.2)
    mean_label = np.mean(ratios_label)
    std_label = np.std(ratios_label)
    plt.title(f"Labeled Gene Ratio (Mean > 5 TPM)", fontsize=12, fontweight='bold')
    plt.xlabel("Transcript-to-Gene Ratio", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.axvline(mean_label, color=color_label, linestyle='dashed', linewidth=1.5)
    plt.text(mean_label + 0.1, plt.ylim()[1] * 0.9, f'Mean: {mean_label:.3f}\nSTD: {std_label:.3f}', color=color_label, fontsize=9)
    # Adjust layout
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_path}")




# esta hay que moverla al explainer folder bro
def plot_distributions_and_roc_with_thresholds(df_current_rbp, rbp_id, optimal_threshold, fpr, tpr, optimal_idx, auc_score, path_save):
    """
    Plots the distributions of scores and the ROC curve with the optimal threshold.

    Args:
        df_current_rbp (DataFrame): DataFrame containing the absolute scores and POSTAR labels.
        rbp_id (str): The ID of the RNA Binding Protein (RBP).
        optimal_threshold (float): The optimal threshold for classification.
        fpr (array-like): False positive rates for the ROC curve.
        tpr (array-like): True positive rates for the ROC curve.
        optimal_idx (int): Index of the optimal threshold in the fpr and tpr arrays.
        auc_score (float): Calculated auc score between a particular RBP Postar_Score and RBP explainability score.
        path_save (str): Path where the figure will be saved.
    """
    # Validate input DataFrame
    if not {'Score', 'Postar_Score'}.issubset(df_current_rbp.columns):
        raise ValueError("DataFrame must contain 'Score' and 'Postar_Score' columns.")
    plt.figure(figsize=(12, 4))
    sns.set(style="whitegrid")
    color_group1 = '#7fc97f'   
    color_group0 = '#beaed4'  
    # Plot distribution of 0s and 1s
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=df_current_rbp, x='Score', hue='Postar_Score', fill=True, 
                palette={1: color_group1, 0: color_group0}, common_norm=False)
    plt.axvline(x=optimal_threshold, color='red', linestyle='--', 
                label=f'Threshold = {optimal_threshold:.2f}')
    plt.title(f'Distribution of 0s and 1s for RBP: {rbp_id}')
    plt.xlabel('Scores')
    plt.ylabel('Density')
    plt.legend(title='Postar', labels=['Class-1', 'Class-0'])
    plt.grid(False)
    # Plot ROC curve
    plt.subplot(1, 2, 2)
    plt.plot(fpr, tpr, label=f'AUC = {auc_score:.2f}')
    plt.scatter(fpr[optimal_idx], tpr[optimal_idx], marker='o', color='red', 
                label=f'Threshold = {optimal_threshold:.2f}')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random')
    plt.title(f'ROC Curve for RBP: {rbp_id}')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    # Save the figure
    plt.savefig(f'{path_save}/figure_{rbp_id}.png', transparent=True)
    plt.close()