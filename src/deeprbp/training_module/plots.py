# src/deeprbp/training_module/plots.py

import os 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Dict, List, Union
#from sklearn.metrics import auc
from ..util.utils import ensure_directory_exists

def scatter_real_vs_pred(
    category: str,
    source_name: str,
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
    source_name (str): Name of the data source, used for saving results.
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
    plt.figure(figsize=(12, 12))
    plt.xlabel('Predicted Values', fontsize=16, fontweight='bold')
    plt.ylabel('Real Values', fontsize=16, fontweight='bold')
    plt.title(f'Real vs Predicted: {category} ({source_name})', fontsize=18, fontweight='bold')
    sns.regplot(
        x=pred, y=labels,
        scatter_kws={'alpha': 0.3, 'color': 'blue'},
        line_kws={'color': 'red', 'lw': 2}
    )
    legend_text = (
        f"Spearman Corr: {metrics['spearman_corr']:.2f}\n"
        f"Pearson Corr: {metrics['pearson_corr']:.2f}\n"
        f"MSE: {metrics['mse']:.2f}\n"
        f"R²: {metrics['r2']:.2f}"
    )
    plt.text(
        0.05, 0.95, legend_text, transform=plt.gca().transAxes,
        fontsize=14, verticalalignment='top',
        bbox=dict(boxstyle="round", edgecolor="black", facecolor="white")
    )
    sns.set_style("whitegrid")
    ensure_directory_exists(output_dir)

    output_path = os.path.join(output_dir, f"{category}-{source_name}.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved at: {output_path}")

def plot_loss_curve(train_history: List[float], val_history: List[float], 
                    title: str = 'Training and Validation Loss', 
                    output_dir: str = None) -> None:
    """
    Plots the training and validation loss curves to visualize the performance of the model over time.

    Parameters:
    - train_history (List[float]): List of training loss values for each epoch.
    - val_history (List[float]): List of validation loss values for each epoch.
    - title (str): Title of the plot (default is 'Training and Validation Loss').
    - output_dir (str): Path to save the plot as a .png file.

    Returns:
    - None: The function will save the plot.
    """
    ensure_directory_exists(output_dir)
    plt.style.use('seaborn-v0_8-muted')
    plt.figure(figsize=(12, 8))
    # Plot both training and validation loss
    plt.plot(train_history, label='Training Loss', color='royalblue', linestyle='-', linewidth=2.5)
    plt.plot(val_history, label='Validation Loss', color='darkorange', linestyle='--', linewidth=2.5)
    plt.title(title, fontsize=20, fontweight='bold', pad=15)
    plt.xlabel('Epoch', fontsize=16, labelpad=10)
    plt.ylabel('Loss', fontsize=16, labelpad=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=14, loc='upper right', frameon=True, shadow=True, fancybox=True)
    plt.tight_layout(pad=2)
    plt.savefig(os.path.join(output_dir, 'loss_curve.png'), dpi=300)
    plt.close()
    print(f"Plot saved to {output_dir}")

def plot_transcript_to_gene_ratio_distributions(
    ratios_pred: np.ndarray,
    ratios_label: np.ndarray,
    source_name: str,
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
    plt.title(f"Predicted Gene Ratio (Mean > 5 TPM) - {source_name}", fontsize=12, fontweight='bold')
    plt.xlabel("Transcript-to-Gene Ratio", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.axvline(mean_pred, color=color_pred, linestyle='dashed', linewidth=1.5)
    plt.text(mean_pred + 0.1, plt.ylim()[1] * 0.9, f'Mean: {mean_pred:.3f}\nSTD: {std_pred:.3f}', color=color_pred, fontsize=9)
    # Labeled ratio histogram
    plt.subplot(1, 2, 2)
    plt.hist(ratios_label, bins=40, color=color_label, edgecolor='black', linewidth=1.2)
    mean_label = np.mean(ratios_label)
    std_label = np.std(ratios_label)
    plt.title(f"Labeled Gene Ratio (Mean > 5 TPM) - {source_name}", fontsize=12, fontweight='bold')
    plt.xlabel("Transcript-to-Gene Ratio", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.axvline(mean_label, color=color_label, linestyle='dashed', linewidth=1.5)
    plt.text(mean_label + 0.1, plt.ylim()[1] * 0.9, f'Mean: {mean_label:.3f}\nSTD: {std_label:.3f}', color=color_label, fontsize=9)
    # Adjust layout
    plt.tight_layout()
    ensure_directory_exists(os.path.dirname(output_path))
    plt.savefig(output_path, dpi=300)
    plt.close()