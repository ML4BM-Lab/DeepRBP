
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os


############################ Plot matriz correlaciones ########################

def plot_upper_triangle_corr_matrix(corr_matrix, complex_idx=None, output_dir=None, subset_size=100): # muy mejorable
    """
    Plot a simple heatmap showing only the upper triangle (including diagonal) of a subset
    of a large correlation matrix, without labels or clustering.

    Parameters:
    - corr_matrix: pd.DataFrame, square correlation matrix.
    - complex_idx: str or int or None (default None)
        Identifier of the protein complex to show in the plot title.
    - output_dir: str or None, path to save the figure (if None, no save)
    - subset_size: int, number of genes to keep for plotting (default 100)
    """
    # Subset matrix
    corr_small = corr_matrix.iloc[:subset_size, :subset_size].copy()
    # Create mask for lower triangle
    mask = np.tril(np.ones_like(corr_small, dtype=bool), k=-1)
    plt.figure(figsize=(8,8))
    sns.set_theme(style="white")
    ax = sns.heatmap(
        corr_small,
        mask=mask,
        cmap="RdYlBu_r",
        square=True,
        cbar_kws={"label": "Correlation"},
        xticklabels=False,
        yticklabels=False,
        linewidths=0,
        vmin=0, vmax=1
    )
    plt.title("Upper Triangle Correlation Heatmap (subset)", fontsize=14, fontweight='bold')
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, f"upper_triangle_heatmap_{complex_idx}.png")
        plt.savefig(path, dpi=300, bbox_inches='tight')
        plt.close()

def plot_correlation_boxplot(flat_complex, flat_outside, complex_idx=None, output_dir=None):
    """
    Plot a boxplot comparing correlation distributions between two groups:
    'Complex' and 'Outside', including complex ID info in the title.

    Parameters:
    - flat_complex: array-like
        Flattened correlation values within the protein complex.
    - flat_outside: array-like
        Flattened correlation values outside the protein complex.
    - complex_idx: str or int or None (default None)
        Identifier of the protein complex to show in the plot title.
    - output_dir: str or None (default None)
        Directory path to save the plot image. If None, the plot is not saved.

    The function creates and displays a boxplot with proper styling.
    """
    data = pd.DataFrame({
        'Correlation': list(flat_complex) + list(flat_outside),
        'Group': ['Complex'] * len(flat_complex) + ['Outside'] * len(flat_outside)
    })
    plt.figure(figsize=(8,6))
    sns.set_theme(style="whitegrid")
    ax = sns.boxplot(x='Group', y='Correlation', data=data, hue='Group',
                     palette=['#4c72b0', '#55a868'], width=0.6, legend=False)
    title = 'Correlation Distribution: Complex vs Outside'
    plt.title(title, fontsize=16, fontweight='bold')
    plt.ylabel('Correlation values', fontsize=14)
    plt.xlabel('')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        file_path = os.path.join(output_dir, f'correlation_boxplot_complex_{complex_idx}.png')
        plt.savefig(file_path)
        plt.close()