
# explain
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
    plt.show()
    plt.close()
    