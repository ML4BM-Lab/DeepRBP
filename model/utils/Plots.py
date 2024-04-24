import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as cm
from collections import namedtuple
from sklearn.decomposition import PCA
import datetime
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from scipy import stats
from scipy.cluster import hierarchy
from scipy.stats import linregress, spearmanr, pearsonr, mannwhitneyu
from sklearn.metrics import auc
#from utils import Utils as utils
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator

def check_create_new_directory(path):
    """ 
    Function that checks if a directory is created and if not create the new directory
    Parameters
    ----------
    path: string with the path
    """
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print(f'[Utils] The new directory {path} is created!')

def plot_real_vs_pred(config, i, pred,labels, spear_cor, pear_cor, path, source, source_pred):
    """
    This function creates a scatter plot of the predicted values vs the real values 
    and also shows the correlation value in the title and in a legend inside a box 
    on the upper left corner of the graph. The title of the plot is also customized 
    with information such as the model used, the number of samples, the learning rate, 
    the optimizer and the number of epochs. The function also saves the plot as a png 
    file in a specific location. At the end, it displays the plot on the screen.
    -------
    """ 
    now = datetime.datetime.now()
    timestamp_str = now.strftime("%H_%M_%S")
    # Plot
    plt.figure(figsize=(10,10))
    plt.xlabel('Predicted Values', fontsize = 14)
    plt.ylabel('Real Values', fontsize = 14)
    plt.title(f'{i} {source_pred} trained with {source}',  fontsize = 16)
    plt.xticks(fontsize=14, rotation=0)
    plt.yticks(fontsize=14, rotation=0)
    sns.regplot(x=pred, y=labels, scatter_kws = {'alpha':0.1})
    plt.legend([f'spear_cor: {spear_cor:.3f}, pear_cor: {pear_cor:.3f}'],
               loc='center', borderaxespad=0.,fontsize=10)
    plt.savefig(f'{path}/real_vs_pred_{i}_in_{source_pred}_trained_with_{source}_{timestamp_str}.png')
    plt.show()
    plt.close()

def plot_expression_ratio_histogram(path_getBM, df_trans_pred, df_trans_label, df_gn, tumor_name, path_save, source, source_pred): # modified 31/08
    """
    This function performs several steps to calculate the ratio of predicted transcript expression to labeled transcript 
    expression for a specific dataset of lung cancer samples (LUAD) and then saves a histogram of the ratio values in a png 
    file. It starts by reading in three CSV files, one containing gene labels, one containing predicted transcript expression 
    values and one containing the getBM reduced values. The function then performs some data cleaning and processing, such 
    as removing duplicate columns, calculating the mean expression of each gene, and filtering out genes with a mean expression 
    less than 5. Next, the function calculates the ratio of predicted expression to labeled expression by grouping the predicted 
    transcript expression values by gene and taking the sum, then dividing that by the corresponding labeled expression values. 
    Then, the function plots a histogram of the ratio values.
    """
    now = datetime.datetime.now()
    timestamp_str = now.strftime("%H_%M_%S")
    getBM = pd.read_csv(path_getBM, index_col=0) # getBM reduced to relate transcripts with genes
    print('getBM shape:', getBM.shape)
    getBM = getBM.sort_values(by='Transcript_ID').reset_index(drop=True)
    
    # 1) Process the Gene expression dataframe (in TPM)
    # - Put the names of the genes to the gene expression dataframe 
    df_gn.columns = getBM.Gene_name
    # - Remove duplicate columns
    df_gn = df_gn.loc[:, ~df_gn.columns.duplicated()]
    # - Calculate the mean for each gene
    gn_mean_exp = df_gn.mean()
    # - Keep just the Genes with a mean expression greather than 5 TPM
    df_gn_filtered = df_gn[gn_mean_exp[gn_mean_exp > 5].index]

    # 2) Aggrupate the Predicted transcript expression by gene and transform to TPM
    # - pred in TPM
    df_trans_pred = np.power(2, df_trans_pred) - 1 
    # - aggrupate the pred transcript expression by gene
    df_trans_pred.columns = getBM.Gene_name
    df_gn_agg_trans_pred = df_trans_pred.groupby(df_trans_pred.columns, axis=1).sum() # Sum the transcript values per gene
    df_gn_agg_trans_pred = df_gn_agg_trans_pred.loc[:,df_gn_filtered.columns]

    # 3) Aggrupate the Theoretical transcript expression by gene and transform to TPM
    df_trans_label = np.power(2, df_trans_label) - 1 
    df_trans_label.columns = getBM.Gene_name
    df_gn_agg_trans_label = df_trans_label.groupby(df_trans_label.columns, axis=1).sum() # Sum the transcript values per gene
    df_gn_agg_trans_label = df_gn_agg_trans_label.loc[:,df_gn_filtered.columns]

    # 4) Calculate the sum(trans_pred)/gene_real & sum(trans_real)/gene_real ratio and plot
    ratio_pred = df_gn_agg_trans_pred/df_gn_filtered
    ratio_pred_values = ratio_pred.values.flatten()
    print('len de ratio_pred_values:', len(ratio_pred_values))
    ratio_pred_values = ratio_pred_values[~np.isnan(ratio_pred_values)]
    print('len de ratio_pred_values after removing nan:', len(ratio_pred_values))
    ratio_pred_values = ratio_pred_values[np.isfinite(ratio_pred_values)]
    print('len de ratio_pred_values after removing inf:', len(ratio_pred_values))

    ratio_label = df_gn_agg_trans_label/df_gn_filtered
    ratio_label_values = ratio_label.values.flatten()
    print('len de ratio_label_values:', len(ratio_label_values))
    ratio_label_values = ratio_label_values[~np.isnan(ratio_label_values)]
    print('len de ratio_label_values after removing nan:', len(ratio_label_values))
    ratio_label_values = ratio_label_values[np.isfinite(ratio_label_values)]
    print('len de ratio_label_values after removing inf:', len(ratio_label_values))

    plt.figure(figsize=(12, 6))
    # 4.1) Pred plot
    plt.subplot(1, 2, 1)
    plt.hist(ratio_pred_values, bins=40, linewidth=1.2)
    plt.title(f'Pred Gene ratio in genes mean > 5TPM in {source_pred}', fontsize=8)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    mean_ratio_pred = np.mean(ratio_pred_values)
    std_ratio_pred = np.std(ratio_pred_values)
    plt.legend(['Mean: ' + str(round(mean_ratio_pred, 3)) + '\nSTD: ' + str(round(std_ratio_pred, 3))])
    # 4.2) Label Plot
    plt.subplot(1, 2, 2)
    plt.hist(ratio_label_values, bins=40, linewidth=1.2)
    plt.title(f'Label Gene ratio in genes mean > 5TPM in {source_pred}', fontsize=8)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    mean_ratio_label = np.mean(ratio_label_values)
    std_ratio_label = np.std(ratio_label_values)
    plt.legend(['Mean: ' + str(round(mean_ratio_label, 3)) + '\nSTD: ' + str(round(std_ratio_label, 3))])
    plt.tight_layout()
    plt.savefig(f'{path_save}/histogram_ratio_{tumor_name}_in{source_pred}_trained_with_{source}_{timestamp_str}.png')
    plt.show()
    plt.close()

def create_custom_legend_analyze_postar(df_melted):
    # Calculate median, minimum, and maximum
    medians = df_melted.groupby('Postar')['Counts'].median()
    mins = df_melted.groupby('Postar')['Counts'].min()
    maxs = df_melted.groupby('Postar')['Counts'].max()
    # Create a custom legend
    custom_legend = [
            f"Median {group}: {int(median)}" for group, median in medians.items()
        ] + [
            f"Min {group}: {int(min_val)}" for group, min_val in mins.items()
        ] + [
            f"Max {group}: {int(max_val)}" for group, max_val in maxs.items()
        ]
    return custom_legend

def plot_analyze_postar_matrix(df_count_rbps_per_gen, df_count_genes_per_rbp, path_save):
    """
    Plot and analyze the POSTAR matrix.
    Parameters:
        df_count_rbps_per_gen (DataFrame): DataFrame containing the number of RBPs per gene.
        df_count_genes_per_rbp (DataFrame): DataFrame containing the number of genes per RBP.
        path_save (str): The path to save the plots.
    """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    df_count_rbps_per_gen_melted = df_count_rbps_per_gen.melt(id_vars=['Genes'], var_name='Postar', value_name='Counts')
    df_count_genes_per_rbp_melted = df_count_genes_per_rbp.melt(id_vars=['RBPs'], var_name='Postar', value_name='Counts')
    custom_palette = {'Class 1': '#536270', 'Class 0': '#CBAC88'}
    # Filtrar los datos para eliminar NaNs y valores 0
    df_count_genes_per_rbp_sin_nan = df_count_genes_per_rbp[(df_count_genes_per_rbp[['Class 0', 'Class 1']] > 0).any(axis=1)]
    print(df_count_genes_per_rbp_sin_nan.mean())
    df_count_rbps_per_gen_sin_nan = df_count_rbps_per_gen[(df_count_rbps_per_gen[['Class 0', 'Class 1']] > 0).any(axis=1)]
    # Subplot 1: RBP en POSTAR3 vs. genes de clase 0 y clase 1
    sns.boxplot(data=df_count_genes_per_rbp_sin_nan[['Class 1', 'Class 0']], palette=custom_palette, ax=ax_row[0], linewidth=1.5)
    ax_row[0].set_title(f'Nº of Genes regulated by RBPs in {tissue.replace("_", " ")}', fontsize=14)
    ax_row[0].set_xlabel('POSTAR3 Class', fontsize=12)
    ax_row[0].set_ylabel('Number of Genes', fontsize=12)
    ax_row[0].set_xticks([0, 1])
    ax_row[0].set_xticklabels(['Class 1', 'Class 0'])
    ax_row[0].tick_params(axis='both', which='major', labelsize=10)
    # Subplot 2: Genes vs. cantidad de genes que regulan
    sns.boxplot(data=df_count_rbps_per_gen_sin_nan[['Class 1', 'Class 0']], palette=custom_palette, ax=ax_row[1], linewidth=1.5)
    ax_row[1].set_title(f'Number of RBPs regulating Genes in {tissue.replace("_", " ")}', fontsize=14)
    ax_row[1].set_xlabel('POSTAR3 Class', fontsize=12)
    ax_row[1].set_ylabel('Number of RBPs', fontsize=12)
    ax_row[1].set_xticks([0, 1])
    ax_row[1].set_xticklabels(['Class 1', 'Class 0'])
    ax_row[1].tick_params(axis='both', which='major', labelsize=10)
    fig.tight_layout()
    plt.savefig(f'{path_save}/details_postar.png')
    plt.show()
    plt.close()

def plot_ktop_rbp_genes(df_list, actual_values_list, x_col, y_col_list, hue, title_list, path_save, plot_type='boxplot'):
    """
    Plot GxRBP scores across genes and across RBPs based on provided data.

    Parameters:
        df_list (list of DataFrame): List of DataFrames containing plot data.
        actual_values_list (list): List of actual values.
        x_col (str): Column name for the x-axis.
        y_col_list (list of str): List of column names for the y-axis.
        hue (str): Column name for hue (e.g., class labels).
        title_list (list of str): List of titles for each plot.
        path_save (str): Path to save the plots.
        plot_type (str, optional): Type of plot ('boxplot' or 'stripplot'). Defaults to 'boxplot'.
    """
    check_create_new_directory(path_save)
    num_plots = len(df_list)
    fig, axes = plt.subplots(1, num_plots, figsize=(11 * num_plots, 8), dpi=300)
    axes = np.ravel(axes)
    custom_palette = {1: '#536270',  # Elegant red shade
                      0: '#CBAC88',  # Elegant blue shade
                      -1: '#8fc38e'}  # Elegant green shade
    for i in range(num_plots):
        # Define plot params
        df = df_list[i].copy()
        df['Postar'].fillna(-1, inplace=True) 
        plot_params = {
            'data': df,
            'y': x_col,
            'x': y_col_list[i],
            "order": actual_values_list[i],
            "hue": hue,
            "hue_order": [1, 0, -1],
            "palette": custom_palette
        }
        # Define the ax & the plot
        ax = axes[i]
        if plot_type == 'boxplot':
            ax = sns.boxplot(ax=ax, **plot_params)
        else:
            ax = sns.stripplot(ax=ax, **plot_params)
        ax.grid(False)
        ax.set_ylabel(f"GxRBPs {x_col}", fontsize=11)
        ax.set_xlabel(f'{y_col_list[i].replace("_", " ")}', fontsize=11)
        ax.set_title(title_list[i], fontsize=16)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_facecolor('white')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_color('lavender')
        # Add the annotations with statannotations 
        if len(df.Postar.unique()) > 1:
            #To do the test we just need the 0s and 1s values
            df = df[df.Postar.isin([0,1])]
            filtered_values = df[df[y_col_list[i]].isin(actual_values_list[i])].groupby(y_col_list[i])['Postar'].nunique()
            actual_values_list[i] = list(filtered_values[filtered_values > 1].index)
            pairs = [
                [(value, 0), (value, 1)] 
                for value in actual_values_list[i]
                if df[df.Gene_ID == value]['Scores'].mean() != 0.0
            ]
            print(pairs)
            annotator = Annotator(ax=ax, pairs=pairs, **plot_params)
            annotator.configure(test="Mann-Whitney", comparisons_correction="bonferroni")
            _, corrected_results = annotator.apply_and_annotate()
        else:
            print('You dont have two Postar groups to perform statistical analysis') 
        # Customize legend for each subplot
        handles, labels = ax.get_legend_handles_labels()
        custom_labels = ['Class-1', 'Class-0', 'Class-NaN']
        if i==0:
            ax.legend(handles, custom_labels, title="Postar", fontsize=11, loc='upper right', bbox_to_anchor=(1.155, 1))
        else:
            ax.legend_.remove()
    plt.tight_layout()
    plt.savefig(f"{path_save}/{title_list[0].replace(' ', '_').replace('-', '_')}.png", transparent=False)
    plt.show()
    plt.close()


def plot_score_results(list_rbps_postar, list_genes_postar, df_combined1, df_combined2, path_save):
    """
    Plot the GxRBP scross across genes or across RBPs

    Parameters:
        list_rbps_postar (list): List of RBPs from POSTAR.
        list_genes_postar (list): List of genes from POSTAR.
        df_combined1 (DataFrame): Combined DataFrame with RBP scores.
        df_combined2 (DataFrame): Combined DataFrame with gene transcripts.
        path_save (str): The path to save the plots.
    """
    i=-1
    while i < len(list_rbps_postar):
        if len(list_rbps_postar) <= 5:
            actual_rbps = list_rbps_postar
            actual_genes = list_genes_postar[i+1:i+5]
        else:
            actual_rbps = list_rbps_postar[i+1:i+5]
            actual_genes = list_genes_postar[i+1:i+5]
        actual_df_combined1 = df_combined1[df_combined1.RBP_name.isin(actual_rbps)]
        actual_df_combined2 = df_combined2[df_combined2.Gene_ID.isin(actual_genes)]
        
        # Figure 1)
        if len(actual_df_combined1.Postar.unique()) > 1 and len(actual_df_combined2.Postar.unique()) > 1: # if there are different classes to compare plot:
            plot_ktop_rbp_genes(df_list=[actual_df_combined1, actual_df_combined2], 
                                actual_values_list=[actual_rbps, actual_genes], 
                                x_col="Scores", 
                                y_col_list=["RBP_name", "Gene_ID"], 
                                hue="Postar", 
                                title_list=[f"Scores by K-{i+1}-{i+5} RBP with more 1-Class Genes", f"Scores by K-{i+1}-{i+5} Genes with more 1-Class RBPs"], 
                                path_save=path_save+'/boxplot_scores_per_rbp_gene', 
                                plot_type='boxplot')

            # Figure 1.1)
            plot_ktop_rbp_genes(df_list=[actual_df_combined1], 
                                actual_values_list=[actual_rbps], 
                                x_col="Num_Trans_Per_Gene", 
                                y_col_list=["RBP_name"], 
                                hue="Postar", 
                                title_list=[f"Nº of trans per gene of Genes with Postar 0 & 1 by K-{i+1}-{i+5} RBP"], 
                                path_save=path_save+'/boxplot_num_trans_per_rbp_gene', 
                                plot_type='boxplot')
        else:
            print(f"There weren't different Postar classes to compare in {i}")
        # Next iteration
        i+=5


def plot_one_ktop_rbp_with_thresholds(df_scores, rbp_name, optimal_thresholds_df, path_save, getBM): 
    """
    Plot the scores of a specific RNA binding protein (RBP) across different classes with thresholds.

    Parameters:
        df_scores (DataFrame): DataFrame containing scores data.
        rbp_name (str): Name of the RNA binding protein to plot.
        optimal_thresholds_df (DataFrame): DataFrame containing optimal thresholds for RBPs.
        path_save (str): Path to save the plot.
        getBM (DataFrame, optional): DataFrame containing gene annotation data. Defaults to None.
    """
    num_plots = 1
    fig, axes = plt.subplots(1, num_plots, figsize=(8 * num_plots, 6))
    axes = np.ravel(axes)
    custom_palette = {1: '#536270',  # Elegant red shade
                    0: '#CBAC88',  # Elegant blue shade
                    -1: '#8fc38e'}  # Elegant green shade
    df_scores = df_scores.copy()
    df_scores = df_scores[df_scores.RBP_name == rbp_name]
    df_scores['Postar'].fillna(-1, inplace=True)
    # Replace 0 values with a small non-zero value
    df_scores['Scores'] = df_scores['Scores'] + 1
    df_scores['Scores'] = np.log2(df_scores['Scores'])
    plot_params = {
            'data': df_scores,
            'y': 'Scores',
            'x': 'RBP_name',
            'hue': 'Postar',
            'hue_order': [1, 0, -1],
            'palette': custom_palette
        }
    # Add calculated threshold
    unique_thres = optimal_thresholds_df[optimal_thresholds_df.RBP_name == rbp_name]['Optimal_Score_Threshold'].values[0]
    # Define the ax & the plot
    ax = axes[0]
    # Overlay strip plot on top of the boxplot (only for Postar == -1)
    sns.stripplot(ax=ax, **plot_params, size=8, jitter=True, dodge=True)
    # Overlay boxplot
    sns.boxplot(ax=ax, **plot_params, fliersize=5)
    # Add a horizontal line at the threshold
    ax.axhline(y=unique_thres, color='r', linestyle='--', label=f'Threshold ({unique_thres:.2f})')
    # Set labels, title, etc.
    ax.set_ylabel("GxRBPs Scores in log$_{2}$-scale", fontsize=12)
    # ax.set_title(f'Scores of {rbp_name} for each of Postar classes', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='right')
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_facecolor('white')
    ax.grid(False)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_color('lavender')
    # Set y-axis limits and manual ticks
    # ax.set_yscale('log', base=2)
    # ax.set_yticks([0, 0.01, 0.1, 1, 2, 4, 6, 8])  # Set manual ticks
    # ax.set_ylim(0, df_scores['Scores'].max())
    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    custom_labels = ['Class-1', 'Class-0', 'Class-NaN']
    ax.legend(handles, custom_labels, title="Postar", fontsize=11, loc='upper right', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    sns.despine() 
    plt.savefig(f"{path_save}/scores_of_{rbp_name}_for_postar_classes_with_thresholds.png", transparent=True)
    plt.show()
    plt.close()


def plot_distributions_and_roc(filtered_df, rbp, optimal_threshold, fpr, tpr, optimal_idx, path_save):
    """
    Plot distributions of scores and ROC curve for a specific RNA binding protein (RBP) with thresholds.

    Parameters:
        filtered_df (DataFrame): DataFrame containing filtered scores data.
        rbp (str): Name of the RNA binding protein.
        optimal_threshold (float): Optimal threshold for the RBP.
        fpr (array): Array of false positive rates.
        tpr (array): Array of true positive rates.
        optimal_idx (int): Index of the optimal threshold in fpr and tpr arrays.
        path_save (str): Path to save the plot.
    """
    path_save = path_save+'/calculate_optimal_thresholds_per_rbp'
    check_create_new_directory(path_save)
    plt.figure(figsize=(12, 4))
    sns.set(style="whitegrid")
    color_group1 = '#536270'  # Ejemplo: Naranja
    color_group0 = '#CBAC88'  # Ejemplo: Azul # '#33FF66' Verde
    # Plot distribution of 0s and 1s
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=filtered_df, x='Scores', hue='Postar', fill=True, palette={1: color_group1, 0: color_group0}, common_norm=False)
    plt.axvline(x=optimal_threshold, color='red', linestyle='--', label=f'Threshold = {optimal_threshold:.2f}')
    plt.title(f'Distribution of 0s and 1s for RBP: {rbp}')
    plt.xlabel('Scores')
    plt.ylabel('Density')
    plt.legend(title='Postar', labels=['Class-1', 'Class-0'])
    plt.grid(False)
    # Plot ROC curve
    plt.subplot(1, 2, 2)
    plt.plot(fpr, tpr, label=f'AUC = {auc(fpr, tpr):.2f}')
    plt.scatter(fpr[optimal_idx], tpr[optimal_idx], marker='o', color='red', label=f'Threshold = {optimal_threshold:.2f}')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random')
    plt.title(f'ROC Curve for RBP: {rbp}')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f'{path_save}/figure_{rbp}.png', transparent=True)
    plt.show()
    plt.close()


def plot_boxplot_with_annotations(df, significant_samples, rbp_interest, experiment, data_type, path_save):
    """
    Plot boxplots with annotations.

    Parameters:
    - df (pandas.DataFrame): DataFrame containing the data to plot.
    - significant_samples (list): List of significant samples.
    - rbp_interest (str): The RNA-binding protein of interest.
    - experiment (str): The name of the experiment.
    - data_type (str): The type of data being plotted (e.g., 'Transcripts', 'Genes').
    - path_save (str): The path where the plot will be saved.

    Returns:
    None
    """
    fig, ax = plt.subplots(figsize=(11, 10), dpi=300)
    # Identify the differential expressed transcripts
    scores_methods = {rbp_interest: df.loc[:, rbp_interest].abs()}
    df_scores = pd.DataFrame(scores_methods)
    df_scores['DE_Limma'] = 'No'
    df_scores.loc[significant_samples, 'DE_Limma'] = 'Yes'
    df_melted = pd.melt(df_scores, id_vars='DE_Limma', var_name=experiment, value_name='Score')
    # Make a copy of the original DataFrame
    df_original = df_melted.copy()
    df_melted["Score"] = df_melted["Score"]+1
    df_melted["Score"] = np.log10(df_melted["Score"])
    # Calculate and display the median of each group
    median_group_yes = df_original[df_original['DE_Limma'] == 'Yes']['Score'].median()
    median_group_no = df_original[df_original['DE_Limma'] == 'No']['Score'].median()
    median_group_yes = round(median_group_yes, 4)
    median_group_no = round(median_group_no, 4)
    palette = {'Yes': '#536270', 'No': '#CBAC88'}
    plot_params = {
        'x': experiment,
        'y': 'Score',
        'hue': 'DE_Limma',
        'palette': palette,
        'hue_order': ['Yes', 'No'] 
    }
    annotator_parameters = {
        'test': 'Mann-Whitney', 
        'comparisons_correction': 'Bonferroni',
        'text_format': 'star'
    }
    sns.stripplot(ax=ax, data=df_melted, **plot_params, jitter=True, dodge=True)
    sns.boxplot(ax=ax, data=df_melted, **plot_params)
    ax.set_title(experiment, fontsize=15)
    #ax.set_ylabel(f'{data_type} Scores in log$_{{10}}$-scale', fontsize=16)
    ax.set_xticklabels([]) 
    ax.set_ylabel(f'Scores in log$_{{10}}$-scale', fontsize=15)
    ax.set_xlabel(rbp_interest, fontsize=15)
    #ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='right', fontsize=12)  # Ajusta el tamaño de la fuente de los ticks
    #ax.tick_params(axis='both', which='major', labelsize=12)  # Ajusta el tamaño de la fuente de los números
    ax.grid(False)
    unique_methods = df_melted[experiment].unique()
    pairs = [[(method, 'No'), (method, 'Yes')] for method in unique_methods]
    annotator = Annotator(ax=ax, pairs=pairs, data=df_original, **plot_params)
    annotator.configure(**annotator_parameters)
    _, corrected_results = annotator.apply_and_annotate()
    handles, labels = ax.get_legend_handles_labels()
    custom_labels = ['DE', 'non-DE']
    ax.legend(handles, custom_labels, title=data_type, fontsize=11, loc='upper right', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    sns.despine() 
    plt.savefig(f'{path_save}/boxplot_methods_scores_{experiment}_{data_type}_{rbp_interest}_y{median_group_yes}_n{median_group_no}.png', transparent=False)
    plt.show()
    plt.close()


