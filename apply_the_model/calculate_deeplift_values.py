
# Load libraries
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from captum.attr import DeepLift, LayerDeepLift
from collections import namedtuple
import os
import sys
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

#ruta_utils = '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/utils'
#os.environ['PYTHONPATH'] = ruta_utils + ':' + os.environ.get('PYTHONPATH', '')
from utils import Utils as utils 
from utils import Plots as plots

def calculate_batch_scores_deeplift(out_node, explainer, rbps_test, base_rbps_test, gn_test):
    """
    This function calculates scores for a batch of input data using an explainer object (DeepLIFT).
    
    Parameters:
        out_node (object): The target node or output for which scores are calculated.
        explainer (object): An explainer object used to attribute importance to input features.
        rbps_test (array-like): The input data for which scores are calculated.
        base_rbps_test (array-like): Baseline input data used for comparison in attribution.
        gn_test (array-like): Additional forward arguments used during attribution.
        
    Returns:
        scores (array-like): Importance scores attributed to input features.
    """
    scores = explainer.attribute(
                    inputs = rbps_test, 
                    baselines = base_rbps_test,
                    target=out_node,
                    additional_forward_args = gn_test)
    return scores

def calculate_rbp_reference(df_scaled_test, select_reference):  
    """
    Calculate the reference for RBP (RNA-binding protein) expression data.
    
    Parameters:
        df_scaled_test (DataFrame): The scaled test data containing RBP expression values.
        select_reference (str): The method to select the reference. Options: 'median', 'knockout'.
        
    Returns:
        base_rbps_test (Tensor): The reference RBP expression data.
    """
    if select_reference == 'median': 
        print('[calculate_rbp_reference] You selected the median of the RBP expression as reference')
        base_rbps_test = torch.tensor(np.median(df_scaled_test.astype(np.float64), axis=0)).float()
        base_rbps_test = torch.reshape(base_rbps_test, (1, base_rbps_test.size()[0]))
        print('reference of the RBP expression:', base_rbps_test)
        print('[calculate_rbp_reference] Calculate the reference of the RBP expression data: median of the test data ... -> DONE')

    if select_reference == 'knockout':
        print('[calculate_rbp_reference] You selected the knockout RBP expression as reference')
        base_rbps_test = torch.zeros(1, len(df_scaled_test.columns))
        print('reference of the RBP expression:', base_rbps_test)
        print('[calculate_rbp_reference] Calculate the reference of the RBP expression data: knockout of the test data ... -> DONE')
    return base_rbps_test

def calculate_deeplift_model_output_scores(model, test_labels, rbps_test, base_rbps_test, test_gn=None):
    """
    Calculate DeepLIFT scores for the model's output layer nodes.

    Parameters:
        model (object): The trained model.
        test_labels (array-like): The labels of the test data.
        rbps_test (array-like): The input data for which scores are calculated.
        base_rbps_test (array-like): Baseline input data used for comparison in attribution.
        test_gn (DataFrame, optional): Additional forward arguments used during attribution.

    Returns:
        list_scores_batch_wise (list): A list of DeepLIFT scores calculated for each output node.
    """
    if (test_gn != None).any().any():
        gn_test = torch.tensor(test_gn.values.astype('float64')).float()

    # Step 1: Create an instance of the Deeplift class
    print('[create_deeplift_scores_dataframe] DeepLIFT in the output layer nodes')
    explainer = DeepLift(model)
    n_out_features = test_labels.shape[1]
   
    # Suppress specific user warnings
    warnings.filterwarnings("ignore", message="Input Tensor .* did not already require gradients")
    warnings.filterwarnings("ignore", message="Setting forward, backward hooks and attributes on non-linear activations")

    # Step 2: Calculate scores batch-wise RBP x S x T
    list_scores_batch = []
    for node in tqdm(range(n_out_features), desc='Calculating scores', unit='node'):
        scores_batch = calculate_batch_scores_deeplift(node, explainer, rbps_test, base_rbps_test, gn_test)
        list_scores_batch.append(scores_batch)

    print('[create_deeplift_scores_dataframe] Calculate deeplift scores for all nodes ... -> DONE')
    return list_scores_batch

def reduce_batch_dimension(list_batch, name_trans, name_rbps, method): 
    """
    Reduce the batch dimension of DeepLIFT scores to generate final TxRBP scores using a specified method.
    
    Parameters:
        list_batch (list): A list of batched DeepLIFT scores.
        name_trans (list): List of target names.
        name_rbps (list): List of feature names (RBP names).
        method (str): The method used for reduction. Options: 'tstat', 'sum_scores'.
        
    Returns:
        df_deeplift_TxRBP (DataFrame): DataFrame containing the reduced DeepLIFT scores.
    """

    print(f'[reduce_batch_dimension] The reduce method selected is {method}')
    
    if method == 'tstat':
        df_deeplift_TxRBP = pd.DataFrame(0, index=name_trans, columns=name_rbps)
        for i, tensor in enumerate(tqdm(list_batch, desc='Calculating t-stat', unit='rbp index')):
            mean = np.mean(tensor.detach().numpy(), axis=0)
            std = np.std(tensor.detach().numpy(), axis=0)
            num_samples = tensor.size()[0] 
            with np.errstate(divide='ignore', invalid='ignore'):
                # Calculate t-statistic and assign zero to values where std is zero.
                #scores = np.where(std != 0, mean / std, 0) OLD VERSION
                scores = np.where(std != 0, mean / (std / np.sqrt(num_samples)), 0) 
            df_deeplift_TxRBP.iloc[i,:] = scores

    if method == 'sum_scores':
        data_scores_batch = torch.stack(list_batch)
        data_sum = torch.sum(data_scores_batch, dim=1).detach().numpy()
        df_deeplift_TxRBP = pd.DataFrame(data_sum, index=name_trans, columns=name_rbps)
    return df_deeplift_TxRBP

def collapse_transcript_scores_to_genes(df_score_TxRBP, getBM):
    """
    Collapse transcript-level scores to gene-level scores.

    Parameters:
        df_score_TxRBP (DataFrame): DataFrame containing scores at the transcript-RBP level.
        getBM (DataFrame): DataFrame containing transcript-to-gene mappings.

    Returns:
        df_score_GxRBP (DataFrame): DataFrame containing collapsed scores at the gene-RBP level.
    """
    getBM = getBM.sort_values(by='Transcript_ID').reset_index(drop=True)
    # From Transcript scores results collapse the scores to genes
    df_score_GxRBP = df_score_TxRBP.copy()
    # Do the absolute value of the Transcript scores (this step it's just in case you didnt perform the abs value of your transcripts before)
    df_score_GxRBP = df_score_GxRBP.abs() # new
    df_score_GxRBP.index = getBM[getBM.Transcript_ID.isin(df_score_TxRBP.index)]['Gene_ID'].tolist()
    #df_score_GxRBP = df_score_GxRBP.groupby(level=0).agg(lambda x: max(x, key=abs)) -> esto sería si necesitamos que tenga en cuenta el signo
    df_score_GxRBP = df_score_GxRBP.groupby(level=0).max()
    return df_score_GxRBP

def perform_deeplift_pipeline(df_scaled_test, test_labels, test_gn, model, path_save, getBM, select_reference='median', method='tstat'):
    """
    Perform a pipeline for DeepLIFT analysis.

    Parameters:
        df_scaled_test (DataFrame): The scaled test data containing RBP expression values.
        test_labels (DataFrame): The labels of the test data.
        test_gn (DataFrame): Additional forward arguments used during attribution.
        model (object): The trained model.
        path_save (str): The path to save the results.
        getBM (DataFrame): DataFrame containing transcript-to-gene mappings.
        select_reference (str, optional): The method to select the reference. Options: 'median', 'knockout'. Default is 'median'.
        method (str, optional): The method used for reduction. Options: 'tstat', 'sum_scores'. Default is 'tstat'.

    Returns:
        df_deeplift_scores_GxRBP (DataFrame): DataFrame containing DeepLIFT scores at the gene-RBP level.
    """
    # 1) Calculate DeepLIFT scores (RBP x T x S)
    rbps_test = torch.tensor(df_scaled_test.values.astype('float64')).float()
    rbps_names =  list(df_scaled_test.columns)
    trans_id = list(test_labels.columns)

    ### Calculate the reference of the RBP expression data
    base_rbps_test = calculate_rbp_reference(df_scaled_test, select_reference)
    print('[perform_deeplift_pipeline] Calculate the reference of the RBP expression data ... -> DONE')
    print('[perform_deeplift_pipeline] reference of the RBP expression:', base_rbps_test)
    print('[perform_deeplift_pipeline] Calculate shapley scores with DeepLIFT')
    
    scores_batch_wise = calculate_deeplift_model_output_scores(model=model, test_labels=test_labels, rbps_test=rbps_test, base_rbps_test=base_rbps_test, test_gn=test_gn)
    print('[perform_deeplift_pipeline] Calculate shapley scores with DeepLIFT ... -> DONE')

    # 2) Reduce batch dimension (RBP x T)
    df_deeplift_scores_TxRBP = reduce_batch_dimension(list_batch=scores_batch_wise, name_trans=trans_id, name_rbps=rbps_names, method=method)
    print('[perform_deeplift_pipeline] Reduce batch dimension ... -> DONE')
    # Set low-expressed genes's scores (mean < 1TPM) to 0.
    df_deeplift_scores_TxRBP.loc[test_gn.mean() < 1, :] = 0 
    print('[perform_deeplift_pipeline] Set low-expressed genes (mean < 1TPM) to 1 of the TxRBP to 0 ... -> DONE')
    # Save TxRBP scores
    utils.check_create_new_directory(path_save)
    df_deeplift_scores_TxRBP.to_csv(f'{path_save}/df_DeepLIFT_{select_reference}_{method}_TxRBPs.csv')
    print(f'[perform_deeplift_pipeline] Saved df_DeepLIFT_TxRBPs.csv in {path_save}')

    # 3) Collapse scores to genes (RBP x G)
    df_deeplift_scores_GxRBP = collapse_transcript_scores_to_genes(df_deeplift_scores_TxRBP, getBM)
    print('[get_deeplift_scores_genes] Group values by gene - take the max absolute value of the transcripts per gene ... -> DONE')

    # Save GxRBP scores
    df_deeplift_scores_GxRBP.to_csv(f'{path_save}/df_DeepLIFT_{select_reference}_{method}_GxRBPs.csv')
    print(f'[perform_deeplift_pipeline] Saved df_DeepLIFT_GxRBPs.csv in {path_save}')
    return df_deeplift_scores_TxRBP, df_deeplift_scores_GxRBP

def match_scores_and_validation_data(df_val_GxRBP, df_score_GxRBP, score_method, path_save=None): 
    """
    Match scores and validation data shapes, and optionally save the results.

    Parameters:
        df_val_GxRBP (DataFrame): Validation data DataFrame containing gene-RBP scores.
        df_score_GxRBP (DataFrame): Scores DataFrame containing gene-RBP scores.
        score_method (str): Method used for scoring.
        path_save (str, optional): Path to save the results. Default is None.

    Returns:
        df_score_GxRBP (DataFrame): Updated Scores DataFrame.
        df_val_nan_included_GxRBP (DataFrame): Modified validation data DataFrame.
    """
    ### Find the matching and non-matching genes and RBPs between the Scores dataframe and Postar
    genes_match = [x for x in df_score_GxRBP.index if x in df_val_GxRBP.index]
    genes_not_match = [x for x in df_score_GxRBP.index if x not in df_val_GxRBP.index] # contains the genes that are present in the Scores dataframe but not in the Postar dataframe.
    rbps_match = [x for x in df_score_GxRBP.columns if x in df_val_GxRBP.columns]
    rbps_not_match = [x for x in df_score_GxRBP.columns if x not in df_val_GxRBP.columns] # contains the RBPs that are present in the Scores dataframe but not in the Postar dataframe.
    print('[match_scores_and_validation_data] Find the matching and non-matching genes and RBPs between the Scores dataframe and Postar ... -> DONE')
    ### Create the df_val_nan_included_GxRBP dataframe with NaN values for the RBPs and genes where we have Scores scores and reorder rows and columns
    df_val_nan_included_GxRBP = pd.concat([df_val_GxRBP, pd.DataFrame(index=genes_not_match, columns=rbps_not_match).fillna(np.nan)], axis=1)
    df_score_GxRBP = df_score_GxRBP.loc[genes_match+genes_not_match, rbps_match+rbps_not_match]
    df_val_nan_included_GxRBP = df_val_nan_included_GxRBP.loc[genes_match+genes_not_match, rbps_match+rbps_not_match]
    print('[match_scores_and_validation_data] Create the df_val_nan_included_GxRBP dataframe with NaN values for the RBPs and genes where we have Scores scores and reorder rows and columns ... -> DONE')
    ## Reorder the rows and cols of the dataframes
    df_score_GxRBP = df_score_GxRBP.sort_index(axis=0).sort_index(axis=1)
    df_val_nan_included_GxRBP = df_val_nan_included_GxRBP.sort_index(axis=0).sort_index(axis=1)
    print('[match_scores_and_validation_data] Reorder the rows and cols of the dataframes ... -> DONE')
    ## Save the score files and return score dataframe and modified postar
    if path_save:
        print(f"[match_scores_and_validation_data] You're saving the files")
        ### Save the final dataframes
        df_score_GxRBP.to_csv(f'{path_save}/df_{score_method}_GxRBPs.csv')
        print(f'[match_scores_and_validation_data] Save the final dataframes in {path_save}\n')
    else:
        print(f"[match_scores_and_validation_data] You're NOT saving any file\n")  
    return df_score_GxRBP, df_val_nan_included_GxRBP
  
def main(df_scaled_test, test_labels, test_gn, model, path_save, df_val_GxRBP, path_data):
    getBM = pd.read_csv(path_data+'/../extra/getBM_reduced.csv', index_col=0)
    # Reorder getBM (in the last model version the data is returned with the transcript values sorted)
    getBM = getBM.sort_values(by='Transcript_ID').reset_index(drop=True)

    # General path
    path_save2 = path_save+'/DeepLIFT'
    utils.check_create_new_directory(path_save2)

    deep_studies_list = ['median', 'knockout'] # different deeplift studies based on the reference
    df_tstat_scores_GxRBP_list = []
    df_sum_scores_GxRBP_list = []
    
    for study in deep_studies_list:
        # path of the study
        path = path_save2+f'/{study}_reference'
        utils.check_create_new_directory(path)

        # 1) Obtain GxRBP score matrix
        print('[main] 1) Obtain GxRBP score matrix')
        print(path+f'/df_{study}_tstat_GxRBPs.csv')
        print(path+f'/df_{study}_sum_scores_GxRBPs.csv')
        print(os.path.isfile(path+f'/df_{study}_tstat_GxRBPs.csv'))
        print(os.path.isfile(path+f'/df_{study}_sum_scores_GxRBPs.csv'))

        if os.path.isfile(path+f'/df_{study}_tstat_GxRBPs.csv') and os.path.isfile(path+f'/df_{study}_sum_scores_GxRBPs.csv'):
            print('[get_deeplift_scores_genes] WARNING: GxRBP scores have been already calculated ... loading the saved GxRBP matrices')
            df_tstat_scores_GxRBP = pd.read_csv(path+f'/df_{study}_tstat_GxRBPs.csv', index_col=0)
            df_sum_scores_GxRBP = pd.read_csv(path+f'/df_{study}_sum_scores_GxRBPs.csv', index_col=0)
            print('[get_deeplift_scores_genes] Loading the saved GxRBP matrices ... -> DONE')

        else:
            print('[get_deeplift_scores_genes] WARNING: Calculating GxRBP scores ...')
            _, df_tstat_scores_GxRBP = perform_deeplift_pipeline(df_scaled_test, test_labels, test_gn, model, path, getBM, select_reference=study, method = 'tstat')
            _, df_sum_scores_GxRBP = perform_deeplift_pipeline(df_scaled_test, test_labels, test_gn, model, path, getBM, select_reference=study, method = 'sum_scores')
            print('[calculate_deeplift][main] Obtain df_score_GxRBP from Deeplift scores ... -> DONE')
        print('[main] 1) Obtain GxRBP score matrix ... -> DONE')

        # 2) Force the matching of the shapes of df_val_GxRBP with the DeepLIFT GxRBP
        print('[main] Force the matching of the shapes of df_val_GxRBP with the DeepLIFT GxRBP')
        df_tstat_scores_GxRBP, df_val_nan_included_GxRBP = match_scores_and_validation_data(df_val_GxRBP=df_val_GxRBP, df_score_GxRBP=df_tstat_scores_GxRBP, path_save=path, score_method=f"{study}_tstat")
        df_sum_scores_GxRBP, _ = match_scores_and_validation_data(df_val_GxRBP=df_val_GxRBP, df_score_GxRBP=df_sum_scores_GxRBP, path_save=path, score_method=f"{study}_sum_scores")
        print('[calculate_deeplift][main] Force the matching of the shapes of df_val_GxRBP with the DeepLIFT GxRBP boolean dataframe (done)')
        print(f'[calculate_deeplift][main] Matched shapes: {df_tstat_scores_GxRBP.shape}, {df_val_nan_included_GxRBP.shape}')
        utils.check_create_new_directory(path+'/Postar')
        print(f'[calculate_deeplift] Save the modified Postar dataframe in {path}/Postar')
        df_val_nan_included_GxRBP.to_csv(f'{path}/Postar/df_Postar_nan_included_GxRBPs.csv')
        print('[main] Force the matching of the shapes of df_val_GxRBP with the DeepLIFT GxRBP ... -> DONE')

        # 2.1) Analyze the matched Postar matrix for this cell line
        utils.analyze_postar_matrix(df_val_nan_included_GxRBP, path_save) 

        # 3) Plot DeepLIFT scores vs Postar
        print('[main] Plotting results')
        path_plot = path+'/plot_genes_results'
        
        folder_exist = utils.check_if_files_exist(path_plot)
        if folder_exist:
            print('[calculate_deeplift][main] The Plot results folder already exists ...')
            print('[calculate_deeplift][main] Plot DeepLIFT scores vs Postar ... -> DONE')
        else:
            print('[calculate_deeplift][main] Creating the files ...')
            utils.check_create_new_directory(path_plot)
            utils.analyze_results_per_rbp_or_gene(df_tstat_scores_GxRBP, df_val_nan_included_GxRBP, path_plot+'/tstat_scores', path_data, getBM)
            utils.analyze_results_per_rbp_or_gene(df_sum_scores_GxRBP, df_val_nan_included_GxRBP, path_plot+'/sum_scores', path_data, getBM)
            print(f'[calculate_deeplift][main] Do plot_results ... -> DONE')
            
        path = path_save2 # Reinitialize the path variable and fill the return score list of ALL genes
        df_tstat_scores_GxRBP_list.append(df_tstat_scores_GxRBP)
        df_sum_scores_GxRBP_list.append(df_sum_scores_GxRBP)
        print('[main] Plotting results ... -> DONE')

    Result = namedtuple('Result', ['df_tstat_scores_GxRBP_list', 'df_sum_scores_GxRBP_list', 'deep_studies_list', 'df_val_nan_included_GxRBP'])
    return Result(df_tstat_scores_GxRBP_list, df_sum_scores_GxRBP_list, deep_studies_list, df_val_nan_included_GxRBP)

