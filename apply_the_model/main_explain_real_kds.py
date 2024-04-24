import argparse
import pandas as pd
import os
import json
import torch
from utils.Utils import check_data_exists, get_data, do_predictions, generate_data, scale_inputs_real_knockout_explainability, evaluate_model, get_significant_samples
from utils.Plots import plot_boxplot_with_annotations
from utils.Config import Config
from DLModelClass.modelsNN import DeepRBP
from calculate_deeplift_values import perform_deeplift_pipeline

def main(path_exp, experiment, rbp_interest, path_train_files, path_data, path_save):  
    # 1. Read and process raw samples to generate the input data for DeepRBP
    
    # Check if you have already processed the raw data
    generate_data_flag = check_data_exists(path_exp, experiment)
    # Define the two experiment conditions
    condition1 = experiment+'_control'
    condition2 = f'{experiment}_{rbp_interest}_kd'
    print(f'condition 1: {condition1}')
    print(f'condition 2: {condition2}')
    description = pd.read_csv(path_exp+f'/info_samples.txt', delimiter='\t')
    print(f'Samples description: {description}')

    # Generate data from raw abundance samples
    if generate_data_flag:
        print('[main] Generate data from raw abundance samples ...')
        generate_data(description, condition1, condition2, path_exp, path_data)
    print('[main] Generate data from raw abundance samples ... -> DONE')

    # 2. Scale control and knockout data
    data_control = scale_inputs_real_knockout_explainability(path_train_files, f'{path_exp}/datasets', condition=condition1)
    data_kd = scale_inputs_real_knockout_explainability(path_train_files, f'{path_exp}/datasets', condition=condition2)
    
    df_rbps_control = data_control.df_scaled_rbps
    df_labels_control = data_control.df_labels
    df_gns_control = data_control.df_gns
    df_rbps_kd = data_kd.df_scaled_rbps
    df_labels_kd = data_kd.df_labels
    df_gns_kd = data_kd.df_gns
    print('[main] Scale control and knockout data ... -> DONE')

    # 3. Evaluate DeepRBP predictor on control and knockout data
    # Load the configuration and the trained model
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    with open(f'{path_train_files}/config.json', 'r') as file:
        config_dict = json.load(file)
    config_obj = Config(**config_dict)
    model = DeepRBP(n_inputs=df_rbps_control.shape[1], n_outputs=df_labels_control.shape[1], config=config_obj, device=device)
    model.load_state_dict(torch.load(path_train_files+'/model.pt', map_location=device))      
    model.eval()
    print('[main] Evaluate DeepRBP predictor on control and knockout data ... -> DONE')

    # Make predictions and calculate pearson correlation, spearman's, and MSE.
    evaluate_model(model, df_rbps_control, df_gns_control, df_labels_control, df_rbps_kd, df_gns_kd, df_labels_kd, device, condition1, condition2)
    print("[main] Make predictions and calculate pearson correlation, spearman's, and MSE ... -> DONE")
    
    # DeepRBP explainability module: calcualte the RBPxTranscripts and RBPxGenes matrices.
    df_deeplift_scores_TxRBP, df_deeplift_scores_GxRBP = perform_deeplift_pipeline(
                                                        df_scaled_test = df_rbps_control, 
                                                        test_labels = df_labels_control, 
                                                        test_gn = df_gns_control, 
                                                        model = model, 
                                                        path_save = path_save, 
                                                        getBM = getBM, 
                                                        select_reference='knockout', 
                                                        method='tstat'
                                                        )
    print("[main] DeepRBP explainability module: calcualte the RBPxTranscripts and RBPxGenes matrices ... -> DONE")
    
    # Plots using the scores agains DE-transcripts and DE-genes
    path_save = path_save+f'/{experiment}'
    os.makedirs(path_save, exist_ok=True)
    significant_samples, significant_samples_gns = get_significant_samples(rbp_interest, experiment, path_exp, df_deeplift_scores_TxRBP)   
    plot_boxplot_with_annotations(df=df_deeplift_scores_TxRBP, significant_samples=significant_samples, rbp_interest=rbp_interest, experiment=experiment, path_save=path_save, data_type='Transcripts')
    plot_boxplot_with_annotations(df=df_deeplift_scores_GxRBP, significant_samples=significant_samples_gns, rbp_interest=rbp_interest, experiment=experiment, path_save=path_save, data_type='Genes')
    print("[main] Plots using the scores agains DE-transcripts and DE-genes ... -> DONE")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DeepRBP explain in real kds parameters')
    parser.add_argument('--path_exp', type=str, help='Path to the experiment directory', required=True)
    parser.add_argument('--experiment', type=str, help='Name of the experiment', required=True)
    parser.add_argument('--rbp_interest', type=str, help='Name of the RBP of interest', required=True)
    parser.add_argument('--path_train_files', type=str, help='Path to the directory containing training files', required=True)
    parser.add_argument('--path_data', type=str, help='Path to the data directory', required=True)
    parser.add_argument('--path_save', type=str, help='Path to save results', required=True)
    args = parser.parse_args()
    
    getBM = pd.read_csv(f'{args.path_data}/extra/getBM_reduced.csv', index_col=0)
    getBM = getBM.sort_values(by='Transcript_ID').reset_index(drop=True)
    
    main(args.path_exp, args.experiment, args.rbp_interest, args.path_train_files, args.path_data, args.path_save)
        
 