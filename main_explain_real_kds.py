import argparse
import pandas as pd
import os
import torch
from utils.Utils import load_model, generate_data, evaluate_model, process_data, check_data_exists, get_significant_samples
from utils.Plots import plot_boxplot_with_annotations
from calculate_deeplift_values import perform_deeplift_pipeline

def main(path_data, experiment, rbp_interest, getBM, path_model, path_result):  
    print('[main] 1) Create Data \n')
    print('[main] 1.1) Check if Data is already created ... \n')
    generate_data_flag = check_data_exists(path_data, experiment)
    print('[main] 1.1) Check if Data is already created ... -> DONE')
    
    condition1 = experiment+'_control'
    condition2 = f'{experiment}_{rbp_interest}_kd'
    df_filereport = pd.read_csv(path_data+f'/info_samples.txt', delimiter='\t')
    
    if generate_data_flag:
        generate_data(df_filereport, condition1, condition2, path_data, path_model)
    
    print('[main] 2) Get data processed (real control and kd data) ... -> DONE \n')
    data_control, data_kd = process_data(path_model, f'{path_data}/datasets', condition1, condition2)
    df_rbps_control = data_control.df_scaled_rbps
    df_labels_control = data_control.df_labels
    df_gns_control = data_control.df_gns
    df_rbps_kd = data_kd.df_scaled_rbps
    df_labels_kd = data_kd.df_labels
    df_gns_kd = data_kd.df_gns
    
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    print('[main] The device selected is:', device)
    model = load_model(path_model,df_rbps_control, df_labels_control, device)
    print('[main] 3) Load the Trained Model ... -> DONE \n')
    print('[main] 3.1) Evaluate the model prediction ...')
    evaluate_model(model, df_rbps_control, df_gns_control, df_labels_control, df_rbps_kd, df_gns_kd, df_labels_kd, device, condition1, condition2)
    print('[main] 3.1) Evaluate the model prediction in (real control and kd data)... -> DONE')
  
    print('[main] 4) Perform explainability from control data with DeepLIFT ...')
    path_save = path_result+f'/explainability_cell_line_kds/{path_model.rsplit("/", 1)[-1]}/{experiment}/DeepLIFT/knockout_reference'
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
    print('[main] 4) Perform explainability from control data ... -> DONE')

    print('[main] 5) Identify which are the transcripts and genes that are differentially expressed in condition 1 and condition 2 after using Limma ...')
    significant_samples, significant_samples_gns = get_significant_samples(rbp_interest, experiment, path_data, df_deeplift_scores_TxRBP)
    print('[main] 5) Identify which are the transcripts and genes that are differentially expressed in condition 1 and condition 2 after using Limma ... -> DONE')
    print('[main] 6) Plot results')
    plot_boxplot_with_annotations(df=df_deeplift_scores_TxRBP, significant_samples=significant_samples, rbp_interest=rbp_interest, experiment=experiment, path_save=path_save, data_type='Transcripts')
    plot_boxplot_with_annotations(df=df_deeplift_scores_GxRBP, significant_samples=significant_samples_gns, rbp_interest=rbp_interest, experiment=experiment, path_save=path_save, data_type='Genes')
    print('[main] 6) Plot results ... -> DONE')
  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DeepRBP explain in real kds parameters')
    parser.add_argument('--path_data', type=str, help='Path to data')
    parser.add_argument('--experiment', type=str, help='Name of the experiment')
    parser.add_argument('--rbp_interest', type=str, help='Name of the RBP of interest')
    parser.add_argument('--path_model', type=str, help='Path to model')
    parser.add_argument('--path_result', type=str, help='Path to result')
    args = parser.parse_args()

    # Load the getBM set
    getBM = pd.read_csv(f"{args.path_model.rsplit('/', 2)[0]}/processed/getBM_reduced.csv", index_col=0)
    # Reorder getBM (in the last model version the data is returned with the transcript values sorted)
    getBM = getBM.sort_values(by='Transcript_ID').reset_index(drop=True)
    # Call main
    main(path_data=args.path_data, experiment=args.experiment, rbp_interest=args.rbp_interest, getBM=getBM, path_model=args.path_model, path_result=args.path_result)

 