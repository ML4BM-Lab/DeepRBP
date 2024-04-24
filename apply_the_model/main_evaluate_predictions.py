import os
import sys
import argparse
import json
import joblib
import numpy as np
import torch
from collections import namedtuple

from utils.Utils import get_data, perform_predictions_on_test_data
from utils.Plots import plot_boxplot_with_annotations
from utils.Config import Config
from DLModelClass.modelsNN import DeepRBP

def main_predictor(args):
    path_data = args.path_data
    
    ## 1. Get and process the test data¶
    ## 1.1. Load the config used in training DeepRBP
    with open(f'{args.path_train_files}/config.json', 'r') as file:
        config_dict = json.load(file)
    config_dict['test_tumor_types'] = 'Liver_Hepatocellular_Carcinoma'
    config = Config(**config_dict)
    print(config.get_config())

    ## 1.2 Get test raw data
    data_test = get_data(config, path_data, set_mode='test')

    ## 1.3. Get the scaler used in training
    scaler_sfs = joblib.load(args.path_train_files+'/scaler_sfs.joblib')
    with open(args.path_train_files+'/sigma_sfs.txt', 'r') as f:
        sigma_sfs = f.readline().strip()
    sigma_sfs = np.float128(sigma_sfs)
    Data_Scale = namedtuple('Data_Scale', ['scaler_sfs', 'sigma_sfs'])
    data_scale = Data_Scale(scaler_sfs, sigma_sfs)

    ## 2. Load the trained model
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    model = DeepRBP(config=config)
    model.load_state_dict(torch.load(f'{args.path_train_files}/model.pt', map_location=device))      
    model.eval()

    ## 3. Make predictions in TCGA and GTeX
    results_tcga, results_gtex = perform_predictions_on_test_data(
                                    config, 
                                    data_test, 
                                    data_scale, 
                                    model, 
                                    device, 
                                    path_data, 
                                    path_save=args.path_save, 
                                    plot_results=args.plot_results)
    results_tcga.to_csv("results_tcga.csv", index=False)
    results_gtex.to_csv("results_gtex.csv", index=False)
    print('[main] Do predictions on test data TCGA and GTEX ... -> DONE\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Configuration parameters for the DeepRBP predictor.')
    parser.add_argument('--path_data', 
                        type=str, 
                        help='Path to the test data.')
    parser.add_argument('--path_train_files', 
                        type=str, 
                        help='Path to the directory containing training files.')
    parser.add_argument('--path_save', 
                        type=str, 
                        help='Path to save the prediction results.')
    parser.add_argument('--plot_results', 
                        type=lambda x: x.lower() == 'true', 
                        default=False, 
                        help='Whether to plot the isoform abundance percentage and the plot real vs predicted scatter plot. Options are "True" or "False".')
    args = parser.parse_args()
    main_predictor(args)
    