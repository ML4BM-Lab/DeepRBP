import argparse
import pandas as pd
import os
import json
import torch
import sys
import uuid
import joblib
import numpy as np

from utils.Utils import get_data, scale_rbp_input_data, create_data_loader, build_optimizer, perform_predictions_on_test_data
from utils.Config import Config
from DLModelClass.modelsNN import DeepRBP


def main_predictor(args):
    path_data = args.path_data
    
    # 1. Load the config of the experiment
    config = Config(
                args.max_node, 
                args.num_hidden_layers, 
                args.learning_rate, 
                args.batch_size, 
                args.epochs, 
                args.same_num_nodes, 
                args.node_divition_factor,
                args.activation_layer, 
                args.sel_optimizer, 
                args.source_train, 
                args.train_tumor_types, 
                args.test_tumor_types)
    print('[main] The config of this experiment was loaded... -> DONE\n')
    print(f'configuration: {config.get_config()}')

    # 2. Create the results folder
    unique_id = uuid.uuid4().hex[:8]
    path_save_files = f'./output/{unique_id}'
    
    if not os.path.exists(path_save_files):
        os.makedirs(path_save_files)
    print('[main] Create the results folder ... -> DONE\n')

    # 3. Get and process the data
    # 3.1. Get the data
    data_raw = get_data(config=config, path=path_data)

    df_train = data_raw.df_train
    df_val = data_raw.df_val
    print(df_train)
    df_train = df_train.drop(['source', 'tumor_type'], axis=1)
    df_val = df_val.drop(['source','tumor_type'], axis=1)
    
    train_labels = data_raw.train_labels
    valid_labels = data_raw.valid_labels
    
    train_gn = data_raw.train_gn
    valid_gn = data_raw.valid_gn

    print(f'Number of samples for training: {df_train.shape[0]}')
    print(f'Number of samples for validation: {df_val.shape[0]}')
    print('[main] Get the data ... -> DONE\n')

    # 3.2. Scale the RBP input data
    data_scale = scale_rbp_input_data(df_train, df_val)
    print('[main] Scale the RBP input data ... -> DONE\n')

    # 3.3. Create training and validation data loaders
    train_loader = create_data_loader(scaled_rbps_df=data_scale.scaledTrain_df, labels_df=train_labels, gn_df=train_gn, config=config)
    val_loader = create_data_loader(scaled_rbps_df=data_scale.scaledValidation_df, labels_df=valid_labels, gn_df=valid_gn, config=config, set_mode='test')
    print('[main] Create training and validation data loaders ... -> DONE\n')

    # 4. Create the DeepRBP model object and get the optimizer
    model = DeepRBP(n_inputs=df_train.shape[1], n_outputs=train_labels.shape[1], config=config)
    optimizer = build_optimizer(model, config.optimizer, config.learning_rate)
    print(model)
    print('[main] Create the DeepRBP model object and get the optimizer ... -> DONE\n')

    # 5. Train the DeepRBP model
    model.fit(
                epochs=config.epochs, 
                train_loader=train_loader, 
                val_loader=val_loader, 
                optimizer=optimizer, 
                path=path_save_files, 
                model_name='model.pt')
    print('[main] Train the DeepRBP model ... -> DONE\n')

    # 6. Save the trained DeepRBP model parameters, the scaler, sigma and the configuration used in results folder
    filename_scaler = os.path.join(path_save_files, 'scaler_sfs.joblib')
    filename_sigma = os.path.join(path_save_files, 'sigma_sfs.txt')
    joblib.dump(data_scale.scaler_sfs, filename_scaler)
    np.savetxt(filename_sigma, [data_scale.sigma_sfs])

    with open(path_save_files+'/config.json', 'w') as file:
        json.dump(config.get_config(), file)
    print('[main] Save the trained DeepRBP model parameters, the scaler, sigma and the configuration used in results folder ... -> DONE\n')

    # 7. Do predictions on test data TCGA and GTEX
    data_test = get_data(config, path_data, set_mode='test')

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    results_tcga, results_gtex = perform_predictions_on_test_data(
                                    config, 
                                    data_test, 
                                    data_scale, 
                                    model, 
                                    device, 
                                    path_data, 
                                    path_save_files=path_save_files, 
                                    plot_results=plot_results)
    results_tcga.to_csv("results_tcga.csv", index=False)
    results_gtex.to_csv("results_gtex.csv", index=False)
    print('[main] Do predictions on test data TCGA and GTEX ... -> DONE\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Configuration parameters for the DeepRBP predictor.')

    parser.add_argument('--sel_optimizer', 
                        type=str, 
                        default='adamW', 
                        help='The optimizer to be used for model training. Options include "adamW", "sgd", etc.')

    parser.add_argument('--source_train', 
                        type=str, 
                        default='TCGA', 
                        help='The source of the training data. Options include "TCGA", "GTEX", etc.')

    parser.add_argument('--max_node', 
                        type=int, 
                        default=1024, 
                        help='The number of nodes in the first hidden layer of the neural network.')

    parser.add_argument('--num_hidden_layers', 
                        type=int, 
                        default=2, 
                        help='The total number of hidden layers in the neural network architecture.')

    parser.add_argument('--learning_rate', 
                        type=float, 
                        default=1e-3, 
                        help='The learning rate used by the optimizer during training.')

    parser.add_argument('--activation_layer', 
                        type=str, 
                        default='relu', 
                        help='The activation function used in the hidden layers. Options include "relu", "sigmoid", etc.')

    parser.add_argument('--node_divition_factor', 
                        type=int, 
                        default=8, 
                        help='The factor by which the number of nodes is divided from one hidden layer to the next.')

    parser.add_argument('--same_num_nodes', 
                        type=lambda x: (str(x).lower() == 'true'), 
                        default=True, 
                        help='Whether to use the same number of nodes in each hidden layer. Options are "True" or "False".')

    parser.add_argument('--batch_size', 
                        type=int, 
                        default=8, 
                        help='The number of samples per batch during training.')

    parser.add_argument('--epochs', 
                        type=int, 
                        default=1000, 
                        help='The total number of training epochs.')

    parser.add_argument('--train_tumor_types', 
                        type=str, 
                        default='all', 
                        help='A comma-separated list of tumor types to train. Use "all" for all types.')

    parser.add_argument('--test_tumor_types',  
                        type=str, 
                        default='all', 
                        help='A comma-separated list of tumor types to test. Use "all" for all types.')

    parser.add_argument('--path_data', 
                        type=str, 
                        default='/scratch/jsanchoz/data', 
                        help='The path to the data directory.')

    parser.add_argument('--plot_results', 
                        type=lambda x: x.lower() == 'true', 
                        default=False, 
                        help='Whether to plot the isoform abundance percentage and the plot real vs predicted scatter plot. Options are "True" or "False".')

    args = parser.parse_args()
    main_predictor(args)

    
      