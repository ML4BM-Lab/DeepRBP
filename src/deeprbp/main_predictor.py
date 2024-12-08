
# main_predictor.py (ESTUDIAR METER ESTE SCRIPT FUERA DEL PAQUETE Y EJECUTA con ayuda del __init__ de deeprbp)

# from datetime import datetime
# from deeprbp.config_loader import load_config
# from deeprbp.data_importer import DataImporter
# from deeprbp.data_splitter import DataSplitter
# from deeprbp.utils import filter_data_by_sample_ids
import argparse
import os
import torch
from torch.utils.data import DataLoader
#import logging
from config_loader import load_config
from processing import DataImporter, DataSplitter, Scaler
from model import PredictorModel
from utils import *
from plots import scatter_real_vs_pred, plot_expression_ratio_histogram 
from evaluation_utils import *
from train_predictor import TrainPredictor 

### Decirle a ChatGPT que use buenas prácticas de eficiencia y reduccion de codigo en este main.

# Script to manage data and configuration for training

# main_predictor.py 
def main(args):
    #logging.basicConfig(level=logging.INFO) crea tu clase (LUIS)
    # Step 1: Load configuration from YAML file
    #config = load_config('/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml')
    config = load_config(args.config_path)
    
    # Step 2: Initialize DataImporter with paths from the config
    data_importer = DataImporter(config['data_paths'])
    
    # Step 3: Load data
    data = data_importer.load()
    
    # Step 4: Select sample IDs based on the provided categories in the config
    sel_sample_ids = select_sample_ids_by_type(
            metadata_df = data['metadata_df'],
            sample_category_col = config['sample_category'], 
            sample_types = config['select_samples']
    )
    # Step 5: Filter the data using the selected sample IDs
    data = filter_data_by_sample_ids(data, sel_sample_ids)
    
    # Step 6: Data splitting based on the configuration settings
    splitter = DataSplitter(data, config)
    train_data, valid_data, test_data = splitter.split_data_sets()
    
    # Step 7: Scaler. No existing scaler and sigma provided, fit the scaler externally
    scaler = Scaler()  # Initialize Scaler with no existing scaler/sigma
    # Fit the scaler using training data
    scaler.fit(train_data['rbp_expr_df'])
    
    # Now use the fitted scaler to transform train, validation, and test sets
    train_data['scaled_rbp_expr_df'] = scaler.transform(train_data['rbp_expr_df'])
    valid_data['scaled_rbp_expr_df'] = scaler.transform(valid_data['rbp_expr_df'])
    test_data['scaled_rbp_expr_df'] = scaler.transform(test_data['rbp_expr_df'])
    
    # Step 8: Save the scaler for later use (optional) and processed datasets
    path_save_data = os.path.join(config['output_dir'], 'data')
    scaler_path = os.path.join(path_save_data, 'scaler_trained')
    scaler.save(scaler_path)
    
    save_processed_data(train_data, os.path.join(path_save_data, 'train_data'))
    save_processed_data(valid_data, os.path.join(path_save_data, 'valid_data'))
    save_processed_data(test_data, os.path.join(path_save_data, 'test_data'))
    #loaded_train_data = load_processed_data(path)

    # Step 9: Create TensorDataset and DataLoaders
    train_dataset = CustomTensorDataset(train_data)  
    valid_dataset = CustomTensorDataset(valid_data) 
    test_dataset = CustomTensorDataset(test_data)   

    train_loader = DataLoader(
            train_dataset, 
            adjust_batch_size(train_dataset, config['training']['batch_size']), 
            drop_last=True, 
            shuffle=True
    )
    
    val_loader = DataLoader(
            valid_dataset, 
            adjust_batch_size(valid_dataset, config['training']['batch_size']*2), 
            shuffle=False
    )  

    test_loader = DataLoader(
         test_dataset, 
         adjust_batch_size(test_dataset, config['training']['batch_size']*2), 
         shuffle=False
    )

    # Step 10: Initialize Model
    model = PredictorModel(
        input_size=len(next(iter(train_dataset))['rbp_expr']), 
        output_size=len(next(iter(train_dataset))['trans_expr']), 
        config=config
    )

    # Step 11: Initialize TrainPredictor
    trainer = TrainPredictor(model=model, config=config)
    # Step 12: Train the model
    train_history, val_history = trainer.fit(train_loader, val_loader, epochs=config['epochs'])
    # Step 13: Save the model and history
    path_save_results = os.path.join(config['output_dir'], 'results')
    trainer.save_model_and_training_history(path_save_results, 'model.pt', train_history, val_history)

    # Step 14 make Predictions and calculate metrics on Training, Validation and Test (general)
    pred_train, label_train, _ = trainer.generate_predictions(train_loader)
    pred_val, label_val, _ = trainer.generate_predictions(val_loader)
    pred_test, label_test, _ = trainer.generate_predictions(test_loader)

    metrics_train = calculate_metrics(pred_train, label_train)
    metrics_val = calculate_metrics(pred_val, label_val)
    metrics_test = calculate_metrics(pred_test, label_test)

    save_metrics_summary(set_names_list = ['train', 'val', 'test'], 
                        metrics_list = [metrics_train, metrics_val, metrics_test], 
                        output_path = f'{path_save_results}/metrics_summary_global.csv')
                   
    # Step 15 make Predictions and calculate metrics on Test (per category)
    # TRAIN
    metrics_list, set_names_list = calculate_metrics_per_category(train_data, 
                                                                  config, 
                                                                  trainer,
                                                                  output_dir=path_save_results,
                                                                  set_name='train')
    save_metrics_summary(set_names_list, metrics_list, 
                         output_path = f'{path_save_results}/train_metrics_summary_per_category.csv')
    # VAL
    metrics_list, set_names_list = calculate_metrics_per_category(valid_data, 
                                                                  config, 
                                                                  trainer,
                                                                  output_dir=path_save_results,
                                                                  set_name='val')
    save_metrics_summary(set_names_list, metrics_list, 
                         output_path = f'{path_save_results}/val_metrics_summary_per_category.csv')

    # TEST
    metrics_list, set_names_list = calculate_metrics_per_category(test_data, 
                                                                  config, 
                                                                  trainer,
                                                                  output_dir=path_save_results,
                                                                  set_name='test')
    save_metrics_summary(set_names_list, metrics_list, 
                         output_path = f'{path_save_results}/test_metrics_summary_per_category.csv')
    ###
    
    # Step 16: Predictions on GTEX - ME QUEDA REVISAR EL TRAIN Y EL VAL EN EL ANTERIOR CODIGO Y HACER EN GTEX Y
   # PEDIR CONSEJO A CHAT PREMIUM DE COMO HE ORGANIZADO TODO EL CODE EN LA PIPELINE A PARTIR DE TODO EL CODE. 

    ###  ###  ###  ###  ###  ###  ###  HASTA AQUI YA REVISADO.

    ### ... WORKING IN PROGRESS ... 





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    args = parser.parse_args()
    main(args)






### basura: 
##################################################################
##################################################################
##################################################################
########################################################################################################

path_configs = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs'
config_path_tcga_train = os.path.join(path_configs, 'config_tcga_train.yaml')
config_path_gtex = os.path.join(path_configs, 'config_gtex.yaml')
config_path_tcga_test = os.path.join(path_configs, 'config_tcga_test.yaml')

# Data
paths_TCGA = {
    "rbp_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/RBPs_log2p_tpm.csv",
    "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/trans_log2p_tpm.csv",
    "metadata_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/phenotype_metadata.csv",
    "gene_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/gn_expr_each_iso_tpm.csv"
    }

paths_TCGA_rep = {
    "rbp_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/rbp_expr_log2p_tpm.csv",
    "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/trans_expr_log2p_tpm.csv",
    "metadata_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/metadata_df.csv",
    "gene_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/gene_expr_tpm.csv"
    }

paths_GTEX = {
    "rbp_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/RBPs_log2p_tpm.csv",
    "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/trans_log2p_tpm.csv",
    "metadata_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/phenotype_metadata.csv",
    "gene_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/gn_expr_each_iso_tpm.csv"
    }

# Caso 1) config train TCGA
config = load_config(config_path_tcga_train)
output_dir = config['output_dir']
batch_size = config['batch_size']

dataset_1 = CustomDataset( 
                    paths=paths_TCGA, 
                    config=config, 
                    output_dir=output_dir,
                    save_files=True
                    )
                    
# Caso 2) config uso repetido de TCGA
config = load_config(config_path_tcga_test)
path_saved_files = '/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA'
scaler, sigma, train_idx, valid_idx, test_idx = CustomDataset.load_scaler_and_idx(path_saved_files)

dataset_2 = CustomDataset(
        paths=paths_TCGA_rep, 
        config=config, 
        train_idx=train_idx,
        valid_idx=valid_idx,
        test_idx=test_idx,
        scaler=scaler,
        sigma=sigma,
        save_files=False
        )

# config GTEX
config = load_config(config_path_gtex)
path_saved_files = '/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-06/train_prediction_model/data/TCGA'
scaler, sigma, _, _, _ = CustomDataset.load_scaler_and_idx(path_saved_files)

dataset = CustomDataset(
        paths=paths_GTEX, 
        config=config, 
        scaler=scaler,
        sigma=sigma,
        save_files=True
        )

####





dataset_2.rbp_names
dataset_2.trans_names
