import os
import numpy as np
import load_config
from data_loader import CustomDataset, create_dataloaders
from deeprbp.module_training.models import TranscriptExpressionPredictor
#import DeepRBP>>> 
from deeprbp.data_loading.config_loader import load_config
from data_loader import CustomDataset, create_dataloaders
from deeprbp.module_training.models import TranscriptExpressionPredictor

paths_TCGA_rep = {
     "rbp_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/rbp_expr_log2p_tpm.csv",
    "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/trans_expr_log2p_tpm.csv",
     "metadata_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/metadata_df.csv",
     "gene_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA/pre-scaling/gene_expr_tpm.csv"
    }

path_configs = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs'
config_path_tcga_train = os.path.join(path_configs, 'config_tcga_train.yaml')
config_path_tcga_test = os.path.join(path_configs, 'config_tcga_test.yaml')
config = load_config(config_path_tcga_train)

config
{'source_name': 'TCGA', 'select_samples': ['Lung_Adenocarcinoma', 'Breast_Invasive_Carcinoma'], 'sample_category': 'detailed_category', 'train_test_split': True, 'test_frac': 0.2, 'train_val_split': True, 'val_frac': 0.15, 'fit_scaler': True, 'optimizer_name': 'adamW', 'max_node': 1024, 'num_hidden_layers': 2, 'node_shrink_factor': 8, 'uniform_nodes': True, 'activation_func': 'relu', 'batch_size': 128, 'learning_rate': '1e-3', 'epochs': 1000, 'cuda': True, 'seed': 0, 'plot_results': False, 'output_dir': '/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model'}
config = load_config(config_path_tcga_test)
config
{'source_name': 'TCGA', 'train_test_split': False, 'train_val_split': False, 'seed': 0, 'batch_size': 128, 'fit_scaler': False, 'cuda': True, 'epochs': 1000, 'select_samples': ['Lung_Adenocarcinoma', 'Breast_Invasive_Carcinoma'], 'sample_category': 'detailed_category', 'output_dir': '/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model'}

config = load_config(config_path_tcga_test)
path_saved_files = '/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-12/train_prediction_model/data/TCGA'
scaler, sigma, train_idx, valid_idx, test_idx = CustomDataset.load_scaler_and_idx(path_saved_files)

dataset = CustomDataset(
         paths=paths_TCGA_rep, 
        config=config, 
        train_idx=train_idx,
         valid_idx=valid_idx,
        test_idx=test_idx,
         scaler=scaler,
         sigma=sigma,
         save_files=False
         )

input_size = len(dataset.rbp_names)
output_size = len(dataset.trans_names)
dataloaders = create_dataloaders(
             dataset=dataset, 
            batch_size=batch_size, 
            )

output_dir = config['output_dir']
batch_size = config['batch_size']
epochs = config['epochs']
dataloaders = create_dataloaders(
             dataset=dataset, 
           batch_size=batch_size, 
            )

train_loader = dataloaders['train_loader']
val_loader = dataloaders['validation_loader']
test_loader = dataloaders['test_loader'] 

config_path_tcga_train = os.path.join(path_configs, 'config_tcga_train.yaml')
config = load_config(config_path_tcga_train)
model = TranscriptExpressionPredictor(input_size=input_size, output_size=output_size, config=config)

model.train_model(
                 epochs=epochs, 
                 train_loader=train_loader, 
                 val_loader=val_loader, 
                path=output_dir, 
                 model_name='model.pt'
               )