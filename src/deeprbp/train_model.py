from data_loader import CustomDataset
from config_loader import load_config
import os
import numpy as np

## version carlosizada
#import DeepRBP

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
    "rbp_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-04/train_prediction_model/data/TCGA/pre-scaling/rbp_expr_log2p_tpm.csv",
    "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-04/train_prediction_model/data/TCGA/pre-scaling/trans_expr_log2p_tpm.csv",
    "metadata_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-04/train_prediction_model/data/TCGA/pre-scaling/metadata_df.csv",
    "gene_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-04/train_prediction_model/data/TCGA/pre-scaling/gene_expr_tpm.csv"
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
path_saved_files = '/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-11-06/train_prediction_model/data/TCGA'
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







# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# comprobaciones que he hecho en este codigo: 
def compare_datasets_with_tolerance(dataset_1, dataset_2):
    dataset_names = [
        ('train_rbps', dataset_1.rbp_expr[dataset_1.train_idx, :], dataset_2.rbp_expr[dataset_2.train_idx, :]),
        ('train_gn', dataset_1.gene_expr[dataset_1.train_idx, :], dataset_2.gene_expr[dataset_2.train_idx, :]),
        ('train_trans', dataset_1.trans_expr[dataset_1.train_idx, :], dataset_2.trans_expr[dataset_2.train_idx, :]),
        ('val_rbps', dataset_1.rbp_expr[dataset_1.valid_idx, :], dataset_2.rbp_expr[dataset_2.valid_idx, :]),
        ('val_gn', dataset_1.gene_expr[dataset_1.valid_idx, :], dataset_2.gene_expr[dataset_2.valid_idx, :]),
        ('val_trans', dataset_1.trans_expr[dataset_1.valid_idx, :], dataset_2.trans_expr[dataset_2.valid_idx, :]),
        ('test_rbps', dataset_1.rbp_expr[dataset_1.test_idx, :], dataset_2.rbp_expr[dataset_2.test_idx, :]),
        ('test_gn', dataset_1.gene_expr[dataset_1.test_idx, :], dataset_2.gene_expr[dataset_2.test_idx, :]),
        ('test_trans', dataset_1.trans_expr[dataset_1.test_idx, :], dataset_2.trans_expr[dataset_2.test_idx, :])
    ]
    for name, data1, data2 in dataset_names:
        if np.allclose(data1, data2, atol=1e-6):  # Usamos una tolerancia de 1e-6
            print(f"Los datasets {name} son iguales (con tolerancia).")
        else:
            print(f"Los datasets {name} NO son iguales (con tolerancia).")

# Llamar a la función para comparar los datasets con tolerancia
compare_datasets_with_tolerance(dataset_1, dataset_2)

def compare_and_show_differences(dataset_1, dataset_2):
    dataset_names = [
        ('train_rbps', dataset_1.rbp_expr[dataset_1.train_idx, :], dataset_2.rbp_expr[dataset_2.train_idx, :]),
        ('train_gn', dataset_1.gene_expr[dataset_1.train_idx, :], dataset_2.gene_expr[dataset_2.train_idx, :]),
        ('train_trans', dataset_1.trans_expr[dataset_1.train_idx, :], dataset_2.trans_expr[dataset_2.train_idx, :]),
        ('val_rbps', dataset_1.rbp_expr[dataset_1.valid_idx, :], dataset_2.rbp_expr[dataset_2.valid_idx, :]),
        ('val_gn', dataset_1.gene_expr[dataset_1.valid_idx, :], dataset_2.gene_expr[dataset_2.valid_idx, :]),
        ('val_trans', dataset_1.trans_expr[dataset_1.valid_idx, :], dataset_2.trans_expr[dataset_2.valid_idx, :]),
        ('test_rbps', dataset_1.rbp_expr[dataset_1.test_idx, :], dataset_2.rbp_expr[dataset_2.test_idx, :]),
        ('test_gn', dataset_1.gene_expr[dataset_1.test_idx, :], dataset_2.gene_expr[dataset_2.test_idx, :]),
        ('test_trans', dataset_1.trans_expr[dataset_1.test_idx, :], dataset_2.trans_expr[dataset_2.test_idx, :])
    ]
    for name, data1, data2 in dataset_names:
        diff = np.abs(data1 - data2)  # Obtener las diferencias absolutas
        if np.all(diff == 0):
            print(f"Los datasets {name} son idénticos.")
        else:
            print(f"Los datasets {name} NO son iguales. Diferencias encontradas:")
            print(diff)  # Imprimir las diferencias

# Llamar a la función para mostrar las diferencias
compare_and_show_differences(dataset_1, dataset_2)
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

dataloaders = create_dataloaders(dataset, batch_size=batch_size)

# Acceder a los loaders que se hayan creado
train_loader = dataloaders.get('train_loader', None)
valid_loader = dataloaders.get('validation_loader', None)
test_loader = dataloaders.get('test_loader', None)

# esto ya es parte del train.py
num_epochs = 10
for epoch in range(num_epochs):
    for batch_index, (rbp_exp, gene_exp, trans_exp) in enumerate(train_loader):
        print(type(rbp_exp), type(gene_exp), type(trans_exp))  # Verificar tipos
        print(rbp_exp, gene_exp, trans_exp)  # Inspeccionar contenido
        print(rbp_exp.shape, gene_exp.shape, trans_exp.shape) 
        break
    break

    
