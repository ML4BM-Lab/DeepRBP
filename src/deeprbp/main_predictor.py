
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
from train_predictor import TrainPredictor 

### (((((POR AQUI BROTHER VAS!!! COMPARA CON EL UNTITLED-1))) ejecuta interactivamente.
#config['source_name']
    #config['data_paths']
    #config['model']
    #config['training']
    #config['output_dir']
    #config['training']['batch_size'] -> antes: config['batch_size']

# Script to manage data and configuration for training
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

    # Initialize TrainPredictor
    trainer = TrainPredictor(model=model, config=config)

    # Train the model
    train_history, val_history = trainer.fit(train_loader, val_loader, epochs=config['epochs'])
    
    # Save the model and metrics
    trainer.save_model_and_metrics(config['output_dir'], 'model.pt', train_history, val_history)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    
    args = parser.parse_args()
    main(args)



for batch in train_loader:
    print(batch)  # <<< Para inspeccionar qué contiene el batch
    break

    inputs, targets, gen_expr = (x.to(self.device).float() for x in batch)




### basura: 
# import pandas as pd
# import os
# import numpy as np
# from tqdm import tqdm  # Asegúrate de que tqdm esté instalado
# import torch  # Necesario para guardar el modelo
# from config_loader import load_config
# from data_loader import CustomDataset, create_dataloaders
# from model import TranscriptExpressionPredictor
# from config_loader import load_config
# from utils import get_sample_ids_by_type, filter_data_by_sample_ids, save_processed_data, load_processed_data, adjust_batch_size



# # -------------------
# # CONFIGURACIÓN
# # -------------------

## ESTO HAY QUE METERLO EN UN MAIN 
# d) clase orquestradora: esta clase orquestaría todo el flujo de trabajo, desde cargar los datos hasta escalar y 
# guardar. No sé aun si iría en el train_predictor.py o donde (probablmemente no la definamos ahí pero la llamaremos ahí seguro)
paths_TCGA = {
    "rbp_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/RBPs_log2p_tpm.csv",
    "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/trans_log2p_tpm.csv",
    "metadata_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/phenotype_metadata.csv",
    "gene_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/gn_expr_each_iso_tpm.csv"
    }
    
# Path a los archivos de configuración
path_configs = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs'
path_config = os.path.join(path_configs, 'config_tcga_train.yaml')

# Cargar configuración [Caso 1) config train TCGA]
config = load_config(path_config)
output_dir = config['output_dir']
#source_name = config['source_name']
batch_size = config['batch_size']
epochs = config['epochs']

# # -------------------
# # CARGA DE DATOS
# # -------------------
# -----------------
data_importer = DataImporter(paths_TCGA)
data = data_importer.load()

# Usar las funciones auxiliares directamente
selected_sample_ids = get_sample_ids_by_type(
        data["metadata_df"], config['sample_category'], config['select_samples']
    )
data = filter_data_by_sample_ids(data, selected_sample_ids)

# Verificamos si se necesita hacer el split
if config["train_test_split"] or config["train_val_split"]:
    # Si es necesario dividir los datos, instanciamos el DataSplitter
    splitter = DataSplitter(data, config)
    train_data, valid_data, test_data = splitter.split_data_sets()  # Llamamos a split_data_sets
else:
    # Si no se necesita división, simplemente usamos los datos tal cual
    train_data, valid_data, test_data = data, None, None  # No se hace ninguna división

# Ahora puedes trabajar directamente con los datasets
print(train_data['metadata_df'])
print(valid_data['metadata_df'])
print(test_data['metadata_df'])
##
# Case 1: No existing scaler and sigma provided, fit the scaler externally
scaler = Scaler()  # Initialize Scaler with no existing scaler/sigma

# Fit the scaler using training data (e.g., train_data['rbp'])
scaler.fit(train_data['rbp_expr_df'])

# Now use the fitted scaler to transform train, validation, and test sets
train_data['scaled_rbp_expr_df'] = scaler.transform(train_data['rbp_expr_df'])
valid_data['scaled_rbp_expr_df'] = scaler.transform(valid_data['rbp_expr_df'])
test_data['scaled_rbp_expr_df'] = scaler.transform(test_data['rbp_expr_df'])

# Guardar el scaler, sigma
scaler_path = os.path.join(output_dir, 'scaler')
scaler.save(folder_path=f"{scaler_path}/scaler_data")
# # Case 2: Existing scaler and sigma provided, no need to fit again
# Cargar el Scaler existente
#loaded_scaler = Scaler.load(folder_path=f"{scaler_path}/scaler")
#train_data['scaled_rbp_expr_df'] = loaded_scaler.transform(train_data['rbp_expr_df'])

# Guardar datasets procesados
save_processed_data(train_data, os.path.join(output_dir, 'train_data'))
save_processed_data(valid_data, os.path.join(output_dir, 'valid_data'))
save_processed_data(test_data, os.path.join(output_dir, 'test_data'))
#loaded_train_data = load_processed_data(path)

# ahora crear el Dataset y Loader
# Paso 1: Crear instancia de CustomTensorDataset
train_dataset = CustomTensorDataset(train_data)  
valid_dataset = CustomTensorDataset(valid_data) 
test_dataset = CustomTensorDataset(test_data)    

# Paso 2: Crear DataLoader
train_loader = DataLoader(
            train_dataset, adjust_batch_size(train_dataset, batch_size), 
            drop_last=True, shuffle=True)
val_loader = DataLoader(
            valid_dataset, adjust_batch_size(valid_dataset, batch_size*2), shuffle=False)
test_loader = DataLoader(test_dataset, 
            adjust_batch_size(test_dataset, batch_size*2), shuffle=False)

model = Model(input_size=1348, output_size=11459, config=config)
# Semilla para reproducibilidad
torch.manual_seed(42)

output_dir = f'{dataset.output_dir}/results_training'
device = model.device
epochs = 100

train_isoform_predictor(model=model, 
            epochs=epochs, 
            train_loader=train_loader, 
            val_loader=val_loader, 
            save_results=False, 
            #output_dir=output_dir, 
            model_name='model.pt', 
            print_every=1, 
            device=device)

##################################################################
##################################################################
##################################################################

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
