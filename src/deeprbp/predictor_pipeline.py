
# main_predictor.py (ESTUDIAR METER ESTE SCRIPT FUERA DEL PAQUETE Y EJECUTA con ayuda del __init__ de deeprbp)

# from datetime import datetime
# from deeprbp.config_loader import load_config
# from deeprbp.data_importer import DataImporter
# from deeprbp.data_splitter import DataSplitter
# from deeprbp.utils import filter_data_by_sample_ids

# main_predictor.py -> predictor_pipeline.py
import argparse
from dataclasses import asdict
import os
from torch.utils.data import DataLoader
#import logging
from config_loader import ConfigParser
from processing import DataImporter, DataSplitter, Scaler
from model import PredictorModel
from utils import *
from plots import plot_loss_curve
from evaluation_utils import *
from train_predictor import TrainPredictor 

class DeepRBPredictorPipeline:
    def __init__(self, config_path, external_config_path=None):
        self.config_parser = ConfigParser(config_path)
        self.base_config = self.config_parser.get_base_config()
        self.training_config = self.config_parser.get_model_training_config()
        # Initialize DataImporter for primary data
        self.data_importer = DataImporter(self.base_config['data_paths'])
        self.scaler = Scaler()
        # Define paths for saving data and results
        self.path_save_data = os.path.join(self.base_config['output_dir'], 'data')
        self.path_save_results = os.path.join(self.base_config['output_dir'], 'results')
        ensure_directory_exists(self.path_save_data)
        ensure_directory_exists(self.path_save_results)
        # Optional: Load external dataset configuration
        if external_config_path:
            self.external_config_parser = ConfigParser(external_config_path)
            self.external_base_config = self.external_config_parser.get_base_config()
            # Modify the output directory to match the base config
            self.external_base_config['output_dir'] = self.base_config['output_dir']
            self.data_importer_external = DataImporter(self.external_base_config['data_paths'])
    def load_data(self, data_importer, base_config):
        data = data_importer.load()
        sel_sample_ids = select_sample_ids_by_type(
            metadata_df = data['metadata_df'],
            sample_category_col = base_config['sample_category'],
            sample_types = base_config['select_samples']
        )
        return filter_data_by_sample_ids(data, sel_sample_ids)
    def split_data(self, data):
        splitter = DataSplitter(data, self.training_config)
        return splitter.split_data_sets()
    def scale_data(self, train_data, valid_data, test_data, external_data=None):
        self.scaler.fit(train_data['rbp_expr_df'])
        for dataset in [train_data, valid_data, test_data]:
            dataset['scaled_rbp_expr_df'] = self.scaler.transform(dataset['rbp_expr_df'])
        # If external data is provided, scale it as well
        if external_data is not None:
            external_data['scaled_rbp_expr_df'] = self.scaler.transform(external_data['rbp_expr_df'])
        return train_data, valid_data, test_data, external_data
    def save_data(self, train_data, valid_data, test_data, external_data=None):
        self.scaler.save(os.path.join(self.path_save_data, 'scaler_trained'))
        for name, dataset in zip(['train', 'valid', 'test'], [train_data, valid_data, test_data]):
            save_processed_data(dataset, os.path.join(self.path_save_data, f'{name}_data'))
        if external_data is not None:
            save_processed_data(external_data, os.path.join(self.path_save_data, 'external_data'))
    def create_data_loaders(self, train_data, valid_data, test_data, external_data=None):
        datasets = [CustomTensorDataset(data) for data in [train_data, valid_data, test_data]]
        loaders = []
        for idx, dataset in enumerate(datasets):
            batch_size = adjust_batch_size(dataset, self.training_config['batch_size'] * (2 if idx > 0 else 1))
            loaders.append(DataLoader(dataset, batch_size=batch_size, shuffle=(idx == 0), drop_last=(idx == 0)))
        # Create DataLoader for external data if it exists
        if external_data is not None:
            external_dataset = CustomTensorDataset(external_data)
            external_loader = DataLoader(external_dataset, batch_size=self.training_config['batch_size'], shuffle=False)
            loaders.append(external_loader)  # Add external loader to the list of loaders
        return loaders 
    def train_model(self, train_loader, val_loader):
        model = PredictorModel(
            input_size=next(iter(train_loader))['rbp_expr'].shape[1],
            output_size=next(iter(train_loader))['trans_expr'].shape[1],
            config=self.training_config
        )
        trainer = TrainPredictor(model=model, config=self.training_config)
        train_history, val_history = trainer.fit(train_loader, val_loader, epochs=self.training_config['epochs'])
        return train_history, val_history, trainer
    def save_model_and_history(self, trainer, train_history, val_history):
        trainer.save_model(self.path_save_results, 'model.pt')
        plot_loss_curve(train_history, val_history, output_dir=self.path_save_results)
    def evaluate_model(self, trainer, train_loader, val_loader, test_loader, external_loader=None):
        preds_labels = [trainer.generate_predictions(loader) for loader in [train_loader, val_loader, test_loader]]
        metrics = [calculate_metrics(pred, label) for pred, label, _ in preds_labels]
        # Evaluate external data if present
        if external_loader is not None:
            external_preds_labels = trainer.generate_predictions(external_loader)
            metrics_external = calculate_metrics(*external_preds_labels)
            metrics.append(metrics_external)  
        save_metrics_summary(
            set_names_list=['train', 'val', 'test'] + (['external'] if external_loader is not None else []),
            metrics_list=metrics,
            output_path=f'{self.path_save_results}/metrics_summary_global.csv'
        )
    def evaluate_metrics_per_category(self, data, trainer, set_name, base_config):
        metrics_list, set_names_list = calculate_metrics_per_category(
            data,
            trainer,
            output_dir=self.path_save_results,
            set_name=set_name,
            sample_category=base_config['sample_category'],
            batch_size=self.training_config['batch_size'],
            source_name=base_config['source_name'],
            plot_results=base_config['plot_results'],
            getBM_path=base_config['data_paths'].get('getBM_path', None)
        )
        save_metrics_summary(
            set_names_list,
            metrics_list,
            output_path=f'{self.path_save_results}/{set_name}_metrics_summary_per_category.csv')
    def run(self):
        # Load and process data
        data = self.load_data(self.data_importer, self.base_config)
        train_data, valid_data, test_data = self.split_data(data)
        external_data = self.load_data(self.data_importer_external, self.external_base_config) if hasattr(self, 'data_importer_external') else None
        # Scale data
        train_data, valid_data, test_data, external_data = self.scale_data(train_data, valid_data, test_data, external_data)
        self.save_data(train_data, valid_data, test_data, external_data)
        # Create data loaders
        train_loader, val_loader, test_loader, external_loader = self.create_data_loaders(train_data, valid_data, test_data, external_data)
        # Train the model
        train_history, val_history, trainer = self.train_model(train_loader, val_loader)
        self.save_model_and_history(trainer, train_history, val_history)
        # Evaluate on training, validation, and test sets
        self.evaluate_model(trainer, train_loader, val_loader, test_loader, external_loader)
        # Evaluate metrics per category
        for dataset, name in zip([train_data, valid_data, test_data], ['train', 'val', 'test']):
            self.evaluate_metrics_per_category(dataset, trainer, name, self.base_config)
        if external_data is not None:
            self.evaluate_metrics_per_category(external_data, trainer, 'external', self.external_base_config)

pipeline = DeepRBPredictorPipeline("/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml")
pipeline.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the DeepRBP predictor training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--external_config_path', type=str, help='Path to the external dataset configuration file.')
    args = parser.parse_args()

    pipeline = DeepRBPredictorPipeline(args.config_path, args.external_config_path)
    pipeline.run()



##### OLD VERSION
# Script to manage data and configuration for training
def main(args):
    #logging.basicConfig(level=logging.INFO) crea tu clase (LUIS)
    # Step 1: Load configuration from YAML file
    # parser = ConfigParser("/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml")
    # base_config = parser.get_base_config()
    
    # Step 2: Initialize DataImporter with paths from the config
    # data_importer = DataImporter(base_config['data_paths'])
    
    # Step 3: Load data
    # data = data_importer.load()
    
    # Step 4: Select sample IDs based on the provided categories in the config
    # sel_sample_ids = select_sample_ids_by_type(
    #         metadata_df = data['metadata_df'],
    #         sample_category_col = base_config['sample_category'], 
    #         sample_types = base_config['select_samples']
    # )
    # Step 5: Filter the data using the selected sample IDs
    # data = filter_data_by_sample_ids(data, sel_sample_ids)
    
    # Obtener la configuración del modelo y del entrenamiento
    # training_config = parser.get_model_training_config()
   
    # # Step 6: Data splitting based on the configuration settings
    # splitter = DataSplitter(data, training_config)
    # train_data, valid_data, test_data = splitter.split_data_sets()

    # Step 7: Scaler. No existing scaler and sigma provided, fit the scaler externally
    # scaler = Scaler()  # Initialize Scaler with no existing scaler/sigma
    # # Fit the scaler using training data
    # scaler.fit(train_data['rbp_expr_df'])
    
    # # Now use the fitted scaler to transform train, validation, and test sets
    # train_data['scaled_rbp_expr_df'] = scaler.transform(train_data['rbp_expr_df'])
    # valid_data['scaled_rbp_expr_df'] = scaler.transform(valid_data['rbp_expr_df'])
    # test_data['scaled_rbp_expr_df'] = scaler.transform(test_data['rbp_expr_df'])
    
    # # Step 8: Save the scaler for later use (optional) and processed datasets
    # path_save_data = os.path.join(base_config['output_dir'], 'data')
    # scaler_path = os.path.join(path_save_data, 'scaler_trained')
    # scaler.save(scaler_path)
    
    # save_processed_data(train_data, os.path.join(path_save_data, 'train_data'))
    # save_processed_data(valid_data, os.path.join(path_save_data, 'valid_data'))
    # save_processed_data(test_data, os.path.join(path_save_data, 'test_data'))
    # #loaded_train_data = load_processed_data(path)

    # Step 9: Create TensorDataset and DataLoaders
    # train_dataset = CustomTensorDataset(train_data)  
    # valid_dataset = CustomTensorDataset(valid_data) 
    # test_dataset = CustomTensorDataset(test_data)   

    # train_loader = DataLoader(
    #         train_dataset, 
    #         adjust_batch_size(train_dataset, training_config['batch_size']), 
    #         drop_last=True, 
    #         shuffle=True
    # )
    
    # val_loader = DataLoader(
    #         valid_dataset, 
    #         adjust_batch_size(valid_dataset, training_config['batch_size']*2), 
    #         shuffle=False
    # )  

    # test_loader = DataLoader(
    #      test_dataset, 
    #      adjust_batch_size(test_dataset, training_config['batch_size']*2), 
    #      shuffle=False
    # )

    # # Step 10: Initialize Model
    # model = PredictorModel(
    #     input_size=len(next(iter(train_dataset))['rbp_expr']), 
    #     output_size=len(next(iter(train_dataset))['trans_expr']), 
    #     config=training_config
    # )

    # # # Step 11: Initialize TrainPredictor
    # trainer = TrainPredictor(model=model, config=training_config)
    # # # Step 12: Train the model
    # train_history, val_history = trainer.fit(train_loader, val_loader, epochs=training_config['epochs'])

    # # Step 13: Save the model and history
    # path_save_results = os.path.join(base_config['output_dir'], 'results')
    # trainer.save_model(path_save_results, 'model.pt')
    # plot_loss_curve(train_history, val_history, output_dir=path_save_results)

    # # Step 14 make Predictions and calculate metrics on Training, Validation and Test (general)
    # pred_train, label_train, _ = trainer.generate_predictions(train_loader)
    # pred_val, label_val, _ = trainer.generate_predictions(val_loader)
    # pred_test, label_test, _ = trainer.generate_predictions(test_loader)

    # metrics_train = calculate_metrics(pred_train, label_train)
    # metrics_val = calculate_metrics(pred_val, label_val)
    # metrics_test = calculate_metrics(pred_test, label_test)

    # save_metrics_summary(set_names_list = ['train', 'val', 'test'], 
    #                     metrics_list = [metrics_train, metrics_val, metrics_test], 
    #                     output_path = f'{path_save_results}/metrics_summary_global.csv')
                   
    # # Step 15: Make predictions and calculate metrics (per category)
    # # TRAIN
    # metrics_list, set_names_list = calculate_metrics_per_category(
    #     train_data,
    #     trainer,
    #     output_dir=path_save_results,
    #     set_name='train',
    #     sample_category=parser.config.sample_category,
    #     batch_size=parser.config.training['batch_size'],
    #     source_name=parser.config.source_name,
    #     plot_results=parser.config.plot_results,
    #     getBM_path=parser.config.data_paths.get('getBM_path', None)
    # )
    # save_metrics_summary(
    #     set_names_list, 
    #     metrics_list, 
    #     output_path=f'{path_save_results}/train_metrics_summary_per_category.csv'
    # )

    # # VAL
    # metrics_list, set_names_list = calculate_metrics_per_category(
    #     valid_data,
    #     trainer,
    #     output_dir=path_save_results,
    #     set_name='val',
    #     sample_category=parser.config.sample_category,
    #     batch_size=parser.config.training['batch_size'],
    #     source_name=parser.config.source_name,
    #     plot_results=parser.config.plot_results,
    #     getBM_path=parser.config.data_paths.get('getBM_path', None)
    # )
    # save_metrics_summary(
    #     set_names_list, 
    #     metrics_list, 
    #     output_path=f'{path_save_results}/val_metrics_summary_per_category.csv'
    # )

    # # TEST
    # metrics_list, set_names_list = calculate_metrics_per_category(
    #     test_data,
    #     trainer,
    #     output_dir=path_save_results,
    #     set_name='test',
    #     sample_category=parser.config.sample_category,
    #     batch_size=parser.config.training['batch_size'],
    #     source_name=parser.config.source_name,
    #     plot_results=parser.config.plot_results,
    #     getBM_path=parser.config.data_paths.get('getBM_path', None)
    # )
    # save_metrics_summary(
    #     set_names_list, 
    #     metrics_list, 
    #     output_path=f'{path_save_results}/test_metrics_summary_per_category.csv'
    # )

    # # Step 16: Predictions on GTEX - ME QUEDA REVISAR EL TRAIN Y EL VAL EN EL ANTERIOR CODIGO Y HACER EN GTEX Y
    # # PEDIR CONSEJO A CHAT PREMIUM DE COMO HE ORGANIZADO TODO EL CODE EN LA PIPELINE A PARTIR DE TODO EL CODE. 
    # # Coger el GTEX y modificar el output_path para que matchee el anterior 
    # parser_gtex = ConfigParser("/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_gtex.yaml")
    # config_gtex = parser_gtex.get_base_config()
    # # modificar el output_path para que apunte al que nosotros queremos en tcga HACER ESTO!!!
    # config_gtex['output_dir'] = base_config['output_dir']

    # # Step 17: Initialize DataImporter with paths from the config
    # data_importer_gtex = DataImporter(config_gtex['data_paths'])
    
    # # Step 18: Load data y scaling
    # data_gtex = data_importer_gtex.load()
    # # Now use the fitted scaler to transform train, validation, and test sets
    # data_gtex['scaled_rbp_expr_df'] = scaler.transform(data_gtex['rbp_expr_df'])
  
    # # Step 19:    
    # save_processed_data(data_gtex, os.path.join(path_save_data, 'data_gtex'))
    # gtex_dataset = CustomTensorDataset(data_gtex)   
    # gtex_loader = DataLoader(
    #      gtex_dataset, 
    #      adjust_batch_size(gtex_dataset, training_config['batch_size']*2), 
    #      shuffle=False
    # )

    # # Step 20 make Predictions (general)
    # pred_gtex, label_gtex, _ = trainer.generate_predictions(gtex_loader)
    # metrics_gtex = calculate_metrics(pred_gtex, label_gtex)

    # save_metrics_summary(set_names_list = ['gtex'], 
    #                     metrics_list = [metrics_gtex], 
    #                     output_path = f'{path_save_results}/metrics_summary_global_gtex.csv')
                   
    # # Step 21: Make predictions and calculate metrics (per category)
    # # TRAIN
    # metrics_list, set_names_list = calculate_metrics_per_category(
    #     data_gtex,
    #     trainer,
    #     output_dir=path_save_results,
    #     set_name='gtex',
    #     sample_category=parser_gtex.config.sample_category,
    #     batch_size=parser_gtex.config.training['batch_size'],
    #     source_name=parser_gtex.config.source_name,
    #     plot_results=parser_gtex.config.plot_results,
    #     getBM_path=parser_gtex.config.data_paths.get('getBM_path', None)
    # )








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
