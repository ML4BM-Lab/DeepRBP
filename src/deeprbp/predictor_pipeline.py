
# predictor_pipeline.py

import argparse
import os
from torch.utils.data import DataLoader
from .config_loader import ConfigParser
from .processing import DataImporter, DatasetLoader, DataSplitter, Scaler
from .models import PredictorModel
from .utils import *
from .plots import plot_loss_curve
from .evaluation_utils import *
from .train_predictor import TrainPredictor 

class DeepRBPredictorPipeline:
    def __init__(self, config_path, external_config_path=None):
        
        self.config_parser = ConfigParser(config_path)
        self.base_config = self.config_parser.get_base_config()
        self.training_config = self.config_parser.get_model_training_config()
        
        # Define paths for saving data and results
        self.path_save_data = os.path.join(self.base_config['output_dir'], 'data')
        self.path_save_results = os.path.join(self.base_config['output_dir'], 'results')
        ensure_directory_exists(self.path_save_data)
        ensure_directory_exists(self.path_save_results)
        
        # Initialize DataImporter and DatasetLoader for primary data
        self.data_importer = DataImporter(self.base_config['data_paths'])
        self.data_loader = DatasetLoader(self.data_importer, self.base_config) 
        self.scaler = Scaler()
        
        # Optional: Load external dataset configuration
        if external_config_path:
            self.external_config_parser = ConfigParser(external_config_path)
            self.external_base_config = self.external_config_parser.get_base_config()
            # Modify the output directory to match the base config
            self.external_base_config['output_dir'] = self.base_config['output_dir']
            self.data_importer_external = DataImporter(self.external_base_config['data_paths'])
            self.data_loader_external = DatasetLoader(self.data_importer_external, self.external_base_config)  # New DataLoader for external

    def load_and_process_data(self):
        data = self.data_loader.load_data()
        external_data = (
            self.data_loader_external.load_data() if hasattr(self, 'data_loader_external') else None
        )
        return data, external_data

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
        loaders = [
            DataLoader(
                dataset,
                batch_size=adjust_batch_size(dataset, self.training_config['batch_size'] * (2 if idx > 0 else 1)),
                shuffle=(idx == 0),
                drop_last=(idx == 0),
            )
            for idx, dataset in enumerate(datasets)
        ]
        external_loader = (
            DataLoader(CustomTensorDataset(external_data), batch_size=self.training_config['batch_size'], shuffle=False)
            if external_data is not None
            else None
        )
        loaders.append(external_loader)
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
        trainer.model.save_model(self.path_save_results, 'model.pt')
        plot_loss_curve(train_history, val_history, output_dir=self.path_save_results)

    def evaluate_model(self, trainer, train_loader, val_loader, test_loader, external_loader=None):
        preds_labels = [trainer.generate_predictions(loader) for loader in [train_loader, val_loader, test_loader]]
        metrics = [calculate_metrics(pred, label) for pred, label, _ in preds_labels]

        # Evaluate external data if present
        if external_loader is not None:
            external_preds_labels = trainer.generate_predictions(external_loader)
            metrics_external = calculate_metrics(*external_preds_labels[:2])
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
        data, external_data = self.load_and_process_data()

        # Split training data
        train_data, valid_data, test_data = self.split_data(data)
       
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

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP predictor training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--external_config_path', type=str, help='Path to the external dataset configuration file.')
    return parser.parse_args()
    
def main():
    args = parse_args()
    pipeline = DeepRBPredictorPipeline(args.config_path, args.external_config_path)
    pipeline.run()

if __name__ == "__main__":
    main()

#python /scratch/jsanchoz/DeepRBP/src/deeprbp/predictor_pipeline.py --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml" --external_config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_gtex.yaml"
