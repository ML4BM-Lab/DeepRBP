# src/deeprbp/training_module/pipeline.py

import os
from torch.utils.data import DataLoader

from ..util.logger import Logger
from ..util.utils import ensure_directory_exists, save_processed_data, adjust_batch_size, CustomTensorDataset
from .evaluation import calculate_metrics, save_metrics_summary, calculate_metrics_per_category
from .plots import plot_loss_curve
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DatasetLoader, DataSplitter, Scaler
from .model import PredictorModel
from .train_model import TrainPredictor

class DeepRBPredictorPipeline:
    """
    A pipeline class for training and evaluating the DeepRBP predictor model.

    This class handles the entire workflow including loading data, 
    processing it, training the model, evaluating its performance, 
    and saving the results. It supports both primary and external datasets.

    Attributes:
        config_parser (ConfigParser): Parser for the main configuration file.
        base_config (dict): Base configuration settings.
        training_config (dict): Configuration settings for training.
        path_save_data (str): Path for saving processed data.
        path_save_results (str): Path for saving results.
        data_importer (DataImporter): Data importer for primary data.
        data_loader (DatasetLoader): Data loader for primary data.
        data_loader_external (DatasetLoader, optional): Data loader for external data.
        scaler (Scaler): Scaler for normalizing data.
        logger (Logger): Logger for logging messages during execution.

    Methods:
        load_and_process_data(): Loads and processes primary and optional external data.
        split_data(data): Splits the dataset into training, validation, and test sets.
        scale_data(train_data, valid_data, test_data, external_data=None): Scales the datasets.
        save_data(train_data, valid_data, test_data, external_data=None): Saves the processed datasets.
        create_data_loaders(train_data, valid_data, test_data, external_data=None): Creates data loaders for training and evaluation.
        train_model(train_loader, val_loader): Trains the predictor model.
        save_model_and_history(trainer, train_history, val_history): Saves the trained model and training history.
        evaluate_model(trainer, train_loader, val_loader, test_loader, external_loader=None): Evaluates model performance on various datasets.
        evaluate_metrics_per_category(data, trainer, set_name, base_config): Evaluates metrics for different categories in the dataset.
        run(): Executes the entire pipeline from data loading to model evaluation.
    """
    def __init__(self, config_path, external_config_path=None, verbose=1):
        self.logger = Logger(verbose)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing the DeepRBPredictorPipeline...", level=1)

        self.config_parser = ConfigParser(config_path)
        self.base_config = self.config_parser.get_base_config()
        self.training_config = self.config_parser.get_model_training_config()
        
        # Define paths for saving data and results
        self.path_save_data = os.path.join(self.base_config['output_dir'], 'data')
        self.path_save_results = os.path.join(self.base_config['output_dir'], 'results')
        ensure_directory_exists(self.path_save_data)
        ensure_directory_exists(self.path_save_results)

        self.logger.log("✅ Directories for saving data and results are ready.", level=1)
        
        # Initialize DataImporter and DatasetLoader for primary data
        self.data_importer = DataImporter(self.base_config['data_paths'])
        self.data_loader = DatasetLoader(self.data_importer, self.base_config) 
        self.scaler = Scaler()
        
        # Optional: Load external dataset configuration
        if external_config_path:
            self.logger.log("📁 Loading external configuration...", level=1)
            self.external_config_parser = ConfigParser(external_config_path)
            self.external_base_config = self.external_config_parser.get_base_config()
            # Modify the output directory to match the base config
            self.external_base_config['output_dir'] = self.base_config['output_dir']
            self.data_importer_external = DataImporter(self.external_base_config['data_paths'])
            self.data_loader_external = DatasetLoader(self.data_importer_external, self.external_base_config) 

    def load_and_process_data(self):
        self.logger.log("🔄 Loading and processing data...", level=1)
        data = self.data_loader.load_data()
        external_data = (
            self.data_loader_external.load_data() if hasattr(self, 'data_loader_external') else None
        )
        self.logger.log("✅ Data loading complete.", level=1)
        return data, external_data

    def split_data(self, data):
        self.logger.log("✂️ Splitting data into train, validation, and test sets...", level=1)
        splitter = DataSplitter(data, self.training_config)
        return splitter.split_data_sets()
    
    def scale_data(self, train_data, valid_data, test_data, external_data=None):
        self.logger.log("📏 Scaling data...", level=1)
        self.scaler.fit(train_data['rbp_expr_log2p_tpm_df'])
        for dataset in [train_data, valid_data, test_data]:
            dataset['scaled_rbp_expr_log2p_tpm_df'] = self.scaler.transform(dataset['rbp_expr_log2p_tpm_df'])

        # If external data is provided, scale it as well
        if external_data is not None:
            external_data['scaled_rbp_expr_log2p_tpm_df'] = self.scaler.transform(external_data['rbp_expr_log2p_tpm_df'])
        self.logger.log("✅ Data scaling complete.", level=1)
        return train_data, valid_data, test_data, external_data
    
    def save_data(self, train_data, valid_data, test_data, external_data=None):
        self.logger.log("💾 Saving processed data...", level=1)
        self.scaler.save(os.path.join(self.path_save_data, 'scaler_trained'))
        for name, dataset in zip(['train', 'valid', 'test'], [train_data, valid_data, test_data]):
            save_processed_data(dataset, os.path.join(self.path_save_data, f'{name}_data'))
            self.logger.log(f"✅ {name.capitalize()} data saved.", level=1)
        if external_data is not None:
            save_processed_data(external_data, os.path.join(self.path_save_data, 'external_data'))
            self.logger.log("✅ External data saved.", level=1)

    def create_data_loaders(self, train_data, valid_data, test_data, external_data=None):
        self.logger.log("🔄 Creating data loaders...", level=1)
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
        self.logger.log("✅ Data loaders created.", level=1)
        return loaders
    
    def train_model(self, train_loader, val_loader):
        self.logger.log("🛠️ Training the model...", level=1)
        model = PredictorModel(
            input_size=next(iter(train_loader))['scaled_rbp_expr_log2p_tpm'].shape[1],
            output_size=next(iter(train_loader))['trans_expr_log2p_tpm'].shape[1],
            config=self.training_config
        )
        trainer = TrainPredictor(model=model, config=self.training_config)
        train_history, val_history = trainer.fit(train_loader, val_loader, epochs=self.training_config['epochs'])
        self.logger.log("✅ Model training complete.", level=1)
        return train_history, val_history, trainer
    
    def save_model_and_history(self, trainer, train_history, val_history):
        self.logger.log("💾 Saving model and training history...", level=1)
        trainer.model.save_model(self.path_save_results, 'model.pt')
        plot_loss_curve(train_history, val_history, output_dir=self.path_save_results)
        self.logger.log("✅ Model and history saved.", level=1)

    def evaluate_model(self, trainer, train_loader, val_loader, test_loader, external_loader=None):
        self.logger.log("📊 Evaluating model performance...", level=1)
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
        self.logger.log("✅ Model evaluation complete.", level=1)

    def evaluate_metrics_per_category(self, data, trainer, set_name, base_config):
        self.logger.log(f"📊 Evaluating metrics per category for {set_name}...", level=1)
        metrics_list, set_names_list = calculate_metrics_per_category(
            data,
            trainer,
            output_dir=self.path_save_results,
            set_name=set_name,
            sample_category=base_config['sample_category'],
            batch_size=self.training_config['batch_size'],
            source_name=base_config['source_name'],
            plot_results=base_config['plot_results'],
            getBM=self.data_importer.getBM if self.data_importer.getBM is not None and not self.data_importer.getBM.empty else None
        )
        save_metrics_summary(
            set_names_list,
            metrics_list,
            output_path=f'{self.path_save_results}/{set_name}_metrics_summary_per_category.csv')
        self.logger.log(f"✅ Metrics for {set_name} evaluated.", level=1)

    def run(self):
        self.logger.log("🚀 Starting the pipeline run...", level=1)

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

        self.logger.log("✅ Pipeline run completed successfully.", level=1)