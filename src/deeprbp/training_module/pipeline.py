# src/deeprbp/training_module/pipeline.py

import os
from torch.utils.data import DataLoader
import pandas as pd
 
from ..util.logger import Logger
from ..util.utils import (
    filter_data_by_sample_ids,
    save_data,
    adjust_batch_size,
    summarize_model,
    CustomTensorDataset,
    format_output
)

from .evaluation import calculate_metrics, calculate_metrics_per_category, calculate_spearman_corr_per_gene
from ..util.plots import plot_loss_curve
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DataSplitter, Scaler
from .model import PredictorModel
from .train_model import TrainPredictor

class DeepRBPredictorPipeline:
    """
    A pipeline class for training and evaluating the DeepRBP predictor model.

    This class handles the entire workflow including loading data, 
    processing it, training the model, evaluating its performance, 
    and saving the results.  

    Args:
        config_path (str): Path to the configuration file of this run.
        output_dir (str): Directory to save the data and results.
        verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during the pipeline execution.
                                 - 0: No logging.
                                 - 1: Basic logging (show progress and essential logs).
                                 - 2: Detailed logging (show additional information).
                                 Default is 1.
    """
    def __init__(self, config_path, output_dir, verbose=1):
        self.logger = Logger(verbose=verbose)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing the DeepRBPredictorPipeline...", level=1)
        
        # Define paths for saving data and results
        self.path_save_data = os.path.join(output_dir, 'data')
        self.path_save_results = os.path.join(output_dir, 'results')
        self.logger.log("✅ Directories for saving data and results are ready.", level=1)
        
        # Set the configuration based on yaml file path
        self.config = ConfigParser(config_path)
        self.data_paths = self.config.get('data_paths')
        self.getBM = pd.read_csv(self.config.get('getBM_path'))
        self.rbp_names = None
        self.gene_names = None
        self.trans_names = None
        
        # Log the retrieved paths
        self.logger.log(f"Data paths: {self.data_paths}", level=1)
        # Initialize DataImporter for importing data
        self.data_importer = DataImporter(self.data_paths)
        # Initialize Scaler object
        self.scaler = None
        # Initialize Trainer object
        self.trainer = None

    def import_data(self):
        """Imports training data and return."""
        # Import training data
        data = self.data_importer.load()    
        self.logger.log("✅ Data import complete.", level=1)
        return data
    
    def filter_samples(self, data):
        """
        Filter a portion of the samples to optimize time and computational resources.

        Args:
            data: DataFrame to filter based on the sample fraction.

        Returns:
            Filtered DataFrame if sample_fraction is defined in the config, else original data.
        """
        sample_fraction = self.config.get('sample_fraction', default=None)
        if sample_fraction is not None:
            self.logger.log(f"[*] Filtering samples to optimize time and resources with fraction: {sample_fraction}...", level=1)
            subset_idx, _ = DataSplitter.split_data_class(
                data=data,
                config=self.config,
                sample_category=self.config.get('sample_category'),
                test_size=sample_fraction
            )
            data_subset = filter_data_by_sample_ids(data, subset_idx)
            self.logger.log("[*] Samples filtered successfully.", level=1)
            return data_subset
        else:
            self.logger.log("No sample_fraction defined in config; returning original data.", level=1)
            return data
        
    def split_data(self, data, test_name='validation'):
        self.logger.log("✂️ Splitting data into train and validation (or test) sets...", level=1)
        splitter = DataSplitter(data, self.config)
        return splitter.split_data_sets(test_name)
    
    def save_split_data(self, train_data, valid_data):
        """
        Saves the training and validation data to specified paths.

        Args:
            train_data: DataFrame containing the training data.
            valid_data: DataFrame containing the validation data.
        """
        # Define paths and filenames for saving data
        data_to_save = [
            (train_data, os.path.join(self.path_save_data, 'Train'), {
                'rbp_df': 'train_RBPs_log2p_tpm.csv',
                'isoform_df': 'train_trans_log2p_tpm.csv',
                'gene_df': 'train_gn_tpm.csv',
                'metadata_df': 'train_phenotype_metadata.csv'
            }),
            (valid_data, os.path.join(self.path_save_data, 'Validation'), {
                'rbp_df': 'val_RBPs_log2p_tpm.csv',
                'isoform_df': 'val_trans_log2p_tpm.csv',
                'gene_df': 'val_gn_tpm.csv',
                'metadata_df': 'val_phenotype_metadata.csv'
            })
        ]
        for data, save_path, custom_names in data_to_save:
            self.logger.log(f"[*] Saving data to: {save_path}...", level=1) 
            save_data(data, save_path, custom_names)
            self.logger.log(f"[*] Data saved successfully.", level=1) 

    def fit_scaler(self, train_data):
        """Fits the scaler to the training data.

        Args:
            train_data (DataFrame): The training data used to fit the scaler.
        """
        self.scaler = Scaler()
        self.scaler.fit(train_data['rbp_df'])
        self.logger.log("✅ Scaler has been fitted to training data.", level=1)
        self.scaler.save(self.path_save_results)
        self.logger.log(f"✅ Scaler has been saved in {self.path_save_results}.", level=1)

    def load_scaler(self, folder_path):
        """
        Load an already trained scaler from the specified directory and assign it to self.scaler.

        Args:
            folder_path (str): The directory where the scaler and sigma are stored.

        Returns:
            None: The scaler is loaded and assigned to self.scaler.
        """
        try:
            self.scaler = Scaler.load(folder_path)
            self.logger.log("✅ Scaler successfully loaded from the specified path.", level=1)
        except FileNotFoundError as e:
            self.logger.error(f"❌ [Scaler:load_scaler] Failed to load scaler: {e}", level=1)
        except Exception as e:
            self.logger.error(f"❌ [Scaler:load_scaler] An unexpected error occurred: {e}", level=1)

    def scale_data(self, *datasets):  
        """Transforms the provided DataFrames using the fitted scaler.

        Args:
            *datasets (DataFrame): The data to be transformed.

        Returns:
             dict or List[dict]: A dictionary or list of dictionaries containing the scaled DataFrames.
        """
        if self.scaler is None:
            self.logger.error("❌ Scaler is not fitted. Please fit the scaler before scaling data.", ValueError)
            raise ValueError("Scaler is not fitted. Please fit the scaler before scaling data.")
        # Transform the provided data using the fitted scaler
        scaled_datasets = []
        for data in datasets:
            # Transform the provided data using the fitted scaler
            scaled_data = data.copy()   
            scaled_data['scaled_rbp_df'] = self.scaler.transform(data['rbp_df'])  # Scale the RBP data 
            scaled_datasets.append(scaled_data)   
        return format_output(scaled_datasets)  
    
    def create_datasets(self, *dataframes):
        """
        Creates CustomTensorDataset instances for provided datasets.

        Initializes CustomTensorDataset objects for the provided datasets.

        Args:
            *dataframes (DataFrame): The DataFrames to create datasets from.

        Returns:
        List[CustomTensorDataset]: A list of CustomTensorDataset instances for multiple datasets.
        """
        datasets = []
        for data in dataframes:
            dataset = CustomTensorDataset(
                data=data.copy(),
                getBM=self.getBM,
                rbp_data_key='scaled_rbp_df', 
                gene_data_key='gene_df', 
                transcript_data_key='isoform_df',
                trans_col_name=self.config.get('trans_col_name'),
                gene_col_name=self.config.get('gene_col_name')
            )
            datasets.append(dataset)
            if len(datasets) == 1:  # Just assign for the first dataset
                self.rbp_names = dataset.rbp_names
                self.gene_names = dataset.gene_names
                self.trans_names = dataset.trans_names
        return datasets
    
    def create_data_loaders(self, *datasets, batch_sizes=None):
        """
        Creates DataLoader instances for provided datasets with specified batch sizes.

        Args:
            *datasets (CustomTensorDataset): Datasets for which to create DataLoaders.
            batch_sizes (list, optional): A list of batch sizes for each dataset. 
                                        If not provided, the first dataset will use config.train_batch_size 
                                        and subsequent datasets will use config.val_batch_size.

        Returns:
            List[DataLoader]: A list of DataLoader instances for multiple datasets.
        Raises:
        ValueError: If train_batch_size or val_batch_size is not defined in the configuration.
        """
        self.logger.log("🔄 Creating data loaders...", level=1)
        data_loaders = []
        # Check if train_batch_size and val_batch_size are defined
        train_batch_size = self.config.get('train_batch_size')
        val_batch_size = self.config.get('val_batch_size')
        # Log error and raise exception if batch sizes are not defined
        if train_batch_size is None:
            error_message = "❌ train_batch_size is not defined in the configuration."
            self.logger.error(error_message, level=1)
            raise ValueError(error_message)
        if val_batch_size is None:
            error_message = "❌ val_batch_size is not defined in the configuration."
            self.logger.error(error_message, level=1)
            raise ValueError(error_message)
        # Default batch sizes if not provided
        if batch_sizes is None:
            batch_sizes = [train_batch_size] + [val_batch_size] * (len(datasets) - 1)
        for idx, dataset in enumerate(datasets):
            current_batch_size = batch_sizes[idx]
            loader = DataLoader(
                dataset,
                batch_size=adjust_batch_size(dataset, current_batch_size),
                shuffle=(current_batch_size == train_batch_size),  # Shuffle only for the dataset with train batch size
                drop_last=(current_batch_size == train_batch_size)  # Drop last batch only for train batch size
            )
            data_loaders.append(loader)
        self.logger.log("✅ Data loaders created.", level=1)
        return data_loaders
    
    def get_loaders(self, *dataframes):
        """
        Encapsulates the creation of datasets and data loaders.

        Args:
            *dataframes (DataFrame): The DataFrames to create datasets and loaders from.

        Returns:
        List[DataLoader]: A list of DataLoader instances for multiple DataFrames.
        """
        datasets = self.create_datasets(*dataframes)  
        return self.create_data_loaders(*datasets) 
    
    def _get_model(self, loader, path_to_weights=None):
        """
        Creates an instance of the PredictorModel. If a path to weights is provided, it loads the model with those weights.

        Args:
            loader (DataLoader): DataLoader for the dataset.
            path_to_weights (str, optional): Path to the file containing the pre-trained weights. Default is None.

        Returns:
            PredictorModel: An instance of the PredictorModel initialized with the configuration or loaded with weights.
        """
        input_size = next(iter(loader))['scaled_rbp_df'].shape[1]
        output_size = next(iter(loader))['isoform_df'].shape[1]
        if path_to_weights:
            # Load the model using the loaded weights
            model = PredictorModel.load_model(path_to_weights, self.config, input_size=input_size, output_size=output_size)
        else:
            model = PredictorModel(input_size=input_size, output_size=output_size, config=self.config)
        # Reassign the model's configuration directly to the pipeline (updating if there were 'unused' variables)
        self.config = model.config
        print("\n[*] Model Summary:")
        summarize_model(model, loader, self.config.get("train_batch_size"))
        return model
    
    def _get_trainer(self, model):
        """
        Creates an instance of the TrainPredictor.

        Args:
            model (PredictorModel): The model instance to be used by the trainer.

        Returns:
            TrainPredictor: An instance of the TrainPredictor initialized with the model.
        """
        trainer = TrainPredictor(
            model=model,
            config=self.config,
            input_features=('scaled_rbp_df', 'gene_df'), 
            output_features=('isoform_df',),
            verbose=self.logger.verbose # Pass the verbosity level from the logger
        )
        return trainer
    
    def setup_model_trainer(self, loader, path_to_weights=None):
        """
        Creates both the model and trainer, returning the trainer. Allows for loading a pre-trained model if a path is provided.

        Args:
            loader (DataLoader): DataLoader for the dataset.
            path_to_weights (str, optional): Path to the file containing the pre-trained weights. Default is None.

        Returns:
            TrainPredictor: An instance of the TrainPredictor initialized with the model.
        """
        model = self._get_model(loader, path_to_weights) # Create the model
        self.trainer = self._get_trainer(model) # Create and store the trainer
        # Mark the trainer as trained if a pre-trained model was loaded
        if path_to_weights:
            self.trainer.is_trained = True  # Mark the trainer as trained

    def train_model(self, train_loader, val_loader, trial=None): 
        """
        Trains the model using the provided data loaders.

        Parameters:
            train_loader (DataLoader): DataLoader for the training dataset.
            val_loader (DataLoader): DataLoader for the validation dataset.
            trial (optuna.Trial, optional): The Optuna trial object for reporting metrics and pruning.

        Returns:
            tuple: A tuple containing:
                - list: History of training losses for each epoch.
                - list: History of validation losses for each epoch.
        """
        self.logger.log("🛠️ Training the model...", level=1)
        # The setup_model_trainer is already called in run, so we assume self.trainer is set up. 
        # Prepare parameters for the fit method
        fit_kwargs = {
            'epochs': self.config.get("num_epochs"),
            'path_save_results': self.path_save_results  
        }
        if trial is not None:
            fit_kwargs['optuna_trial'] = trial   
        # Call trainer.fit with unpacked kwargs
        train_history, val_history = self.trainer.fit(train_loader, val_loader, **fit_kwargs)
        self.logger.log("✅ Model training complete.", level=1)
        return train_history, val_history
    
    def save_history(self, train_history, val_history, plot_name='loss_curve'):
        self.logger.log("💾 Saving training and validation history...", level=1)
        plot_loss_curve(train_history, val_history, output_dir=self.path_save_results, plot_name=plot_name)
        self.logger.log("✅ Saving training and validation history.", level=1)

    def evaluate_model(self, loader, set_name='validation'):
        """
        Evaluate the model's performance on the given dataset.

        This method generates predictions for the dataset, calculates various performance metrics, 
        computes Spearman correlations for each gene, and saves the metrics summary to a CSV file 
        if required.

        Parameters:
        loader (DataLoader): DataLoader for the dataset, which provides batches of input data.
        set_name (str): The name of the dataset being evaluated (e.g., 'train' or 'val').

        Returns:
        tuple: A tuple containing:
            - metrics (dict): A dictionary with general metrics for the dataset.
            - gene_corr_df (DataFrame): DataFrame containing gene IDs, their Spearman correlations per gene, 
            and the maximum transcript correlations for each gene.
        """
        self.logger.log(f"📊 Evaluating model performance on {set_name}...", level=1)
        # Generate predictions
        preds, labels = self.trainer.generate_predictions(loader)
        # Calculate metrics
        metrics_general = calculate_metrics(preds.flatten(), labels.flatten())  
        # Calculate correlations per gene
        corr_per_gene = calculate_spearman_corr_per_gene(self.gene_names, self.trans_names, preds, labels, self.getBM)  
        # Prepare metrics for saving
        metrics_dict = {
            'spearman_corr': metrics_general['spearman_corr'],
            'pearson_corr': metrics_general['pearson_corr'],
            'mse': metrics_general['mse'],
            'r2': metrics_general['r2'],
            'mean_corr_per_gene': corr_per_gene['mean_corr_per_gene'],
            'mean_corr_max_trans_per_gene': corr_per_gene['mean_corr_max_trans_per_gene']
        }
        metrics = pd.DataFrame(metrics_dict, index=[set_name])
        metrics.index.name = "Set"
        # Create DataFrame for gene correlations
        gene_corr_df = pd.DataFrame({
            'gene_id': list(corr_per_gene['gene_corrs_dict'].keys()),  # Gene IDs
            'spearman_corr': list(corr_per_gene['gene_corrs_dict'].values()),  # Spearman correlations
            'gene_id-trans_id': list(corr_per_gene['gene_corrs_dict_max_trans'].keys()),  # Gene IDs
            'max_trans_spearman_corr': list(corr_per_gene['gene_corrs_dict_max_trans'].values())  # Max transcript correlations
        })
        # Save metrics summary to CSV if needed
        if self.config.get('save_results'):
            self.logger.log(f"📂 Saving metrics summary and gene correlation results for the {set_name} dataset...")
            metrics.to_csv(os.path.join(self.path_save_results, f'{set_name}_metrics_summary_global.csv'), index=True)
            gene_corr_df.to_csv(os.path.join(self.path_save_results, f'{set_name}_per_gene_correlations.csv'), index=False)
            self.logger.log(f"✅ Results for the {set_name} dataset have been successfully saved to files.")
        self.logger.log(f"✅ Model evaluation on {set_name} complete.", level=1)
        return metrics, gene_corr_df
    
    def eval_model_per_category(self, data, set_name='validation'): 
        self.logger.log(f"📊 Evaluating metrics per category ...", level=1)
        calculate_metrics_per_category(
              data,
              self.trainer,
              self.path_save_results,
              set_name,
              self.config,
              source_name='TCGA',
              getBM=self.getBM)
        
    def run(self):  
        self.logger.log("🚀 Starting the DeepRBPredictorPipeline run...", level=1)
        # Import data
        self.logger.log("📥 Importing training and test data...", level=1)
        data = self.import_data()
        # Split the filtered training data into train and validation sets
        train_data, valid_data = self.split_data(data)
        # Save the split datasets
        self.save_split_data(train_data, valid_data)   
        # Scale data
        self.fit_scaler(train_data) 
        scaled_train_data, scaled_valid_data = self.scale_data(train_data, valid_data)  
        # Create data loaders
        train_loader, valid_loader = self.get_loaders(scaled_train_data, scaled_valid_data)
        # Setup model and trainer before training
        self.setup_model_trainer(train_loader) # Create the model and trainer, which sets self.trainer
        # Train the model
        train_history, val_history = self.train_model(train_loader, valid_loader)
        self.save_history(train_history, val_history)
        # Evaluate the model on training and validation sets
        for loader, name in zip([train_loader, valid_loader], ['train', 'val']):
            metrics, _ = self.evaluate_model(loader, set_name=name)
            # Optionally, you can log or save metrics here if needed
            self.logger.log(f"Metrics for {name}: {metrics}", level=1)   
        # Evaluate metrics per category
        for dataset, name in zip([scaled_train_data, scaled_valid_data], ['train', 'val']):
            self.eval_model_per_category(dataset, name)  
        self.logger.log("✅ Pipeline run completed successfully.", level=1)



# ######### PROBAR OTRA VEZ TODO JOSEBA, LOS NOMBRES DE LOS METODOS NO ME ACABAN DE CONVENCER.
# config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_train.yaml'
# output_dir = '/scratch/jsanchoz/DeepRBP/output/results'

#pipeline = DeepRBPredictorPipeline(config_path, output_dir)

# # prueba 1: step-by-step
# data = pipeline.import_data()
# train_data, valid_data = pipeline.split_data(data)
 
# pipeline.fit_scaler(train_data) 
# train_data, valid_data = pipeline.scale_data(train_data, valid_data) 

# train_loader, valid_loader = pipeline.get_loaders(train_data, valid_data)
# # ##
# # train_dataset, valid_dataset = pipeline.create_datasets(train_data, valid_data)
# # train_loader, valid_loader = pipeline.create_data_loaders(train_dataset, valid_dataset)
# # ##

# train_history, val_history = pipeline.train_model(train_loader, valid_loader)
# pipeline.save_history(train_history, val_history)

# metrics, gene_corr_df = pipeline.evaluate_model(valid_loader, set_name='val')

# FOR TEST DATA! PRUEBA ESTA JOSEBA!
# config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_test.yaml'
# output_dir = '/scratch/jsanchoz/DeepRBP/output/results'
# trained_files_path = '/scratch/jsanchoz/DeepRBP/output/results/results'

# pipeline2 = DeepRBPredictorPipeline(config_path, output_dir)

# test_data = pipeline2.import_data()
# pipeline2.load_scaler(folder_path=trained_files_path)
# test_data = pipeline2.scale_data(test_data)
# test_loader = pipeline2.get_loaders(test_data)[0]
# # ##
# # test_dataset = pipeline2.create_datasets(test_data)[0]
# # test_loader = pipeline2.create_data_loaders(test_dataset)[0]
# # ##

# pipeline2.setup_model_trainer(test_loader, path_to_weights=f'{trained_files_path}/best_model.pt')
# metrics, gene_corr_df = pipeline2.evaluate_model(test_loader, set_name='test')

