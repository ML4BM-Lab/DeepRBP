# src/deeprbp/training_module/pipeline.py

import os
from torch.utils.data import DataLoader
import pandas as pd
from collections import namedtuple

from ..util.logger import Logger
from ..util.utils import filter_data_by_sample_ids, save_data, adjust_batch_size, summarize_model, CustomTensorDataset
from .evaluation import calculate_metrics, save_metrics_summary, calculate_metrics_per_category, calculate_spearman_corr_per_gene
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
        config_path (str): Path to the configuration file of this run
        output_dir (str): Directory to save the data and results
    """
    def __init__(self, config_path, output_dir):
        self.logger = Logger(verbose=1)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing the DeepRBPredictorPipeline...", level=1)
        # Define paths for saving data and results
        self.path_save_data = os.path.join(output_dir, 'data')
        self.path_save_results = os.path.join(output_dir, 'results')
        self.logger.log("✅ Directories for saving data and results are ready.", level=1)
        # Set the configuration based on yaml file path
        self.config = ConfigParser(config_path)
        self.train_data_paths = self.config.get('data_paths')
        self.getBM = pd.read_csv(self.config.get('getBM_path'))
        # Log the retrieved paths
        self.logger.log(f"Training data paths: {self.train_data_paths}", level=1)
        # Initialize DataImporter for training data
        self.train_data_importer = DataImporter(self.train_data_paths)
    ###
    def import_data(self):
        """Imports training data and return."""
        # Import training data
        data = self.train_data_importer.load()    
        self.logger.log("✅ Data import complete.", level=1)
        return data
    ###
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
    ###
    def split_data(self, data, test_name='validation'):
        self.logger.log("✂️ Splitting data into train, validation, and test sets...", level=1)
        splitter = DataSplitter(data, self.config)
        return splitter.split_data_sets(test_name)
    ###
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
            print(f"[*] Saving data to: {save_path}...")
            save_data(data, save_path, custom_names)
            self.logger.log(f"[*] Data saved successfully.", level=1)   
    ###
    def scale_data(self, train_data, valid_data):
        """
        """
        self.scaler = Scaler()
        train_data['scaled_rbp_df'] = self.scaler.fit_transform(train_data['rbp_df'])
        valid_data['scaled_rbp_df'] = self.scaler.transform(valid_data['rbp_df'])
        return train_data, valid_data
    ###
    def create_datasets(self, train_data, valid_data):
        """
        Creates CustomTensorDataset instances for training and validation datasets.

        Initializes CustomTensorDataset objects for the provided training and validation,

        Args:
            train_data (DataFrame): The DataFrame containing the training data.
            valid_data (DataFrame): The DataFrame containing the validation data.

        Returns:
            tuple: A tuple containing:
                - train_dataset (CustomTensorDataset): The dataset for training data.
                - valid_dataset (CustomTensorDataset): The dataset for validation data.
        """
        # Create CustomTensorDataset for train and valid datasets
        train_dataset, valid_dataset = [
            CustomTensorDataset(
                data,
                self.getBM,
                rbp_data_key='scaled_rbp_df', 
                gene_data_key='gene_df', 
                transcript_data_key='isoform_df',
                trans_col_name=self.config.get('trans_col_name'),
                gene_col_name=self.config.get('gene_col_name')
            ) for data in [train_data, valid_data]
        ]
        return train_dataset, valid_dataset
    ###
    def create_data_loaders(self, train_dataset, valid_dataset):
        """
        Creates DataLoader instances for training and validation datasets.

        Args:
            train_dataset (CustomTensorDataset): Dataset for training data.
            valid_dataset (CustomTensorDataset): Dataset for validation data.

        Returns:
            tuple: A tuple containing:
                - train_loader (DataLoader): DataLoader for training data.
                - val_loader (DataLoader): DataLoader for validation data.
        """
        self.logger.log("🔄 Creating data loaders...", level=1)
        train_loader, val_loader = [
            DataLoader(
                dataset,
                batch_size=adjust_batch_size(dataset, (self.config.get('train_batch_size') if idx == 0 else self.config.get('val_batch_size'))),
                shuffle=(idx == 0),  # Shuffle only for the training set
                drop_last=(idx == 0)  # Drop last batch for training set only
            )
            for idx, dataset in enumerate([train_dataset, valid_dataset])
        ]
        self.logger.log("✅ Data loaders created.", level=1)
        return train_loader, val_loader
    ###
    def get_loaders(self, train_data, valid_data):
        """
        Encapsulates the creation of datasets and data loaders.

        Args:
            train_data (DataFrame): The DataFrame containing the training data.
            valid_data (DataFrame): The DataFrame containing the validation data.

        Returns:
            tuple: A tuple containing:
                - train_loader (DataLoader): DataLoader for training data.
                - val_loader (DataLoader): DataLoader for validation data.
        """
        self.train_dataset, valid_dataset = self.create_datasets(train_data, valid_data)
        return self.create_data_loaders(self.train_dataset, valid_dataset)
    ###
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
        model = PredictorModel(
            input_size=next(iter(train_loader))['scaled_rbp_df'].shape[1],
            output_size=next(iter(train_loader))['isoform_df'].shape[1],
            config=self.config
        ) # model instance doesn't need to make an attribute of the class as we can access to the model calling self.trainer.model
        # Print the model summary
        print("\n[*] Model Summary:")
        summarize_model(model, train_loader, self.config.get("train_batch_size")) # vigila que no haya problemas con el device
        self.trainer = TrainPredictor(
                model=model,
                config=self.config,
                input_features=('scaled_rbp_df', 'gene_df'), 
                output_features=('isoform_df',)
        )
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
    ###
    def save_history(self, train_history, val_history):
        self.logger.log("💾 Saving training and validation history...", level=1)
        #self.trainer.model.save_model(self.path_save_results, 'trained_model.pt')
        plot_loss_curve(train_history, val_history, output_dir=self.path_save_results)
        self.logger.log("✅ Saving training and validation history.", level=1)
    ###
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
        corr_per_gene = calculate_spearman_corr_per_gene(self.train_dataset, preds, labels, self.getBM)  
        # Create DataFrame for gene correlations
        gene_corr_df = pd.DataFrame({
            'gene_id': list(corr_per_gene['gene_corrs_dict'].keys()),  # Gene IDs
            'spearman_corr': list(corr_per_gene['gene_corrs_dict'].values()),  # Spearman correlations
            'max_trans_spearman_corr': list(corr_per_gene['gene_corrs_dict_max_trans'].values())  # Max transcript correlations
        })
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
        # Save metrics summary to CSV if needed
        if self.config.get('save_results'):
            metrics.to_csv(os.path.join(self.path_save_results, f'{set_name}_metrics_summary_global.csv'), index=True)
            gene_corr_df.to_csv(os.path.join(self.path_save_results, f'{set_name}_gene_correlations.csv'), index=False)
        self.logger.log(f"✅ Model evaluation on {set_name} complete.", level=1)
        return metrics, gene_corr_df
    ###
    # def eval_model_per_category(self, data, set_name): 
    #     self.logger.log(f"📊 Evaluating metrics per category ...", level=1)
    #     metrics_list, set_names_list = calculate_metrics_per_category(
    #          data,
    #          self.trainer,
    #          output_dir=self.path_save_results,
    #          set_name=set_name,
    #          sample_category=self.config.get('sample_category'),
    #          batch_size=self.config.get('train_batch_size'),
    #          source_name='TCGA',
    #          plot_results=self.config.get('plot_results'),
    #          getBM=self.getBM  
    #      )
    #     # Save metrics summary
    #     if self.config.get('save_results'):
    #         save_metrics_summary(
    #             set_names_list,
    #             metrics_list,
    #             output_path=f'{self.path_save_results}/{set_name}_metrics_summary_per_category.csv')
    #         self.logger.log(f"✅ Metrics for {set_name} evaluated.", level=1)
    ###
    def run(self): # obsolote, revisar
        self.logger.log("🚀 Starting the DeepRBPredictorPipeline run...", level=1)
        # Import data
        self.logger.log("📥 Importing training and test data...", level=1)
        data = self.import_data()
        # Split the filtered training data into train and validation sets
        train_data, valid_data = self.split_data(data)
        # Save the split datasets
        self.save_split_data(train_data, valid_data)   
        # Scale data
        train_data, valid_data = self.scale_data(train_data, valid_data)  
        # Create data loaders
        train_loader, valid_loader = self.get_loaders(train_data, valid_data)
        # Train the model
        train_history, val_history = self.train_model(train_loader, valid_loader)
        self.save_history(train_history, val_history)
        # Evaluate the model on training and validation sets
        for loader, name in zip([train_loader, valid_loader], ['train', 'val']):
            metrics, _ = self.evaluate_model(loader, set_name=name)
            # Optionally, you can log or save metrics here if needed
            self.logger.log(f"Metrics for {name}: {metrics}", level=1)   
        # Evaluate metrics per category (LA FUNCION INTERNA HAY QUE ACTUALIZAR)
        for dataset, name in zip([train_data, valid_data], ['train', 'val']):
            self.eval_model_per_category(dataset, name)  
        self.logger.log("✅ Pipeline run completed successfully.", level=1)


## PRUEBA AHORA EL CODIGO BRO Y MODIFICA LUEGO EL GRID SEARCH DE OPTUNA!!
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_train.yaml'
output_dir = '/scratch/jsanchoz/DeepRBP/output/results'

pipeline = DeepRBPredictorPipeline(config_path, output_dir)
# prueba 1: step-by-step
data = pipeline.import_data()
train_data, valid_data = pipeline.split_data(data)
#pipeline.save_split_data(train_data, valid_data) 
train_data, valid_data = pipeline.scale_data(train_data, valid_data) 
train_loader, valid_loader = pipeline.get_loaders(train_data, valid_data)
train_history, val_history = pipeline.train_model(train_loader, valid_loader)
pipeline.save_history(train_history, val_history)

# Evaluate the model on training and validation sets
for loader, name in zip([train_loader, valid_loader], ['train', 'val']):
    metrics, _ = pipeline.evaluate_model(loader, set_name=name)
    # Optionally, you can log or save metrics here if needed
    pipeline.logger.log(f"Metrics for {name}: {metrics}", level=1) 

# # Evaluate metrics per category
# for dataset, name in zip([train_data, valid_data], ['train', 'val']):
#     self.eval_model_per_category(dataset, name) #OLD VERSION
    
# self.eval_model_per_category(train_data, valid_data) #NEW
    



# prueba 2: direct run