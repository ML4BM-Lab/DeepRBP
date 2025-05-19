
from ..data_loading.data_loader import DataImporter, DataSplitter, Scaler

import os
import torch
from torch.utils.data import Dataset, DataLoader
import pandas as pd

from typing import Any
from ..util.logger import Logger
from ..util.utils import (
    filter_data_by_sample_ids,
    save_data,
    adjust_batch_size
)

class myTensorDataset(Dataset):
    def __init__(self, 
                 data: dict, 
                 getBM: pd.DataFrame, 
                 rbp_data_key: str = 'scaled_rbp_df', 
                 gene_data_key: str = 'gene_df', 
                 transcript_data_key: str = 'isoform_df',
                 trans_col_name: str = 'Transcript_ID',
                 gene_col_name: str = 'Gene_ID',
                 verbose: int = 0):
        """
        Custom dataset for loading RBP, gene, and transcript expressions.

        Parameters:
        - data: Dictionary containing the DataFrames for RBP, gene, and transcript expressions.
        - getBM: DataFrame containing Transcript_IDs and associated Gene_IDs.
        - rbp_data_key: The name of the DataFrame column for RBP features.
        - gene_data_key: The name of the DataFrame column for gene expression features.
        - transcript_data_key: The name of the DataFrame column for transcript expression features.
        - gene_col_name: The name of the column in getBM that contains the Gene IDs.
        - trans_col_name: The name of the column in getBM that contains the Transcript IDs.
        - verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during dataset initialization.
            - 0: No logging (suppress all output).
            - 1: Basic logging (show feature names and their shapes).
        """
        self.logger = Logger(verbose=verbose) 
        for key in [rbp_data_key, gene_data_key, transcript_data_key]:
            if key not in data:
                raise KeyError(f"{key} not found in data.")
        self.original_data = {key: df.copy() for key, df in data.items()} # Original data
        data_copy = {key: df.copy() for key, df in data.items()}
        # Automatically detect the feature names from the DataFrames
        self.rbp_names = data_copy[rbp_data_key].columns.tolist()
        self.gene_names = data_copy[gene_data_key].columns.tolist()
        self.trans_names = data_copy[transcript_data_key].columns.tolist()
        # Expand the gene matrix data to match the transcript shape data 
        genes_names_each_trans = getBM[getBM[trans_col_name].isin(self.trans_names)][gene_col_name]
        data_copy[gene_data_key] = data_copy[gene_data_key].loc[:, genes_names_each_trans]
        # Store tensors for each feature
        features_data = (rbp_data_key, gene_data_key, transcript_data_key)
        self.features = {feature: torch.tensor(data_copy[f"{feature}"].values, dtype=torch.float32) for feature in features_data} # CHANGED THIS TO FLOAT32 JOSEBA!
        # Log the features and their shapes for debugging
        self.logger.log("Features stored in the dataset:", level=1)
        for feature, tensor in self.features.items():
            self.logger.log(f"{feature}: {tensor} | Shape: {tensor.shape}", level=1)  # Log each feature tensor and its shape
    def __len__(self) -> int:
        return len(next(iter(self.features.values())))
    def __getitem__(self, idx: int):
        return {feature: values[idx] for feature, values in self.features.items()}
    

class PrepareData:
    def __init__(self, config: Any, getBM: pd.DataFrame, output_dir: str, verbose: int = 1, num_workers: int = 0):
        """
        Initializes the PrepareData class with configuration and output directory.

        Args:
            config (object): Class object with the data configuration.
            getBM (pd.DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.
            output_dir (str): Directory where output data will be saved.
            verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during initialization.
                - 0: No logging (suppress all output).
                - 1: Basic logging (show feature names and their shapes).
            num_workers (int, optional): Number of subprocesses to use for data loading. Default is 0.
        Example:
            # Example of initializing with a configuration file and output directory
            prep_data = PrepareData(config_path='path/to/config.ini', output_dir='output_directory/')
        """
        self.logger = Logger(verbose=verbose)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing Data Preparation...", level=1)
        ###
        self.config = config
        self.getBM = getBM
        self.data_paths = self.config.get('data_paths')
        self.logger.log(f"data_paths: {self.data_paths}...", level=1)
        self.sample_category = self.config.get('sample_category')
        self.sample_fraction = self.config.get('sample_fraction', default=None)
        self.trans_col_name=self.config.get('trans_col_name', default="Transcript_ID")
        self.gene_col_name=self.config.get('gene_col_name', default="Gene_ID")
        self.train_batch_size = self.config.get('train_batch_size', default=64)
        self.val_batch_size = self.config.get('val_batch_size', default=128)
        self.num_workers = num_workers
        # Initilizations
        self.path_save_data = os.path.join(output_dir, 'data')
        self.scaler = None
        self.logger.log("Initialization done ✅", level=1)
    ###
    def load_data(self) -> dict[pd.DataFrame]:
        """Loads data from specified paths using the DataImporter class.
        Returns:
            dict[pd.DataFrame]: The loaded data.

        Raises:
            Exception: If loading data fails.

        Example:
            # Load data and print the first few rows
            data = prep_data.load_data()
            print(data.head())
        """
        try:
            self.data_importer = DataImporter(self.data_paths)
            data = self.data_importer.load()    
            self.logger.log("✅ Data import complete.", level=1)
            return data
        except Exception as e:
            self.logger.error(f"❌ Failed to load data: {e}", level=1)
            raise
    ###
    # Alternative method: if you wanna reduce your dataset size you can use this
    def filter_samples(self, data):
        """ Filters a portion of the samples to optimize time and computational resources.
        Args:
            data (dict[pd.DataFrame]): Dictionary of DataFrames to filter based on the sample fraction.

        Returns:
            dict[pd.DataFrame]: Filtered dictionary of dataFrames if sample_fraction is defined in the config, else original data.

        Example:
            # Filter samples if sample fraction is defined
            filtered_data = prep_data.filter_samples(data)
        """
        if self.sample_fraction is not None:
            self.logger.log(f"[*] Filtering samples to optimize time and resources with fraction: {self.sample_fraction}...", level=1)
            subset_idx, _ = DataSplitter.split_data_class(
                data=data,
                config=self.config,
                sample_category=self.sample_category,
                test_size=self.sample_fraction
            )
            data_subset = filter_data_by_sample_ids(data, subset_idx)
            self.logger.log("[*] Samples filtered successfully.", level=1)
            return data_subset
        else:
            self.logger.log("[*] No filtering applied. sample_fraction is not defined.", level=1)
            return data
    ###
    def split_data(self, data: dict[pd.DataFrame], test_name: str = 'validation') -> tuple:
        """Splits data into training and validation (or test) sets.
            
        Args:
            data (dict[pd.DataFrame]): The data to be split.
            test_name (str, optional): Name for the test dataset.

        Returns:
            tuple: A tuple containing the training and validation dict[pd.DataFrame].

        Example:
            # Split the data into training and validation sets
            train_data, valid_data = prep_data.split_data(data)
        """
        self.logger.log("✂️ Splitting data into train and validation (or test) sets...", level=1)
        splitter = DataSplitter(data, self.config)
        return splitter.split_data_sets(test_name)
    ###
    def save_split_data(self, train_data: pd.DataFrame, valid_data: pd.DataFrame):
        """Saves the training and validation data to specified paths.
        Args:
            train_data (dict[pd.DataFrame]): Dictionary of DataFrames containing the training data.
            valid_data (dict[pd.DataFrame]): Dictionary of DataFrames containing the validation data.

        Raises:
            Exception: If saving data fails.

        Example:
            # Save the split data
            prep_data.save_split_data(train_data, valid_data)
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
    ###
    def fit_scaler(self, train_data: dict[pd.DataFrame]):
        """Fits the scaler to the training data.
        Args:
            train_data (dict[pd.DataFrame]): The training data used to fit the scaler.

        Example:
            # Fit the scaler on the training data
            prep_data.fit_scaler(train_data)
        """
        self.scaler = Scaler()
        self.scaler.fit(train_data['rbp_df'])
        self.logger.log("✅ Scaler has been fitted to training data.", level=1)
        self.scaler.save(self.path_save_data)
        self.logger.log(f"✅ Scaler has been saved in {self.path_save_data}.", level=1)
    ###
    # Alternative method: if you already have an scaler you can load it
    def load_scaler(self, folder_path):
        """Load an already trained scaler from the specified directory and assign it to self.scaler.
        Args:
            folder_path (str): The directory where the scaler and sigma are stored.

        Raises:
            FileNotFoundError: If the scaler cannot be loaded from the specified path.

        Example:
            # Load an existing scaler
            prep_data.load_scaler('path/to/scaler/')
        """
        try:
            self.scaler = Scaler.load(folder_path)
            self.logger.log("✅ Scaler successfully loaded from the specified path.", level=1)
        except FileNotFoundError as e:
            self.logger.error(f"❌ [Scaler:load_scaler] Failed to load scaler: {e}", level=1)
        except Exception as e:
            self.logger.error(f"❌ [Scaler:load_scaler] An unexpected error occurred: {e}", level=1)
    ###
    def scale_data(self, data: dict[pd.DataFrame]):
        """Transforms the provided dict[pd.DataFrame] in the RBP expression dataframe using the fitted scaler.

        Args:
            data (dict): Dictionary containing the DataFrame to be transformed.

        Returns:
            dict: A dictionary containing the original DataFrames plus the scaled DataFrame.

        Raises:
            ValueError: If the scaler is not fitted.

        Example:
            # Scale the training dataset
            scaled_train_data = prep_data.scale_data(train_data)
        """
        if self.scaler is None:
            self.logger.error("❌ Scaler is not fitted. Please fit the scaler before scaling data.")
            raise ValueError("Scaler is not fitted. Please fit the scaler before scaling data.")
        # Transform the provided data using the fitted scaler
        scaled_data = data.copy()  # Make a copy of the original data to avoid modifying it
        # Scale the RBP data and add it to the scaled_data dictionary
        scaled_data['scaled_rbp_df'] = self.scaler.transform(data['rbp_df'])  # Scale the RBP data
        return scaled_data
    ###
    def scale_train_val_data(self, train_data: dict[pd.DataFrame], val_data: dict[pd.DataFrame]) -> tuple:
        """Scales both training and validation datasets using the fitted scaler.

        Args:
            train_data (dict): Dictionary containing the training DataFrames to be transformed.
            val_data (dict): Dictionary containing the validation DataFrames to be transformed.

        Returns:
            tuple: A tuple containing the scaled training and validation DataFrames.

        Raises:
            ValueError: If the scaler is not fitted.

        Example:
            # Scale the training and validation datasets
            scaled_train_data, scaled_val_data = prep_data.scale_train_val_data(train_data, val_data)
        """
        scaled_train_data = self.scale_data(train_data)
        scaled_val_data = self.scale_data(val_data)
        return scaled_train_data, scaled_val_data
    ###
    def create_tensor_dataset(self, data: dict[pd.DataFrame]) -> myTensorDataset:
        """Creates a myTensorDataset instance from a provided dictionary of DataFrames.

        Args:
            data (dict): Dictionary containing the scaled DataFrames for RBP, gene, and transcript expressions.

        Returns:
            myTensorDataset: An instance of myTensorDataset.
        """
        dataset = myTensorDataset(
            data=data.copy(),
            getBM=self.getBM,
            rbp_data_key='scaled_rbp_df', 
            gene_data_key='gene_df', 
            transcript_data_key='isoform_df',
            trans_col_name=self.trans_col_name,
            gene_col_name=self.gene_col_name
        )
        return dataset
    ###
    def create_data_loader(self, dataset: myTensorDataset, batch_size: int, shuffle: bool = False, drop_last: bool = False) -> DataLoader:
        """Creates a DataLoader instance for the provided dataset.

        Args:
            dataset (myTensorDataset): The dataset for which to create a DataLoader.
            batch_size (int): The batch size for the DataLoader.
            shuffle (bool, optional): Whether to shuffle the data. Default is False.
            drop_last (bool, optional): Whether to drop the last batch if it is not full. Default is False.

        Returns:
            DataLoader: An instance of DataLoader.
        """
        self.logger.log("🔄 Creating data loader...", level=1)
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=shuffle,  # Shuffle the training data
            drop_last=drop_last,  # Drop the last batch if it is not full
            num_workers=self.num_workers
        )
        self.logger.log("✅ Data loader created.", level=1)
        return loader
    ###
    def create_train_val_datasets(self, scaled_train_data: dict[pd.DataFrame], scaled_val_data: dict[pd.DataFrame]) -> tuple:
        """Creates training and validation tensor datasets from scaled data.

        Args:
            scaled_train_data (dict): Dictionary containing the DataFrames for RBP, gene, and transcript expressions.
            scaled_val_data (dict): Dictionary containing the DataFrames for RBP, gene, and transcript expressions.

        Returns:
            tuple: A tuple containing the training and validation datasets.
        """
        train_dataset = self.create_tensor_dataset(scaled_train_data)
        valid_dataset = self.create_tensor_dataset(scaled_val_data)
        return train_dataset, valid_dataset
    def create_train_val_loaders(self, train_dataset: myTensorDataset, valid_dataset: myTensorDataset) -> tuple:
        """Creates DataLoaders for training and validation datasets.

        Args:
        train_dataset (myTensorDataset): The tensor dataset used for training,
        valid_dataset (myTensorDataset): The tensor dataset used for validation

        Returns:
            tuple: A tuple containing two elements:
                - train_loader (DataLoader): A DataLoader for the training dataset, configured with specified batch size, shuffling enabled, and drop_last set to True.
                - valid_loader (DataLoader): A DataLoader for the validation dataset, configured with specified batch size, shuffling disabled, and drop_last set to False.

        Example:
            # Create DataLoaders for the training and validation datasets
            train_loader, valid_loader = prep_data.create_train_val_loaders(train_dataset, valid_dataset)

        """
        train_loader = self.create_data_loader(train_dataset, batch_size = self.train_batch_size, shuffle = True, drop_last = True)
        valid_loader = self.create_data_loader(valid_dataset, batch_size = self.val_batch_size, shuffle = False, drop_last = False)
        return train_loader, valid_loader





 