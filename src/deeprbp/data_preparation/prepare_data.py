
from ..data_loading.data_loader import DataImporter, DataSplitter, Scaler

import os
import pandas as pd
from typing import Dict, Tuple
from torch.utils.data import DataLoader 

from .tensor_dataset import DeepRBPExpressionDataset
from ..util.logger import Logger
from ..util.utils import save_data

class PrepareData:
    def __init__(self, 
                 getBM: pd.DataFrame, 
                 trans_col_name: str, 
                 gene_col_name: str, 
                 output_dir: str = None, 
                 verbose: int = 1):
        """
        Initializes the PrepareData class for preparing and managing datasets.

        Args:
            getBM (pd.DataFrame): DataFrame containing the mapping of Gene_ID to Transcript_ID.
            trans_col_name (str): The name of the column in getBM representing Transcript_ID.
            gene_col_name (str): The name of the column in getBM representing Gene_ID.
            output_dir (str, optional): Directory where output data will be saved. Default is None.
            verbose (int, optional): Verbosity level for logging. 
                                     - 0: No logging (suppress all output).
                                     - 1: Basic logging (show feature names and their shapes).
            
        Example:
            # Example of initializing the PrepareData class
            prep_data = PrepareData(getBM=your_dataframe, trans_col_name='Transcript_ID', gene_col_name='Gene_ID', output_dir='output_directory/')
        """
        self.verbose = verbose
        self.logger = Logger(self.verbose)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing Data Preparation...", level=1)
        ###
        self.getBM = getBM
        self.trans_col_name = trans_col_name
        self.gene_col_name = gene_col_name
        ###
        # Optional output directory
        if output_dir:
            self.path_save_data = os.path.join(output_dir, 'data')
        else:
            self.path_save_data = None  # Handle case where output_dir is not provided
        self.scaler = None
        self.logger.log("Initialization done ✅", level=1)
    ###
    def load_data(self, path: str) -> Dict[str, pd.DataFrame]:  
        """Loads data from specified paths using the DataImporter class.
        
        Args:
            path (str): The path where data files are located.

        Returns:
            Dict[str, pd.DataFrame]: The loaded data as a dictionary of DataFrames.

        Raises:
            Exception: If loading data fails.

        Example:
            # Load data and print the first few rows
            data = prep_data.load_data(path='path/to/data')
            print(data['metadata_df'].head())
        """
        self.data_importer = DataImporter(path, self.verbose)
        data = self.data_importer.load()    
        self.logger.log("✅ Data import complete.", level=self.verbose)
        return data
    ###
    # Alternative method: if you wanna reduce your dataset size you can use this.
    def filter_samples(self, data: Dict[str, pd.DataFrame], sample_category: str, sample_fraction: float) -> Dict[str, pd.DataFrame]:
        """ Filters a portion of the samples to optimize time and computational resources.
        
        Args:
            data (Dict[str, pd.DataFrame]): Dictionary of DataFrames to filter based on the sample fraction.
            sample_category (str): The column name in metadata used for stratification.
            sample_fraction (float): The fraction of samples to retain.

        Returns:
            Dict[str, pd.DataFrame]: Filtered dictionary of DataFrames if sample_fraction is defined, else original data.

        Example:
            # Filter samples if sample fraction is defined
            filtered_data = prep_data.filter_samples(data, sample_category='detailed_category', sample_fraction=0.1)
        """
        self.logger.log(f"[*] Filtering samples to optimize time and resources with fraction: {sample_fraction}...", level=1)
        data_subset = DataSplitter.split_data_class(data, sample_category, sample_fraction, self.verbose)
        self.logger.log("[*] Samples filtered successfully.", level=self.verbose)
        return data_subset
    ###
    def split_data(self, data: Dict[str, pd.DataFrame], sample_category: str, 
                   test_fraction: float, test_name: str = 'validation') -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """Splits data into training and validation (or test) sets.
            
        Args:
            data (Dict[str, pd.DataFrame]): The data to be split.
            sample_category (str): The column name in metadata used for stratification.
            test_fraction (float): The fraction of the dataset to be used as the validation set.
            test_name (str, optional): Name for the validation dataset. Default is 'validation'.

        Returns:
            Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]: A tuple containing the training and validation datasets.

        Example:
            # Split the data into training and validation sets
            train_data, valid_data = prep_data.split_data(data, sample_category='detailed_category', test_fraction=0.2)
        """
        self.logger.log("✂️ Splitting data into train and validation (or test) sets...", level=self.verbose)
        splitter = DataSplitter(data, sample_category)
        return splitter.split_data_sets(test_fraction, test_name)
    ###
    def save_split_data(self, train_data: pd.DataFrame, valid_data: pd.DataFrame):
        """Saves the training and validation data to specified paths.
        Args:
            train_data (dict[pd.DataFrame]): Dictionary of DataFrames containing the training data.
            valid_data (dict[pd.DataFrame]): Dictionary of DataFrames containing the validation data.

        Raises:
            Exception: If output_dir is not set.

        Example:
            # Save the split data
            prep_data.save_split_data(train_data, valid_data)
        """
        if self.path_save_data is None:
            self.logger.error("❌ Output directory is not set. Cannot save data.")
            raise ValueError("Output directory is not defined. Set output_dir during initialization.")
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
            self.logger.log(f"[*] Saving data to: {save_path}...", level=self.verbose) 
            save_data(data, save_path, custom_names)
            self.logger.log(f"[*] Data saved successfully.", level=self.verbose) 
    ###
    def fit_scaler(self, train_data: dict[pd.DataFrame]):
        """Fits the scaler to the training data.
        Args:
            train_data (dict[pd.DataFrame]): The training data used to fit the scaler.

        Raises:
            ValueError: If output_dir is not defined.

        Example:
            # Fit the scaler on the training data
            prep_data.fit_scaler(train_data)
        """
        if self.path_save_data is None:
            self.logger.error("❌ Output directory is not set. Cannot save the scaler.")
            raise ValueError("Output directory is not defined. Set output_dir during initialization.")
        self.scaler = Scaler()
        self.scaler.fit(train_data['rbp_df'])
        self.logger.log("✅ Scaler has been fitted to training data.", level=self.verbose)
        self.scaler.save(self.path_save_data)
        self.logger.log(f"✅ Scaler has been saved in {self.path_save_data}.", level=self.verbose)
    ###
    # Alternative method: if you already have an scaler you can load it
    def load_scaler(self, folder_path: str):
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
            self.logger.log("✅ Scaler successfully loaded from the specified path.", level=self.verbose)
        except FileNotFoundError as e:
            self.logger.error(f"❌ [Scaler:load_scaler] Failed to load scaler: {e}", level=self.verbose)
        except Exception as e:
            self.logger.error(f"❌ [Scaler:load_scaler] An unexpected error occurred: {e}", level=self.verbose)
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
    def create_tensor_dataset(self, data: dict[pd.DataFrame]) -> DeepRBPExpressionDataset:
        """Creates a DeepRBPExpressionDataset instance from a provided dictionary of DataFrames.

        Args:
            data (dict): Dictionary containing the scaled DataFrames for RBP, gene, and transcript expressions.

        Returns:
            DeepRBPExpressionDataset: An instance of DeepRBPExpressionDataset.
        """
        self.logger.log("🔄 Creating tensor dataset ...", level=self.verbose)
        dataset = DeepRBPExpressionDataset(
            data=data.copy(),
            getBM=self.getBM,
            rbp_data_key='scaled_rbp_df', 
            gene_data_key='gene_df', 
            transcript_data_key='isoform_df',
            trans_col_name=self.trans_col_name,
            gene_col_name=self.gene_col_name,
            verbose=self.verbose
        )
        self.logger.log("✅ Tensor dataset created.", level=self.verbose)
        return dataset
    ###
    def create_data_loader(self, dataset: DeepRBPExpressionDataset, batch_size: int, shuffle: bool = False, drop_last: bool = False, 
                           num_workers: int = 0) -> DataLoader:
        """Creates a DataLoader instance for the provided dataset.

        Args:
            dataset (DeepRBPExpressionDataset): The dataset for which to create a DataLoader.
            batch_size (int): The batch size for the DataLoader.
            shuffle (bool, optional): Whether to shuffle the data. Default is False.
            drop_last (bool, optional): Whether to drop the last batch if it is not full. Default is False.
            num_workers (int, optional): Number of subprocesses to use for data loading. Default is 0.

        Returns:
            DataLoader: An instance of DataLoader.
        """
        self.logger.log("🔄 Creating data loader...", level=self.verbose)
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=shuffle,  # Shuffle the training data
            drop_last=drop_last,  # Drop the last batch if it is not full
            num_workers=num_workers
        )
        self.logger.log("✅ Data loader created.", level=self.verbose)
        return loader
