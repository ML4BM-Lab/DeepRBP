# src/deeprbp/data_loading/data_loader.py

import numpy as np
import pandas as pd
import joblib
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split as sk_train_test_split
from typing import Dict, Tuple

from ..util.logger import Logger
from ..util.utils import print_section_separator, filter_data_by_sample_ids

class DataImporter:
    def __init__(self, base_path: str, verbose: int = 1, 
                 sample_category: str = None, select_category: str = None, 
                 disease_condition: str = None, select_condition: list = None,
                 sample_fraction: float = None):
        """
        This class is designed to facilitate the loading of multiple datasets related to RNA-binding proteins (RBP),
        gene expression, isoform expression, and associated metadata from a specified base path.

        Parameters:
        - base_path (str): Base path where the data files are located (e.g., train or test).
        - verbose (int): Verbosity level for logging. Default is 1, which enables logging.
                         Higher values may provide more detailed logging output.
        - sample_category (str): Column name used for filtering samples based on category. 
        - select_category (str): Specific category of samples to retain.
        - disease_condition (str): Column name used for filtering samples based on disease conditions.
        - select_condition (list): List of conditions to retain in the filtered dataset.
        - sample_fraction (float): Fraction of samples to retain for dataset size reduction.
        """
        self.logger = Logger(verbose)
        self.base_path = base_path
        self.paths = self._generate_data_patterns()
        self.sample_category = sample_category
        self.select_category = select_category
        self.disease_condition = disease_condition
        self.select_condition = select_condition
        self.sample_fraction = sample_fraction
    ###
    def _generate_data_patterns(self) -> Dict[str, str]:
        """Generates the file paths based on the base path provided."""
        return {
            'rbp_path': os.path.join(self.base_path, f"RBPs_log2p_tpm.csv"),
            'isoform_expr_path': os.path.join(self.base_path, f"trans_log2p_tpm.csv"),
            'gene_expr_path': os.path.join(self.base_path, f"gn_tpm.csv"),
            'metadata_path': os.path.join(self.base_path, f"phenotype_metadata.csv")
        }
    ###
    def filter_data_by_category_and_condition(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        """Filters the dataset based on the specified sample category and disease condition(s).

        Args:
            data (Dict[str, pd.DataFrame]): The dictionary containing expression data and metadata dataframes
            
        Returns:
            Dict[str, pd.DataFrame]: The filtered dataset based on the specified category and condition.
        """
        self.logger.log(f"\n[*] Filtering samples by category {self.select_category} and disease conditon {self.select_condition}", level=1)
        selected_sample_ids = data['metadata_df'][
            (data['metadata_df'][self.sample_category] == self.select_category) &
            (data['metadata_df'][self.disease_condition].isin(self.select_condition))
        ].index.tolist()
        data = filter_data_by_sample_ids(data, selected_sample_ids)
        self.logger.log("Data filtered by sample category and disease condition successfully.", level=1)
        return data
    ###
    def reduce_dataset_size(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        """Filters a portion of the samples to optimize time and computational resources.
        
        Args:
            data (Dict[str, pd.DataFrame]): The dictionary containing expression data and metadata dataframes.

        Returns:
            Dict[str, pd.DataFrame]: Filtered dictionary of DataFrames.

        """
        self.logger.log(f"\n[*] 🧩 Filtering samples to optimize time and resources with fraction: {self.sample_fraction}...", level=1)
        data_subset = DataSplitter.split_data_class(data, self.sample_category, self.sample_fraction, self.logger.verbose)
        self.logger.log("[*] Samples filtered successfully.\n", level=self.logger.verbose)
        return data_subset
    ###
    def load(self) -> Dict[str, pd.DataFrame]:
        """Load data dictionary from the specified paths and store it in the object.
            If sample_category and select_category are provided, the data will be filtered accordingly.
            If sample_fraction is provided, the dataset size will be reduced.
        """
        try:
            self.logger.log("📥 Loading data from specified paths...")
            data = {}
            for key, path in self.paths.items():
                if os.path.exists(path):  # Check if the file exists before loading
                    data[f"{key.split('_')[0]}_df"] = pd.read_csv(path, index_col=0)
                    self.logger.log(f"✅ Loaded {key} data from {path}.")
                else:
                    self.logger.error(f"❌ File not found: {path}")
            # Optionally filter the data if filtering parameters are provided
            if self.sample_category and self.select_category and self.disease_condition and self.select_condition:
                data = self.filter_data_by_category_and_condition(data)
            # Optionally reduce dataset size if sample_fraction is provided
            if self.sample_fraction is not None:
                data = self.reduce_dataset_size(data)
            #print_section_separator()
            return data  # Return raw loaded data
        except Exception as e:
            self.logger.error(f"❌ Error loading files: {e}", level=1)
            raise
    
###
# path = '/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Train'

# # Ejemplo 1: Carga Completa Sin Filtrado
# data_importer = DataImporter(base_path=path, verbose=1)
# data = data_importer.load()
# print("Datos cargados:", data.keys())

# # Ejemplo 2: Carga Con Filtrado
# data_importer = DataImporter(
#     base_path=path,
#     verbose=1,
#     sample_category="detailed_category",
#     select_category="Liver_Hepatocellular_Carcinoma",
#     disease_condition="sample_type",
#     select_condition=["Primary_Tumor", "Solid_Tissue_Normal"]
# )
# filtered_data = data_importer.load()
# print("Datos filtrados:", filtered_data.keys())

# # Ejemplo 3: Carga Completa Con Reducción de Tamaño
# data_importer = DataImporter(
#     base_path=path,
#     verbose=1,
#     sample_category="detailed_category",
#     sample_fraction=0.1  # Reducir el tamaño al 10%
# )
# reduced_data = data_importer.load()
# print("Datos cargados y reducidos:", reduced_data.keys())


class DataSplitter:
    def __init__(self, data: Dict[str, pd.DataFrame], sample_category: str = "detailed_category", verbose: int = 1):
        """
        Initialize the DataSplitter class for splitting datasets into training, validation/test sets.

        This class is designed to facilitate the stratified splitting of a dataset based on specified categories in the metadata.
        The main functionalities include:
        - Performing stratified splits to ensure that each subset (training, validation, and test) 
          maintains the same distribution of categories as in the full dataset.
        - Assigning labels to samples to indicate their respective sets (training, validation or testing).
    
        Parameters:
        - data (Dict[str, pd.DataFrame]): A dictionary containing the datasets, where 'metadata_df' is a 
          DataFrame that holds metadata relevant for stratification during data splitting.
        - sample_category (str, optional): The column name in 'metadata_df' used for stratification. 
          Default is 'detailed_category'. This determines how samples are divided to ensure proportional 
          representation of categories in each subset.
        - verbose (int): Verbosity level for logging. Default is 1, which enables logging.
                         Higher values may provide more detailed logging output.
        """
        self.verbose = verbose
        self.logger = Logger(self.verbose)
        self.data = data
        self.sample_ids_index = self.data['metadata_df'].index # Index of sample IDs
        self.sample_category_labels = self.data['metadata_df'][sample_category] # Series of sample categories
    def split_data(self, test_size: float) -> Tuple[list, list]:
        """Perform a stratified train-test split based on the detailed_category in metadata.
        Parameters:
        - test_size (float): The fraction of the dataset to be used as the test set.

        Returns:
        Tuple[list, list]: Two lists containing the training and testing sample IDs.
        """
        self.logger.log(f'📊 [split_data] Performing a stratified data split with fraction division equal to {test_size}...')
        # Perform the stratified split based on patient IDs (still strings at this point) and sample category
        train_id, test_id = sk_train_test_split(self.sample_ids_index, test_size=test_size, stratify=self.sample_category_labels, random_state=42)
        return train_id, test_id
    @classmethod
    def split_data_class(cls, data: Dict[str, pd.DataFrame], sample_category: str, 
                         test_size: float, verbose: int) -> Dict[str, pd.DataFrame]:
        """Class method to split data data in two and just return the remained portion, instantiating the class and calling the instance method.
        It is intended to be used for creating a toy dataset to optimize the computational resources.

        Parameters:
        - data (Dict[str, pd.DataFrame]): A dictionary containing the datasets.
        - sample_category (str): The column name in 'metadata_df' used for stratification.
        - test_size (float): The fraction of the dataset to be used as the test set.
        - verbose (int): Verbosity level for logging.

        Returns:
        Dict[str, pd.DataFrame]: A dictionary containing the subset of the original data corresponding to the training samples for toy dataset.
        """
        instance = cls(data, sample_category, verbose)
        _, subset_idx2 = instance.split_data(test_size)
        data_subset = filter_data_by_sample_ids(data, subset_idx2)
        return data_subset
    def add_sample_set_label(self, set_name: str, samples: list):
        """Add a column to the metadata DataFrame indicating whether samples belong to the training, validation or test set.

        Parameters:
        - set_name (str): The name of the set (e.g., 'training', 'validation' or 'testing').
        - samples (list): The list of sample IDs to label.
        """
        sample_set = set(samples)
        if 'set_type' not in self.data['metadata_df'].columns:
            self.data['metadata_df']['set_type'] = 'unknown'
        self.data['metadata_df'].loc[self.sample_ids_index.isin(sample_set), 'set_type'] = set_name
        self.logger.log(f"🏷️  Added sample set label '{set_name}' for {len(samples)} samples.")
        #print_section_separator()
    def split_data_sets(self, test_fraction, test_name='testing') -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        This function handles the splitting of data into training and testing (or 'validation') sets based on the 
        provided test fraction.
        
        Parameters:
        - test_fraction (float): The fraction of the dataset to be used for the test set.
        - test_name (str, optional): The name to use for the test set label (default is 'testing').

        Returns:
        Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]: Two dictionaries containing the training and testing datasets.
        """
        train_data, test_data = {}, {}
        # Determine the type of split and log the corresponding message
        if test_name == 'testing':
            self.logger.log("🔄 Performing train/test split...")
            log_message = "📝 Writing set_type in metadata after train/test split"
        else:
            self.logger.log("🔄 Performing train/val split...")
            log_message = "📝 Writing set_type in metadata after train/val split"
        # Perform the data splitting
        self.train_id, self.test_id = self.split_data(test_size=test_fraction)
        # Log the corresponding message in a single line
        self.logger.log(log_message)
        # Add labels to the sample sets
        self.add_sample_set_label(set_name='training', samples=self.train_id)
        self.add_sample_set_label(set_name=test_name, samples=self.test_id)
        # Filter the data by sample IDs
        train_data = filter_data_by_sample_ids(self.data, self.train_id)
        test_data = filter_data_by_sample_ids(self.data, self.test_id)
        self.logger.log("✅ Data splitting completed.", level=self.verbose)
        print_section_separator()
        return train_data, test_data
    
class Scaler:
    def __init__(self, existing_scaler=None, existing_sigma=None):
        """
        This class provides functionality for fitting a StandardScaler to a training dataset,
        transforming datasets by normalizing their values, and clipping them to specified bounds.

        Parameters:
        existing_scaler (object, optional): An existing scaler to be used if fit_scaler is False.
        existing_sigma (float, optional): An existing sigma value to be used if fit_scaler is False.
        """
        self.logger = Logger(verbose=1)
        self.scaler = existing_scaler
        self.sigma = existing_sigma
        if self.scaler is not None and self.sigma is not None:
            self.logger.log("✅ [Scaler] Existing scaler and sigma loaded. Ready for transformation.")
        #print_section_separator()
    def fit(self, train_set):
        """
        Fit a StandardScaler to the training dataset and compute the standard deviation (sigma) for clipping.

        This method calculates the mean and standard deviation of the training dataset, which is used 
        for normalizing and clipping the data during transformation.

        Parameters:
        train_set (DataFrame): The training dataset used to fit the scaler.
        """
        self.logger.log("🔄 [Scaler:fit] Fitting a StandardScaler to the training set...")
        self.scaler = StandardScaler()
        self.scaler.fit(train_set)
        self.sigma = np.std(self.scaler.transform(train_set).flatten().astype(np.float64))
        self.logger.log(f"✅ [Scaler:fit] Sigma computed: {self.sigma:.4f}\n")
    def transform(self, transform_set):
        """
        Normalize and clip the data using the fitted scaler and computed sigma.

        This method scales the input dataset using the fitted scaler, then clips the scaled values to 
        the range defined by +/- 2 sigma, and normalizes them.

        Parameters:
        - transform_set (DataFrame): The dataset (pandas DataFrame) to be transformed.

        Returns:
        - DataFrame: The normalized and clipped data, returned as a pandas DataFrame.
        """
        if self.scaler is None or self.sigma is None:
            self.logger.error("❌ [Scaler:transform] Scaler and sigma must be fitted or loaded before transformation.", ValueError)
        # Scale the data
        scaled_set = pd.DataFrame(
            self.scaler.transform(transform_set),
            index=transform_set.index,
            columns=transform_set.columns
        )
        print('\n')
        self.logger.log(f"[Scaler:transform] Mean after scaling: {scaled_set.mean().mean():.4f}")
        self.logger.log(f"[Scaler:transform] Std after scaling: {scaled_set.std().mean():.4f}")
        # Clip and normalize
        self.logger.log("🔄 [Scaler:transform] Clipping the data...")
        scaled_set = np.clip(scaled_set, -2 * self.sigma, 2 * self.sigma, axis=1)
        scaled_set += 2 * self.sigma
        scaled_set /= 4 * self.sigma
        return scaled_set.astype(np.float64)
    def save(self, folder_path):
        """
        Save the fitted scaler and the computed sigma to the specified directory.

        This method stores the scaler as a joblib file and the sigma value as a NumPy array in the specified folder.

        Parameters:
        - folder_path (str): The directory where the scaler and sigma will be saved.
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            self.logger.log(f"📁 [Scaler:save] Created directory: {folder_path}")
        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy')
        if self.scaler is not None:
            joblib.dump(self.scaler, scaler_file)
            self.logger.log(f"✅ [Scaler:save] Scaler saved to: {scaler_file}")
        else:
            self.logger.error("❌ [Scaler:save] No scaler available to save.", ValueError)
        if self.sigma is not None:
            np.save(sigma_file, np.array(self.sigma, dtype=np.float64))
            self.logger.log(f"✅ [Scaler:save] Sigma saved to: {sigma_file}")
        else:
            self.logger.error("❌ [Scaler:save] No sigma available to save.", ValueError)
        #print_section_separator()
    @classmethod
    def load(cls, folder_path):
        """
        Load the scaler and sigma from the specified directory.

        This method retrieves the scaler and sigma values from the specified folder, 
        allowing for reuse of previously fitted scaling parameters.

        Parameters:
        - folder_path (str): The directory where the scaler and sigma are stored.

        Returns:
        - Scaler: A new instance of the Scaler class with the loaded scaler and sigma.
        """
        logger = Logger(verbose=1) 
        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy')
        if not os.path.exists(scaler_file):
            logger.error(f"❌ [Scaler:load] Scaler file not found at: {scaler_file}", FileNotFoundError)
        if not os.path.exists(sigma_file):
            logger.error(f"❌ [Scaler:load] Sigma file not found at: {sigma_file}", FileNotFoundError)
        # Load scaler and sigma
        scaler = joblib.load(scaler_file)
        sigma = np.load(sigma_file)
        logger.log(f"✅ [Scaler:load] Scaler loaded from: {scaler_file}")
        logger.log(f"✅ [Scaler:load] Sigma loaded from: {sigma_file}")
        #print_section_separator()
        return cls(existing_scaler=scaler, existing_sigma=sigma)
    
    # You can load the scaler later from disk if needed
    # scaler = Scaler.load(folder_path=f"{scaler_path}/scaler")
