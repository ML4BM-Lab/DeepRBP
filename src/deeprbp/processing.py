
# processing.py
import numpy as np
import pandas as pd
import joblib
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split as sk_train_test_split
from typing import Dict, Tuple, Optional, List
from .utils import *
from .logger import Logger

# Ejemplos de uso:
# nosotros asumimos que el usuario su data lo ha transformado de la siguiente manera antes de importarlo en esta clase:
    # rbp --> log2(TPM+1)
    # transcripts --> log2(TPM+1)
    # genes --> TPM

# a) Clase de configuración de datos: Para separar la carga de datos y la preparación de un conjunto de datos:
class DataImporter:
    def __init__(self, paths):
        """
        Initialize the DataImporter with file paths.

        Parameters:
        - paths (dict): Dictionary with paths to the RBP, gene, isoform expression datasets, and metadata.
        """
        self.logger = Logger(verbose=1)
        if "metadata_path" not in paths or not paths["metadata_path"]:
            self.logger.error("The 'metadata_path' is required and must be provided.", ValueError)

        self.paths = paths
        self.data = None

    def load(self):
        """
        Load data dictionary from the specified paths and store it in the object.
        """
        try:
            self.logger.log("Loading data from specified paths...")
            self.data = {
                    "rbp_expr_df": pd.read_csv(self.paths["rbp_path"], index_col=0), # Ensure user knows expected format
                    "gene_expr_df": pd.read_csv(self.paths["gene_expr_path"], index_col=0),
                    "trans_expr_df": pd.read_csv(self.paths["isoform_expr_path"], index_col=0),
                    "metadata_df": pd.read_csv(self.paths["metadata_path"], index_col=0),
            }
        
        except FileNotFoundError as e:
            self.logger.error(f"Error loading file: {e}", FileNotFoundError)

        self._check_consistency()  # Verify the indices are consistent
        return self.data
    
    def _check_consistency(self):
        """
        Ensure the indices are consistent across all datasets.
        """
        rbp_ids = self.data["rbp_expr_df"].index
        gene_ids = self.data["gene_expr_df"].index
        trans_ids = self.data["trans_expr_df"].index
        metadata_ids = self.data["metadata_df"].index

        if not (rbp_ids.equals(gene_ids) and gene_ids.equals(trans_ids) and trans_ids.equals(metadata_ids)):
            self.logger.error("Patient IDs are not aligned across datasets.", ValueError)

class DatasetLoader:
    def __init__(self, data_importer, base_config):
        self.data_importer = data_importer
        self.base_config = base_config

    def load_data(self):
        # Cargar los datos utilizando DataImporter
        data = self.data_importer.load()
        
        # Filtrar las muestras basadas en la configuración proporcionada
        sel_sample_ids = select_sample_ids_by_type(
            metadata_df=data['metadata_df'],
            sample_category_col=self.base_config['sample_category'],
            sample_types=self.base_config['select_samples']
        )
        
        # Filtrar los datos por las muestras seleccionadas
        return filter_data_by_sample_ids(data, sel_sample_ids)

# b)  Clase de división de datos: Para separar la lógica de la división de los datos de la clase principal
class DataSplitter:
    def __init__(self, data: Dict[str, pd.DataFrame], config: Dict, sample_category: str = "detailed_category"):
        """
        Initialize the DataSplitter with the dataset and configuration.

        Parameters:
        - data: Dictionary containing the dataset, including metadata_df.
        - config: Dictionary with configuration settings, including test/train split parameters.
        - sample_category: The category used for stratification (default 'detailed_category').
        """
        self.logger = Logger(verbose=1)
        self.data = data
        self.config = config
        self.sample_category = sample_category
        self.sample_ids = self.data['metadata_df'].index
        self.id2index_mapping = self._generate_id2index_mapping()

    def _generate_id2index_mapping(self) -> Dict[str, int]:
        """Generate the mapping from patient IDs to indices."""
        self.logger.log("Generating ID to index mapping.")
        return {patient_id: idx for idx, patient_id in enumerate(self.sample_ids)}
    
    def split_data(self, data, test_size):
        """Perform a stratified train-test split based on the detailed_category in metadata."""
        sample_category = data['metadata_df'][self.sample_category]
        self.logger.log(f'[split_data] Performing a stratified data split with fraction division equal to {test_size}')
        
        # Step 1: Perform the stratified split based on patient IDs (still strings at this point) and sample category
        train_idx, test_idx = sk_train_test_split(
            data['metadata_df'].index,
            test_size=test_size,
            stratify=sample_category,
            random_state=self.config['seed']
        )
        # Step 2: Convert string indices (train_idx, test_idx) to integer indices
        train_idx = [self.id2index_mapping[patient_id] for patient_id in train_idx]
        test_idx = [self.id2index_mapping[patient_id] for patient_id in test_idx]
        return train_idx, test_idx
    
    def add_sample_set_label(self, set_name: str, samples: list):
        """ Add a column to the metadata DataFrame indicating whether samples belong to the training or test set."""
        sample_set = set(samples)
        if 'set_type' not in self.data['metadata_df'].columns:
            self.data['metadata_df']['set_type'] = 'unknown'
        self.data['metadata_df'].loc[self.data['metadata_df'].index.isin(sample_set), 'set_type'] = set_name
        self.logger.log(f"Added sample set label '{set_name}' for {len(samples)} samples.")

    def split_data_sets(self) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        This function handles the splitting of data into train, validation, and test sets based on the config.
        """
        train_data, valid_data, test_data = {}, {}, {}
        if self.config["train_test_split"]:
            self.logger.log("Performing train/test split...")
            self.train_idx, self.test_idx = self.split_data(data=self.data, test_size=self.config['test_fraction'])
            self.logger.log("Writing set_type in metadata after train/test split")
            self.add_sample_set_label(set_name='training', samples=index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='testing', samples=index2id(self.id2index_mapping, self.test_idx))
            train_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.train_idx))
            test_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.test_idx))

        if self.config["train_val_split"]:
            self.logger.log("Performing train/val split...")
            training_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.train_idx)) 
            self.train_idx, self.valid_idx = self.split_data(data=training_data, test_size=self.config['val_fraction'])
            self.logger.log("Writing set_type in metadata after train/val split")
            self.add_sample_set_label(set_name='training', samples=index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='validation', samples=index2id(self.id2index_mapping, self.valid_idx))
            valid_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.valid_idx))
        
        if not train_data and not valid_data and not test_data:
            self.logger.warn("No splits performed, returning only test data.")
            test_data = self.data  # Return all data as test data
        return train_data, valid_data, test_data

# c) Clase de escalado
class Scaler:
    def __init__(self, existing_scaler=None, existing_sigma=None):
        """
        Initialize the Scaler class.

        Parameters:
        existing_scaler (object, optional): An existing scaler to be used if fit_scaler is False.
        existing_sigma (float, optional): An existing sigma value to be used if fit_scaler is False.
        """
        self.logger = Logger(verbose=1)
        self.scaler = existing_scaler
        self.sigma = existing_sigma
        if self.scaler is not None and self.sigma is not None:
            self.logger.log("[Scaler] Existing scaler and sigma loaded. Ready for transformation.")
        else:
            self.logger.warn("[Scaler] No existing scaler or sigma provided. Please fit before using.")

    def fit(self, train_set):
        """
        Fit a StandardScaler to the training set and compute sigma for clipping.

        Parameters:
        train_set (DataFrame): The training dataset used to fit the scaler.
        """
        self.logger.log("[Scaler:fit] Fitting a StandardScaler to the training set...")
        self.scaler = StandardScaler()
        self.scaler.fit(train_set)
        self.sigma = np.std(self.scaler.transform(train_set).flatten().astype(np.float64))
        self.logger.log(f"[Scaler:fit] Sigma computed: {self.sigma}")

    def transform(self, transform_set):
        """
        Normalize and clip the data using the fitted scaler and sigma.

        Parameters:
        transform_set (DataFrame): The dataset to be transformed.

        Returns:
        DataFrame: The normalized and clipped data.
        """
        if self.scaler is None or self.sigma is None:
            self.logger.error("[Scaler:transform] Scaler and sigma must be fitted or loaded before transformation.", ValueError)
        
        # Scale the data
        scaled_set = pd.DataFrame(
            self.scaler.transform(transform_set),
            index=transform_set.index,
            columns=transform_set.columns
        )
        
        self.logger.log(f"[Scaler:transform] Mean after scaling: {scaled_set.mean().mean():.4f}")
        self.logger.log(f"[Scaler:transform] Std after scaling: {scaled_set.std().mean():.4f}")
        
        # Clip and normalize
        self.logger.log("[Scaler:transform] Clipping the data...")
        scaled_set = np.clip(scaled_set, -2 * self.sigma, 2 * self.sigma, axis=1)
        scaled_set += 2 * self.sigma
        scaled_set /= 4 * self.sigma
        return scaled_set.astype(np.float64)
    
    def save(self, folder_path):
        """
        Save the scaler and sigma to specified folder.

        Parameters:
        folder_path (str): The directory where the scaler and sigma will be saved.
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            self.logger.log(f"[Scaler:save] Created directory: {folder_path}")

        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy')
        if self.scaler is not None:
            joblib.dump(self.scaler, scaler_file)
            self.logger.log(f"[Scaler:save] Scaler saved to: {scaler_file}")
        else:
            self.logger.error("[Scaler:save] No scaler available to save.", ValueError)
        if self.sigma is not None:
            np.save(sigma_file, np.array(self.sigma, dtype=np.float64))
            self.logger.log(f"[Scaler:save] Sigma saved to: {sigma_file}")
        else:
            self.logger.error("[Scaler:save] No sigma available to save.", ValueError)
        
    @classmethod
    def load(cls, folder_path):
        """
        Load the scaler and sigma from a specified folder.

        Parameters:
        folder_path (str): The directory where the scaler and sigma are stored.

        Returns:
        Scaler: A new instance of the Scaler class with loaded scaler and sigma.
        """
        logger = Logger(verbose=1) 
        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy')
        if not os.path.exists(scaler_file):
            logger.error(f"[Scaler:load] Scaler file not found at: {scaler_file}", FileNotFoundError)
        if not os.path.exists(sigma_file):
            logger.error(f"[Scaler:load] Sigma file not found at: {sigma_file}", FileNotFoundError)
        # Load scaler and sigma
        scaler = joblib.load(scaler_file)
        sigma = np.load(sigma_file)
        logger.log(f"[Scaler:load] Scaler loaded from: {scaler_file}")
        logger.log(f"[Scaler:load] Sigma loaded from: {sigma_file}")
        return cls(existing_scaler=scaler, existing_sigma=sigma)
    
    # You can load the scaler later from disk if needed
    # scaler = Scaler.load(folder_path=f"{scaler_path}/scaler")
    
   