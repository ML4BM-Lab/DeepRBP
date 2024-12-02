
# processing.py
import numpy as np
import pandas as pd
import joblib
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split as sk_train_test_split
from typing import Dict, Tuple, Optional, List

from utils import index2id, filter_data_by_sample_ids

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
        # Ensure metadata_path is provided
        if "metadata_path" not in paths or not paths["metadata_path"]:
            raise ValueError("The 'metadata_path' is required and must be provided.")

        self.paths = paths
        self.data = None

    def load(self):
        """
        Load data dictionary from the specified paths and store it in the object.
        """
        try:
            self.data = {
                    "rbp_expr_df": pd.read_csv(self.paths["rbp_path"], index_col=0), # Ensure user knows expected format
                    "gene_expr_df": pd.read_csv(self.paths["gene_expr_path"], index_col=0),
                    "trans_expr_df": pd.read_csv(self.paths["isoform_expr_path"], index_col=0),
                    "metadata_df": pd.read_csv(self.paths["metadata_path"], index_col=0),
            }
        
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Error loading file: {e}")
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
            raise ValueError(
                "Patient IDs are not aligned across datasets "
                "(RBP expression, gene expression, isoform expression, and metadata). "
                "Please ensure all datasets have the same index and order."
            )

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
        self.data = data
        self.config = config
        self.sample_category = sample_category
        self.sample_ids = self.data['metadata_df'].index
        self.id2index_mapping = self._generate_id2index_mapping()

    def _generate_id2index_mapping(self) -> Dict[str, int]:
        """Generate the mapping from patient IDs to indices."""
        return {patient_id: idx for idx, patient_id in enumerate(self.sample_ids)}
    
    def split_data(self, data, test_size):
        """Perform a stratified train-test split based on the detailed_category in metadata."""
        sample_category = data['metadata_df'][self.sample_category]
        print(f'[split_data] Performing a stratified data split with fraction division equal to {test_size}')
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

    def split_data_sets(self) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        This function handles the splitting of data into train, validation, and test sets based on the config.
        """
        train_data, valid_data, test_data = {}, {}, {}
        if self.config["training"]["train_test_split"]:
            print("Performing train/test split...")
            self.train_idx, self.test_idx = self.split_data(data=self.data, test_size=self.config['test_fraction'])
            print("Writing set_type in metadata after train/test split")
            self.add_sample_set_label(set_name='training', samples=index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='testing', samples=index2id(self.id2index_mapping, self.test_idx))
            train_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.train_idx))
            test_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.test_idx))

        if self.config["training"]["train_val_split"]:
            print("Performing train/val split...")
            training_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.train_idx)) 
            self.train_idx, self.valid_idx = self.split_data(data=training_data, test_size=self.config['val_fraction'])
            print("Writing set_type in metadata after train/val split")
            self.add_sample_set_label(set_name='training', samples=index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='validation', samples=index2id(self.id2index_mapping, self.valid_idx))
            valid_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.valid_idx))
        
        if not train_data and not valid_data and not test_data:
            print("No splits performed, returning only test data.")
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
        self.scaler = existing_scaler
        self.sigma = existing_sigma
        if self.scaler is not None and self.sigma is not None:
            print("[Scaler] Existing scaler and sigma loaded. Ready for transformation.")
        else:
            print("[Scaler] No existing scaler or sigma provided. Please fit before using.")

    def fit(self, train_set):
        """
        Fit a StandardScaler to the training set and compute sigma for clipping.

        Parameters:
        train_set (DataFrame): The training dataset used to fit the scaler.
        """
        print("[Scaler:fit] Fitting a StandardScaler to the training set...")
        self.scaler = StandardScaler()
        self.scaler.fit(train_set)
        self.sigma = np.std(self.scaler.transform(train_set).flatten().astype(np.float64))
        print(f"[Scaler:fit] Sigma computed: {self.sigma}")

    def transform(self, transform_set):
        """
        Normalize and clip the data using the fitted scaler and sigma.

        Parameters:
        transform_set (DataFrame): The dataset to be transformed.

        Returns:
        DataFrame: The normalized and clipped data.
        """
        if self.scaler is None or self.sigma is None:
            raise ValueError("[Scaler:transform] Scaler and sigma must be fitted or loaded before transformation.")
        # Scale the data
        scaled_set = pd.DataFrame(
            self.scaler.transform(transform_set),
            index=transform_set.index,
            columns=transform_set.columns
        )
        print(f"[Scaler:transform] Mean after scaling: {scaled_set.mean().mean():.4f}")
        print(f"[Scaler:transform] Std after scaling: {scaled_set.std().mean():.4f}")
        # Clip and normalize
        print("[Scaler:transform] Clipping the data...")
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
        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy')
        if self.scaler is not None:
            joblib.dump(self.scaler, scaler_file)
            print(f"[Scaler:save] Scaler saved to: {scaler_file}")
        else:
            raise ValueError("[Scaler:save] No scaler available to save.")
        if self.sigma is not None:
            np.save(sigma_file, np.array(self.sigma, dtype=np.float64))
            print(f"[Scaler:save] Sigma saved to: {sigma_file}")
        else:
            raise ValueError("[Scaler:save] No sigma available to save.")
        
    @classmethod
    def load(cls, folder_path):
        """
        Load the scaler and sigma from a specified folder.

        Parameters:
        folder_path (str): The directory where the scaler and sigma are stored.

        Returns:
        Scaler: A new instance of the Scaler class with loaded scaler and sigma.
        """
        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy')
        if not os.path.exists(scaler_file):
            raise FileNotFoundError(f"[Scaler:load] Scaler file not found at: {scaler_file}")
        if not os.path.exists(sigma_file):
            raise FileNotFoundError(f"[Scaler:load] Sigma file not found at: {sigma_file}")
        # Load scaler and sigma
        scaler = joblib.load(scaler_file)
        sigma = np.load(sigma_file)
        print(f"[Scaler:load] Scaler loaded from: {scaler_file}")
        print(f"[Scaler:load] Sigma loaded from: {sigma_file}")
        return cls(existing_scaler=scaler, existing_sigma=sigma)
    
    # You can load the scaler later from disk if needed
    # scaler = Scaler.load(folder_path=f"{scaler_path}/scaler")
    
   