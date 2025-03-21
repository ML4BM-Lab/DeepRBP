# src/deeprbp/data_loading/data_loader.py

import numpy as np
import pandas as pd
import joblib
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split as sk_train_test_split
from typing import Dict, Tuple

from ..util.logger import Logger
from ..util.utils import filter_data_by_sample_ids, index2id

# new
class DataImporter:
    def __init__(self, paths):
        """
        This class is designed to facilitate the loading of multiple datasets related to RNA-binding proteins (RBP),
        gene expression, isoform expression, and associated metadata from specified file paths.

        Parameters:
        - paths (dict): A dictionary containing the file paths for the datasets. The expected keys include:
          - 'rbp_path': Path to the RBP expression data file.
          - 'gene_expr_path': Path to the gene expression data file.
          - 'isoform_expr_path': Path to the isoform expression data file.
          - 'counts_path': Path to the gene counts data file.
          - 'metadata_path': Path to the metadata file, which is required.
        """
        self.logger = Logger(verbose=1)

        if "metadata_path" not in paths or not paths["metadata_path"]:
            self.logger.error("❌ The 'metadata_path' is required and must be provided.", ValueError)
        self.paths = paths
        self.data = {}

    def load(self):
        """
        Load data dictionary from the specified paths and store it in the object.
        """
        try:
            self.logger.log("📥 Loading data from specified paths...")
            for key, path in self.paths.items():
                if path:  # Only load if the path is provided
                    self.data[f"{key.split('_')[0]}_df"] = pd.read_csv(path, index_col=0)
                    self.logger.log(f"✅ Loaded {key} data.")
        except FileNotFoundError as e:
            self.logger.error(f"❌ Error loading file: {e}", FileNotFoundError)
        return self.data  # Return raw loaded data

# old 
# class DataImporter:
#     def __init__(self, paths):
#         """
#         This class is designed to facilitate the loading of multiple datasets related to RNA-binding proteins (RBP),
#         gene expression, isoform expression, and associated metadata from specified file paths.
#         The main functionalities include:
#         - Loading data from CSV files into pandas DataFrames.
#         - Transforming expression data for further analysis.
#         - Ensuring consistency across different datasets by filtering based on common indices.

#         Parameters:
#         - paths (dict): A dictionary containing the file paths for the datasets. The expected keys include:
#           - 'rbp_path': Path to the RBP expression data file.
#           - 'gene_expr_path': Path to the gene expression data file.
#           - 'isoform_expr_path': Path to the isoform expression data file.
#           - 'counts_path': Path to the gene counts data file.
#           - 'metadata_path': Path to the metadata file, which is required.
#           - 'getBM_path': Path to the getBM data file (if applicable).
#         """
#         self.logger = Logger(verbose=1)
#         if "metadata_path" not in paths or not paths["metadata_path"]:
#             self.logger.error("❌ The 'metadata_path' is required and must be provided.", ValueError)

#         self.paths = paths
#         self.data = {}

#     def load(self):
#         """
#         Load data dictionary from the specified paths and store it in the object.
#         """
#         try:
#             self.logger.log("📥 Loading data from specified paths...")

#             # Load data conditionally based on provided paths
#             if "rbp_path" in self.paths:
#                 self.data["rbp_expr_tpm_df"] = pd.read_csv(self.paths["rbp_path"], index_col=0)
#                 self.logger.log("✅ Loaded RBP expression data.")

#             if "gene_expr_path" in self.paths:
#                 self.data["gene_expr_tpm_df"] = pd.read_csv(self.paths["gene_expr_path"], index_col=0)
#                 self.logger.log("✅ Loaded gene expression data.")

#             if "isoform_expr_path" in self.paths:
#                 self.data["trans_expr_tpm_df"] = pd.read_csv(self.paths["isoform_expr_path"], index_col=0)
#                 self.logger.log("✅ Loaded isoform expression data.")

#             if "counts_path" in self.paths:
#                 self.data["gn_counts_df"] = pd.read_csv(self.paths["counts_path"], index_col=0)
#                 self.logger.log("✅ Loaded gene counts data.")

#             if "metadata_path" in self.paths:
#                 self.data["metadata_df"] = pd.read_csv(self.paths["metadata_path"], index_col=0)
#                 self.logger.log("✅ Loaded metadata.")

#             if "getBM_path" in self.paths:
#                 self.getBM = pd.read_csv(self.paths["getBM_path"])
#                 self.logger.log("✅ Loaded getBM data.")

#         except FileNotFoundError as e:
#             self.logger.error(f"❌ Error loading file: {e}", FileNotFoundError)

#         self._transform_data()  # Apply transformations
#         self._expand_gene_matrix() # Expand gene matrix
#         self._check_consistency() # Verify the indices are consistent
#         return self.data
    
#     def _transform_data(self):
#         """
#         Transform RBP and isoform expression data to log2(TPM + 1).
#         Store the results in new keys in the data dictionary.
#         """
#         self.logger.log("🔄 Starting transformation of RBP and isoform expression data...")
#         if "rbp_expr_tpm_df" in self.data:
#             self.data["rbp_expr_log2p_tpm_df"] = np.log2(self.data["rbp_expr_tpm_df"] + 1)
#             self.logger.log("✅ Applied log2 transformation to RBP expression data.")

#         if "trans_expr_tpm_df" in self.data:
#             self.data["trans_expr_log2p_tpm_df"] = np.log2(self.data["trans_expr_tpm_df"] + 1)
#             self.logger.log("✅ Applied log2 transformation to isoform expression data.")

#     def _expand_gene_matrix(self):
#         """
#         Expand the gene expression matrix to match each isoform expression.
#         Store the results in 'gn_expr_each_iso_tpm_df'.
#         """
#         self.logger.log("🔄 Starting expansion of gene expression matrix to match isoform expression...")
#         if "gene_expr_tpm_df" in self.data and "trans_expr_tpm_df" in self.data:
#             list_transcripts = self.data["trans_expr_tpm_df"].columns.tolist()

#             # Get the Gene_IDs list mapped to each Transcript_ID
#             list_genes_mapped_to_trans = [
#                 self.getBM.loc[self.getBM['Transcript_ID'] == trans_id, 'Gene_ID'].values[0]
#                 for trans_id in list_transcripts
#             ]

#             # Expand the gene expression matrix
#             self.data["gn_expr_each_iso_tpm_df"] = self.data["gene_expr_tpm_df"].loc[:, list_genes_mapped_to_trans]
#             # Set column names for genes mapped to transcripts.
#             self.data["gn_expr_each_iso_tpm_df"].columns = list_transcripts
#             self.logger.log("✅ Expanded gene expression matrix to match isoform expression.")

#     def _check_consistency(self):
#         """
#         Ensure the indices are consistent across all datasets.
#         If not, retain only the common samples by intersecting the indices.
#         """
#         indices = {key: df.index for key, df in self.data.items()}
#         common_indices = set.intersection(*[set(index) for index in indices.values()])

#         for key, df in self.data.items():
#             if not df.index.isin(common_indices).all():
#                 self.data[key] = df.loc[common_indices]
#                 self.logger.log(f"Filtered {key} to keep only common samples.")

#         # Verify that all DataFrames now have the same index
#         first_index = next(iter(indices.values())).tolist()
#         if not all(df.index.tolist() == first_index for df in self.data.values()):
#             self.logger.error("❌ Error: Not all DataFrames have the same indices after filtering.", ValueError)

# class DatasetLoader:
#     """
#     This class is designed to facilitate the loading of datasets using a DataImporter instance,
#     and subsequently filter the loaded data based on specified configuration settings. The main functionalities include:
#     - Loading data from various sources, such as expression and metadata files.
#     - Filtering the loaded data samples based on user-defined criteria, ensuring that only relevant samples are retained 
#         for analysis.

#     Parameters:
#     - data_importer (DataImporter): An instance of the DataImporter class responsible for loading datasets.
#         This instance should be initialized with the appropriate file paths to the datasets.
    
#     - base_config (dict): A configuration dictionary that specifies the criteria for filtering the loaded datasets.
#         It should include:
#         - 'sample_category': The column name in the metadata DataFrame used for categorizing samples.
#         - 'select_samples': A list of sample types to be selected for analysis, based on the sample category.
#     """
#     def __init__(self, data_importer, base_config):
#         self.data_importer = data_importer
#         self.base_config = base_config

#     def load_data(self):
#         # Load the data using DataImporter
#         data = self.data_importer.load()
#         # Filter the samples based on the provided configuration
#         sel_sample_ids = select_sample_ids_by_type(
#             metadata_df=data['metadata_df'],
#             sample_category_col=self.base_config['sample_category'],
#             sample_types=self.base_config['select_samples']
#         )
#         # Filter the data by the selected samples
#         return filter_data_by_sample_ids(data, sel_sample_ids)

##### de aqui para arriba está sujeto a cambios.
class DataSplitter:
    def __init__(self, data: Dict[str, pd.DataFrame], config: Dict, sample_category: str = "detailed_category"):
        """
        Initialize the DataSplitter class for splitting datasets into training, validation, and test sets.

        This class is designed to facilitate the stratified splitting of a dataset based on specified categories in the metadata.
        The main functionalities include:
        - Performing stratified splits to ensure that each subset (training, validation, and test) 
          maintains the same distribution of categories as in the full dataset.
        - Assigning labels to samples to indicate their respective sets (training, validation, testing).
        - Providing mappings from patient IDs to indices for easier data handling during the splitting process.

        Parameters:
        - data (Dict[str, pd.DataFrame]): A dictionary containing the datasets, with 'metadata_df' being a 
          DataFrame that holds metadata relevant for stratification during data splitting.
          
        - config (Dict): A configuration dictionary that contains settings for data splitting. This includes 
          parameters for train/test and train/validation splitting, such as:
          - 'train_test_split': A boolean indicating if a train/test split should be performed.
          - 'train_val_split': A boolean indicating if a train/validation split should be performed.
          - 'test_fraction': The fraction of the dataset to be used as the test set.
          - 'val_fraction': The fraction of the training set to be used as the validation set.
          
        - sample_category (str, optional): The column name in 'metadata_df' used for stratification. 
          Default is 'detailed_category'. This determines how samples are divided to ensure proportional 
          representation of categories in each subset.
        """
        self.logger = Logger(verbose=1)
        self.data = data
        self.config = config
        self.sample_category = sample_category
        self.sample_ids = self.data['metadata_df'].index
        self.id2index_mapping = self._generate_id2index_mapping()

    def _generate_id2index_mapping(self) -> Dict[str, int]:
        """Generate the mapping from patient IDs to indices."""
        self.logger.log("🔄 Generating ID to index mapping...")
        return {patient_id: idx for idx, patient_id in enumerate(self.sample_ids)}
    
    def split_data(self, data, test_size):
        """Perform a stratified train-test split based on the detailed_category in metadata."""
        sample_category = data['metadata_df'][self.sample_category]
        self.logger.log(f'📊 [split_data] Performing a stratified data split with fraction division equal to {test_size}...')
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
        self.logger.log(f"🏷️ Added sample set label '{set_name}' for {len(samples)} samples.")

    def split_data_sets(self) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        This function handles the splitting of data into train, validation, and test sets based on the config.
        """
        train_data, valid_data, test_data = {}, {}, {}
        if self.config.get("train_test_split", False):
            self.logger.log("🔄 Performing train/test split...")
            self.train_idx, self.test_idx = self.split_data(data=self.data, test_size=self.config['test_fraction'])
            self.logger.log("📝 Writing set_type in metadata after train/test split")
            self.add_sample_set_label(set_name='training', samples=index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='testing', samples=index2id(self.id2index_mapping, self.test_idx))
            train_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.train_idx))
            test_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.test_idx))

        if self.config.get("train_val_split", False): 
            self.logger.log("🔄 Performing train/val split...")
            training_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.train_idx)) 
            self.train_idx, self.valid_idx = self.split_data(data=training_data, test_size=self.config['val_fraction'])
            self.logger.log("📝 Writing set_type in metadata after train/val split")
            self.add_sample_set_label(set_name='training', samples=index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='validation', samples=index2id(self.id2index_mapping, self.valid_idx))
            valid_data = filter_data_by_sample_ids(self.data, index2id(self.id2index_mapping, self.valid_idx))
            
        if not train_data and not valid_data and not test_data:
            self.logger.warn("⚠️ No splits performed, returning only test data.")
            test_data = self.data  # Return all data as test data
        return train_data, valid_data, test_data

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
        else:
            self.logger.warn("⚠️ [Scaler] No existing scaler or sigma provided. Please fit before using.")

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
        self.logger.log(f"✅ [Scaler:fit] Sigma computed: {self.sigma:.4f}")

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
        return cls(existing_scaler=scaler, existing_sigma=sigma)
    
    # You can load the scaler later from disk if needed
    # scaler = Scaler.load(folder_path=f"{scaler_path}/scaler")
