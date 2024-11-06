# El data está en : /scratch/jsanchoz/DeepRBP/data/training_module/processed
# La idea es que este loader sea capaz de cargarte las matrices de rbp, trans, y gene, que te divida en train, test y val
# usando un porcentaje de los tipos tumorales. Con las de TCGA puedes hacer tb toy version para optimización, GTEX que vaya todo para el test.
# El scaling tb sea uno de los procesos

# Faltaría como mucho ahora ya tener en cuenta esto de aquí del toy (14/10):
# #toy_set = False # with frac = 0.4 # esta funcionalidad podría no tener mucho sentido tenerla aqui y seria mejor 
# # en un arch_optimization.py usar simplemente la funcion de split_data 0.4%.

############################################################################################################
import os
import joblib
import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader, SubsetRandomSampler
from sklearn.model_selection import train_test_split as sk_train_test_split
from sklearn.preprocessing import StandardScaler
from typing import Dict, Optional, List

#from deeprbp.config_loader import get_config 

# Ejemplos de uso:
# nosotros asumimos que el usuario su data lo ha transformado de la siguiente manera antes de importarlo en esta clase:
    # rbp --> log2(TPM+1)
    # transcripts --> log2(TPM+1)
    # genes --> TPM

def create_dataloaders(dataset, batch_size: int, valid_batch_size: Optional[int] = None, pin_memory: bool = True):
    """
    Creates DataLoaders based on the partitions of the custom dataset. If the training, validation, 
    or test partitions are None, the corresponding DataLoader will not be created.

    Args:
        dataset (CustomDataset): The custom dataset that contains the indices of the partitions.
        batch_size (int): Batch size for the training and test loaders.
        valid_batch_size (int, optional): Batch size for the validation loader. If not specified,
                                          it defaults to double the batch_size.
        pin_memory (bool): Whether to use pin_memory for the DataLoader. Default is True.
    Returns:
        dict: A dictionary with the available loaders: 'train_loader', 'validation_loader', 'test_loader'.
    """
    dataloaders = {}
    # Set validation batch size (default is twice the training batch size)
    valid_batch_size = valid_batch_size or batch_size * 2
    # Create DataLoader for the training set if train_idx is not None
    if dataset.train_idx is not None:
        print('[create_dataloaders] Creating loader for training data ...')
        train_sampler = SubsetRandomSampler(dataset.train_idx)
        dataloaders['train_loader'] = DataLoader(
            dataset, 
            batch_size=batch_size, 
            sampler=train_sampler, 
            drop_last=True, 
            pin_memory=pin_memory
        )
    # Create DataLoader for the validation set if valid_idx is not None
    if dataset.valid_idx is not None:
        print('[create_dataloaders] Creating loader for validation data ...')
        valid_sampler = SubsetRandomSampler(dataset.valid_idx)
        dataloaders['validation_loader'] = DataLoader(
            dataset, 
            batch_size=valid_batch_size, 
            sampler=valid_sampler, 
            drop_last=True, 
            pin_memory=pin_memory
        )
    # Create DataLoader for the test set if test_idx is not None
    if dataset.test_idx is not None:
        print('[create_dataloaders] Creating loader for test data ...')
        test_sampler = SubsetRandomSampler(dataset.test_idx)
        dataloaders['test_loader'] = DataLoader(
            dataset, 
            batch_size=valid_batch_size,  # Can use the same size as validation
            sampler=test_sampler, 
            drop_last=True, 
            pin_memory=pin_memory
        )
    return dataloaders

#### version carlosizada ###############################################
class CustomDataset(Dataset):
    def __init__(self, 
                 paths: Dict[str, str], 
                 config: Optional[Dict] = None,
                 train_idx: Optional[List[int]] = None, 
                 valid_idx: Optional[List[int]] = None, 
                 test_idx: Optional[List[int]] = None,
                 scaler=None,
                 sigma=None,
                 sample_category: str = 'detailed_category',
                 save_files: bool = True,
                 output_dir: Optional[str] = None):
        
        self.config = config
        self.paths = paths
        self.save_files = save_files
        self.output_dir = output_dir
        self.fit_scaler = self.config['fit_scaler']
        self.scaler = scaler
        self.sigma = sigma

        if scaler is not None and sigma is not None:
            self._init_with_scaler(train_idx, valid_idx, test_idx)
        else:
            self._init_without_scaler(config, sample_category)

        # Perform scaling
        self._perform_scaling()
        # Save post-scaling data if needed
        self._save_post_scaling_data()
        # Cache the data as numpy arrays for fast access
        self._cache_data()

    def _init_without_scaler(self, config, sample_category): 
        """Inicialización cuando scaler y sigma no están definidos."""
        print('[init] Initialization without scaler')
        self.sample_category = sample_category
        self.source_name = self.config['source_name']
        self.sample_types = self.config['select_samples']
        
        # Load and process data
        self.data = self._load_data()
        self._check_consistency()
        
        # Process samples 
        self.data = self.filter_samples(self.data, self._select_samples())
        
        # Assign sample_ids and id2index_mapping
        self._assign_sample_ids_and_mapping()  # Llamada al nuevo método
        self._perform_splits()
        
        # Save pre-scaling data if needed
        self._save_pre_scaling_data()
       
    def _init_with_scaler(self, train_idx, valid_idx, test_idx):
        """Inicialización cuando scaler y sigma están definidos."""
        print('[init] Initialization with scaler')
        self.data = self._load_data()
        self._check_consistency()

        # Assign sample_ids and id2index_mapping
        self._assign_sample_ids_and_mapping()  

        if train_idx is None and valid_idx is None and test_idx is None:
            self.test_idx = [self.id2index_mapping[patient_id] for patient_id in self.sample_ids]
            self.train_idx = train_idx  # No training set
            self.valid_idx = valid_idx  # No validation set
            self.add_sample_set_label(set_name='testing', samples=self.index2id(self.id2index_mapping, self.test_idx))
        
        else:
            self.train_idx, self.valid_idx, self.test_idx = train_idx, valid_idx, test_idx

    def _assign_sample_ids_and_mapping(self):
        """
        Assigns sample_ids and id2index_mapping, ensuring that it only happens once.
        This step is critical as it stores the mapping from the sample id to the selected index.
        """
        self.sample_ids = self.data['metadata_df'].index
        self.id2index_mapping = self.id_to_index

    @property
    def id_to_index(self):
        """Mapping of patient IDs to their respective indices."""
        return {patient_id: idx for idx, patient_id in enumerate(self.sample_ids)}
    
    def _load_data(self):
        """
        Load data from the specified paths.
        """
        try:
            return {
                "rbp_expr_df": pd.read_csv(self.paths["rbp_path"], index_col=0),
                "gene_expr_df": pd.read_csv(self.paths["gene_expr_path"], index_col=0),
                "trans_expr_df": pd.read_csv(self.paths["isoform_expr_path"], index_col=0),
                "metadata_df": pd.read_csv(self.paths["metadata_path"], index_col=0)
            }
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Error loading file: {e}")
        
    def _check_consistency(self):
        """Ensure the indices are consistent across all datasets."""
        rbp_ids = self.data["rbp_expr_df"].index
        gene_ids = self.data["gene_expr_df"].index
        trans_ids = self.data["trans_expr_df"].index
        metadata_ids = self.data["metadata_df"].index
        if not (rbp_ids.equals(gene_ids) and gene_ids.equals(trans_ids) and trans_ids.equals(metadata_ids)):
            raise ValueError("Patient IDs are not aligned across matrices (RBP expression, gene expression, isoform expression, and metadata). Please ensure all matrices have the same index and order.")
        
    def _select_samples(self):
            """Select samples based on the config."""
            if self.sample_types == 'all':
                unique_types = self.data["metadata_df"][self.sample_category].unique()
                return self.data["metadata_df"].index[self.data["metadata_df"][self.sample_category].isin(unique_types)]
            elif isinstance(self.sample_types, list):
                return self.data["metadata_df"].index[self.data["metadata_df"][self.sample_category].isin(self.sample_types)]
            else:
                raise ValueError("select_samples should be 'all' or a list of strings.")
            
    def filter_samples(self, data: Dict[str, pd.DataFrame], selected_sample_ids: pd.Index): # esta igual es de utils
        """Filter the data based on the selected sample IDs."""
        return {
            key: df.loc[selected_sample_ids] for key, df in data.items()
        }
    
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
    
    def index2id(self, id_to_index_mapping: Dict[str, int], index_list: List[int]) -> pd.Index:
        """ Converts a list of indices into their corresponding sample IDs based on the id_to_index_mapping."""
        return pd.Index([key for key, idx in id_to_index_mapping.items() if idx in index_list])
    
    def add_sample_set_label(self, set_name: str, samples: list):
        """ Add a column to the metadata DataFrame indicating whether samples belong to the training or test set."""
        sample_set = set(samples)
        if 'set_type' not in self.data['metadata_df'].columns:
            self.data['metadata_df']['set_type'] = 'unknown'
        self.data['metadata_df'].loc[self.data['metadata_df'].index.isin(sample_set), 'set_type'] = set_name
    
    def _perform_splits(self) -> None:
        """
        This function handles the splitting of data into train, validation, and test sets based on the config.
        """
        if self.config["train_test_split"]:
            print("Performing train/test split...")
            self.train_idx, self.test_idx = self.split_data(
                data=self.data, 
                test_size=self.config['test_frac']
            )
            print("Writing set_type in metadata after train/test split")
            self.add_sample_set_label(set_name='training', samples=self.index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='testing', samples=self.index2id(self.id2index_mapping, self.test_idx))
        if self.config["train_val_split"]:
            print("Performing train/val split...")
            training_data = self.filter_samples(self.data, self.index2id(self.id2index_mapping, self.train_idx)) 
            self.train_idx, self.valid_idx = self.split_data(
                data=training_data, 
                test_size=self.config['val_frac']
            )
            print("Writing set_type in metadata after train/test split")
            self.add_sample_set_label(set_name='training', samples=self.index2id(self.id2index_mapping, self.train_idx))
            self.add_sample_set_label(set_name='validation', samples=self.index2id(self.id2index_mapping, self.valid_idx))
        if not self.config["train_test_split"] and not self.config["train_val_split"]:
            self.test_idx = [self.id2index_mapping[patient_id] for patient_id in self.sample_ids]
            self.train_idx = None  # No training set
            self.valid_idx = None  # No validation set
            self.add_sample_set_label(set_name='testing', samples=self.index2id(self.id2index_mapping, self.test_idx))
    
    def fit_rbp_scaler(self, train_set):
        """Fit a StandardScaler to the training set."""
        print('[fit_rbp_scaler] Fitting a StandardScaler to the training set...')
        scaler = StandardScaler()
        scaler.fit(train_set)
        sigma = np.std(scaler.transform(train_set).flatten().astype(np.float64))
        return scaler, sigma
    
    def normalize_and_clip_rbp_expression(self, transform_set, scaler, sigma):
        """Scale the RBP expression with clipping and normalization."""
        scaled_set = pd.DataFrame(scaler.transform(transform_set), index=transform_set.index, columns=transform_set.columns)
        print(f'[normalize_and_clip_rbp_expression] Mean after scaling: {scaled_set.mean().mean()}')
        print(f'[normalize_and_clip_rbp_expression] Std after scaling: {scaled_set.std().mean()}')
        print(f'[normalize_and_clip_rbp_expression] Clipping the data...')
        scaled_set = np.clip(scaled_set, -2*sigma, 2*sigma, axis=1)
        scaled_set += 2*sigma
        scaled_set /= 4*sigma
        return scaled_set.astype(np.float64)
    
    def apply_transform_to_samples(self, sample_ids, data_type):
        """Transform the data for the given sample IDs."""
        self.data['rbp_expr_df'].loc[sample_ids] = self.normalize_and_clip_rbp_expression(
            transform_set=self.data['rbp_expr_df'].loc[sample_ids], 
            scaler=self.scaler, 
            sigma=self.sigma
        )
    
    def _perform_scaling(self):
        """ Perform scaling on the training data and set scaler and sigma. """
        print("[_perform_scaling] Performing the scaling...")
        
        # Check if we have training samples and retrieve their IDs if needed
        if self.train_idx is not None:
            selected_train_sample_ids = self.index2id(self.id2index_mapping, self.train_idx)
        print(self.fit_scaler)
        if self.fit_scaler: # We need to fit the scaler with this data
            if self.train_idx is None:
                raise ValueError("[_perform_scaling] Error: train_idx is None. No training samples available for scaling.")
            print("[_perform_scaling] Fitting the scaler with the training data...")
            # Fit the scaler using the training data
            scaler, sigma = self.fit_rbp_scaler(train_set=self.filter_samples(self.data, selected_train_sample_ids)['rbp_expr_df'])
            self.scaler = scaler
            self.sigma = sigma

        else: # Use an existing scaler and sigma
            if self.scaler is None or self.sigma is None:
                raise ValueError("[_perform_scaling] Error: scaler or sigma is not set. Pre-trained scaler and sigma are required.")
            print("[_perform_scaling] Using pre-trained scaler and sigma.")
        # Transform data
        if self.train_idx is not None:
            print("[_perform_scaling] Transforming the training data...")
            self.apply_transform_to_samples(selected_train_sample_ids, 'train')
        if self.test_idx is not None:
            print("[_perform_scaling] Transforming the test data...")
            selected_test_sample_ids = self.index2id(self.id2index_mapping, self.test_idx)
            self.apply_transform_to_samples(selected_test_sample_ids, 'test')
        if self.valid_idx is not None:
            print("[_perform_scaling] Transforming the validation data...")
            selected_valid_sample_ids = self.index2id(self.id2index_mapping, self.valid_idx)
            self.apply_transform_to_samples(selected_valid_sample_ids, 'validation')
    
    def _save_pre_scaling_data(self):
        """ Save pre-scaling data if required. """
        if self.save_files and self.output_dir is not None:
            self.save_data(
                output_dir=self.output_dir, 
                folder_name=f"data/{self.source_name}/pre-scaling", 
                file_names=['rbp_expr_log2p_tpm', 'gene_expr_tpm', 'trans_expr_log2p_tpm', 'metadata_df']
            )
    
    def _save_post_scaling_data(self):
        """ Save post-scaling data if required. """
        if self.save_files and self.output_dir is not None:
            self.save_data(
                output_dir=self.output_dir,
                folder_name=f"data/{self.source_name}/post-scaling",
                file_names=['scaled_rbp_expr_log2p_tpm'],
                save_first_only=True
            )
            if self.fit_scaler:
                self._save_scaler_and_idx(folder_name=f"data/{self.source_name}")
    
    def save_data(self, output_dir: str, folder_name: str, file_names: list, save_first_only: bool = False) -> None:
        """ Save the RBP, gene, trans expression DataFrames and metadata to CSV files."""
        folder_path = os.path.join(output_dir, folder_name)
        os.makedirs(folder_path, exist_ok=True)
        keys = ['rbp_expr_df', 'gene_expr_df', 'trans_expr_df', 'metadata_df']
        # Determine the range of keys to save based on save_first_only
        range_to_save = 1 if save_first_only else len(keys)
        for i in range(range_to_save):
            key = keys[i]
            file_name = file_names[i]
            df = self.data[key]
            file_path = os.path.join(folder_path, f"{file_name}.csv")
            df.to_csv(file_path, index=True)
            print(f"[save_data] DataFrame '{key}' has been successfully exported to '{file_path}'.")
    
    def _save_scaler_and_idx(self, folder_name) -> None:
        """
        Save the scaler, sigma, and the split indices (train, valid, test) to files.
        """
        # Define file paths
        folder_path = os.path.join(self.output_dir, folder_name)
        scaler_file = os.path.join(folder_path, 'scaler.joblib')
        sigma_file = os.path.join(folder_path, 'sigma.npy') 
        train_idx_file = os.path.join(folder_path, 'train_idx.txt')
        valid_idx_file = os.path.join(folder_path, 'valid_idx.txt')
        test_idx_file = os.path.join(folder_path, 'test_idx.txt')
        # Save the scaler using joblib
        joblib.dump(self.scaler, scaler_file)
        print(f"[_save_scaler_and_idx] Scaler saved to {scaler_file}")
        # Save the sigma as a binary file (.npy for high precision)
        np.save(sigma_file, np.array(self.sigma, dtype=np.float64))
        print(f"[_save_scaler_and_idx] Sigma saved to {sigma_file}")
        # Save train, valid, and test indices to text files
        np.savetxt(train_idx_file, self.train_idx, fmt='%d')
        print(f"[_save_scaler_and_idx] Train indices saved to {train_idx_file}")
        np.savetxt(test_idx_file, self.test_idx, fmt='%d')
        print(f"[_save_scaler_and_idx] Test indices saved to {test_idx_file}")
        if self.valid_idx is not None:
            np.savetxt(valid_idx_file, self.valid_idx, fmt='%d')
            print(f"[_save_scaler_and_idx] Validation indices saved to {valid_idx_file}")
    
    @staticmethod
    def load_scaler_and_idx(path_to_saved_model_data: str):
        """Loads the scaler, sigma, and dataset indices (train, valid, and test) from saved files."""
        print('[load_scaler_and_idx] Loading the scaler, sigma and dataset indices')
        try:
            scaler = joblib.load(f'{path_to_saved_model_data}/scaler.joblib')
            sigma = np.load(f'{path_to_saved_model_data}/sigma.npy')
            try:
                train_idx = np.loadtxt(f'{path_to_saved_model_data}/train_idx.txt', dtype=int).tolist()
            except OSError:
                train_idx = None
            try:
                valid_idx = np.loadtxt(f'{path_to_saved_model_data}/valid_idx.txt', dtype=int).tolist()
            except OSError:
                valid_idx = None
            try:
                test_idx = np.loadtxt(f'{path_to_saved_model_data}/test_idx.txt', dtype=int).tolist()
            except OSError:
                test_idx = None
            return scaler, sigma, train_idx, valid_idx, test_idx
        except FileNotFoundError as e:
            print(f"[load_scaler_and_idx] File not found: {e}")
            return None, None, None, None, None
        except Exception as e:
            print(f"[load_scaler_and_idx] An error occurred: {e}")
            return None, None, None, None, None
    
    def _cache_data(self):
        """ Cache data as numpy arrays for fast access. """
        self.rbp_expr = self.data['rbp_expr_df'].values.astype(np.float64)
        self.gene_expr = self.data['gene_expr_df'].values.astype(np.float64)
        self.trans_expr = self.data['trans_expr_df'].values.astype(np.float64)
        self.metadata = self.data['metadata_df']

        self.rbp_names = self.data['rbp_expr_df'].columns.tolist()  # List of RBP gene names
        #self.gene_names = self.data['gene_expr_df'].columns.tolist()  # List of genes names
        self.trans_names = self.data['trans_expr_df'].columns.tolist()  # List of transcript names
        # self.getBM
        # Importante: here put also rbp_names, gene_names, trans_names (mierda no tengo ni el self.getBM ni los nombres de los genes reales)
        # voy a tener que en preprocessing dejar la matriz de genes como estaba y que sea aqui donde la transformamos en el formato de trans!
    
    def __getitem__(self, idx: int):
        """Return a sample from the dataset (RBP, gene, and trans expressions)."""
        rbp_exp = torch.tensor(self.rbp_expr[idx], dtype=torch.float64)
        gene_exp = torch.tensor(self.gene_expr[idx], dtype=torch.float64)
        trans_exp = torch.tensor(self.trans_expr[idx], dtype=torch.float64)
        return rbp_exp, gene_exp, trans_exp  # Solo devolvemos las expresiones
    
    def __len__(self) -> int:
        """Return the total number of samples."""
        return len(self.rbp_expr)
    
    def get_metadata(self, idx: int):
        """Return metadata for a specific sample."""
        return self.metadata.iloc[idx]

###
   
   



# #### TESTS DE CODIGO (DECIRLE A XABIER GARROTE QUE LOS PROGRAME, tiene que hacer tests para meter todas las funcionalidades)
# ############################################################################################################################################################
# ### meter tests a este código:
# total_samples = len(dataset)  # O el número total de muestras en tu dataset
# num_train = len(train_idx)
# num_val = len(valid_idx)
# num_test = len(test_idx)

# print(f'Train: {num_train}, Validation: {num_val}, Test: {num_test}')
# print(f'Total: {num_train + num_val + num_test} / {total_samples}')

# assert len(set(train_idx) & set(valid_idx)) == 0, "Hay índices duplicados entre entrenamiento y validación"
# assert len(set(train_idx) & set(test_idx)) == 0, "Hay índices duplicados entre entrenamiento y prueba"
# assert len(set(valid_idx) & set(test_idx)) == 0, "Hay índices duplicados entre validación y prueba"

# assert all(0 <= idx < total_samples for idx in train_idx), "Índice fuera de rango en entrenamiento"
# assert all(0 <= idx < total_samples for idx in valid_idx), "Índice fuera de rango en validación"
# assert all(0 <= idx < total_samples for idx in test_idx), "Índice fuera de rango en prueba"

# from collections import Counter

# train_class_distribution = Counter(your_labels[train_idx])
# val_class_distribution = Counter(your_labels[valid_idx])
# test_class_distribution = Counter(your_labels[test_idx])

# print(f'Distribución de clases en entrenamiento: {train_class_distribution}')
# print(f'Distribución de clases en validación: {val_class_distribution}')
# print(f'Distribución de clases en prueba: {test_class_distribution}')

# # distribuicon de clases de tipo tumoral
# dataset.data['metadata_df'].loc[selected_train_sample_ids].detailed_category.value_counts()
# dataset.data['metadata_df'].loc[selected_valid_sample_ids].detailed_category.value_counts()
# dataset.data['metadata_df'].loc[selected_test_sample_ids].detailed_category.value_counts()

###
# Verificar que el input está en las unidades que debería estar:

# Verificar que las expresiones son en TPM
# def verify_tpm(df_gene_tpm, df_trans_tpm):
#     gene_sums = df_gene_tpm.sum(axis=0)
#     trans_sums = df_trans_tpm.sum(axis=0)
    
#     # Las sumas deben estar cerca de 1 millón
#     print("Sumas por columna (genes):", gene_sums.describe())
#     print("Sumas por columna (transcritos):", trans_sums.describe())

# # Verificar que la expresión de los transcritos tiende a la de los genes
# def verify_transcript_gene_expression(df_gene_tpm, df_trans_tpm, getBM):
#     # Unir los DataFrames para agrupar transcritos por gen
#     merged = getBM.merge(df_trans_tpm, left_on='Transcript_ID', right_index=True)
    
#     # Agrupar por Gene_ID y sumar la expresión de transcritos por gen
#     transcript_sum_per_gene = merged.groupby('Gene_ID').sum()

#     # Reindexar el DataFrame de genes según el índice de transcript_sum_per_gene
#     gene_expression = df_gene_tpm.loc[transcript_sum_per_gene.index]

#     # Comparar la expresión de los genes con la suma de sus transcritos
#     correlation = gene_expression.corrwith(transcript_sum_per_gene)
    
#     print("Correlación entre expresión de genes y la suma de sus transcritos:")
#     print(correlation.describe())

# # Cargar los DataFrames y realizar las verificaciones
# verify_tpm(df_gene_tpm, df_trans_tpm)
# verify_transcript_gene_expression(df_gene_tpm, df_trans_tpm, getBM)
