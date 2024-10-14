# El data está en : /scratch/jsanchoz/DeepRBP/data/training_module/processed
# La idea es que este loader sea capaz de cargarte las matrices de rbp, trans, y gene, que te divida en train, test y val
# usando un porcentaje de los tipos tumorales. Con las de TCGA puedes hacer tb toy version para optimización, GTEX que vaya todo para el test.
# El scaling tb sea uno de los procesos

# Faltaría como mucho ahora ya tener en cuenta esto de aquí del toy (14/10):
# #toy_set = False # with frac = 0.4 # esta funcionalidad podría no tener mucho sentido tenerla aqui y seria mejor 
# # en un arch_optimization.py usar simplemente la funcion de split_data 0.4%.

############################################################################################################
import os
import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader, SubsetRandomSampler
from sklearn.model_selection import train_test_split as sk_train_test_split
from sklearn.preprocessing import StandardScaler
from typing import Dict, Optional, List

from config_loader import get_config 

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

class CustomDataset(Dataset):
    def __init__(self, 
                 config: Dict, 
                 source_train: Optional[str] = None, 
                 paths: Optional[Dict[str, str]] = None, 
                 sample_category: str = 'detailed_category',
                 train_idx: Optional[List[int]] = None, 
                 valid_idx: Optional[List[int]] = None, 
                 test_idx: Optional[List[int]] = None,
                 output_dir: Optional[str] = None,
                 save_files: bool = False,
                 scaler=None,
                 sigma=None) -> None:
        """
        Custom Dataset that loads RBP, gene expression, isoform expression data 
        along with metadata to perform sample selection
        Args:
            config (dict): Configuration dictionary.
            source_train (str, optional): Dataset name in config to load paths from, e.g., 'TCGA' or 'GTEX'.
            paths (dict, optional): Dictionary with paths to the data files. If not provided, it uses the config paths for source_train.
            sample_category (str): Sample category column name in metadata.
            train_idx (list, optional): Predefined list of training indices. If provided, it skips the splitting process.
            valid_idx (list, optional): Predefined list of validation indices. If provided, it skips the splitting process.
            test_idx (list, optional): Predefined list of test indices. If provided, it skips the splitting process.
            output_dir (str, optional): Directory where processed data will be saved, including pre-scaling and post-scaling data files.
            save_files (bool): Flag to determine whether to save data files; defaults to False.
        """
        # Initialize
        self.config = config
        self.source_train = source_train
        self.paths = paths
        self.sample_category = sample_category
        self.train_idx, self.valid_idx, self.test_idx = train_idx, valid_idx, test_idx
        self.output_dir = output_dir
        self.save_files = save_files
        self.scaler = scaler
        self.sigma = sigma
        # Load data from paths
        self.data = self._load_data()
        # Check if IDs are consistent
        self._check_consistency()
        # Select samples based on the selected tumor types in config
        self.data = self.filter_samples(self.data, self._select_samples())
        # Store the sample (patient) IDs
        self.sample_ids = self.data['metadata_df'].index
        self.id2index_mapping = self.id_to_index
        # Splitting of data into train, validation, and test sets based on the config and provided indices.
        self._perform_splits()
        # Save pre-scaling data if save_files is True and output_dir is provided
        if self.save_files and self.output_dir is not None:
            self.save_data(
                        output_dir = self.output_dir, 
                        folder_name = "data/pre-scaling", 
                        file_names = ['rbp_expr_log2p_tpm', 'gene_expr_tpm', 'trans_expr_log2p_tpm', 'metadata_df']
                        )
        # Perform scaling
        self._perform_scaling()
        # Save post-scaling data if save_files is True and output_dir is provided
        # Achtung: Here attention, you are saving all the rbp data scaled together but scaling is performed taking into account
        # training, validation and test samples individually.
        if self.save_files and self.output_dir is not None:
            self.save_data(
                output_dir=self.output_dir,
                folder_name="data/post-scaling",
                file_names=['scaled_rbp_expr_log2p_tpm'],
                save_first_only=True
            )
        # Cache the data as numpy arrays for fast access
        self.rbp_expr = self.data['rbp_expr_df'].values
        self.gene_expr = self.data['gene_expr_df'].values
        self.trans_expr = self.data['trans_expr_df'].values
        self.metadata = self.data['metadata_df']
    #
    @property
    def id_to_index(self):
        """Mapping of patient IDs to their respective indices."""
        return {patient_id: idx for idx, patient_id in enumerate(self.sample_ids)}
    def _load_data(self):
        """
        Internal method to load data either from config or custom paths.
        """
        # Use paths from config if not explicitly provided
        if self.paths is None and self.source_train is not None:
            paths = self.config['data']['training_module'][self.source_train]
        else:
            paths = self.paths
        try:
            return {
                "rbp_expr_df": pd.read_csv(paths["rbp_path"], index_col=0),
                "gene_expr_df": pd.read_csv(paths["gene_expr_path"], index_col=0),
                "trans_expr_df": pd.read_csv(paths["isoform_expr_path"], index_col=0),
                "metadata_df": pd.read_csv(paths["metadata_path"], index_col=0)
            }
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Error loading file: {e}")
    ###
    def _check_consistency(self):
        """Ensure the indices are consistent across all datasets."""
        rbp_ids = self.data["rbp_expr_df"].index
        gene_ids = self.data["gene_expr_df"].index
        trans_ids = self.data["trans_expr_df"].index
        metadata_ids = self.data["metadata_df"].index
        if not (rbp_ids.equals(gene_ids) and gene_ids.equals(trans_ids) and trans_ids.equals(metadata_ids)):
            raise ValueError("Patient IDs are not aligned across matrices (RBP expression, gene expression, isoform expression, and metadata). Please ensure all matrices have the same index and order.")
    ###
    def _select_samples(self):
        """Select samples based on the tumor types specified in the config."""
        selected_tumor_types = self.config["train_prediction_model"]["train_tumor_types"]
        if selected_tumor_types == 'all':
            return self.data["metadata_df"].index
        else:
            return self.data["metadata_df"].index[self.data["metadata_df"][self.sample_category].isin(selected_tumor_types)]
    ###
    def filter_samples(self, data: Dict[str, pd.DataFrame], selected_sample_ids: pd.Index): # esta igual es de utils
        """Filter the data based on the selected sample IDs."""
        return {
            key: df.loc[selected_sample_ids] for key, df in data.items()
        }
    ###
    def split_data(self, data, test_size):
        """Perform a stratified train-test split based on the detailed_category in metadata."""
        sample_category = data['metadata_df'][self.sample_category]
        print(f'[split_data] Performing a stratified data split with fraction division equal to {test_size}')
        # Step 1: Perform the stratified split based on patient IDs (still strings at this point) and sample category
        train_idx, test_idx = sk_train_test_split(
            data['metadata_df'].index,
            test_size=test_size,
            stratify=sample_category,
            random_state=self.config['train_prediction_model']['seed']
        )
        # Step 2: Convert string indices (train_idx, test_idx) to integer indices
        train_idx = [self.id2index_mapping[patient_id] for patient_id in train_idx]
        test_idx = [self.id2index_mapping[patient_id] for patient_id in test_idx]
        return train_idx, test_idx
    ###
    def index2id(self, id_to_index_mapping: Dict[str, int], index_list: List[int]) -> pd.Index:
        """ Converts a list of indices into their corresponding sample IDs based on the id_to_index_mapping."""
        return pd.Index([key for key, idx in id_to_index_mapping.items() if idx in index_list])
    ###
    def add_sample_set_label(self, set_name: str, samples: list):
        """ Add a column to the metadata DataFrame indicating whether samples belong to the training or test set."""
        sample_set = set(samples)
        if 'set_type' not in self.data['metadata_df'].columns:
            self.data['metadata_df']['set_type'] = 'unknown'
        self.data['metadata_df'].loc[self.data['metadata_df'].index.isin(sample_set), 'set_type'] = set_name
    ###
    def _perform_splits(self) -> None:
        """
        This function handles the splitting of data into train, validation, and test sets based on the config and provided indices.
        """
        # Case 1: If train_idx or test_idx is provided, respect those values
        if self.train_idx is not None and self.test_idx is not None:
            if self.config['train_prediction_model']['train_test_split']:
                print("Warning: Both train_idx and test_idx are provided. Skipping train/test split.")
        #
        # Case 2: If both train_test_split and train_val_split are False, assign all samples to test_idx
        elif not self.config['train_prediction_model']['train_test_split'] and not self.config['train_prediction_model']['train_val_split']:
            print("Assigning all samples to test set as both train_test_split and train_val_split are False.")
            self.test_idx = [self.id2index_mapping[patient_id] for patient_id in self.sample_ids]
            self.train_idx = None  # No training set
            self.valid_idx = None  # No validation set
        #
        # Case 3: Perform train/test split if necessary and if indices are not provided
        elif self.config['train_prediction_model']['train_test_split']:
            print("Performing train/test split...")
            self.train_idx, self.test_idx = self.split_data(
                data=self.data, 
                test_size=self.config['train_prediction_model']['test_frac']
            )
        # Write set_type in metadata After train/test split
        if self.train_idx is not None:
            self.add_sample_set_label(set_name='training', samples=self.index2id(self.id2index_mapping, self.train_idx))
        if self.test_idx is not None:
            self.add_sample_set_label(set_name='testing', samples=self.index2id(self.id2index_mapping, self.test_idx))
        #
        # Case 4: Handle validation split if required
        if self.train_idx is not None and self.valid_idx is None:
            if self.config['train_prediction_model']['train_val_split']:
                print("Performing train/val split...")
                training_data = self.filter_samples(self.data, self.index2id(self.id2index_mapping, self.train_idx)) 
                self.train_idx, self.valid_idx = self.split_data(
                    data=training_data, 
                    test_size=self.config['train_prediction_model']['val_frac']
                )
        elif self.valid_idx is not None and self.config['train_prediction_model']['train_val_split']:
            print("Warning: valid_idx is provided. Skipping train/val split.")
        #
        # Write set_type in metadata After train/val split
        if self.train_idx is not None:
            self.add_sample_set_label(set_name='training', samples=self.index2id(self.id2index_mapping, self.train_idx))
        if self.valid_idx is not None:
            self.add_sample_set_label(set_name='validation', samples=self.index2id(self.id2index_mapping, self.valid_idx))
    ###
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
    ###
    def fit_scaler(self, train_set):
        """Fit a StandardScaler to the training set."""
        print('[fit_scaler] Fitting a StandardScaler to the training set...')
        scaler = StandardScaler()
        scaler.fit(train_set)
        sigma = np.std(scaler.transform(train_set).flatten())
        return scaler, sigma
    ###
    def scale_rbp_expression(self, transform_set, scaler, sigma):
        """Scale the RBP expression with clipping and normalization."""
        print(f'[scale_rbp_expression] Scaling the data...')
        scaled_set = pd.DataFrame(scaler.transform(transform_set), index=transform_set.index, columns=transform_set.columns)
        print(f'[scale_rbp_expression] Mean after scaling: {scaled_set.mean().mean()}')
        print(f'[scale_rbp_expression] Std after scaling: {scaled_set.std().mean()}')
        print(f'[scale_rbp_expression] Clipping the data...')
        scaled_set = np.clip(scaled_set, -2*sigma, 2*sigma, axis=1)
        scaled_set += 2*sigma
        scaled_set /= 4*sigma
        return scaled_set
    ###
    def _perform_scaling(self):
        """ Perform scaling on the training data and set scaler and sigma. """
        print("[_perform_scaling] Performing the scaling...")
        # Check if the scaler and sigma are already set
        if self.scaler is not None and self.sigma is not None:
            print("[_perform_scaling] Warning: Scaler and sigma are already set. Skipping fitting the scaler with the training data.")
        # If the scaler is not set, we need to fit it
        elif self.scaler is None and self.sigma is None:
            if self.train_idx is None:
                raise ValueError("[_perform_scaling] Error: train_idx is None. No training samples available for scaling.")
            else:
                selected_train_sample_ids = self.index2id(self.id2index_mapping, self.train_idx)   
                if len(selected_train_sample_ids) == 0:
                    raise ValueError("[_perform_scaling] Error: No training samples available for scaling.")
                print("[_perform_scaling] Fitting the scaler with the training data...")
                scaler, sigma = self.fit_scaler(train_set=self.filter_samples(self.data, selected_train_sample_ids)['rbp_expr_df'])
                self.scaler = scaler
                self.sigma = sigma
        # Transform the training data
        print("[_perform_scaling] Scaling the training data...")
        if len(selected_train_sample_ids) > 0:
            self.data['rbp_expr_df'].loc[selected_train_sample_ids] = self.scale_rbp_expression(
                                    transform_set=self.data['rbp_expr_df'].loc[selected_train_sample_ids], 
                                    scaler=self.scaler, 
                                    sigma=self.sigma
                                    )
        # Transform the test data
        print("[_perform_scaling] Scaling the test data...")
        selected_test_sample_ids = self.index2id(self.id2index_mapping, self.test_idx)
        if len(selected_test_sample_ids) > 0:
            self.data['rbp_expr_df'].loc[selected_test_sample_ids] = self.scale_rbp_expression(
                                    transform_set=self.data['rbp_expr_df'].loc[selected_test_sample_ids], 
                                    scaler=self.scaler, 
                                    sigma=self.sigma)
        # Transform the validation data
        if self.valid_idx is not None:
            selected_valid_sample_ids = self.index2id(self.id2index_mapping, self.valid_idx)
            if len(selected_valid_sample_ids) == 0:
                print("[_perform_scaling] Warning: No validation samples available for scaling.")
            else:
                print("[_perform_scaling] Scaling the validation data...")
                self.data['rbp_expr_df'].loc[selected_valid_sample_ids] = self.scale_rbp_expression(
                                    transform_set=self.data['rbp_expr_df'].loc[selected_valid_sample_ids], 
                                    scaler=self.scaler, 
                                    sigma=self.sigma)
    ###
    def __getitem__(self, idx: int):
        """Return a sample from the dataset (RBP, gene, and trans expressions)."""
        rbp_exp = torch.tensor(self.rbp_expr[idx], dtype=torch.float32)
        gene_exp = torch.tensor(self.gene_expr[idx], dtype=torch.float32)
        trans_exp = torch.tensor(self.trans_expr[idx], dtype=torch.float32)
        return rbp_exp, gene_exp, trans_exp  # Solo devolvemos las expresiones
    ###
    def __len__(self) -> int:
        """Return the total number of samples."""
        return len(self.rbp_expr)
    def get_metadata(self, idx: int):
        """Return metadata for a specific sample."""
        return self.metadata.iloc[idx]
   
################################# main code

# hola joseba vas por aquí: 

### Save the scaler and sigma used in Training to a file
        # filename_scaler = os.path.join(path_save_files, 'scaler_sfs.joblib')
        # filename_sigma = os.path.join(path_save_files, 'sigma_sfs.txt')
        # joblib.dump(data_scale.scaler_sfs, filename_scaler)
        # np.savetxt(filename_sigma, [data_scale.sigma_sfs])
 # Load or save the scaler and sigma used in the model training
        # scaler_sfs = joblib.load(path_save_files+'/scaler_sfs.joblib')
        # with open(path_save_files+'/sigma_sfs.txt', 'r') as f:
        #     sigma_sfs = f.readline().strip()
        # sigma_sfs = np.float128(sigma_sfs)

## luego:
# 1. Usando source_train desde el config:
# dataset = CustomDataset(config=config, source_train='TCGA', output_dir=output_dir, save_files=True)
# 2. Usando rutas personalizadas:
# paths = { # esto sobre todo para volver a ejecutar un pre-scaling data ya utilizado o unos paths q no estén en la config.
#     "rbp_path": "/custom/path/rbp.csv",
#     "isoform_expr_path": "/custom/path/isoform.csv",
#     "metadata_path": "/custom/path/metadata.csv",
#     "gene_expr_path": "/custom/path/gene.csv"
# }
# dataset = CustomDataset(config=config, paths=paths, output_dir=output_dir, paths=paths)

config = get_config()
output_dir = config['train_prediction_model']['output_dir']
batch_size = config['train_prediction_model']['batch_size']
dataset = CustomDataset(config=config, source_train='TCGA', output_dir=output_dir, save_files=True)
#dataset = CustomDataset(config=config, source_train='TCGA')

dataloaders = create_dataloaders(dataset, batch_size=batch_size)

# Acceder a los loaders que se hayan creado
train_loader = dataloaders.get('train_loader', None)
valid_loader = dataloaders.get('validation_loader', None)
test_loader = dataloaders.get('test_loader', None)

# esto ya es parte del train.py
num_epochs = 10
for epoch in range(num_epochs):
    for batch_index, (rbp_exp, gene_exp, trans_exp) in enumerate(train_loader):
        print(type(rbp_exp), type(gene_exp), type(trans_exp))  # Verificar tipos
        print(rbp_exp, gene_exp, trans_exp)  # Inspeccionar contenido
        print(rbp_exp.shape, gene_exp.shape, trans_exp.shape) 
        break
    break

    





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
