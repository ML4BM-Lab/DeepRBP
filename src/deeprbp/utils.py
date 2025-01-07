# utils.py
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
from typing import List, Dict, Optional, Union
import os

class CustomTensorDataset(Dataset):
    def __init__(self, data, features_data=('rbp_expr', 'gene_expr', 'trans_expr')):
        """
        Custom dataset for loading RBP, gene, and transcript expressions
        """
        self.features = {feature: torch.tensor(data[f"{feature}_df"].values, dtype=torch.float64) for feature in features_data}
    def __len__(self) -> int:
        return len(next(iter(self.features.values())))
    def __getitem__(self, idx: int):
        return {feature: values[idx] for feature, values in self.features.items()}
       
def select_sample_ids_by_type(metadata_df: pd.DataFrame, sample_category_col: str, sample_types):
    """
    Select sample IDs from the metadata based on the provided sample types.

    Args:
        sample_category_col (str): Column in metadata to categorize samples.
        sample_types (str or list): 'all' to select all samples, or a list of specific sample types.

    Returns:
        pd.Index: Index of selected sample IDs.
    """
    # Handle 'all' case
    if 'all' in sample_types:
        return metadata_df.index
    # Handle specific sample types
    return metadata_df.index[metadata_df[sample_category_col].isin(sample_types)]

def filter_data_by_sample_ids(data: Dict[str, pd.DataFrame], selected_sample_ids: pd.Index) -> Dict[str, pd.DataFrame]:
    """
    Filter the data based on the selected sample IDs.

    Args:
        data (Dict[str, pd.DataFrame]): Dictionary containing the data to be filtered.
        selected_sample_ids (pd.Index): List of sample IDs to retain.

    Returns:
        dict: Filtered data dictionary.
    """
    return {key: df.loc[selected_sample_ids] for key, df in data.items()}

def index2id(id_to_index_mapping: Dict[str, int], index_list: List[int]) -> pd.Index:
    """ Converts a list of indices into their corresponding sample IDs based on the id_to_index_mapping."""
    return pd.Index([key for key, idx in id_to_index_mapping.items() if idx in index_list])

def save_processed_data(data, path):
    os.makedirs(path, exist_ok=True)
    for key, df in data.items():
        if isinstance(df, pd.DataFrame):
            df.to_csv(os.path.join(path, f'{key}.csv'))

def load_processed_data(path):
    data = {}
    for file in os.listdir(path):
        if file.endswith('.csv'):
            key = file.replace('.csv', '')
            data[key] = pd.read_csv(os.path.join(path, file), index_col=0)
    return data

def adjust_batch_size(dataset, batch_size):
    """
    Adjust the batch size based on the number of samples in the dataset.

    Args:
        dataset (torch.utils.data.Dataset): The dataset containing the data samples.
        batch_size (int): The desired batch size.

    Returns:
        int: The adjusted batch size, which is the minimum of the number of samples in the dataset and the requested batch size.
    """
    return min(len(dataset), batch_size)

def ensure_directory_exists(directory: str) -> None:
    """
    Ensures the given directory exists. If it doesn't, it creates it.

    Parameters:
    -----------
    directory : str
        Path of the directory to check/create.
    
    Returns:
    --------
    None
    """
    if directory and not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

def get_gene_info(gene_ids_or_names: List[str], getBM: pd.DataFrame, return_type: str = 'names') -> List[Optional[str]]:
    """
    Retrieves Gene_names from Gene_IDs or Gene_IDs from Gene_names using the getBM DataFrame.

    Parameters:
    gene_ids_or_names (List[str]): A list of Gene_IDs or Gene_names.
    getBM (pd.DataFrame): A DataFrame that contains the relationship between Gene_ID and Gene_name.
    return_type (str): Indicates what to return:
                       - 'names': return Gene_names for given Gene_IDs
                       - 'ids': return Gene_IDs for given Gene_names

    Returns:
    List[Optional[str]]: A list of Gene_names or Gene_IDs corresponding to the provided input.
                         If a Gene_ID or Gene_name does not have a corresponding entry, None will be returned.
    """
    getBM_subset = getBM[['Gene_ID', 'Gene_name']].drop_duplicates()
    if return_type == 'names':
        # Convert Gene_IDs to Gene_names
        getBM_subset.set_index('Gene_ID', inplace=True)
        return getBM_subset.loc[gene_ids_or_names]['Gene_name'].values.tolist()
    
    elif return_type == 'ids':
        # Convert Gene_names to Gene_IDs
        getBM_subset.set_index('Gene_name', inplace=True)
        #gene_ids = getBM_subset.loc[gene_ids_or_names]
        # Prepare to retrieve Gene_IDs and print duplicates
        gene_ids = []
        for gene_name in gene_ids_or_names:
            if gene_name in getBM_subset.index:
                occurrences = getBM_subset.loc[gene_name]
                if len(occurrences) > 1:
                    print(f"Duplicate entries found for Gene_name '{gene_name}': {occurrences['Gene_ID'].values.tolist()}")
                # Use the first occurrence as the final result
                if isinstance(occurrences, pd.DataFrame):
                    first_gene_id = occurrences['Gene_ID'].iloc[0]  # Safe access for DataFrame
                else:
                    first_gene_id = occurrences['Gene_ID']  # Direct access for Series
                print(f"Using '{first_gene_id}' as Gene_ID for Gene_name '{gene_name}'")
                gene_ids.append(first_gene_id)
            else:
                gene_ids.append(None)
        return gene_ids
    else:
        raise ValueError("Invalid return_type. Use 'names' or 'ids'.")
