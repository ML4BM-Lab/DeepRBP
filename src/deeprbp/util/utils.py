# src/deeprbp/util/utils.py

import pandas as pd
import torch
from torchinfo import summary
from torch.utils.data import Dataset
from typing import List, Dict, Optional
import os

class CustomTensorDataset(Dataset):
    def __init__(self, data, 
                 getBM, 
                 rbp_data_key='scaled_rbp_df', 
                 gene_data_key='gene_df', 
                 transcript_data_key='isoform_df',
                 trans_col_name='Transcript_ID',
                 gene_col_name='Gene_ID'
                 ):
        """
        Custom dataset for loading RBP, gene, and transcript expressions
        Parameters:
        - data: Dictionary containing the DataFrames for RBP, gene, and transcript expressions.
        - getBM: DataFrame containing Transcript_IDs and associated Gene_IDs.
        - rbp_data_key: The name of the DataFrame column for RBP features.
        - gene_data_key: The name of the DataFrame column for gene expression features.
        - transcript_data_key: The name of the DataFrame column for transcript expression features.
        - gene_col_name: The name of the column in getBM that contains the Gene IDs.
        - trans_col_name: The name of the column in getBM that contains the Transcript IDs.
        """
        for key in [rbp_data_key, gene_data_key, transcript_data_key]:
            if key not in data:
                raise KeyError(f"{key} not found in data.")
        self.data = data # Original data
        # Automatically detect the feature names from the DataFrames
        self.rbp_names = data[rbp_data_key].columns.tolist()
        self.gene_names = data[gene_data_key].columns.tolist()
        self.trans_names = data[transcript_data_key].columns.tolist()
        # Expand the gene matrix data to match the transcript shape data 
        genes_names_each_trans = getBM[getBM[trans_col_name].isin(self.trans_names)][gene_col_name]
        data[gene_data_key] = data[gene_data_key].loc[:, genes_names_each_trans]
        # Store tensors for each feature
        features_data = (rbp_data_key, gene_data_key, transcript_data_key)
        self.features = {feature: torch.tensor(data[f"{feature}"].values, dtype=torch.float64) for feature in features_data}
    def __len__(self) -> int:
        return len(next(iter(self.features.values())))
    def __getitem__(self, idx: int):
        return {feature: values[idx] for feature, values in self.features.items()}
    
def print_section_separator(char="═", width=50):
    """
    Prints a decorative separator line with a minimalist design.
    
    Parameters:
        char (str): The character to use for the separator line.
        width (int): The width of the separator line.
    """
    line = char * width
    print(f"\n{line}\n")

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

def save_data(data, path, custom_names=None):
    os.makedirs(path, exist_ok=True)
    if custom_names is None:
        custom_names = {}
    for key, df in data.items():
        if isinstance(df, pd.DataFrame):
            file_name = custom_names.get(key, f'{key}.csv')
            df.to_csv(os.path.join(path, file_name))

def load_data(path):
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
    
def calculate_category_proportions(data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Calculate the proportions of different tumor types samples based on the 'detailed_category' in the metadata. To 
    be sure that the stratified split is always done correctly.

    Parameters:
    - data (Dict[str, pd.DataFrame]): A dictionary containing the dataset, including 'metadata_df' 
      with tumor type information.

    Returns:
    - pd.DataFrame: A DataFrame containing the proportion of each tumor type.
    """
    # Access the metadata DataFrame
    metadata_df = data['metadata_df']
    # Count the number of samples by tumor type
    category_counts = metadata_df['detailed_category'].value_counts()
    # Calculate the proportion of each tumor type over the total samples
    total_samples = len(metadata_df)
    category_proportions = category_counts / total_samples
    # Prepare the results as a DataFrame
    category_proportions = category_proportions.reset_index()
    category_proportions.columns = ['detailed_category', 'proportion']
    return category_proportions

def summarize_model(model, train_loader, batch_size, device='cpu'):
    """
    Generates dummy inputs based on the training data and summarizes the model.

    Args:
        model: The model to summarize.
        train_loader: The DataLoader for the training dataset.
        batch_size: The batch size to use for the dummy inputs.
        device: The device to which the inputs should be sent ('cpu' or 'cuda').

    Returns:
        None
    """
    # Create dummy inputs
    rbp_input_size = next(iter(train_loader))['scaled_rbp_df'].shape[1]  # Number of features for rbp_df
    gene_input_size = next(iter(train_loader))['gene_df'].shape[1]      # Number of features for gene_df
    # Generate dummy inputs
    dummy_rbp = torch.randn(batch_size, rbp_input_size).to(device)  # Use 'cpu' or 'cuda' based on your setup
    dummy_gene = torch.randn(batch_size, gene_input_size).to(device)  # Use 'cpu' or 'cuda' based on your setup
    # Use torchinfo to summarize the model
    summary(model, input_data=(dummy_rbp, dummy_gene))

