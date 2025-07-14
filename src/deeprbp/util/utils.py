# src/deeprbp/util/utils.py

import pandas as pd
import torch
from typing import List, Dict
import os
import GPUtil
import random
import numpy as np
import torch
import pytorch_lightning as pl
import re
import warnings

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
    adjusted_size = min(len(dataset), batch_size)
    if adjusted_size < batch_size:
        warnings.warn(f"Requested batch size {batch_size} is greater than the number of samples in the dataset ({len(dataset)}). "
                      f"Adjusting to {adjusted_size}.", UserWarning)
    return adjusted_size

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

def log_section_separator(title, max_length=70):
    title = title.upper()
    title_length = len(title)
    line_length = min(max_length, title_length + 4)             
    padding_length = line_length - 2 - title_length
    left_padding = padding_length // 2
    right_padding = padding_length - left_padding
    separator_line = "#" * line_length
    title_line = "#" + " " * left_padding + title + " " * right_padding + "#"
    print(separator_line)
    print(title_line)
    print(separator_line)

def print_gpu_memory_info():
    # Get the list of GPUs
    gpus = GPUtil.getGPUs()
    # Iterate over the GPUs and print their details
    for gpu in gpus:
        memory_free_gb = gpu.memoryFree / 1024  # Convert free memory to GB
        memory_used_gb = gpu.memoryUsed / 1024  # Convert used memory to GB
        gpu_load_percentage = gpu.load * 100  # Convert GPU load to percentage
        print(f"GPU: {gpu.id} {gpu.name} | Memory Free: {memory_free_gb:.2f} GB | Memory Used: {memory_used_gb:.2f} GB | GPU Load: {gpu_load_percentage:.2f}%")
    print()  # Empty line for better readability

def set_random_seed(seed=123):
    """Sets the random seed for reproducibility in training.

    Args:
        seed (int): The seed value to set for random number generators.
    """
    # Ensure that training is as reproducible as possible
    torch.manual_seed(seed)  # Random seed for CPU tensors
    torch.cuda.manual_seed_all(seed)  # Seed for generating random numbers for the current GPU
    torch.backends.cudnn.deterministic = True  # Ensure deterministic computations
    torch.backends.cudnn.benchmark = False
    random.seed(seed)
    np.random.seed(seed)
    os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'  # Set the environment variable for deterministic behavior
    pl.seed_everything(seed, workers=seed)  # Set seed for PyTorch Lightning
    # Set PyTorch print options
    torch.set_printoptions(precision=5, threshold=10_000)

def print_if_main(*args):
    """Print the message if the current process is the main process 
    (LOCAL_RANK == 0) or if running on a CPU."""
    local_rank = int(os.environ.get("LOCAL_RANK", 0))
    if local_rank == 0 or not torch.cuda.is_available():
        print(*args)

def setup_output_directory(output_dir: str) -> str:
    """Determine the output directory based on the execution context.
    
    Args:
        output_dir (str): The base output directory.

    Returns:
        str: The unique output directory based on the execution context.
    """
    local_rank = int(os.environ.get("LOCAL_RANK", 0))
    print(f"Local Rank: {local_rank}")
    if local_rank == 0 or not torch.cuda.is_available():
        # Running on CPU or main process (rank 0)
        unique_output_dir = output_dir
    else:
        # Running on GPU, create unique directory for this GPU
        unique_output_dir = os.path.join(output_dir, f"gpu_{local_rank}")
    if unique_output_dir:
        os.makedirs(unique_output_dir, exist_ok=True)
        print(f"Directory ensured: {unique_output_dir}")
    return unique_output_dir

def find_latest_checkpoint(checkpoint_dir: str) -> str:
    """Function to find the latest checkpoint file in the specified directory"""
    # List all files in the checkpoint directory
    files = os.listdir(checkpoint_dir)
    # Filter for files that end with .ckpt
    ckpt_files = [f for f in files if f.endswith('.ckpt')]
    # Check if there's exactly one checkpoint file
    if len(ckpt_files) == 1:
        return os.path.join(checkpoint_dir, ckpt_files[0])
    else:
        raise ValueError(f"Expected one checkpoint file in {checkpoint_dir}, found: {len(ckpt_files)}")

def find_best_checkpoint(checkpoint_dir: str) -> str:
    """Function to find the checkpoint file with the lowest validation loss in the specified directory."""
    # List all files in the checkpoint directory
    files = os.listdir(checkpoint_dir)
    # Filter for files that end with .ckpt and extract val_loss
    ckpt_files = [f for f in files if f.endswith('.ckpt')]
    # Dictionary to hold filename and corresponding val_loss
    val_loss_dict = {}
    for f in ckpt_files:
        # Use regex to extract the val_loss value from the filename
        match = re.search(r'validation_loss=([0-9.]+)', f)
        if match:
            # Clean the val_loss string to avoid conversion errors
            val_loss_str = match.group(1).rstrip('.')  # Remove any trailing dot
            try:
                val_loss = float(val_loss_str)
                val_loss_dict[f] = val_loss
            except ValueError:
                print(f"Skipping file {f} due to conversion error.")
    # Check if there are any checkpoint files found
    if not val_loss_dict:
        raise ValueError(f"No valid checkpoint files found in {checkpoint_dir}")
    # Find the checkpoint file with the minimum val_loss
    best_ckpt_file = min(val_loss_dict, key=val_loss_dict.get)
    return os.path.join(checkpoint_dir, best_ckpt_file)
