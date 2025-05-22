
import pandas as pd
from tabulate import tabulate 
import torch
from torch.utils.data import Dataset

from ..util.logger import Logger

class DeepRBPExpressionDataset(Dataset):
    def __init__(self, 
                 data: dict, 
                 getBM: pd.DataFrame, 
                 rbp_data_key: str = 'scaled_rbp_df', 
                 gene_data_key: str = 'gene_df', 
                 transcript_data_key: str = 'isoform_df',
                 trans_col_name: str = 'Transcript_ID',
                 gene_col_name: str = 'Gene_ID',
                 verbose: int = 0):
        """
        Custom dataset for loading RBP, gene, and transcript expressions.

        Parameters:
        - data: Dictionary containing the DataFrames for RBP, gene, and transcript expressions.
        - getBM: DataFrame containing Transcript_IDs and associated Gene_IDs.
        - rbp_data_key: The name of the DataFrame column for RBP features.
        - gene_data_key: The name of the DataFrame column for gene expression features.
        - transcript_data_key: The name of the DataFrame column for transcript expression features.
        - gene_col_name: The name of the column in getBM that contains the Gene IDs.
        - trans_col_name: The name of the column in getBM that contains the Transcript IDs.
        - verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during dataset initialization.
            - 0: No logging (suppress all output).
            - 1: Basic logging (show feature names and their shapes).
        """
        self.verbose = verbose
        self.logger = Logger(self.verbose) 
        for key in [rbp_data_key, gene_data_key, transcript_data_key]:
            if key not in data:
                raise KeyError(f"{key} not found in data.")
        self.original_data = {key: df.copy() for key, df in data.items()} # Original data
        data_copy = {key: df.copy() for key, df in data.items()}
        # Automatically detect the feature names from the DataFrames
        self.rbp_names = data_copy[rbp_data_key].columns.tolist()
        self.gene_names = data_copy[gene_data_key].columns.tolist()
        self.trans_names = data_copy[transcript_data_key].columns.tolist()
        # Expand the gene matrix data to match the transcript shape data 
        genes_names_each_trans = getBM[getBM[trans_col_name].isin(self.trans_names)][gene_col_name]
        data_copy[gene_data_key] = data_copy[gene_data_key].loc[:, genes_names_each_trans]
        # Store tensors for each feature
        features_data = (rbp_data_key, gene_data_key, transcript_data_key)
        self.features = {feature: torch.tensor(data_copy[f"{feature}"].values, dtype=torch.float32) for feature in features_data} # CHANGED THIS TO FLOAT32 JOSEBA!
    def print_tensor_shapes(self):
        self.logger.log("Features stored in the dataset:", level=self.verbose)
        table_data = []
        for feature, tensor in self.features.items():
            table_data.append([feature, tensor.shape])
        print(tabulate(table_data, headers=["Feature", "Shape"], tablefmt="grid"))
    def __len__(self) -> int:
        return len(next(iter(self.features.values())))
    def __getitem__(self, idx: int):
        return {feature: values[idx] for feature, values in self.features.items()}