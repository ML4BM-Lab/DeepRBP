# DeepRBP/src/deeprbp/deeplift_handler.py

import torch
import numpy as np
import pandas as pd
from tqdm import tqdm
from .logger import Logger

class DeepLiftHandler:
    ## Solve this: # Suppress specific user warnings
    # warnings.filterwarnings("ignore", message="Input Tensor .* did not already require gradients")
    # warnings.filterwarnings("ignore", message="Setting forward, backward hooks and attributes on non-linear activations")
    
    def __init__(self, model, data, explain_config):
        """
        Initialize the DeepLiftHandler with model, data, and explain configuration.

        Parameters:
            model: The model to be explained.
            data (dict): Dictionary containing the scaled RBP expression DataFrame under 'scaled_rbp_expr_df'
                     and gene expression DataFrame under 'gene_expr_df'.
            explain_config (dict): Configuration for explanations.
        """
        self.deeplift_explainer = DeepLift(model)
        self.data = data
        self.explain_config = explain_config
        self.logger = Logger(verbose=1)
       
        # Define rbps_id and trans_id from the data provided
        self.rbps_id = list(data['scaled_rbp_expr_df'].columns)
        self.trans_id = list(data['gene_expr_df'].columns)
        self.len_out_features = len(self.trans_id)

    def prepare_rbp_tensors(self):
        """
        Prepare scaled RBP tensor and reference RBP tensor based on the configuration.

        Parameters:
            explain_config (dict): Dictionary containing the configuration. 
                               Expected key: 'reference_data' with options:
                               - 'median_reference': Median of the scaled RBP expression.
                               - 'knockdown_reference': Zeroed tensor (knockdown).
                               - 'half_reference': Half-zeroed tensor.
        Returns:
            tuple: A tuple containing:
                - scaled_rbp_tensor (Tensor): Tensor of scaled RBP expression values.
                - gn_tensor (Tensor): Tensor of gene expression values
                - reference_rbp_tensor (Tensor): Tensor used as a reference for RBP expression.
        """
        scaled_rbp_tensor = torch.tensor(data['scaled_rbp_expr_df'].values, dtype=torch.float32)
        gn_tensor = torch.tensor(data['gene_expr_df'].values, dtype=torch.float32)
        self.logger.log("[prepare_rbp_tensors] Converted scaled RBP expression and gene expression to tensors (float32).")

        # Determine reference tensor based on configuration
        reference_type = self.explain_config.get('reference_data')
        if reference_type == 'median_reference':
            self.logger.log("[prepare_rbp_tensors] Using median of RBP expression as reference.")
            reference_rbp_tensor = torch.tensor(np.median(data['scaled_rbp_expr_df'], axis=0), dtype=torch.float32)
            reference_rbp_tensor = torch.reshape(reference_rbp_tensor, (1, reference_rbp_tensor.size(0)))
            self.logger.log("[prepare_rbp_tensors] Median reference tensor created.")
        
        elif reference_type == 'knockdown_reference':
            self.logger.log("[prepare_rbp_tensors] Using knockdown (zeroed) RBP expression as reference.")
            reference_rbp_tensor = torch.zeros(1, data['scaled_rbp_expr_df'].shape[1], dtype=torch.float32)
            self.logger.log("[prepare_rbp_tensors] Knockdown reference tensor created.")
        
        elif reference_type == 'half_reference':
            self.logger.log("[prepare_rbp_tensors] Using half-zeroed RBP expression as reference.")
            reference_rbp_tensor = torch.ones(1, data['scaled_rbp_expr_df'].shape[1], dtype=torch.float32) * 0.5
            self.logger.log("[prepare_rbp_tensors] Half-zeroed reference tensor created.")
        
        else:
            self.logger.error(f"[prepare_rbp_tensors] Unknown reference type: '{reference_type}'", ValueError)
            raise ValueError(f"Unknown reference type: {reference_type}")

        return scaled_rbp_tensor, gn_tensor, reference_rbp_tensor

    def _compute_attribution_scores_for_batch(self, target_node, scaled_rbp_inputs, reference_rbp_inputs, gene_expression_inputs):
        """
        Compute attribution scores for a batch of input data using the DeepLift explainer.

        Parameters:
            target_node (int): Index of the target node or output for which scores are calculated.
            scaled_rbp_inputs (Tensor): Tensor of scaled RBP input data.
            reference_rbp_inputs (Tensor): Baseline tensor (reference) for comparison during attribution.
            gene_expression_inputs (Tensor): Additional input arguments (e.g., gene expression data).

        Returns:
            scores (Tensor): Importance scores attributed to input features.
        """
        scores = self.deeplift_explainer.attribute(
            inputs=scaled_rbp_inputs,
            baselines=reference_rbp_inputs,
            target=target_node,
            additional_forward_args=gene_expression_inputs,
        )
        return scores

    def compute_attribution_scores(self, scaled_rbp_tensor, reference_rbp_tensor, gn_tensor):
        """
        Compute attribution scores batch-wise for RBP x S x T.

        Parameters:
            scaled_rbp_tensor (torch.Tensor): The scaled RBP tensor.
            reference_rbp_tensor (torch.Tensor): The reference RBP tensor.
            gn_tensor (torch.Tensor): The gene expression tensor.

        Returns:
            list: A list of batch attribution scores.
        """
        list_batch_scores = []
        for target_node in tqdm(range(self.len_out_features), desc="Calculating scores", unit="node"):
            scores_batch = self._compute_attribution_scores_for_batch(
                target_node, scaled_rbp_tensor, reference_rbp_tensor, gn_tensor
            )
            list_batch_scores.append(scores_batch)
        return list_batch_scores

    def reduce_batch_dimension(self, list_batch_scores):
        """
        Reduce the batch dimension of DeepLIFT scores to generate final TxRBP scores using a specified method.

        Parameters:
            list_batch_scores (list): List of computed attribution scores for a batch using DeepLift.
            explain_config (dict): Dictionary containing the configuration.
                               Expected key: 'batch_reduction_method' with options:
                               - 't-statistic': Computes the t-statistic across samples.
                               - 'sum_scores': Sums the scores across samples.

        Returns:
            pd.DataFrame: DataFrame containing the reduced DeepLIFT scores.
        """
        batch_reduction_type = self.explain_config.get('batch_reduction_method')
        self.logger.log(f"[reduce_batch_dimension] The reduce method selected is {batch_reduction_type}.")

        if batch_reduction_type not in ['t-statistic', 'sum_scores']:
            self.logger.error(f"[reduce_batch_dimension] Unknown batch reduction method: '{batch_reduction_type}'", ValueError)
            raise ValueError(f"Invalid batch reduction method: {batch_reduction_type}")

        batch_scores_np = [tensor.detach().numpy() for tensor in list_batch_scores]
        if batch_reduction_type == 't-statistic':
            self.logger.log("[reduce_batch_dimension] Calculating t-statistic.")
            scores_stack = np.stack(batch_scores_np, axis=0)
            mean = np.mean(scores_stack, axis=1)
            std = np.std(scores_stack, axis=1)
            num_samples = scores_stack.shape[1]
            with np.errstate(divide='ignore', invalid='ignore'):
                t_stat_scores = np.where(std != 0, mean / (std / np.sqrt(num_samples)), 0)

            zero_std_indices = np.argwhere(std == 0)
            if zero_std_indices.size > 0:
                self.logger.log(f"[reduce_batch_dimension] Found {len(zero_std_indices)} positions with std=0.")
                for idx in zero_std_indices:
                    trans = trans_id[idx[0]]
                    rbp = rbps_id[idx[1]]
                    self.logger.log(f"[reduce_batch_dimension] Zero std at Transcript: {trans}, RBP: {rbp}.")
            
            df_deeplift_TxRBP = pd.DataFrame(t_stat_scores, index=trans_id, columns=rbps_id)
            self.logger.log("[reduce_batch_dimension] T-statistic reduction completed.")

        elif batch_reduction_type == 'sum_scores':
            self.logger.log("[reduce_batch_dimension] Calculating sum of scores.")
            scores_stack = np.stack(batch_scores_np, axis=0)
            sum_scores = np.sum(scores_stack, axis=1)
            df_deeplift_TxRBP = pd.DataFrame(sum_scores, index=trans_id, columns=rbps_id)
            self.logger.log("[reduce_batch_dimension] Sum reduction completed.")

        return df_deeplift_TxRBP

    def filter_scores_for_low_expressed_genes(self, deeplift_scores, gene_expr_df, threshold=1):
        """
        Set low-expressed genes (mean expression < threshold) to 0 in the TxRBP scores DataFrame.

        Parameters:
            deeplift_scores (pd.DataFrame): DataFrame with DeepLIFT scores (TxRBP).
            gene_expr_df (pd.DataFrame): DataFrame with gene expression values.
            threshold (float): Expression threshold (default=1 TPM).

        Returns:
            pd.DataFrame: Updated DeepLIFT scores with low-expressed genes set to 0.
        """
        low_expr_genes = gene_expr_df.mean() < threshold
        deeplift_scores.loc[low_expr_genes, :] = 0
        self.logger.log(f"Low-expressed genes (mean expression < {threshold} TPM) have been excluded from the TxRBP scores.")
        return deeplift_scores
