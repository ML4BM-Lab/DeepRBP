#src/deeprbp/explainability_module/deeplift_handler.py

import torch
from captum.attr import DeepLift
import numpy as np
import pandas as pd
from tqdm import tqdm

from ..util.logger import Logger

import warnings
warnings.filterwarnings('ignore')
#warnings.filterwarnings("ignore", message="Input Tensor .* did not already require gradients")
warnings.filterwarnings("ignore", message="Setting forward, backward hooks and attributes on non-linear activations")
    
class DeepLiftHandler:
    """
    DeepLiftHandler class for managing the DeepLIFT explanation process within the DeepRBP explainability module.

    This class prepares reference data, computes attribution scores using the DeepLIFT method,
    and aggregates results to provide insights into the contributions of RNA-binding proteins (RBPs)
    to gene expression.

    Attributes:
        logger (Logger): Logger instance for tracking the progress and results of the explanation process.
        config_explain (ConfigParser): Configuration parser instance for the explanation process, containing settings and paths.
        dataset (DeepRBPExpressionDataset): An instance of DeepRBPExpressionDataset containing RBP, gene, and transcript expression data.
        deeplift_explainer (DeepLift): Instance of the DeepLift explainer used to compute attribution scores.
        
        
        verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during training.
                                 - 0: No logging (suppress device information and other logs).
                                 - 1: Basic logging (show training progress and essential logs).
                                 - 2: Detailed logging (show device information).
                                 Default is 1.  
    """
    def __init__(self, config_explain, dataset, model, verbose=1):
        """
        Initialize the DeepLiftHandler with the model, dataset, explanation configuration, and gene mappings.

        Parameters:
            model: The trained machine learning model to be explained.
            config_explain (ConfigParser): Configuration parser instance containing settings and paths for the explanation process.
        """
        self.logger = Logger(verbose=verbose)
        self.config_explain = config_explain
        self.dataset = dataset
        self.deeplift_explainer = DeepLift(model)
        self.logger.log("✅ DeepLiftHandler initialized successfully.", level=1)
    ###
    def prepare_rbp_reference_tensor(self):
        """
        Prepare reference RBP tensor based on the configuration.
        Returns:
            reference_rbp_tensor (Tensor): Tensor used as a reference for RBP expression.
        """
        # Determine reference tensor based on configuration
        reference_type = self.config_explain.get('reference_data')
        if reference_type == 'median_reference':
            self.logger.log("[prepare_rbp_tensors] Using median of RBP expression as reference.")
            reference_rbp_tensor = torch.tensor(np.median(self.dataset.features['scaled_rbp_df'], axis=0), dtype=torch.float32)
            reference_rbp_tensor = torch.reshape(reference_rbp_tensor, (1, reference_rbp_tensor.size(0)))
            self.logger.log("[prepare_rbp_tensors] Median reference tensor created.")
        elif reference_type == 'knockout_reference':
            self.logger.log("[prepare_rbp_tensors] Using knockout (zeroed) RBP expression as reference.")
            reference_rbp_tensor = torch.zeros(1, self.dataset.features['scaled_rbp_df'].shape[1], dtype=torch.float32)
            self.logger.log("[prepare_rbp_tensors] knockout reference tensor created.")
        elif reference_type == 'half_reference':
            self.logger.log("[prepare_rbp_tensors] Using half-zeroed RBP expression as reference.")
            reference_rbp_tensor = torch.ones(1, self.dataset.features['scaled_rbp_df'].shape[1], dtype=torch.float32) * 0.5
            self.logger.log("[prepare_rbp_tensors] Half-zeroed reference tensor created.")
        else:
            self.logger.error(f"[prepare_rbp_tensors] Unknown reference type: '{reference_type}'", ValueError)
            raise ValueError(f"Unknown reference type: {reference_type}")
        return reference_rbp_tensor
    ###
    def compute_attribution_scores(self, reference_rbp_tensor):
        """
        Compute attribution scores batch-wise for RBP x S x T.

        Parameters:
            reference_rbp_tensor (torch.Tensor): The reference RBP tensor.
            
        Returns:
            list: A list of batch attribution scores.
        """
        scaled_rbp_tensor, gene_tensor = self.dataset.features['scaled_rbp_df'], self.dataset.features['gene_df']
        list_batch_scores = []
        for target_node in tqdm(range(self.dataset.features['isoform_df'].shape[1]), desc="Calculating scores", unit="node"):
            scores_batch = self._compute_attribution_scores_for_batch(
                target_node, scaled_rbp_tensor, reference_rbp_tensor, gene_tensor
            )
            list_batch_scores.append(scores_batch)
        return list_batch_scores
    ###
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
        self.logger.log(f"scaled_rbp_inputs dtype: {scaled_rbp_inputs.dtype}", level=2)
        self.logger.log(f"reference_rbp_inputs dtype: {reference_rbp_inputs.dtype}", level=2)
        self.logger.log(f"gene_expression_inputs dtype: {gene_expression_inputs.dtype}", level=2)
        scores = self.deeplift_explainer.attribute(
            inputs=scaled_rbp_inputs,
            baselines=reference_rbp_inputs,
            target=target_node,
            additional_forward_args=gene_expression_inputs,
        )
        return scores
    ###
    def reduce_batch_dimension(self, list_batch_scores):
        """
        Reduce the batch dimension of DeepLIFT scores to generate final TxRBP scores using a specified method.

        Parameters:
            list_batch_scores (list): List of computed attribution scores for a batch using DeepLift.
            config_explain (ConfigParser): Configuration parser instance
                               Expected key: 'batch_reduction_method' with options:
                               - 't-statistic': Computes the t-statistic across samples.
                               - 'sum_scores': Sums the scores across samples.

        Returns:
            pd.DataFrame: DataFrame containing the reduced DeepLIFT scores.
        """
        batch_reduction_type = self.config_explain.get('batch_reduction_method')
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
                result_scores = np.where(std != 0, mean / (std / np.sqrt(num_samples)), 0)
            zero_std_indices = np.argwhere(std == 0)
            if zero_std_indices.size > 0:
                self.logger.log(f"[reduce_batch_dimension] Found {len(zero_std_indices)} positions with std=0.")
                for idx in zero_std_indices:
                    trans = self.dataset.trans_names[idx[0]]
                    rbp = self.dataset.rbp_names[idx[1]]
                    self.logger.log(f"[reduce_batch_dimension] Zero std at Transcript: {trans}, RBP: {rbp}.")
            self.logger.log("[reduce_batch_dimension] T-statistic reduction completed.")
        elif batch_reduction_type == 'sum_scores':
            self.logger.log("[reduce_batch_dimension] Calculating sum of scores.")
            scores_stack = np.stack(batch_scores_np, axis=0)
            result_scores = np.sum(scores_stack, axis=1)
            self.logger.log("[reduce_batch_dimension] Sum reduction completed.")
        df_deeplift_TxRBP = pd.DataFrame(result_scores, index=self.dataset.trans_names, columns=self.dataset.rbp_names)
        return df_deeplift_TxRBP
    ###
    def calculate_scores_transcript_level(self):
        """
        Calculate attribution scores at the transcript level using the DeepLIFT method.

        Returns:
            pd.DataFrame: A DataFrame containing the calculated DeepLIFT scores at the transcript level.
        """
        # Prepare RBP reference tensor based on the selected configuration
        self.logger.log("🔧 Preparing RBP reference tensor...", level=1)
        reference_tensor = self.prepare_rbp_reference_tensor()  
        self.logger.log("✅ RBP reference tensor prepared successfully.", level=1)
        # Compute attribution scores for batch
        self.logger.log("🔢 Computing attribution scores for batch...", level=1)
        list_batch_scores = self.compute_attribution_scores(reference_tensor)  
        self.logger.log(f"✅ Computing attribution scores for batch.", level=1)
        # Reduce the batch dimension to obtain final scores (RBP x T)
        self.logger.log("🔽 Reducing batch dimension to aggregate scores...", level=1)
        df_scores_TxRBP = self.reduce_batch_dimension(list_batch_scores) 
        self.logger.log("✅ Batch dimension reduction complete. Scores aggregated.", level=1)
        return df_scores_TxRBP

