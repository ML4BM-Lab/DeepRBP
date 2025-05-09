#src/deeprbp/explainability_module/deeplift_handler.py

import torch
from captum.attr import DeepLift
import numpy as np
import pandas as pd
from tqdm import tqdm

from ..util.logger import Logger
from ..util.utils import get_gene_info

import warnings
warnings.filterwarnings('ignore')
#warnings.filterwarnings("ignore", message="Input Tensor .* did not already require gradients")
warnings.filterwarnings("ignore", message="Setting forward, backward hooks and attributes on non-linear activations")
    
class DeepLiftHandler:
    """
    DeepLiftHandler class for managing the DeepLIFT explanation process for DeepRBP explainability module.

    This class is responsible for preparing the input data, computing attribution scores using the DeepLIFT method,
    and filtering and aggregating the results to provide insights into the contributions of RBPs to gene expression.

    Attributes:
        deeplift_explainer (DeepLift): Instance of the DeepLift explainer used to compute attribution scores.
        data (dict): A dictionary containing the necessary data for the explanation process, including scaled RBP
                     expression data and gene expression data.
        explain_config (dict): Configuration settings for the explanation process, such as reference data and
                               methods for reducing batch dimensions.
        base_config (dict): Configuration settings for loading data and saving results.
        logger (Logger): Logger instance for logging the progress and results of the explanation process.
        rbps_id (list): List of RNA-binding protein IDs.
        trans_id (list): List of transcript IDs.
        getBM (DataFrame): DataFrame containing additional gene information, loaded from specified paths.
        len_out_features (int): Number of output features (transcripts) for which scores will be computed.
    """
    def __init__(self, model, data, base_config, explain_config):
        """
        Initialize the DeepLiftHandler with model, data, and explain configuration.

        Parameters:
            model: The model to be explained. This should be an instance of a trained machine learning model that
                   supports explanations through the DeepLift method.
            data (dict): A dictionary containing necessary data for the explanation process.
                It should include:
                - 'scaled_rbp_expr_log2p_tpm_df': DataFrame containing the scaled RBP (RNA-binding protein) expression data.
                - 'gn_expr_each_iso_tpm_df': DataFrame containing the gene expression data.
            base_config (dict): A dictionary containing various configuration settings and paths that are
                                essential for loading data, saving results, and other operational parameters. 
            explain_config (dict): Configuration for explanations.
        """
        self.deeplift_explainer = DeepLift(model)
        self.data = data
        self.explain_config = explain_config
        self.base_config = base_config
        self.logger = Logger(verbose=1)

        # Define rbps_id and trans_id from the data provided
        self.rbps_id = list(data['scaled_rbp_expr_log2p_tpm_df'].columns)
        self.trans_id = list(data['gn_expr_each_iso_tpm_df'].columns)
        self.getBM = pd.read_csv(self.base_config['data_paths'].get('getBM_path', None)).drop_duplicates()
        self.len_out_features = len(self.trans_id)

    def prepare_data_tensors(self):
        """
        Prepare input tensors and reference RBP tensor based on the configuration.

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

        scaled_rbp_tensor = torch.tensor(self.data['scaled_rbp_expr_log2p_tpm_df'].values, dtype=torch.float32)
        gn_tensor = torch.tensor(self.data['gn_expr_each_iso_tpm_df'].values, dtype=torch.float32)
        self.logger.log("[prepare_rbp_tensors] Converted scaled RBP expression and gene expression to tensors (float32).")
        
        # Determine reference tensor based on configuration
        reference_type = self.explain_config.get('reference_data')

        if reference_type == 'median_reference':
            self.logger.log("[prepare_rbp_tensors] Using median of RBP expression as reference.")
            reference_rbp_tensor = torch.tensor(np.median(self.data['scaled_rbp_expr_log2p_tpm_df'], axis=0), dtype=torch.float32)
            reference_rbp_tensor = torch.reshape(reference_rbp_tensor, (1, reference_rbp_tensor.size(0)))
            self.logger.log("[prepare_rbp_tensors] Median reference tensor created.")
       
        elif reference_type == 'knockdown_reference':
            self.logger.log("[prepare_rbp_tensors] Using knockdown (zeroed) RBP expression as reference.")
            reference_rbp_tensor = torch.zeros(1, self.data['scaled_rbp_expr_log2p_tpm_df'].shape[1], dtype=torch.float32)
            self.logger.log("[prepare_rbp_tensors] Knockdown reference tensor created.")
        
        elif reference_type == 'half_reference':
            self.logger.log("[prepare_rbp_tensors] Using half-zeroed RBP expression as reference.")
            reference_rbp_tensor = torch.ones(1, self.data['scaled_rbp_expr_log2p_tpm_df'].shape[1], dtype=torch.float32) * 0.5
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
                    trans = self.trans_id[idx[0]]
                    rbp = self.rbps_id[idx[1]]
                    self.logger.log(f"[reduce_batch_dimension] Zero std at Transcript: {trans}, RBP: {rbp}.")
            
            df_deeplift_TxRBP = pd.DataFrame(t_stat_scores, index=self.trans_id, columns=self.rbps_id)
            self.logger.log("[reduce_batch_dimension] T-statistic reduction completed.")
        
        elif batch_reduction_type == 'sum_scores':
            self.logger.log("[reduce_batch_dimension] Calculating sum of scores.")
            scores_stack = np.stack(batch_scores_np, axis=0)
            sum_scores = np.sum(scores_stack, axis=1)
            df_deeplift_TxRBP = pd.DataFrame(sum_scores, index=self.trans_id, columns=self.rbps_id)
            self.logger.log("[reduce_batch_dimension] Sum reduction completed.")
        return df_deeplift_TxRBP
    
    def filter_scores_for_low_expressed_transcripts(self, deeplift_scores, trans_expr_df):
        """
        Filters the DeepLIFT scores for transcripts that never express

        Args:
            deeplift_scores (pd.DataFrame): DataFrame containing DeepLIFT RBP scores indexed by transcripts.
            trans_expr_df (pd.DataFrame): DataFrame containing expression levels of transcripts indexed by samples in TPM.

        Returns:
            pd.DataFrame: Updated DeepLIFT scores DataFrame with scores of low-expressed transcripts set to 0.
        """
        # Set scores to 0 for transcripts that never express (constant zero)
        transcripts_never_expressed = (trans_expr_df.sum(axis=0) == 0)  # Transcripts that have a total expression of 0
        deeplift_scores.loc[transcripts_never_expressed, :] = 0
        # Count the transcripts that never express
        num_never_expressed = transcripts_never_expressed.sum()
        self.logger.log(f"Filter results: {num_never_expressed} transcripts never expressed (total expression = 0).")
        return deeplift_scores
    
    def filter_scores_for_low_expressed_genes(self, deeplift_scores, gene_expr_df, threshold=1):
        """
        Set low-expressed genes (mean expression < threshold) to 0 in the TxRBP scores DataFrame.

        Parameters:
            deeplift_scores (pd.DataFrame): DataFrame with DeepLIFT scores (TxRBP).
            gene_expr_df (pd.DataFrame): DataFrame with gene expression values in TPM.
            threshold (float): Expression threshold (default=1 TPM).

        Returns:
            pd.DataFrame: Updated DeepLIFT scores with low-expressed genes set to 0.
        """
        low_expr_genes = gene_expr_df.mean() < threshold
        deeplift_scores.loc[low_expr_genes, :] = 0
        self.logger.log(f"Low-expressed genes (mean expression < {threshold} TPM) have been excluded from the TxRBP scores.")
        return deeplift_scores
    
    def collapse_transcript_scores_to_genes(self, deeplift_scores):
        """
        Collapse DeepLIFT scores from transcript level to gene level by aggregating the scores 
        for each gene-RBP pair. Specifically, this method transforms the given DataFrame of 
        scores into a long format, merges it with gene information, and then aggregates 
        to find the maximum absolute score for each gene-RBP pair. Finally, it creates a 
        pivot table to summarize the results.

        Parameters:
            deeplift_scores (pd.DataFrame): DataFrame containing DeepLIFT scores at the transcript level,
                                                where rows correspond to transcripts and columns 
                                                correspond to RBPs. 
        Returns:
            tuple: A tuple containing:
                - result_table (pd.DataFrame): A DataFrame with RBP and transcript details,
                - deeplift_scores_genes (pd.DataFrame): DataFrame with DeepLIFT scores (GxRBP).
        """
        # Transform the wide format DataFrame into a long format
        deeplift_scores_long = deeplift_scores.stack().reset_index()
        deeplift_scores_long.columns = ['Transcript_ID', 'RBP_ID', 'Score']  
        
        # Get RBP names from their IDs
        deeplift_scores_long['RBP_name'] = get_gene_info(deeplift_scores_long['RBP_ID'], self.getBM, return_type='names')
        
        # Merge with gene information to get Gene_IDs and additional metadata
        deeplift_scores_long = deeplift_scores_long.merge(self.getBM, on='Transcript_ID', how='left')
        
        # Determine collapse method based on configuration
        collapse_type = self.explain_config.get('gene_collapse_method')
        
        if collapse_type == 'max_absolute_value':
            deeplift_scores_long['Score_abs'] = deeplift_scores_long['Score'].abs() 
            
            # Find the index of the maximum score for each Gene-RBP pair
            max_indices = deeplift_scores_long.loc[deeplift_scores_long.groupby(['Gene_ID', 'RBP_ID'])['Score_abs'].idxmax()]  
            result_table = max_indices[['RBP_ID', 'RBP_name', 'Gene_ID', 'Gene_name', 
                                        'Transcript_ID', 'Transcript_name', 
                                        'Transcript_biotype', 'Score']].reset_index(drop=True)
            
            # Count the number of transcripts per Gene_ID
            num_transcripts_per_gene = self.getBM['Gene_ID'].value_counts().reset_index()
            num_transcripts_per_gene.columns = ['Gene_ID', 'Num_trans_per_gene']
            
            # Merge the count of transcripts with result_table
            result_table = result_table.merge(num_transcripts_per_gene, on='Gene_ID', how='left')
            
            # Create a pivot table to summarize scores by Gene_ID and RBP_ID
            deeplift_scores_genes = result_table.pivot_table(
                    index='Gene_ID', 
                    columns='RBP_ID', 
                    values='Score', 
                    aggfunc='first'
                )
        return result_table, deeplift_scores_genes[self.rbps_id]
