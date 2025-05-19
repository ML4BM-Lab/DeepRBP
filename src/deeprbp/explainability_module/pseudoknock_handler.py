# src/deeprbp/explanability_module/pseudoknock_handler.py

import numpy as np
import pandas as pd
import torch
from scipy import stats
#from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from ..util.logger import Logger
from ..util.utils import print_section_separator
 
class PseudoKnockHandler:
    """
    PseudoKnockHandler class for managing the explanation process using the in-silico knock pipeline in 
    DeepRBP explainability module.

    This class computes attribution scores using [insertar aqui]
    and aggregates results to provide insights into the contributions of RNA-binding proteins (RBPs)
    to gene expression.

    Attributes:
        logger (Logger): Logger instance for tracking the progress and results of the explanation process.
        config_explain (ConfigParser): Configuration parser instance for the explanation process, containing settings and paths.
        dataset (CustomTensorDataset): An instance of CustomTensorDataset containing RBP, gene, and transcript expression data.
        model: The model used for generating predictions.
        verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during training.
                                 - 0: No logging (suppress device information and other logs).
                                 - 1: Basic logging (show training progress and essential logs).
                                 - 2: Detailed logging (show device information).
                                 Default is 1.  
    """
    def __init__(self, config_explain, dataset, model, verbose=2):
        self.logger = Logger(verbose=verbose)
        self.config_explain = config_explain
        self.dataset = dataset
        self.model = model
        self.logger.log("✅ PseudoKnockHandler initialized successfully.", level=1)
    
    def modify_rbp_expression(self, condition, col):
        """
        Modifies the expression of a specific RBP for all patients based on the given condition.

        This function applies the following changes to the specified column in the RBP expression tensor:
        - 'Kout': Sets the expression to 0 (knockout).
        - 'half-Kout': Sets the expression to 0.5 (knockdown). AQUI NO SERÍA MEJOR RESTARLE A CADA PACIENTE SU EXPRESION 3 VECES?
        - 'Kup': Sets the expression to 1 (overexpression).
        - 'control': No changes are made.

        Args:
            condition (str): Condition to apply ('Kout', 'half-Kout', 'Kup', or 'control').
            col (int): Column index of the RBP to modify.

        Returns:
            torch.Tensor: The modified RBP expression tensor.
        """
        self.logger.log(f"🔧 Modifying RBP expression for condition '{condition}' on column index {col}.", level=2)
        knock_rbp_tensor = self.dataset.features['scaled_rbp_df'].clone() 
        if condition == 'Kout':
            knock_rbp_tensor[:, col] = 0  # Establecer toda la columna a 0
            self.logger.log(f'Performing a knockout (Kout) on column {col}', level=2)
        elif condition == 'half-Kout':
            knock_rbp_tensor[:, col] = 0.5  # Establecer toda la columna a 0.5
            self.logger.log(f'Performing a half-knockout on column {col}', level=2)
        elif condition == 'Kup':
            knock_rbp_tensor[:, col] = 1  # Establecer toda la columna a 1
            self.logger.log(f'Performing a knockup (Kup) on column {col}', level=2)
        else:
            self.logger.log(f'Condition is control, passing any change', level=2)
        return knock_rbp_tensor
    
    def calculate_log2fold_change(self, pred_condition1, pred_condition2):
        """
        Calculate the log2-fold change between two prediction conditions for a batch of samples (S).

        This function converts prediction values to Transcripts Per Million (TPM) 
        and computes the log2-fold change (log2FC) as the ratio of TPM values.

        Args:
            pred_condition1 (np.ndarray): Predictions for condition 1.
            pred_condition2 (np.ndarray): Predictions for condition 2.

        Returns:
            pd.DataFrame: A DataFrame containing the log2-fold change between condition 1 and condition 2.
        
        Note: 
        Assigning a Zero Value for the Fold Change
        If pred_tpm1 is 0 and pred_tpm2 is 0, it can be considered that there is no expression in either condition, and a fold change of 0 can be assigned.
        If pred_tpm1 is 0 and pred_tpm2 is greater than 0, the fold change can be assigned a value of 0, indicating that there is no expression in the first condition.
        If pred_tpm1 is greater than 0 and pred_tpm2 is 0, the fold change can be considered nan, indicating overexpression in condition 1. ASK ABOUT THIS
        """
        self.logger.log("🔢 Calculating log2-fold change between two prediction conditions...", level=2)
        pred_tpm1 = np.power(2, pred_condition1) - 1
        pred_tpm2 = np.power(2, pred_condition2) - 1
        # Calculate fold change with handling for zero values
        fold_change = np.zeros_like(pred_tpm1, dtype=float)
        fold_change[(pred_tpm1 > 0) & (pred_tpm2 > 0)] = pred_tpm1[(pred_tpm1 > 0) & (pred_tpm2 > 0)] / pred_tpm2[(pred_tpm1 > 0) & (pred_tpm2 > 0)]
        fold_change[(pred_tpm1 == 0) & (pred_tpm2 > 0)] = 0  # Condition 1 is zero, assign 0
        fold_change[(pred_tpm1 > 0) & (pred_tpm2 == 0)] = np.nan  # Condition 2 is zero, assign nan
        fold_change[(pred_tpm1 == 0) & (pred_tpm2 == 0)] = 0  # Both are zero, assign 0 (or you can use np.nan)
        # Create a DataFrame from the fold change results
        dataFC = pd.DataFrame(data=fold_change, columns=self.dataset.trans_names)
        # Count NaN values in the DataFrame
        nan_count = np.isnan(dataFC).sum().sum() 
        self.logger.log(f"Number of NaN values in fold change DataFrame: {nan_count}", level=2)    
        # Find and report the minimum and maximum value in dataFC
        self.logger.log(f"Minimum value in fold change DataFrame: {dataFC.min().min()}", level=2)  
        self.logger.log(f"Maximum value in fold change DataFrame: {dataFC.max().max()}", level=2)  
        # Calculate the log2(FC)
        data_log2FC = np.log2(dataFC)
        self.logger.log("✅ Log2-fold change calculation completed.", level=2)
        return data_log2FC
    
    def perform_ttest_on_transcripts(self, data_log2FC): 
        """
        Perform a one-sample t-test on transcript log2FC sample (S) data 
        (for the mean of ONE group of scores.).

        This function melts the provided DataFrame of log2 fold change values 
        into a long format, performs a one-sample t-test for each transcript 
        to determine if the mean expression differs significantly from zero, 
        and calculates adjusted p-values to account for multiple comparisons.

        Parameters:
        data_log2FC (pd.DataFrame): A DataFrame containing log2 fold change values 
                                    of transcripts with transcripts in columns and 
                                    samples in rows.

        Returns:
        pd.DataFrame: A DataFrame containing the transcript IDs, t-statistics, 
                    p-values, and adjusted p-values from the t-tests.
        """
        self.logger.log("🔍 Performing one-sample t-test on log2 fold change data...", level=2)
        # Melt fold change data into long format with Sample_ID, Transcript_ID, and log2FC
        melted_fc_data = pd.melt(data_log2FC.reset_index(), 
                            id_vars='index'
                        ).rename(columns={'index': 'Sample_ID', 'variable': 'Transcript_ID', 'value': 'log2FC'})
        # Perform t-test    
        result_ttest = (
            melted_fc_data.groupby('Transcript_ID')['log2FC']
            .apply(lambda x: stats.ttest_1samp(x, 0))
            .apply(pd.Series)
            .rename(columns={0: 't_stat', 1: 'p_value'})
            .reset_index()
        )
        # Filter out NaN values from t-test results
        result_ttest = result_ttest.dropna(subset=['p_value', 't_stat']).reset_index(drop=True)
        # Calculate adjusted p-values using Benjamini-Hochberg method
        #result_ttest['adj_p_value'] = multipletests(result_ttest['p_value'], method='fdr_bh')[1] # not used now but could be important
        self.logger.log("✅ T-test completed. Results prepared.", level=2)
        return result_ttest
    
    def compute_attribution_scores(self):
        """
        Compute attribution scores for RBP x T.

        Returns:
        pd.DataFrame: A DataFrame where rows are transcripts and columns are RBPs containing t-statistics.
        """
        self.logger.log("🔢 Initializing computation of attribution scores for RBPs...", level=1)
        # Initialize the DataFrame to hold the scores
        df_scores_TxRBP = pd.DataFrame(0, index=self.dataset.trans_names, columns=self.dataset.rbp_names)
        # Iterate through each RBP
        for index_col, rbp_name in tqdm(enumerate(self.dataset.rbp_names), total=len(self.dataset.rbp_names), desc="Processing RBPs"):
            self.logger.log(f"🔄 Processing RBP {index_col}: {rbp_name}...", level=2)
            df_scores_TxRBP[rbp_name] = self._compute_attribution_scores_for_rbp(index_col)
            self.logger.log(f"✅ Finished processing RBP {rbp_name}.", level=2)
            print(df_scores_TxRBP)
            print_section_separator()
        self.logger.log("✅ Attribution scores computation completed for all RBPs.", level=1)
        return df_scores_TxRBP
    
    def _compute_attribution_scores_for_rbp(self, index_col):
        """
        Computes transcript (T) attribution scores for a specific RBP.

        Parameters:
        index_col: Index of the RBP in the dataset.

        Returns:
        pd.Series: A Series containing t-statistics for each transcript.
        """
        # Get the modified expression tensors for each condition
        tensor_condition1, tensor_condition2 = (
            self.modify_rbp_expression(self.dataset, self.config_explain.get(condition), index_col)
            for condition in ('condition1', 'condition2')
        )
        # Check if tensors are equal
        if np.allclose(tensor_condition1, tensor_condition2, atol=1e-6):
            self.logger.log(f"Tensors for {self.dataset.rbp_names[index_col]} are equal, skipping calculations.", level=1)
            return pd.Series(0, index=self.dataset.trans_names)  # Return zeros if tensors are equal
        # Generate predictions
        with torch.no_grad():
            pred_condition1, pred_condition2 = (
                self.model(tensor, self.dataset.features['gene_df']).detach().numpy() 
                for tensor in (tensor_condition1, tensor_condition2)
            )
        # Calculate the log2-Fold-Change (FC)
        data_log2FC = self.calculate_log2fold_change(pred_condition1, pred_condition2)
        # Perform t-test and get results
        result_ttest = self.perform_ttest_on_transcripts(data_log2FC)
        print(result_ttest)
        # Create a temporary Series from the t-test results
        return pd.Series(result_ttest.set_index('Transcript_ID')['t_stat'])
    
    def calculate_scores_transcript_level(self):
        """
        Calculate attribution scores at the transcript level using the PseudoKnockHandler method.

        Returns:
            pd.DataFrame: A DataFrame containing the calculated PseudoKnockHandler scores at the transcript level.
        """
        # Compute attribution scores for transcripts
        self.logger.log("🔢 Computing attribution scores for transcripts...", level=1)
        df_scores_TxRBP = self.compute_attribution_scores()
        self.logger.log("✅ Attribution scores for transcripts computed successfully.", level=1)
        return df_scores_TxRBP
        
       
# # JOSEBA AQUI!!!
# config_path_train = '/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor/results/config.yaml'
# config_path_explain = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_explain_alternative.yaml'
# config_explain = ConfigParser(config_path_explain)
# output_dir = '/scratch/jsanchoz/DeepRBP/output/results/explainability'

# explainer = ExplainerModel(config_path_explain, config_path_train, output_dir)
# data = explainer.load_process_scale_data()
# dataset = explainer.build_tensor_dataset(data)

# model = explainer.load_trained_predictor_model()

# explainer_handler = explainer.initialize_explainer_handler(model)
 
# df_scores_TxRBP = compute_attribution_scores(explainer_handler.dataset, explainer_handler.model)


# explainer_handler.dataset
# explainer_handler.model
# explainer_handler.config_explain




