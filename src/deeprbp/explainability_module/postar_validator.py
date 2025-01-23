# 
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from typing import Tuple
from sklearn.metrics import roc_curve, roc_auc_score
from logger import Logger
from utils import ensure_directory_exists
from deeprbp.module_training.plots import plot_distributions_and_roc_with_thresholds

class PostarValidator:
    """
    A class to validate the results of the ExplainerModel against POSTAR experimental data.

    This class is responsible for loading POSTAR data that contains GxRBP (Gene to RNA Binding Protein) relationships, 
    matching these scores with those computed by the ExplainerModel, and providing analysis and visualization of the results.

    Attributes:
        explainer_model (ExplainerModel): An instance of the ExplainerModel that has computed scores for GxRBP.
        explain_config (dict): Configuration parameters from the ExplainerModel.
        postar_score_genes (DataFrame): Raw DataFrame containing POSTAR GxRBP scores before alignment.
        postar_score_genes_with_nan (DataFrame): Aligned Postar DataFrame with NaN for non-matching genes/RBPs.
        deeplift_scores_genes (pd.DataFrame): DataFrame containing DeepLIFT scores at the gene level.
        logger (Logger): Logger instance for logging messages and errors.
    """
    def __init__(self, explainer_model, df_results_summary, verbose = 1):
        """
        Initializes the PostarValidator with the given explainer model and verbosity level.

        Args:
            explainer_model (ExplainerModel): An instance of the ExplainerModel to validate.
            verbose (int): The level of verbosity for logging (default is 1).
        """
        self.explainer_model = explainer_model
        self.explain_config = explainer_model.explain_config
        self.df_results_summary = df_results_summary
        self.postar_score_genes_with_nan = pd.DataFrame()
        self.deeplift_scores_genes = pd.DataFrame()
        self.df_count_rbps_per_gen = pd.DataFrame()
        self.df_count_genes_per_rbp = pd.DataFrame()
        self.list_rbps_postar = list()
        self.list_genes_postar = list()
        self.optimal_thresholds_df = pd.DataFrame()
        self.auc_df = pd.DataFrame()
        self.logger = Logger(verbose)
        self.logger.log("Initialized PostarValidator.")
    
    def validate_dataframe(self, df: pd.DataFrame, expected_columns: list):
        """Validate that the DataFrame contains the expected columns."""
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            self.logger.error(f"Missing columns in DataFrame: {missing_columns}")
            raise ValueError(f"DataFrame is missing the following columns: {missing_columns}")
    
    def load_postar_data(self):
        """Load POSTAR experimental data with GxRBP relationships.
        This method reads the POSTAR data from a CSV file specified in the configuration,
        renames the axes to 'RBP_ID' and 'Gene_ID', and stores the data in the postar_score_genes attribute.
        """
        try:
            self.logger.log("Loading POSTAR data...")
            postar_score_genes = pd.read_csv(
                os.path.join(self.explain_config['postar_matrix_path'], self.explain_config['postar_file']),
                index_col=0
            ).rename_axis('RBP_ID', axis=1).rename_axis('Gene_ID')
            self.logger.log("POSTAR data loaded successfully.")
            return postar_score_genes
        
        except FileNotFoundError as e:
            self.logger.error(f"File not found: {self.explain_config['postar_matrix_path']}/{self.explain_config['postar_file']}")
        except Exception as e:
            self.logger.error(f"Failed to load POSTAR data: {e}")
    
    def _match_scores_and_postar_data(self, postar_score_genes: pd.DataFrame, deeplift_scores_genes: pd.DataFrame):
        """Match DeepLIFT scores with POSTAR data and return aligned DataFrames.

        This method matches the provided DeepLIFT scores with the POSTAR scores,
        ensuring both DataFrames have the same shape and corresponding genes and RBPs.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Modified DeepLIFT scores DataFrame and POSTAR scores DataFrame
            with matching genes and RBPs.
        """
        # Find matching and non-matching genes and RBPs
        genes_match = [x for x in deeplift_scores_genes.index if x in postar_score_genes.index]
        genes_not_match = [x for x in deeplift_scores_genes.index if x not in postar_score_genes.index]
        rbps_match = [x for x in deeplift_scores_genes.columns if x in postar_score_genes.columns]
        rbps_not_match = [x for x in deeplift_scores_genes.columns if x not in postar_score_genes.columns]
        
        self.logger.log("Finding matching and non-matching genes and RBPs...")
        self.logger.log(f"Matching genes: {len(genes_match)} matched, {len(genes_not_match)} not matched.")
        self.logger.log(f"Matching RBPs: {len(rbps_match)} matched, {len(rbps_not_match)} not matched.")
        
        # Create a DataFrame for POSTAR scores with NaN values for non-matching RBPs and genes
        nan_df = pd.DataFrame(index=genes_not_match, columns=rbps_not_match, dtype=np.float32).fillna(np.nan)
        postar_score_genes_with_nan = pd.concat([postar_score_genes, nan_df], axis=0)
        
        # Reindex the DataFrames to align their shapes
        self.postar_score_genes_with_nan = postar_score_genes_with_nan.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
        self.deeplift_scores_genes = deeplift_scores_genes.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
        
        # Set the axis names
        self.postar_score_genes_with_nan.index.name = 'Gene_ID'
        self.postar_score_genes_with_nan.columns.name = 'RBP_ID'
        self.logger.log("Scores matched successfully.")
    
    def _add_postar_data_to_results(self):
        """
        Complete the results table with POSTAR information by melting the POSTAR scores DataFrame,
        determining known Gene_IDs and RBP_IDs, and merging the results with df_results_summary.
        """
        melted_postar = self.postar_score_genes_with_nan.reset_index().melt(id_vars='Gene_ID', var_name='RBP_ID', value_name='Postar_Score')
        
        # Validate the melted DataFrame
        self.validate_dataframe(melted_postar, ['Gene_ID', 'RBP_ID', 'Postar_Score'])
        
        # Determine if the Gene_ID is known: A Gene_ID is known if at least one score is not NaN.
        known_genes = melted_postar.groupby('Gene_ID')['Postar_Score'].transform(lambda x: x.notna().any())
        melted_postar['known_Gene'] = known_genes
        
        # Determine if the RBP_ID is known: An RBP_ID is known if at least one score is not NaN.
        known_rbps = melted_postar.groupby('RBP_ID')['Postar_Score'].transform(lambda x: x.notna().any())
        melted_postar['known_RBP'] = known_rbps
        
        # Combine the resulting table with df_results_summary
        self.df_results_summary = self.df_results_summary.merge(melted_postar, on=['Gene_ID', 'RBP_ID'], how='left')
        self.logger.log("Results successfully combined with POSTAR data.")
    
    def process_postar_data(self, postar_score_genes: pd.DataFrame, deeplift_scores_genes: pd.DataFrame):
        # Match DeepLIFT scores with POSTAR data and return aligned DataFrames
        self._match_scores_and_postar_data(postar_score_genes, deeplift_scores_genes)
        # Complete the results table with POSTAR information
        self._add_postar_data_to_results()
    
    def _count_classes(self, data: pd.DataFrame, axis: int) -> pd.DataFrame:
        """
        Count the occurrences of each class (0, 1, NaN) in the given DataFrame along the specified axis.

        Args:
            data (pd.DataFrame): The DataFrame to analyze.
            axis (int): The axis to apply the counting (0 for columns, 1 for rows).

        Returns:
            pd.DataFrame: A DataFrame containing the counts of each class.
        """
        self.logger.log(f"Counting classes along axis {axis}...")
        class_counts = pd.DataFrame()
        class_counts['Class 0'] = data.apply(lambda x: (x == 0).sum(), axis=axis)
        class_counts['Class 1'] = data.apply(lambda x: (x == 1).sum(), axis=axis)
        class_counts['Class NaN'] = data.apply(lambda x: x.isna().sum(), axis=axis)
        self.logger.log("Class counting completed.")
        return class_counts
    
    def count_and_sort_postar_matrix(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Analyze the aligned POSTAR matrix to count the number of RBPs per gene and the number of genes per RBP,
        and sort the results based on the number of Class 1 occurrences.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: DataFrames containing counts of RBPs per gene and genes per RBP,
                                                both sorted by the number of Class 1 occurrences.
        """
        self.logger.log("Counting and sorting POSTAR matrix...")
       
        # Count RBPs per Gene
        self.df_count_rbps_per_gen = self._count_classes(self.postar_score_genes_with_nan, axis=1)
        self.df_count_rbps_per_gen['Genes'] = self.df_count_rbps_per_gen.index
        self.df_count_rbps_per_gen = self.df_count_rbps_per_gen.reset_index(drop=True).sort_values(by='Class 1', ascending=False).reset_index(drop=True)
        
        # Count Genes per RBP
        self.df_count_genes_per_rbp = self._count_classes(self.postar_score_genes_with_nan, axis=0)
        self.df_count_genes_per_rbp['RBPs'] = self.df_count_genes_per_rbp.index
        self.df_count_genes_per_rbp = self.df_count_genes_per_rbp.reset_index(drop=True).sort_values(by='Class 1', ascending=False).reset_index(drop=True)
        self.logger.log("Count and sort completed.")
        
        # Store the ordered list of RBPs and Genes
        self.list_rbps_postar = self.df_count_genes_per_rbp['RBPs'].values.tolist()
        self.list_genes_postar = self.df_count_rbps_per_gen['Genes'].values.tolist()
        return self.df_count_rbps_per_gen, self.df_count_genes_per_rbp
    
    def calculate_rbp_thresholds(self, path_save_results) -> pd.DataFrame:
        """
        Calculate optimal RBP thresholds from the combined results DataFrame.

        Returns:
            pd.DataFrame: A DataFrame containing the optimal score thresholds for each RBP.
        """
        self.logger.log("[calculate_rbp_thresholds] Calculating optimal thresholds per RBP...")
        
        # Validate the results summary DataFrame
        self.validate_dataframe(self.df_results_summary, ['Postar_Score', 'RBP_ID', 'Score'])
        
        # Filter combined results to retain rows with valid Postar_Score values (0 or 1)
        combined_results_filtered = self.df_results_summary[self.df_results_summary['Postar_Score'].isin([0, 1])]
        
        # Calculate absolute scores
        self.logger.log('Using ABSOLUTE scores for calculating the threshold scores.')
        combined_results_filtered['Score'] = combined_results_filtered['Score'].abs()
        
        # Lists to store thresholds and AUCs
        list_thresholds = []
        list_aucs = []
        list_unique_rbps = combined_results_filtered['RBP_ID'].unique().tolist()
        threshold_figures_path = os.path.join(path_save_results, 'Threshold_figures')
        ensure_directory_exists(threshold_figures_path)
        
        for rbp_id in tqdm(list_unique_rbps, desc="Calculating Optimal Thresholds"):
            df_current_rbp = combined_results_filtered[combined_results_filtered['RBP_ID'] == rbp_id]
            
            # Calculate threshold
            fpr, tpr, thresholds = roc_curve(df_current_rbp['Postar_Score'], df_current_rbp['Score'])
            optimal_idx = np.argmax(tpr - fpr)
            optimal_threshold = thresholds[optimal_idx]
            list_thresholds.append({'RBP_ID': rbp_id, 'Optimal_Score_Threshold': optimal_threshold})
            
            # Calculate AUC
            auc_score = roc_auc_score(df_current_rbp['Postar_Score'], df_current_rbp['Score'])
            list_aucs.append({'RBP_ID': rbp_id, 'AUC': auc_score})

            # Verificar si hay ambos 0s y 1s en df_current_rbp antes de graficar
            if df_current_rbp['Postar_Score'].nunique() == 2:  
                plot_distributions_and_roc_with_thresholds(
                    df_current_rbp, 
                    rbp_id, 
                    optimal_threshold, 
                    fpr, 
                    tpr, 
                    optimal_idx, 
                    auc_score,
                    path_save=threshold_figures_path
                )
            else:
                self.logger.log(f"Skipping plot for RBP: {rbp_id} as it does not contain both classes.")
        
        self.optimal_thresholds_df = pd.DataFrame(list_thresholds)
        self.auc_df = pd.DataFrame(list_aucs)
        self.logger.log("Optimal thresholds calculated successfully.")
        return self.optimal_thresholds_df, self.auc_df
    
    def return_summary_results(self) -> pd.DataFrame:
        """Return the complete results summary DataFrame."""
        self.logger.log("Returning the complete results summary.")
        return self.df_results_summary
    
    def save_results(self, path_save_results: str) -> None:
        """
        Save the optimal thresholds and results summary DataFrames to CSV files.

        Args:
            path_save_results (str): The path where the results will be saved.
        """
        ensure_directory_exists(path_save_results)
        results_summary_path = os.path.join(path_save_results, 'df_results_summary.csv')
        thresholds_path = os.path.join(path_save_results, 'optimal_thresholds.csv')
        rbps_per_gen_path = os.path.join(path_save_results, 'count_rbps_per_gen.csv')
        genes_per_rbp_path = os.path.join(path_save_results, 'count_genes_per_rbp.csv')
        list_rbps_postar_ordered_path = os.path.join(path_save_results, 'list_rbps_postar_ordered.csv')
        list_genes_postar_ordered_path = os.path.join(path_save_results, 'list_genes_postar_ordered.csv')
        auc_results_path = os.path.join(path_save_results, 'auc_results.csv')
        
        self.df_results_summary.to_csv(results_summary_path, index=False)
        self.logger.log(f"Results summary saved to {results_summary_path}")
        self.optimal_thresholds_df.to_csv(thresholds_path, index=False)
        self.logger.log(f"Optimal thresholds saved to {thresholds_path}")
        self.df_count_rbps_per_gen.to_csv(rbps_per_gen_path, index=False)
        self.logger.log(f"Counts of RBPs per gene saved to {rbps_per_gen_path}")
        self.df_count_genes_per_rbp.to_csv(genes_per_rbp_path, index=False)
        self.logger.log(f"Counts of genes per RBP saved to {genes_per_rbp_path}")
        pd.DataFrame(self.list_rbps_postar, columns=['RBPs']).to_csv(list_rbps_postar_ordered_path, index=False)
        self.logger.log(f"RBPs list saved to: {list_rbps_postar_ordered_path}")
        pd.DataFrame(self.list_genes_postar, columns=['Genes']).to_csv(list_genes_postar_ordered_path, index=False)
        self.logger.log(f"Genes list saved to: {list_genes_postar_ordered_path}")
        self.auc_df.to_csv(auc_results_path, index=False)
        self.logger.log(f"Mean AUC results saved to {auc_results_path}")

# HAY QUE CALCULAR EL ROC TOTAL
