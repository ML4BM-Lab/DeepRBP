# src/deeprbp/explainability_module/postar_validator.py

import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from typing import Tuple
from sklearn.metrics import roc_curve, roc_auc_score

from ..util.utils import ensure_directory_exists
from .results_visualization.plots import plot_distributions_and_roc_with_thresholds
from ..util.logger import Logger
from ..data_loading.config_loader import ConfigParser

class PostarValidator:
    """
    A class to validate the results of the ExplainerModel against POSTAR experimental data.

    This class is responsible for loading POSTAR data that contains GxRBP (Gene to RNA Binding Protein) relationships, 
    matching these scores with those computed by the ExplainerModel, and providing analysis and visualization of the results.

    Attributes:
        logger (Logger): Logger instance for logging messages and errors.
        config_parser (ConfigParser): Instance for parsing configuration files.
        base_config (dict): Base configuration parameters.
        explain_config (dict): Configuration parameters specific to explainability.
        df_results_summary (pd.DataFrame): DataFrame containing the combined results of validation.
        postar_score_genes_with_nan (pd.DataFrame): DataFrame of POSTAR scores with NaN for non-matching genes/RBPs.
        explain_scores_genes (pd.DataFrame): DataFrame containing explainability scores at the gene level.
        df_count_rbps_per_gen (pd.DataFrame): DataFrame counting the number of RBPs per gene.
        df_count_genes_per_rbp (pd.DataFrame): DataFrame counting the number of genes per RBP.
        list_rbps_postar (list): List of RBPs present in the POSTAR dataset.
        list_genes_postar (list): List of genes present in the POSTAR dataset.
        optimal_thresholds_df (pd.DataFrame): DataFrame containing optimal score thresholds for each RBP.
        auc_df (pd.DataFrame): DataFrame containing AUC results for RBPs.
        path_save_results (str): Path for saving results to CSV files.
    """
    def __init__(self, config_path_explain, verbose = 1):     
        """
        Initializes the PostarValidator with the specified configuration path and verbosity level.

        Args:
            config_path_explain (str): Path to the configuration file for explainability settings.
            verbose (int): Verbosity level for logging (default is 1).
        """
        self.logger = Logger(verbose)
        self.config_parser = ConfigParser(config_path_explain)
        self.base_config = self.config_parser.get_base_config()
        self.explain_config = self.config_parser.get_explainability_config()
        self.logger.log("Initialized PostarValidator. 🌟")

        # Initialize attributes
        self.df_results_summary = None
        self.postar_score_genes_with_nan = pd.DataFrame()
        self.explain_scores_genes = pd.DataFrame()

        self.df_count_rbps_per_gen = pd.DataFrame()
        self.df_count_genes_per_rbp = pd.DataFrame()
        self.list_rbps_postar = list()
        self.list_genes_postar = list()
        self.optimal_thresholds_df = pd.DataFrame()
        self.auc_df = pd.DataFrame()

        # Define paths for saving data and results
        self.path_save_results = os.path.join(
            self.base_config['output_dir'], 
            'results',
            f"{self.explain_config['explanation_method']}_{self.explain_config['reference_data']}_{self.explain_config['batch_reduction_method']}_{self.explain_config['gene_collapse_method']}"
        )
        ensure_directory_exists(self.path_save_results)

    def perform_validation(self, explainability_results):
        """
        Perform validation of the explainability results against POSTAR data.

        This method loads the POSTAR data, matches it with the explainability scores, 
        and calculates optimal thresholds and AUC scores for each RBP.

        Args:
            explainability_results (dict): Dictionary containing the results from the ExplainerModel, 
                                            including 'df_scores_GxRBP' and 'result_table'.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
                - optimal_thresholds_df (pd.DataFrame): DataFrame with the optimal score thresholds for each RBP.
                - auc_df (pd.DataFrame): DataFrame with AUC results for each RBP.
        
        Raises:
        ValueError: If the required keys are missing in explainability_results.
        """
        # Validate input data
        if 'df_scores_GxRBP' not in explainability_results or 'result_table' not in explainability_results:
            self.logger.error("Missing required keys in explainability_results.")
            raise ValueError("explainability_results must contain 'df_scores_GxRBP' and 'result_table'.")

        df_postar_scores = self._load_postar_data()
        self._match_scores_and_postar_data(df_postar_scores, explainability_results['df_scores_GxRBP'])
        self._add_postar_data_to_results(explainability_results['result_table'])
        self._count_and_sort_postar_matrix()
        optimal_thresholds, auc_results = self.calculate_rbp_thresholds()
        return optimal_thresholds, auc_results

    def _load_postar_data(self):
        """Load POSTAR experimental data with GxRBP relationships.

        This method reads the POSTAR data from a CSV file specified in the configuration,
        renames the axes to 'RBP_ID' and 'Gene_ID', and stores the data in the `postar_score_genes_with_nan` attribute.

        Returns:
            pd.DataFrame: DataFrame containing the loaded POSTAR data.
        """
        try:
            self.logger.log("Loading POSTAR data... 🔄")
            postar_score_genes = pd.read_csv(
                os.path.join(self.explain_config['postar_matrix_path'], self.explain_config['postar_file']),
                index_col=0
            ).rename_axis('RBP_ID', axis=1).rename_axis('Gene_ID')
            self.logger.log("POSTAR data loaded successfully. ✅")
            return postar_score_genes
        
        except FileNotFoundError as e:
            self.logger.error(f"File not found: {self.explain_config['postar_matrix_path']}/{self.explain_config['postar_file']}")
        except Exception as e:
            self.logger.error(f"Failed to load POSTAR data: {e}")
    
    def _match_scores_and_postar_data(self, postar_score_genes: pd.DataFrame, explain_scores_genes: pd.DataFrame):
        """Match explainability scores with POSTAR data.

        This method aligns the explainability scores with the POSTAR scores, ensuring both DataFrames
        have the same genes and RBPs, and fills non-matching scores with NaN.

        Args:
            postar_score_genes (pd.DataFrame): DataFrame containing POSTAR GxRBP scores.
            explain_scores_genes (pd.DataFrame): DataFrame containing explainability scores from the ExplainerModel.

        Returns:
            None
        """
        # Find matching and non-matching genes and RBPs
        self.logger.log("Finding matching and non-matching genes and RBPs... 🔍")
        genes_match = [x for x in explain_scores_genes.index if x in postar_score_genes.index]
        genes_not_match = [x for x in explain_scores_genes.index if x not in postar_score_genes.index]
        rbps_match = [x for x in explain_scores_genes.columns if x in postar_score_genes.columns]
        rbps_not_match = [x for x in explain_scores_genes.columns if x not in postar_score_genes.columns]
        
        self.logger.log(f"Matching genes: {len(genes_match)} matched, {len(genes_not_match)} not matched.")
        self.logger.log(f"Matching RBPs: {len(rbps_match)} matched, {len(rbps_not_match)} not matched.")
        
        # Create a DataFrame for POSTAR scores with NaN values for non-matching RBPs and genes
        nan_df = pd.DataFrame(index=genes_not_match, columns=rbps_not_match, dtype=np.float32).fillna(np.nan)
        postar_score_genes_with_nan = pd.concat([postar_score_genes, nan_df], axis=0)
        
        # Reindex the DataFrames to align their shapes
        self.postar_score_genes_with_nan = postar_score_genes_with_nan.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
        self.explain_scores_genes = explain_scores_genes.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
        
        # Set the axis names
        self.postar_score_genes_with_nan.index.name = 'Gene_ID'
        self.postar_score_genes_with_nan.columns.name = 'RBP_ID'
        self.logger.log("Scores matched successfully. 🎯")
    
    def validate_dataframe(self, df: pd.DataFrame, expected_columns: list):
        """Validate that the DataFrame contains the expected columns.

        Args:
            df (pd.DataFrame): The DataFrame to validate.
            expected_columns (list): List of expected column names.

        Raises:
            ValueError: If any expected columns are missing from the DataFrame.
        """
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            self.logger.error(f"Missing columns in DataFrame: {missing_columns}")
            raise ValueError(f"DataFrame is missing the following columns: {missing_columns}")
        
    def _add_postar_data_to_results(self, df_results_summary):
        """
        Complete the results table with POSTAR information.

        This method melts the POSTAR scores DataFrame to identify known Gene_IDs and RBP_IDs,
        and merges the results with the provided summary DataFrame.

        Args:
            df_results_summary (pd.DataFrame): Summary DataFrame containing results to be augmented with POSTAR data.

        Returns:
            None
        """
        melted_postar = self.postar_score_genes_with_nan.reset_index().melt(id_vars='Gene_ID', var_name='RBP_ID', value_name='Postar_Score')
        
        # Validate the melted DataFrame
        self.logger.log("Validating melted DataFrame for POSTAR data... 🔍")
        self.validate_dataframe(melted_postar, ['Gene_ID', 'RBP_ID', 'Postar_Score'])
        
        # Determine if the Gene_ID is known: A Gene_ID is known if at least one score is not NaN.
        known_genes = melted_postar.groupby('Gene_ID')['Postar_Score'].transform(lambda x: x.notna().any())
        melted_postar['known_Gene'] = known_genes
        #self.logger.log("Identified known Gene_IDs. 🧬")
        
        # Determine if the RBP_ID is known: An RBP_ID is known if at least one score is not NaN.
        known_rbps = melted_postar.groupby('RBP_ID')['Postar_Score'].transform(lambda x: x.notna().any())
        melted_postar['known_RBP'] = known_rbps
        self.logger.log("Identified known RBP_IDs. 🔬")
        
        # Combine the resulting table with df_results_summary
        self.df_results_summary = df_results_summary.merge(melted_postar, on=['Gene_ID', 'RBP_ID'], how='left')
        self.logger.log("Results successfully combined with POSTAR data. ✅")
    
    def _count_classes(self, data: pd.DataFrame, axis: int) -> pd.DataFrame:
        """
        Count the occurrences of each class (0, 1, NaN) in the given DataFrame along the specified axis.

        Args:
            data (pd.DataFrame): The DataFrame to analyze.
            axis (int): The axis to apply the counting (0 for columns, 1 for rows).

        Returns:
            pd.DataFrame: A DataFrame containing the counts of each class.
        """
        self.logger.log(f"Counting classes along axis {axis}... 🔍")
        class_counts = pd.DataFrame()
        class_counts['Class 0'] = data.apply(lambda x: (x == 0).sum(), axis=axis)
        class_counts['Class 1'] = data.apply(lambda x: (x == 1).sum(), axis=axis)
        class_counts['Class NaN'] = data.apply(lambda x: x.isna().sum(), axis=axis)
        self.logger.log("Class counting completed. ✅")
        return class_counts
    
    def _count_and_sort_postar_matrix(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Analyze the aligned POSTAR matrix to count the number of RBPs per gene and the number of genes per RBP,
        and sort the results based on the number of Class 1 occurrences.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: DataFrames containing counts of RBPs per gene and genes per RBP,
                                                both sorted by the number of Class 1 occurrences.
        """
        self.logger.log("Counting and sorting POSTAR matrix... ⏳")
       
        # Count RBPs per Gene
        self.df_count_rbps_per_gen = self._count_classes(self.postar_score_genes_with_nan, axis=1)
        self.df_count_rbps_per_gen['Genes'] = self.df_count_rbps_per_gen.index
        self.df_count_rbps_per_gen = self.df_count_rbps_per_gen.reset_index(drop=True).sort_values(by='Class 1', ascending=False).reset_index(drop=True)
        
        # Count Genes per RBP
        self.df_count_genes_per_rbp = self._count_classes(self.postar_score_genes_with_nan, axis=0)
        self.df_count_genes_per_rbp['RBPs'] = self.df_count_genes_per_rbp.index
        self.df_count_genes_per_rbp = self.df_count_genes_per_rbp.reset_index(drop=True).sort_values(by='Class 1', ascending=False).reset_index(drop=True)
        self.logger.log("Count and sort completed. ✅")
        
        # Store the ordered list of RBPs and Genes
        self.list_rbps_postar = self.df_count_genes_per_rbp['RBPs'].values.tolist()
        self.list_genes_postar = self.df_count_rbps_per_gen['Genes'].values.tolist()
        return self.df_count_rbps_per_gen, self.df_count_genes_per_rbp
    
    def calculate_rbp_thresholds(self) -> pd.DataFrame:
        """
        Calculate optimal RBP thresholds from the combined results DataFrame.

        This method computes the optimal score thresholds for each RBP based on the Area Under the Curve (AUC) for each RBP.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
                - optimal_thresholds_df (pd.DataFrame): DataFrame with the optimal score thresholds for each RBP.
                - auc_df (pd.DataFrame): DataFrame with AUC results for each RBP.
        """
        self.logger.log("[calculate_rbp_thresholds] Calculating optimal thresholds per RBP... 🔄")

        # Validate the results summary DataFrame
        self.validate_dataframe(self.df_results_summary, ['Postar_Score', 'RBP_ID', 'Score'])
        
        # Filter combined results to retain rows with valid Postar_Score values (0 or 1)
        combined_results_filtered = self.df_results_summary[self.df_results_summary['Postar_Score'].isin([0, 1])]
        
        # Calculate absolute scores
        self.logger.log('Using ABSOLUTE scores for calculating the threshold scores. 📊')
        #combined_results_filtered['Score'] = combined_results_filtered['Score'].abs()
        combined_results_filtered.loc[:, 'Score'] = combined_results_filtered['Score'].abs()  # Use .loc to avoid SettingWithCopyWarning
    
        # Lists to store thresholds and AUCs
        list_thresholds = []
        list_aucs = []
        list_unique_rbps = combined_results_filtered['RBP_ID'].unique().tolist()
        threshold_figures_path = os.path.join(self.path_save_results, 'Threshold_figures')
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

            if df_current_rbp['Postar_Score'].nunique() == 2:  # # Check if there are both 0s and 1s in df_current_rbp before plotting
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
                self.logger.log(f"Skipping plot for RBP: {rbp_id} as it does not contain both classes. ⚠️")

        self.optimal_thresholds_df = pd.DataFrame(list_thresholds)
        self.auc_df = pd.DataFrame(list_aucs)
        self.logger.log("Optimal thresholds calculated successfully. ✅")
        return self.optimal_thresholds_df, self.auc_df
    
    def return_summary_results(self) -> pd.DataFrame:
        """Return the complete results summary DataFrame.

        Returns:
            pd.DataFrame: The DataFrame containing the combined results summary.
        
        Raises:
            UserWarning: If called before performing validation, warning that results are not available.
        """
        if self.df_results_summary is None:
            self.logger.warn(
                "You have not executed 'perform_validation()'. The summary results are not available.",
                level=1
            )
        self.logger.log("Returning the complete results summary. 📋")
        return self.df_results_summary
    
    def save_results(self) -> None:
        """
        Save the optimal thresholds and results summary DataFrames to CSV files.

        This method saves various results including the summary DataFrame, optimal thresholds, 
        counts of RBPs per gene, counts of genes per RBP, and AUC results to specified paths.

        Returns:
            None
        """
        self.logger.log(f"Saving results to {self.path_save_results}... 📂")
        ensure_directory_exists(self.path_save_results)
        results_summary_path = os.path.join(self.path_save_results, 'df_results_summary.csv')
        thresholds_path = os.path.join(self.path_save_results, 'optimal_thresholds.csv')
        rbps_per_gen_path = os.path.join(self.path_save_results, 'count_rbps_per_gen.csv')
        genes_per_rbp_path = os.path.join(self.path_save_results, 'count_genes_per_rbp.csv')
        list_rbps_postar_ordered_path = os.path.join(self.path_save_results, 'list_rbps_postar_ordered.csv')
        list_genes_postar_ordered_path = os.path.join(self.path_save_results, 'list_genes_postar_ordered.csv')
        auc_results_path = os.path.join(self.path_save_results, 'auc_results.csv')
        
        self.df_results_summary.to_csv(results_summary_path, index=False)
        self.logger.log(f"Results summary saved to {results_summary_path} 📂")
        self.optimal_thresholds_df.to_csv(thresholds_path, index=False)
        self.logger.log(f"Optimal thresholds saved to {thresholds_path} 📂")
        self.df_count_rbps_per_gen.to_csv(rbps_per_gen_path, index=False)
        self.logger.log(f"Counts of RBPs per gene saved to {rbps_per_gen_path} 📂")
        self.df_count_genes_per_rbp.to_csv(genes_per_rbp_path, index=False)
        self.logger.log(f"Counts of genes per RBP saved to {genes_per_rbp_path} 📂")
        pd.DataFrame(self.list_rbps_postar, columns=['RBPs']).to_csv(list_rbps_postar_ordered_path, index=False)
        self.logger.log(f"RBPs list saved to: {list_rbps_postar_ordered_path} 📂")
        pd.DataFrame(self.list_genes_postar, columns=['Genes']).to_csv(list_genes_postar_ordered_path, index=False)
        self.logger.log(f"Genes list saved to: {list_genes_postar_ordered_path} 📂")
        self.auc_df.to_csv(auc_results_path, index=False)
        self.logger.log(f"Mean AUC results saved to {auc_results_path} 📂")


