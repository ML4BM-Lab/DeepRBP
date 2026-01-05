
# src/deeprbp/explanability_module/postar_validation/postar_validator.py
import argparse
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from typing import Optional
from sklearn.metrics import roc_curve, roc_auc_score

from ...util.logger import Logger
from .plot_utils import plot_distributions_and_roc_with_thresholds

class PostarValidator:
    """
    A class to validate the results of the explainability scores against POSTAR experimental data.

    This class is responsible for loading POSTAR data that contains GxRBP (Gene to RNA Binding Protein) relationships, 
    matching these scores with those computed by the explainability method used, and providing analysis and visualization of the results.

    Attributes:
        postar_matrix_dir (str): Directory path to the POSTAR matrix.
        postar_file (str): Filename of the POSTAR binary matrix (genes x RBPs).
        scores_result_dir (str): Directory path for the scores matrices and results table.
        verbose (int): Verbosity level for logging messages.
        logger (Logger): An instance of the Logger class for logging progress and messages.
        postar_data (DataFrame): DataFrame containing the POSTAR data loaded from the specified file.
        output_dir (str): Directory path where validation results will be saved. 
    """
    def __init__(self, postar_matrix_dir: str, postar_file: str, scores_result_dir: str, output_dir: str, verbose=1, pvalues_csv: Optional[str] = None):
        self.postar_matrix_dir = postar_matrix_dir
        self.postar_file = postar_file
        self.scores_result_dir = scores_result_dir
        self.path_save_results = os.path.join(output_dir, 'postar_validation')
        os.makedirs(output_dir, exist_ok=True)
        self.verbose = verbose
        self.logger = Logger(verbose=verbose)  # Initialize the logger
        self.df_count_rbps_per_gen = pd.DataFrame()
        self.df_count_genes_per_rbp = pd.DataFrame()
        self.updated_results_summary = pd.DataFrame()
        self.pvalues_csv = pvalues_csv
        self._pval_map = None   # dict: RBP_ID -> string ya formateada

    def load_postar_data(self):
        """Load the POSTAR data from the specified file."""
        postar_path = os.path.join(self.postar_matrix_dir, self.postar_file)
        try:
            self.logger.log("Loading POSTAR data...", level=1)
            postar_data = pd.read_csv(postar_path, index_col=0)   
            # Rename index and columns for clarity
            postar_data.rename_axis("gene_id", axis=0, inplace=True)  # Rename index to gene_id
            postar_data.rename_axis("rbp_gene_id", axis=1, inplace=True)  # Rename columns to rbp_gene_id
            self.logger.log("POSTAR data loaded successfully.", level=1)
            return postar_data
        except FileNotFoundError as e:
            self.logger.error(f"File not found: {postar_path}")
            raise e  # Re-raise the exception for further handling if necessary
        except Exception as e:
            self.logger.error(f"Failed to load POSTAR data: {e}")
            raise e  # Re-raise the exception for further handling if necessary

    def load_explainability_scores(self):
        """Load the explainability scores and results table from the specified directory."""
        try:
            df_scores_path = os.path.join(self.scores_result_dir, 'df_scores_GxRBP.csv')
            result_table_path = os.path.join(self.scores_result_dir, 'result_table.csv')
            self.logger.log("Loading explainability scores...", level=1)
            # Load the explainability scores
            df_scores_GxRBP = pd.read_csv(df_scores_path, index_col=0)   
            self.logger.log("Explainability scores loaded successfully.", level=1)
            # Load the results table
            results_summary = pd.read_csv(result_table_path, index_col=0)   
            self.logger.log("Results table loaded successfully.", level=1)
            return df_scores_GxRBP, results_summary
        except FileNotFoundError as e:
            self.logger.error(f"File not found: {e.filename}")
            raise e  # Re-raise the exception for further handling if necessary
        except Exception as e:
            self.logger.error(f"Failed to load explainability scores: {e}")
            raise e  # Re-raise the exception for further handling if necessary

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

    def _count_and_sort_postar_matrix(self, matched_postar_data):
        """
        Analyze the aligned POSTAR matrix to count the number of RBPs per gene and the number of genes per RBP,
        and sort the results based on the number of Class 1 occurrences.
        Args:
            matched_postar_data (pd.DataFrame): A DataFrame containing the aligned POSTAR matrix scores
        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: DataFrames containing counts of RBPs per gene and genes per RBP,
                                                both sorted by the number of Class 1 occurrences.
        """
        self.logger.log("Counting and sorting POSTAR matrix... ⏳")
        # Count RBPs per Gene
        self.df_count_rbps_per_gen = self._count_classes(matched_postar_data, axis=1)
        self.df_count_rbps_per_gen['Genes'] = self.df_count_rbps_per_gen.index
        self.df_count_rbps_per_gen = self.df_count_rbps_per_gen.reset_index(drop=True).sort_values(by='Class 1', ascending=False).reset_index(drop=True)
        # Count Genes per RBP
        self.df_count_genes_per_rbp = self._count_classes(matched_postar_data, axis=0)
        self.df_count_genes_per_rbp['RBPs'] = self.df_count_genes_per_rbp.index
        self.df_count_genes_per_rbp = self.df_count_genes_per_rbp.reset_index(drop=True).sort_values(by='Class 1', ascending=False).reset_index(drop=True)
        self.logger.log("Count and sort completed. ✅")

    def _match_scores_and_postar_data(self, postar_data: pd.DataFrame, scores_data: pd.DataFrame):
        """Match explainability scores with POSTAR data.

        This method aligns the explainability scores with the POSTAR scores, ensuring both DataFrames
        have the same genes and RBPs, and fills non-matching scores with NaN.

        Args:
            postar_data (pd.DataFrame): DataFrame containing POSTAR GxRBP scores.
            scores_data (pd.DataFrame): DataFrame containing explainability scores from the ExplainerModel.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Matched POSTAR data and matched scores data.
        """
        # Find matching and non-matching genes and RBPs
        self.logger.log("Finding matching and non-matching genes and RBPs... 🔍")
        genes_match = [x for x in scores_data.index if x in postar_data.index]
        self.genes_not_match = [x for x in scores_data.index if x not in postar_data.index]
        rbps_match = [x for x in scores_data.columns if x in postar_data.columns]
        rbps_not_match = [x for x in scores_data.columns if x not in postar_data.columns]
        self.logger.log(f"Matching genes: {len(genes_match)} matched, {len(self.genes_not_match)} not matched.", level=1)
        self.logger.log(f"Matching RBPs: {len(rbps_match)} matched, {len(rbps_not_match)} not matched.", level=1)
        # Create a DataFrame for POSTAR scores with NaN values for non-matching RBPs and genes
        nan_df = pd.DataFrame(index=self.genes_not_match, columns=rbps_not_match, dtype=np.float32).fillna(np.nan)
        matched_postar_data = pd.concat([postar_data.copy(), nan_df], axis=0)
        # Reindex the DataFrames to align their shapes
        matched_postar_data = matched_postar_data.loc[genes_match, rbps_match + rbps_not_match]
        matched_scores_data = scores_data.copy().loc[genes_match, rbps_match + rbps_not_match]
        # Set the axis names
        matched_postar_data.index.name = 'Gene_ID'
        matched_postar_data.columns.name = 'RBP_ID'
        # Log the genes that were not matched and removed
        if self.genes_not_match:
            self.logger.log(f"Genes that did not match: {', '.join(self.genes_not_match)}. They will be eliminated in _integrate_postar_into_summary internal method", level=1)
        self.logger.log("Scores matched successfully. ")
        return matched_postar_data, matched_scores_data

    def _integrate_postar_into_summary(self, matched_postar_data, results_summary):
        """
        Complete the results summary with POSTAR information.

        This method melts the POSTAR scores DataFrame and merges the results with the provided summary DataFrame.

        Args:
            matched_postar_data (pd.DataFrame): DataFrame containing matched POSTAR GxRBP scores.
            results_summary (pd.DataFrame): Summary DataFrame containing additional information on explainability results

        Returns:
            pd.DataFrame: A DataFrame containing the combined results with POSTAR data.
        """
        melted_postar = matched_postar_data.copy().reset_index().melt(id_vars='Gene_ID', var_name='RBP_ID', value_name='Postar_Score')
        # Determine if the RBP_ID is known
        known_rbps = melted_postar.groupby('RBP_ID')['Postar_Score'].transform(lambda x: x.notna().any())
        melted_postar['known_RBP'] = known_rbps
        self.logger.log("Identified known RBP_IDs. 🔬")
        # Combine the resulting table with df_results_summary
        updated_results_summary = results_summary.copy().merge(melted_postar, on=['Gene_ID', 'RBP_ID'], how='left')
        # Remove genes that did not match and log the action
        if self.genes_not_match:
            initial_count = updated_results_summary.shape[0]
            updated_results_summary = updated_results_summary[~updated_results_summary['Gene_ID'].isin(self.genes_not_match)]
            final_count = updated_results_summary.shape[0]
            self.logger.log(f"Removed {initial_count - final_count} entries from results summary that did not match with POSTAR data.", level=1)
        self.logger.log("Results successfully combined with POSTAR data. ✅")
        return updated_results_summary

    def process_postar_and_scores(self, postar_data: pd.DataFrame, scores_data: pd.DataFrame, results_summary: pd.DataFrame):
        """
        Match the POSTAR data with explainability scores, count and sort the POSTAR matrix,
        and integrate POSTAR data into the summary.

        Args:
            postar_data (pd.DataFrame): DataFrame containing POSTAR GxRBP scores.
            scores_data (pd.DataFrame): DataFrame containing explainability scores.
            results_summary (pd.DataFrame): Summary DataFrame containing additional information on explainability results

        Returns:
            pd.DataFrame: Updated results summary with integrated POSTAR information.
        """
        # Match POSTAR data with explainability scores
        matched_postar_data, _ = self._match_scores_and_postar_data(postar_data, scores_data)
        # Count and sort the POSTAR matrix
        self._count_and_sort_postar_matrix(matched_postar_data)
        # Integrate POSTAR into the summary
        updated_results_summary = self._integrate_postar_into_summary(matched_postar_data, results_summary)
        return updated_results_summary

    def _load_pvalues_map(self):
        """Carga p-values por RBP si se proporciona CSV. Devuelve dict RBP_ID -> texto."""
        if not self.pvalues_csv or not os.path.exists(self.pvalues_csv):
            return {}
        dfp = pd.read_csv(self.pvalues_csv)
        cols = [c for c in dfp.columns]
        key = 'p_adj' if 'p_adj' in cols else ('p_value' if 'p_value' in cols else ('p' if 'p' in cols else None))
        if key is None or 'RBP_ID' not in cols:
            self.logger.log(f"[p-values] CSV provided but required columns not found. "
                            f"Need 'RBP_ID' and one of ['p_adj','p_value','p']. Ignoring.", level=1)
            return {}
        # formatea p en notación científica corta y añade estrellas si quieres
        def _fmt(p):
            if pd.isna(p):
                return None
            try:
                p = float(p)
            except Exception:
                return None
            if p == 0:
                txt = "p<1e-300"
            elif p < 1e-3:
                txt = f"p={p:.1e}"
            else:
                txt = f"p={p:.3f}"
           
            stars = ("ns" if p >= 0.05 else ("*" if p >= 1e-2 else ("**" if p >= 1e-3 else ("***" if p >= 1e-4 else "****"))))
            return f"{'p_adj' if key=='p_adj' else 'p'}: {txt}{(' ' + stars) if stars else ''}"

        mp = {}
        for _, r in dfp.iterrows():
            txt = _fmt(r.get(key))
            if txt:
                mp[str(r['RBP_ID'])] = txt
        self.logger.log(f"[p-values] Loaded p-values for {len(mp)} RBPs.", level=1)
        return mp

    def calculate_rbp_thresholds(self, updated_results_summary: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate optimal RBP thresholds from the combined results DataFrame.

        This method computes the optimal score thresholds for each RBP based on the Area Under the Curve (AUC) for each RBP.

        Args: pd.DataFrame: Updated explainainability results summary with integrated POSTAR information.
        
        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
                - optimal_thresholds_df (pd.DataFrame): DataFrame with the optimal score thresholds for each RBP.
                - auc_df (pd.DataFrame): DataFrame with AUC results for each RBP.
        """
        self.logger.log("[calculate_rbp_thresholds] Calculating optimal thresholds per RBP... 🔄")
        self.updated_results_summary = updated_results_summary
        # Filter combined results to retain rows with valid Postar_Score values (0 or 1)
        combined_results_filtered = self.updated_results_summary.copy()[self.updated_results_summary['Postar_Score'].isin([0, 1])]
        # Calculate absolute scores
        self.logger.log('Using ABSOLUTE scores for calculating the threshold scores. 📊')
        combined_results_filtered.loc[:, 'Score'] = combined_results_filtered['Score'].abs()   
        
        # Lists to store thresholds and AUCs
        list_thresholds = []
        list_aucs = []
        list_unique_rbps = combined_results_filtered['RBP_ID'].unique().tolist()
        threshold_figures_path = os.path.join(self.path_save_results, 'threshold_figures')
        os.makedirs(threshold_figures_path, exist_ok=True)
        
        # Carga p-values si hay
        if self._pval_map is None:
            self._pval_map = self._load_pvalues_map()

        for rbp_id in tqdm(list_unique_rbps, desc="Calculating Optimal Thresholds"):
            df_current_rbp = combined_results_filtered[combined_results_filtered['RBP_ID'] == rbp_id]
            rbp_display_name = (df_current_rbp['RBP_name'].dropna().astype(str).mode().iat[0]
                        if 'RBP_name' in df_current_rbp.columns and not df_current_rbp['RBP_name'].dropna().empty
                        else rbp_id)
            
            # Calculate threshold
            fpr, tpr, thresholds = roc_curve(df_current_rbp['Postar_Score'], df_current_rbp['Score'])
            optimal_idx = np.argmax(tpr - fpr)
            optimal_threshold = thresholds[optimal_idx]
            list_thresholds.append({'RBP_ID': rbp_id, 'RBP_name': rbp_display_name, 'Optimal_Score_Threshold': optimal_threshold})
            
            # Calculate AUC
            auc_score = roc_auc_score(df_current_rbp['Postar_Score'], df_current_rbp['Score'])
            list_aucs.append({'RBP_ID': rbp_id, 'AUC': auc_score})
            
            if df_current_rbp['Postar_Score'].nunique() == 2:  # Check if there are both 0s and 1s in df_current_rbp before plotting
                
                stats_label = self._pval_map.get(str(rbp_id), None) # texto estadístico opcional para la leyenda
                plot_distributions_and_roc_with_thresholds(
                    df_current_rbp=df_current_rbp, 
                    rbp_id=rbp_id, 
                    optimal_threshold=optimal_threshold, 
                    fpr=fpr, 
                    tpr=tpr, 
                    optimal_idx=optimal_idx, 
                    auc_score=auc_score,
                    path_save=threshold_figures_path, 
                    rbp_display_name=rbp_display_name,  # << nombre bonito en el título y archivo
                    stats_label=stats_label
                )

            else:
                self.logger.log(f"Skipping plot for RBP: {rbp_id} as it does not contain both classes. ⚠️")

        self.optimal_thresholds_df = pd.DataFrame(list_thresholds)
        self.auc_df = pd.DataFrame(list_aucs)
        self.logger.log("Optimal thresholds calculated successfully. ✅", level=1)
        self.logger.log(f"Mean AUC results: {self.auc_df.AUC.mean()}", level=1)
        return self.optimal_thresholds_df, self.auc_df

    def save_results(self) -> None:
        """
        Save the optimal thresholds and results summary DataFrames to CSV files.

        This method saves various results including the summary DataFrame, optimal thresholds, 
        counts of RBPs per gene, counts of genes per RBP, and AUC results to specified paths.

        Returns:
            None
        """
        self.logger.log(f"Saving results to {self.path_save_results}... 📂")
        results_summary_path = os.path.join(self.path_save_results, 'result_table_completed.csv')
        thresholds_path = os.path.join(self.path_save_results, 'optimal_thresholds.csv')
        rbps_per_gen_path = os.path.join(self.path_save_results, 'count_rbps_per_gen.csv')
        genes_per_rbp_path = os.path.join(self.path_save_results, 'count_genes_per_rbp.csv')
        auc_results_path = os.path.join(self.path_save_results, 'auc_results.csv')
  
        self.updated_results_summary.to_csv(results_summary_path, index=False)
        self.logger.log(f"Results summary saved to {results_summary_path} 📂")
        self.optimal_thresholds_df.to_csv(thresholds_path, index=False)
        self.logger.log(f"Optimal thresholds saved to {thresholds_path} 📂")
        self.df_count_rbps_per_gen.to_csv(rbps_per_gen_path, index=False)
        self.logger.log(f"Counts of RBPs per gene saved to {rbps_per_gen_path} 📂")
        self.df_count_genes_per_rbp.to_csv(genes_per_rbp_path, index=False)
        self.logger.log(f"Counts of genes per RBP saved to {genes_per_rbp_path} 📂")
        self.auc_df.to_csv(auc_results_path, index=False)
        self.logger.log(f"Mean AUC results saved to {auc_results_path} 📂")

def parse_args():   
    parser = argparse.ArgumentParser(description="Validate explainability scores against POSTAR data.")
    parser.add_argument("--postar_matrix_dir", type=str, required=True, help="Directory path to the POSTAR matrix.")
    parser.add_argument("--postar_file", type=str, required=True, help="Filename of the POSTAR binary matrix (genes x RBPs).")
    parser.add_argument("--scores_result_dir", type=str, required=True, help="Directory path for the scores matrices and results table.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory path where validation results will be saved.")
    parser.add_argument("--verbose", type=int, default=1, help="Verbosity level for logging messages (default is 1).")
    parser.add_argument("--pvalues_csv", type=str, default=None, help="Optional CSV with per-RBP p-values (columns: RBP_ID + p_adj|p_value|p).")
    return parser.parse_args()

def main():
    args = parse_args()
    postar_validator = PostarValidator(
        postar_matrix_dir=args.postar_matrix_dir,
        postar_file=args.postar_file,
        scores_result_dir=args.scores_result_dir,
        output_dir=args.output_dir,
        verbose=args.verbose,
        pvalues_csv=args.pvalues_csv
    )   
    # Load POSTAR data
    postar_data = postar_validator.load_postar_data()
    # Load explainability scores
    scores_data, results_summary = postar_validator.load_explainability_scores()
    # Process POSTAR and scores together
    updated_results_summary = postar_validator.process_postar_and_scores(postar_data, scores_data, results_summary)
    # Calculate optimal thresholds
    postar_validator.calculate_rbp_thresholds(updated_results_summary)  
    # Save results
    postar_validator.save_results()

if __name__ == "__main__":
    main()
