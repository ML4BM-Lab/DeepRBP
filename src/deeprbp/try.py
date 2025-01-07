### TRYING THE NEW CODE!!
# aqui hay que poner los puntos relativos.
import os
import numpy as np
import pandas as pd
from typing import Tuple
from logger import Logger
from models import PredictorModel, ExplainerModel

config_path_explain = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_explain.yaml'
config_path_train = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml'

# en explainer_validatior_postar.py
class ExplainerValidatorPostar:
    """
    A class to validate the results of the ExplainerModel against POSTAR experimental data.

    This class is responsible for loading POSTAR data that contains GxRBP (Gene to RNA Binding Protein) relationships, 
    matching these scores with those computed by the ExplainerModel, and providing analysis and visualization of the results.

    Attributes:
        explainer_model (ExplainerModel): An instance of the ExplainerModel that has computed scores for GxRBP.
        explain_config (dict): Configuration parameters from the ExplainerModel.
        postar_score_genes_raw (DataFrame): Raw DataFrame containing POSTAR  GxRBP scores before alignment.
        postar_score_genes_with_nan_aligned (DataFrame): Aligned Postar DataFrame with NaN for non-matching genes/RBPs.
        deeplift_scores_genes_aligned (DataFrame) : Aligned Explainer Dataframe for comparation with Postar matrix.
        logger (Logger): Logger instance for logging messages and errors.
    """
    def __init__(self, explainer_model, verbose = 1):
        """
        Initializes the ExplainerValidatorPostar with the given explainer model and verbosity level.

        Args:
            explainer_model (ExplainerModel): An instance of the ExplainerModel to validate.
            verbose (int): The level of verbosity for logging (default is 1).
        """
        self.explainer_model = explainer_model
        self.explain_config = explainer_model.explain_config
        self.postar_score_genes_raw = None
        self.postar_score_genes_with_nan_aligned = None
        self.deeplift_scores_genes_aligned = None
        self.logger = Logger(verbose)
    def load_postar_data(self):
        """Load POSTAR experimental data with GxRBP relationships.
        This method reads the POSTAR data from a CSV file specified in the configuration,
        renames the axes to 'RBP_ID' and 'Gene_ID', and stores the data in the postar_score_genes attribute.
        """
        try:
            self.logger.log("Loading POSTAR data...")
            self.postar_score_genes_raw = pd.read_csv(
                os.path.join(self.explain_config['postar_matrix_path'], self.explain_config['postar_file']),
                index_col=0
            ).rename_axis('RBP_ID', axis=1).rename_axis('Gene_ID')
            self.logger.log("POSTAR data loaded successfully.")
        except FileNotFoundError as e:
            self.logger.error(f"File not found: {e}", exception_type=FileNotFoundError)
        except Exception as e:
            self.logger.error(f"Failed to load POSTAR data: {e}", exception_type=RuntimeError)
    def match_scores_and_postar_data(self, deeplift_scores_genes: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Match DeepLIFT scores with POSTAR data and return aligned DataFrames.

        This method matches the provided DeepLIFT scores with the POSTAR scores,
        ensuring both DataFrames have the same shape and corresponding genes and RBPs.

        Args:
            deeplift_scores_genes (pd.DataFrame): DataFrame containing DeepLIFT scores at the gene level.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Modified DeepLIFT scores DataFrame and POSTAR scores DataFrame
            with matching genes and RBPs.
        """
        if self.postar_score_genes_raw is None:
            self.logger.error("POSTAR scores must be loaded before matching.", exception_type=RuntimeError)
        # Find matching and non-matching genes and RBPs
        genes_match = [x for x in deeplift_scores_genes.index if x in self.postar_score_genes_raw.index]
        genes_not_match = [x for x in deeplift_scores_genes.index if x not in self.postar_score_genes_raw.index]
        rbps_match = [x for x in deeplift_scores_genes.columns if x in self.postar_score_genes_raw.columns]
        rbps_not_match = [x for x in deeplift_scores_genes.columns if x not in self.postar_score_genes_raw.columns]
        self.logger.log("Finding matching and non-matching genes and RBPs...")
        # Create a DataFrame for POSTAR scores with NaN values for non-matching RBPs and genes
        nan_df = pd.DataFrame(index=genes_not_match, columns=rbps_not_match, dtype=np.float32).fillna(np.nan)
        postar_score_genes_with_nan = pd.concat([self.postar_score_genes_raw, nan_df], axis=0)
        # Reindex the DataFrames to align their shapes
        self.deeplift_scores_genes_aligned = deeplift_scores_genes.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
        self.postar_score_genes_with_nan_aligned = postar_score_genes_with_nan.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
        return self.deeplift_scores_genes_aligned, self.postar_score_genes_with_nan_aligned


# integrar esto: # 2.1) Analyze the matched Postar matrix for this cell line
df_count_rbps_per_gen, df_count_genes_per_rbp = count_and_sort_postar_matrix(matched_postar_scores)

def count_and_sort_postar_matrix(matched_postar_scores):
    """
    Analyze the POSTAR matrix to count the number of RBPs per gene and the number of genes per RBP,
    and sort the results based on the number of Class 1 occurrences.

    Parameters:
        matched_postar_scores (pd.DataFrame): DataFrame containing POSTAR scores with Gene_ID as index and RBP_ID as columns.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: DataFrames containing counts of RBPs per gene and genes per RBP,
                                            both sorted by the number of Class 1 occurrences.
    """
    # Create a filtered copy of the matched POSTAR scores
    matched_postar_scores_filtered = matched_postar_scores.copy()
    # Number of RBPs per Gene
    df_count_rbps_per_gen = pd.DataFrame()
    df_count_rbps_per_gen['Class 0'] = matched_postar_scores_filtered.apply(lambda x: (x == 0).sum(), axis=1)
    df_count_rbps_per_gen['Class 1'] = matched_postar_scores_filtered.apply(lambda x: (x == 1).sum(), axis=1)
    df_count_rbps_per_gen['Class NaN'] = matched_postar_scores_filtered.apply(lambda x: x.isna().sum(), axis=1)
    df_count_rbps_per_gen['Genes'] = df_count_rbps_per_gen.index
    df_count_rbps_per_gen = df_count_rbps_per_gen.reset_index(drop=True)
    # Sort the RBPs per Gene DataFrame by Class 1 in descending order
    df_count_rbps_per_gen = df_count_rbps_per_gen.sort_values(by='Class 1', ascending=False).reset_index(drop=True)
    # Number of genes per RBP
    df_count_genes_per_rbp = pd.DataFrame()
    df_count_genes_per_rbp['Class 0'] = matched_postar_scores_filtered.apply(lambda x: (x == 0).sum(), axis=0)
    df_count_genes_per_rbp['Class 1'] = matched_postar_scores_filtered.apply(lambda x: (x == 1).sum(), axis=0)
    df_count_genes_per_rbp['Class NaN'] = matched_postar_scores_filtered.apply(lambda x: x.isna().sum(), axis=0)
    df_count_genes_per_rbp['RBPs'] = df_count_genes_per_rbp.index
    df_count_genes_per_rbp = df_count_genes_per_rbp.reset_index(drop=True)
    # Sort the Genes per RBP DataFrame by Class 1 in descending order
    df_count_genes_per_rbp = df_count_genes_per_rbp.sort_values(by='Class 1', ascending=False).reset_index(drop=True)
    return df_count_rbps_per_gen, df_count_genes_per_rbp



############################################## llamada
    
# Uso de la clase
explainer_model = ExplainerModel(config_path_explain=config_path_explain, config_path_train=config_path_train)
data = explainer_model.load_and_process_data()
model = explainer_model.load_model()
outputs = explainer_model.perform_explainer()

# Acceso a los resultados
df_scores_TxRBP = outputs['df_scores_TxRBP']
df_scores_GxRBP = outputs['df_scores_GxRBP']
result_table = outputs['result_table']

# Imprimir resultados
print("Transcript Scores DataFrame (TxRBP):")
print(df_scores_TxRBP)
print("Gene Scores DataFrame (GxRBP):")
print(df_scores_GxRBP)
print("Result Table:")
print(result_table)

postar_validator = ExplainerValidatorPostar(explainer_model)
# Cargar datos de POSTAR
postar_validator.load_postar_data()

df_scores_GxRBP_aligned, df_postar_GxRBP_aligned = postar_validator.match_scores_and_postar_data(df_scores_GxRBP)

# /scratch/jsanchoz/DeepRBP/src/deeprbp/postar_utils.py
                              
# code to yet develop 

# 2.1) Analyze the matched Postar matrix for this cell line
df_count_rbps_per_gen, df_count_genes_per_rbp = count_and_sort_postar_matrix(matched_postar_scores)

# 3) Plot Explainer computed scores vs Postar

# 3) Plot DeepLIFT scores vs Postar (esto ayudarme de un x.py - piensa un nombre guay)
# ME HE QUEDADO AQUI BROTHER !!!
