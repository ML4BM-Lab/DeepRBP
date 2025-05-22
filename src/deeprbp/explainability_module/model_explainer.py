# src/deeprbp/explainability_module/model.py

import os
import pandas as pd
from collections import namedtuple

from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, Scaler
from ..util.logger import Logger
from ..util.utils import DeepRBPExpressionDataset, filter_data_by_sample_ids, get_gene_info

from ..training_module.model import PredictorModel
from .deeplift_handler import DeepLiftHandler
from .pseudoknock_handler import PseudoKnockHandler

class ExplainerModel(): 
    """
    Class for handling the explainability of models using various methods.

    Attributes:
        logger (Logger): An instance of the Logger for logging progress and messages.
        path_save_results (str): The directory path where results will be saved.
        config_explain (ConfigParser): Configuration parser for explainability settings.
        config_train (ConfigParser): Configuration parser used for training the predictor model.
        data_paths (dict): A dictionary containing file paths for loading data
        getBM (DataFrame): DataFrame mapping Transcript_IDs to their associated Gene_IDs.
        sample_category (str): Column name in metada that refers to samples' tumor type.
        select_category (str): The specific tumor type(s) you want to select.
        disease_condition (str): The column in the metadata to stratify on.
        select_condition (str): Conditions for filtering the samples, such as "Primary_Tumor" or potentially including "Solid_Tissue_Normal".
        data_importer (DataImporter): An instance for importing datasets.
        scaler (Scaler): An instance of the Scaler for data preprocessing.
        dataset (DeepRBPExpressionDataset): The dataset used for explainability, created from the loaded data.
        model (PredictorModel): The trained predictor model used for making predictions.
    """
    def __init__(self, config_path_explain, config_path_train, output_dir, verbose=1):
        """
        Initializes the ExplainerModel with configuration paths.

        Parameters:
            config_path_explain (str): Path to the YAML configuration file for explainability.
            config_path_train (str): Path to the YAML configuration file for model training.
            output_dir (str): Path to the directory where results will be saved.
            verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during execution.
                                     - 0: No logging.
                                     - 1: Basic logging (show progress and essential logs).
                                     - 2: Detailed logging (show additional information).
                                     Default is 1.
        """
        self.logger = Logger(verbose=verbose)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing the DeepRBP-Explainer model...", level=1)
        # Define paths for saving results
        self.path_save_results = os.path.join(output_dir, 'results')
        os.makedirs(self.path_save_results, exist_ok=True) 
        self.logger.log("✅ Directory for results is ready.", level=1)
        # Set the configuration based on yaml file path
        self.config_explain = ConfigParser(config_path_explain)
        self.config_train = ConfigParser(config_path_train)
        # Load data paths from the explain configuration
        self.data_paths = self.config_explain.get('data_paths')
        self.logger.log(f"Data paths: {self.data_paths}", level=1)
        # Load additional configurations
        self.getBM = pd.read_csv(self.config_train.get('getBM_path'))
        self.sample_category = self.config_explain.get('sample_category')
        self.select_category = self.config_explain.get('select_category')
        self.disease_condition = self.config_explain.get('disease_condition')
        self.select_condition = self.config_explain.get('select_condition')
        # Initialize DataImporter for importing data
        self.data_importer = DataImporter(self.data_paths)
        # Initialize Scaler object
        self.scaler = Scaler.load(self.config_explain.get('scaler_dir'))
        self.logger.log("✅ Trained scaler is ready.", level=1)
        # Initialize dataset and model attributes
        self.dataset = None
        self.model = None

    def load_process_scale_data(self):
        """Load, process, and scale the training data based on the specified sample category."""
        # Import training data
        data = self.data_importer.load()    
        self.logger.log("✅ Data import complete.", level=1)
        # Print the selected category and condition for user reference
        self.logger.log(f"Selected category for filtering: '{self.select_category}'", level=1)
        self.logger.log(f"Selected conditions for filtering: {self.select_condition}", level=1)
        # Filter just the desired tumor type samples
        selected_sample_ids = data['metadata_df'][
                (data['metadata_df'][self.sample_category] == self.select_category) &
                (data['metadata_df'][self.disease_condition].isin(self.select_condition))
            ].index.tolist()
        data = filter_data_by_sample_ids(data, selected_sample_ids)
        self.logger.log("✅ Data filtered by sample category and disease condition.", level=1)
        # Scale the RBP data
        self.logger.log("🔧 Scaling RBP data...", level=1)
        data['scaled_rbp_df'] = self.scaler.transform(data['rbp_df']) # Scale the RBP data 
        self.logger.log("✅ RBP data scaling complete.", level=1)
        self.logger.log("✅ Data has been load, processed, and scaled successfully.", level=1)
        return data  
    
    def build_tensor_dataset(self, data):
        """Create a DeepRBPExpressionDataset instance"""
        self.logger.log("Starting the creation of the DeepRBPExpressionDataset...", level=1)
        dataset =  DeepRBPExpressionDataset(
                data=data.copy(),
                getBM=self.getBM,
                rbp_data_key='scaled_rbp_df', 
                gene_data_key='gene_df', 
                transcript_data_key='isoform_df',
                trans_col_name=self.config_train.get('trans_col_name'),
                gene_col_name=self.config_train.get('gene_col_name')
            )
        self.logger.log("✅ DeepRBPExpressionDataset instance created successfully.", level=1)
        self.logger.log("Features data stored in DeepRBPExpressionDataset:", level=2)
        for feature_name, tensor in dataset.features.items():
            self.logger.log(f"{feature_name}: Shape: {tensor.shape}", level=2)
        return dataset
    
    def load_trained_predictor_model(self):  
        """Load the trained predictor model."""  
        self.logger.log("📦 Loading the trained predictor model...") 
        model = PredictorModel.load_model(  
            path_to_weights = os.path.join(self.config_explain.get('trained_model_dir'), self.config_explain.get('model_file')),  
            config = self.config_train,
            input_size = self.dataset.features['scaled_rbp_df'].shape[1], 
            output_size = self.dataset.features['isoform_df'].shape[1]
        )  
        self.logger.log("✅ Model loaded successfully.", level=1)  
        return model 
    
    def initialize_explainer_handler(self):
        """Initialize the ExplainerHandler with the loaded model and data."""
        self.logger.log("🔍 Initializing Explainer handler...", level=1)
        explain_method = self.config_explain.get('explanation_method') 
        if explain_method == "DeepLIFT":
            explainer_handler = DeepLiftHandler(self.config_explain, self.dataset, self.model)
        elif explain_method == "Pseudoknockdown":   
            explainer_handler = PseudoKnockHandler(self.config_explain, self.dataset, self.model)
        else:
            self.logger.error(f"❌ Unknown explanation method: {explain_method}", level=1)
            raise ValueError(f"Unknown explanation method: {explain_method}")
        self.logger.log(" ✅ Explainer handler initialized successfully.", level=1)
        return explainer_handler  
    
    def filter_scores_for_low_expressed_transcripts(self, deeplift_scores):
        """
        Filters the DeepLIFT scores for transcripts that never express

        Args:
            deeplift_scores (pd.DataFrame): DataFrame containing DeepLIFT RBP scores indexed by transcripts.
            
            Returns:
                pd.DataFrame: Updated DeepLIFT scores DataFrame with scores of low-expressed transcripts set to 0.
        """
        self.logger.log("🔍 Filtering scores for low-expressed transcripts...", level=1)
        # Convert log2-transcripts per million (log2p(tpm)) expression data to transcripts per million (tpm)
        trans_expr = 2 ** self.dataset.features['isoform_df'] - 1
        trans_expr_df = pd.DataFrame(trans_expr.numpy(), columns=self.dataset.trans_names)
        # Identify transcripts that have a total expression of 0 (never expressed)
        transcripts_never_expressed = (trans_expr_df.sum(axis=0) == 0)  # Transcripts that have a total expression of 0
        # Set DeepLIFT scores to 0 for the transcripts that are never expressed
        deeplift_scores.loc[transcripts_never_expressed, :] = 0
        # Count the number of transcripts that are never expressed
        num_never_expressed = transcripts_never_expressed.sum()
        self.logger.log(f"✅ Filter results: {num_never_expressed} transcripts never expressed (total expression = 0).", level=1)
        return deeplift_scores
    
    def filter_scores_for_low_expressed_genes(self, deeplift_scores, threshold=1):
        """
        Set low-expressed genes (mean expression < threshold) to 0 in the TxRBP scores DataFrame.

        Parameters:
            deeplift_scores (pd.DataFrame): DataFrame with DeepLIFT scores (TxRBP).
            threshold (float): Expression threshold (default=1 TPM).

        Returns:
            pd.DataFrame: Updated DeepLIFT scores with low-expressed genes set to 0.
        """
        self.logger.log("🔍 Filtering scores for low-expressed genes...", level=1)
        gene_expr_df = pd.DataFrame(self.dataset.features['gene_df'].numpy(), columns=self.dataset.trans_names)  # DataFrame with gene expression values in TPM. 
        # colnames are the transcript ids are in this case we are working with the expanded matrix
        low_expr_genes = gene_expr_df.mean() < threshold
        deeplift_scores.loc[low_expr_genes, :] = 0
        self.logger.log(f"✅ Low-expressed genes (mean expression < {threshold} TPM) have been excluded from the TxRBP scores.")
        return deeplift_scores
    
    def collapse_transcript_scores_to_genes(self, deeplift_scores):
        """
        Collapse DeepLIFT scores from transcript level to gene level by aggregating the scores for each gene-RBP pair. 
        Specifically, this method transforms the given DataFrame of scores into a long format, merges it with gene information, 
        and then aggregates to find the maximum absolute score for each gene-RBP pair. Finally, it creates a pivot table to summarize the results.

        Parameters:
        deeplift_scores (pd.DataFrame): DataFrame containing DeepLIFT scores at the transcript level,
                                                  where rows correspond to transcripts and columns 
                                                  correspond to RBPs. 
        Returns:
        GeneScores: A namedtuple containing:
            - result_table (pd.DataFrame): A DataFrame with RBP and transcript details.
            - df_scores_GxRBP (pd.DataFrame): DataFrame with DeepLIFT scores (GxRBP).
        """ 
        self.logger.log("🔍 Collapsing transcripts scores to genes...", level=1)
        GeneScores = namedtuple('GeneScores', ['result_table', 'df_scores_GxRBP']) 
        # Transform the wide format DataFrame into a long format
        deeplift_scores_long = deeplift_scores.stack().reset_index()
        deeplift_scores_long.columns = ['Transcript_ID', 'RBP_ID', 'Score']  
        # Get RBP names from their IDs
        deeplift_scores_long['RBP_name'] = get_gene_info(deeplift_scores_long['RBP_ID'], self.getBM, return_type='names')
        # Merge with gene information to get Gene_IDs and additional metadata
        deeplift_scores_long = deeplift_scores_long.merge(self.getBM, on='Transcript_ID', how='left')
        # Determine collapse method based on configuration
        collapse_type = self.config_explain.get('gene_collapse_method')
        self.logger.log(f'Determining gene collapse method {collapse_type} based on configuration.', level=1)
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
            df_scores_GxRBP = result_table.pivot_table(
                index='Gene_ID', 
                columns='RBP_ID', 
                values='Score', 
                aggfunc='first'
            )
        self.logger.log("✅ Collapsing transcripts scores to genes...", level=1)
        # Return the results as a namedtuple
        return GeneScores(result_table=result_table, df_scores_GxRBP=df_scores_GxRBP[self.dataset.rbp_names])
    
    def prepare_data_and_model(self):
        """Load and process data, and prepare the trained model for explainability."""
        self.logger.log("🚀 Starting data and model preparation...", level=1)
        # Load and process the data for explainability
        data = self.load_process_scale_data()
        # Build the input tensor dataset
        self.dataset = self.build_tensor_dataset(data)
        # Load the trained model
        self.model = self.load_trained_predictor_model()
        self.logger.log(" ✅ Data and model preparation done...", level=1)
        return self.dataset, self.model
    
    def run_explainability_pipeline(self):
        """Perform the explanation method based on the configured explainer handler and return the results."""
        self.logger.log("🚀 Starting the explanation process...")
        # Prepare data and model
        self.prepare_data_and_model()    
        # Initialize explainer
        explainer_handler = self.initialize_explainer_handler()
        # Calculate scores at transcript level
        df_scores_TxRBP = explainer_handler.calculate_scores_transcript_level()
        # Filter scores for low-expressed transcripts
        df_scores_TxRBP = self.filter_scores_for_low_expressed_transcripts(df_scores_TxRBP) #trans_expr_df=data['trans_expr_tpm_df']  
        # Filter scores for low-expressed genes
        df_scores_TxRBP = self.filter_scores_for_low_expressed_genes(df_scores_TxRBP, threshold=1) #gene_expr_df=data['gn_expr_each_iso_tpm_df'], 
        # Collapse scores to genes (RBP x G)
        results = self.collapse_transcript_scores_to_genes(df_scores_TxRBP)  
        # Save results as CSV files
        self.save_results(df_scores_TxRBP, results.df_scores_GxRBP, results.result_table)
        self.logger.log("✅ Explanation process completed successfully.")
        return {
            'df_scores_TxRBP': df_scores_TxRBP,
            'df_scores_GxRBP': results.df_scores_GxRBP,
            'result_table': results.result_table
        }
    
    def save_results(self, df_scores_TxRBP, df_scores_GxRBP, result_table):
        """Save the results as CSV files."""
        try:
            self.logger.log("💾 Saving results...")
            df_scores_TxRBP.to_csv(os.path.join(self.path_save_results, 'df_scores_TxRBP.csv'), index=True)
            df_scores_GxRBP.to_csv(os.path.join(self.path_save_results, 'df_scores_GxRBP.csv'), index=True)
            result_table.to_csv(os.path.join(self.path_save_results, 'result_table.csv'), index=True)
            self.logger.log("✅ Results saved successfully.")
        except Exception as e:
            self.logger.error(f"❌ Error saving results: {e}")


# # (DEBUGGING) Load and process the data for explainability
# config_path_train = '/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor/results/config.yaml'
# config_path_explain = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_explain.yaml'
# alternative: config_path_explain = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_explain_alternative.yaml'

# output_dir = '/scratch/jsanchoz/DeepRBP/output/results/explainability_deeplift'
# explainer = ExplainerModel(config_path_explain, config_path_train, output_dir)

# data = explainer.load_process_scale_data()
# dataset = explainer.build_tensor_dataset(data)

# #scaled_rbp_tensor = dataset.features['scaled_rbp_df'] # TENSORS ALREADY CREATED tensor_dataset.features['scaled_rbp_df'].shape[1]
# # gene_tensor = dataset.features['gene_df']
# # transcript_tensor = dataset.features['isoform_df']                             tensor_dataset.features['isoform_df'].shape[1]
# model = explainer.load_trained_predictor_model()

# explainer_handler = explainer.initialize_explainer_handler(model)
# # Prepare RBP reference tensor
# reference_tensor = explainer_handler.prepare_rbp_reference_tensor(dataset) # esto es solo para deeplift???
# # Compute attribution scores
# list_batch_scores = explainer_handler.compute_attribution_scores(dataset, reference_tensor) # POR AQUI BROTHER!!!    
# # Reduce batch dimension (RBP x T)    
# df_scores_TxRBP = explainer_handler.reduce_batch_dimension(list_batch_scores)    
# # Filter scores for low-expressed transcripts
# df_scores_TxRBP = explainer_handler.filter_scores_for_low_expressed_transcripts( ## AQUI YA ME DA ERROR!
#     deeplift_scores=df_scores_TxRBP)  #trans_expr_df=data['trans_expr_tpm_df']  
# # Filter scores for low-expressed genes
# df_scores_TxRBP = explainer_handler.filter_scores_for_low_expressed_genes(
#     deeplift_scores=df_scores_TxRBP, #gene_expr_df=data['gn_expr_each_iso_tpm_df'], 
#     threshold=1
# )
# # Collapse scores to genes (RBP x G)
# results = explainer_handler.collapse_transcript_scores_to_genes(df_scores_TxRBP) # aqui devolver un elemento?
# # result_table, df_scores_GxRBP
# ########################################################################################
#   # Define paths for saving results
#         # self.path_save_results = os.path.join(
#         #     self.base_config['output_dir'], 
#         #     'results',
#         #     f"{self.explain_config['explanation_method']}_{self.explain_config['reference_data']}_{self.explain_config['batch_reduction_method']}_{self.explain_config['gene_collapse_method']}"
#         # )
