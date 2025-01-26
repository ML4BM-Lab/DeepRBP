# src/deeprbp/explainability_module/model.py

import os

from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DatasetLoader, Scaler
from ..util.logger import Logger
from ..util.utils import ensure_directory_exists
from ..training_module.model import PredictorModel
from .deeplift_handler import DeepLiftHandler

class ExplainerModel(): 
    """
    Class for handling the explainability of models using various methods.

    Attributes:
        logger (Logger): Logger instance for logging progress and messages.
        config_parser (ConfigParser): Parser for configuration files.
        base_config (dict): Base configuration settings.
        explain_config (dict): Explainability configuration settings.
        training_config (dict): Model training configuration settings.
        path_save_results (str): Path to save results.
        data_importer (DataImporter): Data importer for loading datasets.
        data_loader (DatasetLoader): Loader for datasets.
        scaler (Scaler): Scaler for data preprocessing.
    """
    def __init__(self, config_path_explain, config_path_train):
        """
        Initializes the ExplainerModel with configuration paths.

        Parameters:
            config_path_explain (str): Path to the explainability configuration file.
            config_path_train (str): Path to the training configuration file.
        """
        self.logger = Logger(verbose=1)
        
        # Load configurations
        self.logger.log("📁 Loading configurations...")
        self.config_parser = ConfigParser(config_path_explain)
        self.base_config = self.config_parser.get_base_config()
        self.explain_config = self.config_parser.get_explainability_config()
        self.config_train_parser = ConfigParser(config_path_train)
        self.training_config = self.config_train_parser.get_model_training_config()
        
        # Define paths for saving data and results
        self.path_save_results = os.path.join(
            self.base_config['output_dir'], 
            'results',
            f"{self.explain_config['explanation_method']}_{self.explain_config['reference_data']}_{self.explain_config['batch_reduction_method']}_{self.explain_config['gene_collapse_method']}"
        )
        ensure_directory_exists(self.path_save_results)
        
        # Initialize data processing components
        self.logger.log("🔄 Initializing data processing components...")
        self.data_importer = DataImporter(self.base_config['data_paths'])
        self.data_loader = DatasetLoader(self.data_importer, self.base_config)
        self.scaler = Scaler.load(self.explain_config['scaler_path'])
        self.logger.log("✅ Data processing components initialized.")
    
    def load_and_process_data(self):
        """Load and preprocess the data for explainability."""
        self.logger.log("📊 Loading and processing data...")
        data = self.data_loader.load_data()
        data['scaled_rbp_expr_log2p_tpm_df'] = self.scaler.transform(data['rbp_expr_log2p_tpm_df'])
        self.logger.log("✅ Data loaded and processed successfully.")
        return data
   
    def load_trained_predictor_model(self):  
        """Load the trained predictor model."""  
        self.logger.log("📦 Loading the trained predictor model...") 
        model = PredictorModel.load_model(  
            path_to_weights=os.path.join(self.explain_config['trained_model_path'], self.explain_config['model_file']),  
            config=self.training_config  
        )  
        self.logger.log("✅ Model loaded successfully.")  
        return model  
    
    def initialize_explainer_handler(self, model, data):
        """Initialize the ExplainerHandler with the loaded model and data."""
        self.logger.log("🔍 Initializing explainer handler...")
        reference_type = self.explain_config.get('explanation_method')
        if reference_type == "DeepLIFT":
            explainer_handler = DeepLiftHandler(model, data, self.base_config, self.explain_config)
            self.logger.log("✅ DeepLiftHandler initialized successfully.")
            self.logger.log("DeepLiftHandler initialized successfully.")
        elif reference_type == "Pseudoknocking":  # yet to develop (work in progress)
            self.logger.warn("⚠️ Pseudoknocking method is still in development.")
            # Handle other explanation methods as needed
        return explainer_handler
    
    def perform_explainer(self):
        """Perform the explanation method based on the configured explainer handler and return the results."""
        self.logger.log("🚀 Starting the explanation process...")
        
        # Load and process the data for explainability
        data = self.load_and_process_data()
        
        # Load the trained model
        model = self.load_trained_predictor_model()
        
        # Initialize explainer
        explainer_handler = self.initialize_explainer_handler(model, data)
        
        # Prepare RBP tensors
        self.logger.log("📐 Preparing RBP tensors...")
        scaled_rbp_tensor, gn_tensor, reference_rbp_tensor = explainer_handler.prepare_rbp_tensors()
        
        # Compute attribution scores
        self.logger.log("🔢 Computing attribution scores...")
        list_batch_scores = explainer_handler.compute_attribution_scores(scaled_rbp_tensor, reference_rbp_tensor, gn_tensor)
        
        # Reduce batch dimension (RBP x T)
        self.logger.log("🔽 Reducing batch dimension...")
        df_scores_TxRBP = explainer_handler.reduce_batch_dimension(list_batch_scores)
        
        # Filter scores for low-expressed transcripts
        self.logger.log("🔍 Filtering scores for low-expressed transcripts...")
        df_scores_TxRBP = explainer_handler.filter_scores_for_low_expressed_transcripts(
            deeplift_scores=df_scores_TxRBP,
            trans_expr_df=data['trans_expr_tpm_df']
        )
        
        # Filter scores for low-expressed genes
        self.logger.log("🔍 Filtering scores for low-expressed genes...")
        df_scores_TxRBP = explainer_handler.filter_scores_for_low_expressed_genes(
            deeplift_scores=df_scores_TxRBP,
            gene_expr_df=data['gn_expr_each_iso_tpm_df'], 
            threshold=1
        )
        
        # Collapse scores to genes (RBP x G)
        self.logger.log("📊 Collapsing scores to genes...")
        result_table, df_scores_GxRBP = explainer_handler.collapse_transcript_scores_to_genes(df_scores_TxRBP)
        
        # Save results as CSV files
        self.save_results(df_scores_TxRBP, df_scores_GxRBP, result_table)
        self.logger.log("✅ Explanation process completed successfully.")
        return {
            'df_scores_TxRBP': df_scores_TxRBP,
            'df_scores_GxRBP': df_scores_GxRBP,
            'result_table': result_table
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