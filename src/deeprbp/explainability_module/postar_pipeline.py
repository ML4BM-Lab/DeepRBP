# src/deeprbp/explainability_module/postar_pipeline.py

from ..util.logger import Logger
from .model import ExplainerModel
from .postar_validator import PostarValidator

class DeepRBPostarExplainabilityPipeline:
    """
    A pipeline class for calculating explainability scores using the DeepRBP predictor model in conjunction 
    with a tissue-specific Postar binary matrix. This class computes score thresholds for known RNA-binding proteins 
    (RBPs) based on the Postar dataset. It evaluates the area under the curve (AUC) between scores from the 
    0-Postar class and the 1-Postar class. Here, the 0-Postar class indicates that a specific RBP has shown 
    attachment to a particular gene in experimental CLIP (cross-linking immunoprecipitation) studies, while the 
    1-Postar class indicates the absence of such attachment.
    """
    def __init__(self, config_path_explain, config_path_train, verbose=1):
        self.logger = Logger(verbose)  # Initialize the logger with verbosity level
        self.logger.log("📁 Initializing the DeepRBPostarExplainabilityPipeline...", level=1)

        self.config_path_explain = config_path_explain
        self.config_path_train = config_path_train

        # Initialize the Explainer model class & POSTAR validator class
        self.explainer_model = ExplainerModel(self.config_path_explain, self.config_path_train)
        self.postar_validator = PostarValidator(self.config_path_explain)
            
    def display_results(self, optimal_thresholds, auc_results):
        """Display the results of the POSTAR analysis."""
        df_results_summary = self.postar_validator.return_summary_results()
        self.logger.log("Result Table:")
        print(df_results_summary)
        self.logger.log("Optimal Thresholds DataFrame:")
        print(optimal_thresholds)
        self.logger.log("Mean AUC results:")
        print(auc_results.AUC.mean())

    def run(self):
        self.logger.log("Running the explainability pipeline... 🔄")
        # Perform the explainability and get results
        explainability_results = self.explainer_model.perform_explainer()

        # Perform validation and processing of POSTAR data
        optimal_thresholds, auc_results = self.postar_validator.perform_validation(explainability_results)

        # Display summary results
        self.display_results(optimal_thresholds, auc_results)

        # Save results
        self.postar_validator.save_results()
