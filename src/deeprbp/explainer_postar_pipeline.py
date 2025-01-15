
# src/deeprbp/explainer_postar_pipeline.py
import argparse
from .models import ExplainerModel
from .explainer_validator_postar import ExplainerValidatorPostar

def run(config_path_explain, config_path_train):
    # Initialize the explainer model
    explainer_model = ExplainerModel(config_path_explain=config_path_explain, config_path_train=config_path_train)
    
    # Load and process the data for explainability
    data = explainer_model.load_and_process_data()
    print("Data loaded and processed:")
    print(data)
    
    # Load the trained model
    model = explainer_model.load_trained_predictor_model()
    print("Trained model loaded:")
    print(model)

    # Perform the explainability
    outputs = explainer_model.perform_explainer()

    # Show the results
    df_transcript_scores = outputs['df_scores_TxRBP']
    df_gene_scores = outputs['df_scores_GxRBP']
 
    # Print results
    print("Transcript Scores DataFrame (TxRBP):")
    print(df_transcript_scores)
    print("Gene Scores DataFrame (GxRBP):")
    print(df_gene_scores)
 
    # Initialize the POSTAR validator
    validator = ExplainerValidatorPostar(explainer_model, outputs['result_table'])

    # Load POSTAR data using the explainer_model's configuration
    df_postar_scores = validator.load_postar_data()

    # Process POSTAR data
    validator.process_postar_data(df_postar_scores, df_gene_scores)

    # Count and reorder POSTAR data
    df_rbps_per_gene_count, df_genes_per_rbp_count = validator.count_and_sort_postar_matrix()
    print("Count number of RBPs regulating genes in POSTAR:")
    print(df_rbps_per_gene_count)
    print("Count number of genes regulating RBPs in POSTAR:")
    print(df_genes_per_rbp_count)

    # Calculate and return the summary results
    df_results_summary = validator.return_summary_results()
    print("Result Table:")
    print(df_results_summary)

    # Calculate thresholds for RBPs and AUC
    optimal_thresholds, auc_results = validator.calculate_rbp_thresholds(explainer_model.path_save_results)
    print("Optimal Thresholds DataFrame:")
    print(optimal_thresholds)
    print("Mean auc results")
    print(auc_results.mean())

    # Save these results
    validator.save_results(explainer_model.path_save_results)

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP explainer and validate the scores using POSTAR.')
    parser.add_argument('--config_path_explain', type=str, required=True, help='Path to the configuration file for explainability scores.')
    parser.add_argument('--config_path_train', type=str, required=True, help='Path to the configuration file for training the predictor model.')
    return parser.parse_args()
    
def main():
    args = parse_args()
    run(args.config_path_explain, args.config_path_train)

if __name__ == "__main__":
    main()

#python /scratch/jsanchoz/DeepRBP/src/deeprbp/explainer_postar_pipeline.py --config_path_explain "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_explain.yaml" --config_path_train "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml"
# run-deeprbp-explainer-postar --config_path_explain "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_explain.yaml" --config_path_train "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml"