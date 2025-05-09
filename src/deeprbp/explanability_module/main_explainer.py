# src/deeprbp/explainability_module/main_explainer.py

import argparse
from .model_explainer import ExplainerModel  

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP explainer to calculate transcript x RBP and genes x RBP explainability scores.')
    parser.add_argument('--config_path_explain', type=str, required=True, help='Path to the configuration file for explainability.')
    parser.add_argument('--config_path_train', type=str, required=True, help='Path to the configuration file for training the predictor model.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory where the results will be saved.')  
    return parser.parse_args()
    
def main():
    args = parse_args()
    explainer = ExplainerModel(args.config_path_explain, args.config_path_train, args.output_dir)
    explainer.run_explainability_pipeline()

if __name__ == "__main__":
    main()
