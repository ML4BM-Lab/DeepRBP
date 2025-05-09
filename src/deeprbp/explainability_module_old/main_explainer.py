# src/deeprbp/explainability_module/main_explainer.py

import argparse
from .postar_pipeline import DeepRBPostarExplainabilityPipeline  

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP explainer and validate the scores using POSTAR.')
    parser.add_argument('--config_path_explain', type=str, required=True, help='Path to the configuration file for explainability scores.')
    parser.add_argument('--config_path_train', type=str, required=True, help='Path to the configuration file for training the predictor model.')
    return parser.parse_args()
    
def main():
    args = parse_args()
    pipeline = DeepRBPostarExplainabilityPipeline(args.config_path_explain, args.config_path_train)
    pipeline.run()

if __name__ == "__main__":
    main()

