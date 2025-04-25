# src/deeprbp/training_module/main_predictor.py

import argparse
from .pipeline import DeepRBPredictorPipeline   

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP predictor training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results.')
    return parser.parse_args()
    
def main():
    args = parse_args()
    # Create an instance of the DeepRBPredictorPipeline with the provided configuration and output directory
    pipeline = DeepRBPredictorPipeline(args.config_path, args.output_dir)
    # Start the training process by running the pipeline
    pipeline.run()

if __name__ == "__main__":
    main()


