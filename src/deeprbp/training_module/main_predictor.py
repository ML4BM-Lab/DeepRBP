# src/deeprbp/training_module/main_predictor.py

import argparse
from .pipeline import DeepRBPredictorPipeline  # Adjust the import based on your module structure

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP predictor training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--external_config_path', type=str, help='Path to the external dataset configuration file.')
    return parser.parse_args()
    
def main():
    args = parse_args()
    pipeline = DeepRBPredictorPipeline(args.config_path, args.external_config_path)
    pipeline.run()

if __name__ == "__main__":
    main()

#python /scratch/jsanchoz/DeepRBP/src/deeprbp/training_module/main_predictor.py --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml" --external_config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_gtex.yaml"
