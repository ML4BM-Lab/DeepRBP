# src/deeprbp/training_module/main_predictor.py

import argparse
import os
import torch

from ..data_loading.config_loader import ConfigParser
from ..data_preparation.data_module import DeepRBPDataModule
from .utils_training import train_and_evaluate_model

from ..util.utils import (print_if_main, setup_output_directory, 
                          print_gpu_memory_info, set_random_seed)

def main():
    args = parse_args()
    
    # Load configuration and auxiliary file
    print_if_main('\n[main_predictor] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path) 

    # Determine the output directory based on gpu rank or cpu device
    output_dir = setup_output_directory(args.output_dir)
    print_if_main('\n[main_predictor] Output directory for main process: ', output_dir)

    # Load data module and prepare data for training
    print_if_main('\n[main_predictor] 🚀 Initializing DataModule...')
    dm = DeepRBPDataModule(config, output_dir)

    # Starting model training and evaluation
    print_if_main("\n[main_predictor] 🧠 Starting model training and evaluation...")
    train_and_evaluate_model(dm, args, config, output_dir)

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP predictor training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results.')
    parser.add_argument('--epochs', help='Training epochs for training', type=int, default=10)
    parser.add_argument('--num_workers', help='DataLoader number of workers', type=int, default=0)
    parser.add_argument('--min_delta', help='Minimum change to qualify as an improvement (for early stopping)', type=float, default=0.001)
    parser.add_argument('--patience', help='How many epochs to wait after the last improvement (for early stopping)', type=int, default=30)
    parser.add_argument('--save_top_k', 
                        help='The best k models according to the MSE validation will be saved. '
                             'If save_top_k == 0, no models are saved. '
                             'If save_top_k == -1, all models are saved.', 
                        type=int, 
                        default=1)
    parser.add_argument('--verbose', help='Verbosity level (0: no prints, 1: general prints, 2: debug prints)', type=int, default=1)
    return parser.parse_args()

if __name__ == "__main__":
    if os.environ.get("LOCAL_RANK")=="0": 
        print_gpu_memory_info()
    set_random_seed()
    # Set precision
    torch.set_float32_matmul_precision("high")
    main()