# src/deeprbp/training_module/tumor_specific_training/main_specific_vs_general.py

import argparse
import pandas as pd
import os
import torch

from ...data_loading.config_loader import ConfigParser
from ...data_preparation.data_module import DeepRBPDataModule

from ...util.utils import (print_if_main, setup_output_directory,
                          print_gpu_memory_info, set_random_seed) # aqui tengo que decidir, el numero de epochs deberia ser para todos el mismo 

from ..utils_training import train_and_evaluate_model
from .utils_specific import collect_cross_tumor_metrics
from .plot_specific import plot_metric_confusion_matrix
from .constants import LIST_TUMOR_TYPES, BATCH_SIZE_BY_TUMOR

def main(): 
    args = parse_args()

    config = ConfigParser(
        train_path_files="/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Train",
        test_path_files="/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Test",
        getBM_path = "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv",
        gene_col_name = "Gene_ID",
        trans_col_name = "Transcript_ID",
        train_batch_size = 8, 
        val_batch_size = 8,
        sample_category="detailed_category",
        select_category="",
        test_fraction=0.2,
        seed=42,
        cuda = True,
        plot_results = False
        )
    print("Configuration loaded successfully.")

    # Training tumor-specific models
    for ttype in LIST_TUMOR_TYPES:
        print(f'\n🔧 [Tumor Type: {ttype}] Setting configuration and initializing...')
        # Clonamos la config base y actualizamos el tumor
        config.update('select_category', ttype)
        print_if_main(config)

        bs = BATCH_SIZE_BY_TUMOR.get(ttype, 8)
        config.update('train_batch_size', bs)
        config.update('val_batch_size', bs)
        print_if_main(f"[main_specific_vs_general] 🧪 Using batch size {bs} (train/val) for {ttype}")

        # Determine the output directory based on gpu rank or cpu device
        output_dir = setup_output_directory(os.path.join(args.output_base_dir, ttype))
        print_if_main('\n[main_specific_vs_general] Output directory for main process: ', output_dir)

        print_if_main('[main_specific_vs_general] 🚀 Initializing DataModule...')
        dm = DeepRBPDataModule(config, output_dir)

        # Starting model training and evaluation
        print_if_main("\n[main_predictor] 🧠 Starting model training and evaluation...")
        train_and_evaluate_model(dm, args, config, output_dir)
        print_if_main('\n')

    # ── sincroniza a todos antes de recolectar/plotear
    if torch.distributed.is_available() and torch.distributed.is_initialized():
        torch.distributed.barrier()

    def _is_main_process() -> bool:
        # Si estamos en DDP, sólo rank global 0
        if torch.distributed.is_available() and torch.distributed.is_initialized():
            return torch.distributed.get_rank() == 0
        # Si no hay grupo inicializado (ejecución 1 GPU / CPU),
        # usa RANK de torchrun si existe, y si no, asume proceso único.
        return os.environ.get("RANK", "0") == "0"

    # Collect results
    if _is_main_process:
        metric_names = ['spearman_corr', 'pearson_corr', 'mse', 'r2', 'mean_corr_per_gene', 'mean_corr_max_trans_per_gene']
        subset = ['Liver_Hepatocellular_Carcinoma', 'Kidney_Clear_Cell_Carcinoma', 'Acute_Myeloid_Leukemia']

        df_metrics = {
            metric: collect_cross_tumor_metrics(args.output_base_dir, metric, args.all_output_dir)
            for metric in metric_names
        }

        for metric_name, df_metric in df_metrics.items():
            # Plot full matrix
            plot_metric_confusion_matrix(
                df_metric, 
                metric_name=metric_name, 
                output_path=os.path.join(args.output_base_dir, f"{metric_name}.png")
            )
            
            # Plot subset
            df_subset = df_metric.loc[[*subset, 'all'], subset]
            plot_metric_confusion_matrix(
                df_subset, 
                metric_name=metric_name, 
                output_path=os.path.join(args.output_base_dir, f"{metric_name}_txiki.png")
            )
    
def parse_args():
    parser = argparse.ArgumentParser(description='Pipeline for training DeepRBP models on specific tumor types and comparing them against a general model trained on all tumor types.')
    parser.add_argument('--output_base_dir', type=str, default='/scratch/jsanchoz/DeepRBP/output/results/tumor_specific_training',
                        help='Base directory to save the training results.')
    parser.add_argument('--all_output_dir', type=str, help="Directory containing the 'test_tumor_category_results.csv' file for the model trained "
                "using all tumor types together. This model's metrics will be added as an extra row ('all') "
                "in the cross-tumor metric matrix.")
    parser.add_argument('--epochs', type=int, default=100, help='Training epochs for model training.')
    parser.add_argument('--num_workers', type=int, default=0, help='Number of workers for DataLoader.')
    parser.add_argument('--min_delta', type=float, default=0.001, 
                        help='Minimum change to qualify as an improvement (for early stopping).')
    parser.add_argument('--patience', type=int, default=30, help='How many epochs to wait after the last improvement (for early stopping).')
    parser.add_argument('--save_top_k', type=int, default=1,
                        help='The best k models according to validation loss will be saved. '
                             'If save_top_k == 0, no models are saved. '
                             'If save_top_k == -1, all models are saved.')
    parser.add_argument('--verbose', type=int, default=1, help='Verbosity level (0: silent, 1: normal, 2: debug).')
    return parser.parse_args()

if __name__ == "__main__":
    if os.environ.get("LOCAL_RANK")=="0": 
        print_gpu_memory_info()
    set_random_seed()

    # Set precision
    torch.set_float32_matmul_precision("high")
    main()