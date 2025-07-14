# src/deeprbp/training_module/benchmark_methods/train_baselines.py

import argparse
import torch
import pandas as pd
import os
import numpy as np
import time
import traceback

from ...data_loading.config_loader import ConfigParser
from ...data_preparation.data_module import DeepRBPDataModule

from .benchmark_utils import generate_X_y_data, create_multi_output_regressor, set_seed
from ...util.utils import setup_output_directory
from ..evaluation import calculate_general_metrics, spearmanr_per_gene

def parse_args():
    parser = argparse.ArgumentParser(description='Run the isoform abundance predictor using a selected ML algorithm as a benchmark.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results.')
    parser.add_argument('--algorithm', type=str, required=True, choices=[
        'svr', 'decision_tree', 'random_forest', 'gradient_boosting',
        'xgboost', 'lightgbm', 'knn', 'elastic_net', 'ridge'
    ], help='Algorithm to use for benchmarking.')
    return parser.parse_args()

def compute_metrics(y_true, y_pred, gene_names, trans_names, getBM):
    metrics = calculate_general_metrics(y_true.flatten(), y_pred.flatten())
    spearman = spearmanr_per_gene(
        gene_names=gene_names, 
        trans_names=trans_names, 
        getBM=getBM, 
        outputs=y_pred,
        labels=y_true
    )
    return metrics, spearman

def main():
    set_seed(42) 
    args = parse_args()
    
    # Load configuration and auxiliary file
    print('\n[train_baselines] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path)

    # Create method-specific output directory
    method_output_dir = os.path.join(args.output_dir, args.algorithm)
    output_dir = setup_output_directory(method_output_dir)
    print('\n[train_baselines] Output directory: ', output_dir)

    # Load data module and prepare data
    print('\n[train_baselines] 🚀 Initializing DataModule...')
    dm = DeepRBPDataModule(config, output_dir)
    print('[train_baselines] 🚀 Preparing data for training...')
    dm.setup('fit')  

    # Label log2(tpm+1) Values used only for final metric calculation
    y_log2p_tpm_train = dm.train_dataset.to_numpy()['isoform_df']  
    y_log2p_tpm_val = dm.val_dataset.to_numpy()['isoform_df']  
    gene_tpm_train = dm.train_dataset.to_numpy()['gene_df'] 
    gene_tpm_val = dm.val_dataset.to_numpy()['gene_df']  

    print('\n[train_baselines] 🚀 Generate X (input) and abundances (y) for modelling...')
    X_train, y_train = generate_X_y_data(dm.train_dataset.to_numpy(), calculate_abundance=True)  
    X_val, _ = generate_X_y_data(dm.val_dataset.to_numpy(), calculate_abundance=True)  

    gene_names = dm.train_dataset.gene_names
    trans_names = dm.train_dataset.trans_names

    # Benchmark selected method
    model_name = args.algorithm
    result_row = {'algorithm': model_name}
    start_time = time.time()

    try:
        print(f'\n[train_baselines] Evaluating: {model_name}')
        model = create_multi_output_regressor(model_name)

        print('\n[train_baselines] Fitting the model...')
        model.fit(X_train, y_train)

        print('\n[train_baselines] Creating predictions...')
        y_pred_train = model.predict(X_train)
        y_pred_val = model.predict(X_val) 

        print('\n[train_baselines] Scaling predicted abundances by gene TPM and applying log2 transformation...')
        tpm_pred_train = np.clip(y_pred_train * gene_tpm_train, 0, None)
        tpm_pred_val = np.clip(y_pred_val * gene_tpm_val, 0, None)
        y_log2_pred_train = np.log2(tpm_pred_train + 1)
        y_log2_pred_val = np.log2(tpm_pred_val + 1)

        print('\n[train_baselines] Calculating metrics for training and validation data...')
        train_metrics, train_spearman = compute_metrics(
            y_log2p_tpm_train, y_log2_pred_train, gene_names, trans_names, dm.getBM
        )
        val_metrics, val_spearman = compute_metrics(
            y_log2p_tpm_val, y_log2_pred_val, gene_names, trans_names, dm.getBM
        )

        result_row.update({
            'train_spearman': train_metrics['spearman_corr'],
            'train_pearson': train_metrics['pearson_corr'],
            'train_mse': train_metrics['mse'],
            'train_r2': train_metrics['r2'],
            'train_mean_corr_per_gene': train_spearman['mean_corr'],
            'train_mean_corr_max_transcript': train_spearman['mean_corr_max'],
            'val_spearman': val_metrics['spearman_corr'],
            'val_pearson': val_metrics['pearson_corr'],
            'val_mse': val_metrics['mse'],
            'val_r2': val_metrics['r2'],
            'val_mean_corr_per_gene': val_spearman['mean_corr'],
            'val_mean_corr_max_transcript': val_spearman['mean_corr_max'],
            'status': 'OK',
            'runtime_seconds': round(time.time() - start_time, 2)
        })

        print('\n[train_baselines] Algorithm execution done ✅')

    except Exception as e:
        result_row.update({
            'status': 'ERROR',
            'runtime_seconds': round(time.time() - start_time, 2),
            'error_message': traceback.format_exc()
        })
        print(f'\n[train_baselines] ❌ Error with {model_name}: {str(e)}\n')

    # Save results
    results_df = pd.DataFrame([result_row]).set_index("algorithm")
    results_df.to_csv(os.path.join(output_dir, "benchmark_results.csv"))
    print('\n[train_baselines] ✅ Results saved to benchmark_results.csv')

if __name__ == "__main__":
    main()
