# /src/deeprbp/training_module/hyperparameter_optimization/grid_search_optuna.py

import psutil
import os
import logging
import sys
import argparse
import pandas as pd
import optuna
from optuna.samplers import TPESampler
#import optuna.visualization.matplotlib as optuna_plt ESTO PUEDO INTENTAR DESARROLLARLO TB DEL ARCHIVO DE LUIS
import optuna.logging
import time
import torch

from ..pipeline import DeepRBPredictorPipeline
from ...util.utils import print_section_separator

# Define the default configuration path and output directory
default_config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml'
default_output_dir = '/scratch/jsanchoz/DeepRBP/stuff/hyperpameter_optimization'

def parse_args():   
    parser = argparse.ArgumentParser(description='Hyperparameter optimization using Optuna.')
    parser.add_argument('-c', '--config_path_file', type=str, 
                        default=default_config_path,
                        help='Path to the configuration YAML file (default: %(default)s)')   
    parser.add_argument('-n', '--n_trials', type=int, required=True, help='Number of trials for Optuna')
    parser.add_argument('-o', '--output_dir', type=str, 
                        default=default_output_dir, 
                        help='Output directory for results (default: %(default)s)')
    return parser.parse_args()

def main():
    args = parse_args()
    # Load, process, and scale the data
    pipeline, train_dataset, val_dataset = load_and_process_data(args.config_path_file, args.output_dir)
    # pipeline, train_dataset, val_dataset = load_and_process_data(default_config_path, default_output_dir)

    # Optimization with Optuna using TPESampler and MedianPruner
    print("[*] Starting optimization with Optuna...")
    study_name = "DeepRBPredictor-optimization"
    optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))

    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5, n_startup_trials=5)    
    study = optuna.create_study(study_name=study_name, 
                                sampler=TPESampler(seed=pipeline.config.get('seed')), 
                                direction='minimize', 
                                pruner=pruner) # Recommended budgets with this sampler (#trials: 100-1000)

    # Execute the optimization
    print(f"[*] Running {args.n_trials} trials...")
    study.optimize(lambda trial: objective(trial, pipeline, train_dataset, val_dataset), n_trials=args.n_trials) 
    
    # Print the results of the best trial
    print(f'Best trial: {study.best_trial}, with parameters: {study.best_params} '
      f'and objective value: {study.best_value}\n')

    # Save the results to a CSV file
    df_results = study.trials_dataframe()
    print(df_results)
    results_file_path = os.path.join(args.output_dir, f'results_{study_name}.csv')
    df_results.to_csv(results_file_path, index=False)
    print(f"[*] Results saved to: {results_file_path}\n")

def load_and_process_data(config_path, output_dir):
    """
    Load, process, and scale the data, returning the datasets.

    Args:
        config (dict): Configuration object dictionary containing paths and parameters.
        output_dir (str): Output directory for saving processed data.

    Returns:
        tuple: pipeline, train_dataset, valid_dataset
    """
    # Load the data using the data importer and getBM that relates transcript_id with gene_id info
    pipeline = DeepRBPredictorPipeline(config_path, output_dir)
    data = pipeline.import_data() 
    # Filter a portion of the samples to optimize time and computational resources
    data_subset = pipeline.filter_samples(data)
    # Perform a train/validation split with the data subset
    train_data, val_data = pipeline.split_data(data_subset)
    # Save the training and validation samples used in the optimization to the specified output directory
    pipeline.save_split_data(train_data, val_data)
    # Scale data
    pipeline.fit_scaler(train_data) 
    scaled_train_data, scaled_val_data = pipeline.scale_data(train_data, val_data) 
    # Create custom data sets
    train_dataset, val_dataset = pipeline.create_datasets(scaled_train_data, scaled_val_data)
    return pipeline, train_dataset, val_dataset

def consolidate_metrics(train_metrics, val_metrics):
    """
    Consolidates training and validation metrics into a single dictionary.

    Args:
        train_metrics (pd.DataFrame): DataFrame containing training metrics.
        val_metrics (pd.DataFrame): DataFrame containing validation metrics.

    Returns:
        dict: A dictionary containing consolidated metrics with appropriate naming.
    """
    # Ensure the input DataFrames contain the required metrics
    required_columns = ['spearman_corr', 'pearson_corr', 'mse', 'r2', 'mean_corr_per_gene', 'mean_corr_max_trans_per_gene']
    for col in required_columns:
        if col not in train_metrics.columns or col not in val_metrics.columns:
            raise ValueError(f"Missing required metric '{col}' in input DataFrames.")
    # Create a combined dictionary using a dictionary comprehension
    combined_metrics = {
        **{f'train_{col}': train_metrics.loc['train', col] for col in required_columns},
        **{f'val_{col}': val_metrics.loc['val', col] for col in required_columns}
    }
    return combined_metrics

def objective(trial, pipeline, train_dataset, val_dataset):    
    process = psutil.Process()
    mem_before = process.memory_info().rss / (1024 ** 3)  # Memory before in GB

    # Suggest Optuna: Sample hyperparameters for this Trial   
    num_hidden_layers = trial.suggest_int('num_hidden_layers', 0, 4)
    hidden1_nodes = trial.suggest_categorical('hidden1_nodes', [64, 128, 256, 512, 1024, 2048, 4096])  
    uniform_nodes = trial.suggest_categorical('uniform_nodes', [True, False])
    node_shrink_factor = trial.suggest_categorical('node_shrink_factor', [2, 4, 8])
    activation_func = trial.suggest_categorical('activation_func', ["relu", "tanh", "sigmoid"])
    learning_rate = trial.suggest_categorical('learning_rate', [0.0001, 0.001, 0.01])
    optimizer_name = trial.suggest_categorical('optimizer_name', ['sgd90', 'asgd', 'adam', 'adagrad', 'adadelta', 'adamW'])
    train_batch_size = trial.suggest_categorical('batch_size', [32, 64, 128, 256, 512, 1024, 2048]) 
    num_epochs = trial.suggest_categorical('num_epochs', [50, 100, 500, 1000, 2000, 3000])   
    batch_norm_eps = trial.suggest_categorical('batch_norm_eps', [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]) 
    batch_norm_momentum = trial.suggest_categorical('batch_norm_momentum', [0.1, 0.5, 0.9]) 

    # prueba estos valores a ver
    # num_hidden_layers = 2; hidden1_nodes = 1024; uniform_nodes = False; node_shrink_factor = 2; 
    # activation_func = "relu"; optimizer_name = "adamW"; learning_rate = 0.0001; train_batch_size = 128; num_epochs = 50

    pipeline.config.update("num_hidden_layers", num_hidden_layers)
    pipeline.config.update("hidden1_nodes", hidden1_nodes)
    pipeline.config.update("uniform_nodes", uniform_nodes)
    pipeline.config.update("node_shrink_factor", node_shrink_factor)
    pipeline.config.update("activation_func", activation_func)
    pipeline.config.update("learning_rate", learning_rate)
    pipeline.config.update("optimizer_name", optimizer_name)
    pipeline.config.update("train_batch_size", train_batch_size)
    pipeline.config.update("num_epochs", num_epochs) 
    pipeline.config.update("batch_norm_eps", batch_norm_eps) 
    pipeline.config.update("batch_norm_momentum", batch_norm_momentum) 
    print(pipeline.config)
    
    # Log the hyperparameters for the current trial
    print(f"[*] Trial {trial.number}:")
    print(f"    - Number of Hidden Layers: {pipeline.config.get('num_hidden_layers')}")
    print(f"    - Hidden Layer 1 Nodes: {pipeline.config.get('hidden1_nodes')}")
    print(f"    - Use Uniform Nodes: {pipeline.config.get('uniform_nodes')}")
    print(f"    - Node Shrink Factor: {pipeline.config.get('node_shrink_factor')}")
    print(f"    - Activation Function: {pipeline.config.get('activation_func')}")
    print(f"    - Learning Rate: {pipeline.config.get('learning_rate')}")
    print(f"    - Optimizer: {pipeline.config.get('optimizer_name')}")
    print(f"    - Training Batch Size: {pipeline.config.get('train_batch_size')}")
    print(f"    - Number of Epochs: {pipeline.config.get('num_epochs')}")
    print(f"    - Batch Normalization Epsilon (eps): {pipeline.config.get('batch_norm_eps')}")
    print(f"    - Batch Normalization Momentum: {pipeline.config.get('batch_norm_momentum')}")
    print("\n")
    
    try:
        # Create DataLoader instances
        train_loader, val_loader = pipeline.create_data_loaders(train_dataset, val_dataset)
        # Create model instance
        pipeline.setup_model_trainer(train_loader)
        # Train model
        train_history, val_history = pipeline.train_model(train_loader, val_loader, trial)
        timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")   
        plot_name = f'trial_{trial.number}_{timestamp}' 
        pipeline.save_history(train_history, val_history, plot_name)

        # Evaluate the model
        train_metrics, _ = pipeline.evaluate_model(train_loader, set_name='train')
        val_metrics, _ = pipeline.evaluate_model(val_loader, set_name='val')
        metrics = consolidate_metrics(train_metrics, val_metrics)

        mem_after = process.memory_info().rss / (1024 ** 3)  # memory after in GB
        print(f"Memory used during this trial: {mem_after - mem_before:.4f} GB")

        print("\n")
        return val_history[-1]
    
    except Exception as e:
        # Catch any exception raised in the try block
        print(f"Training failed due to an error: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    
    finally:
        # Update trial.params if the corresponding pipeline.config value is 'unused'
        for key in ['hidden1_nodes', 'uniform_nodes', 'node_shrink_factor', 
                    'activation_func', 'batch_norm_eps', 'batch_norm_momentum']:
            # Check if the pipeline.config value is 'unused' and update trial.params
            if pipeline.config.get(key) == 'unused':
                trial.params[key] = pipeline.config.get(key)  # Overwrite trial parameter
        print(trial.params)
        print("\n")
        
        # Log metrics to the trial
        for key in metrics:
            trial.set_user_attr(key, metrics[key])

        # Cleanup: Free resources
        pipeline.trainer = None
        torch.cuda.empty_cache()  # If using GPU, clear the cache
        print_section_separator()

if __name__ == "__main__":
    main()

# The authors should provide the full table of results for the hyperparameter optimization runs
# to be able to validate the claim that more complex models (more hidden layers) are necessary.
# Verify that the suggested trial is elegible
# 1) si num_hidden_layers == 0 -> hidden1_nodes, uniform_nodes, node_shrink_factor y activation_func no aplican "not_used" (bien)
# 2) si num_hidden_layers == 1 -> uniform_nodes y node_shrink_factor no aplican "not_used" (bien)
# 3) si num_hidden_layers == 2 -> uniform_nodes no aplica "not_used" (bien)
# 4) si num_hidden_layers == 3 -> no hay error (con uniform_nodes = False)
# 5) si num_hidden_layers == 4 -> si hidden1_nodes es 64 y node_shrink_factor es 8, 128, 256, 512, 1024, 2048 error! (con uniform_nodes = False)
# return model parameters (esto mejorar luego)


 