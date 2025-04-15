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
from torch.utils.data import DataLoader

from ...data_loading.config_loader import ConfigParser
from ...data_loading.data_loader import DataImporter, DataSplitter, Scaler
from ..train_model import TrainPredictor
from ..model import PredictorModel
from ...util.plots import plot_loss_curve
from ..evaluation import calculate_metrics, calculate_spearman_corr_per_gene
from ...util.utils  import (
    CustomTensorDataset,
    filter_data_by_sample_ids,
    save_data,
    adjust_batch_size,
    summarize_model,
    print_section_separator
)

# Define the default configuration path and output directory
default_config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml'
default_output_dir = '/scratch/jsanchoz/DeepRBP/stuff/'

def parse_args():   
    parser = argparse.ArgumentParser(description='Hyperparameter optimization using Optuna.')
    parser.add_argument('-c', '--config_path_file', type=str, 
                        default=default_config_path,
                        help='Path to the configuration YAML file (default: %(default)s)')   
    parser.add_argument('-n', '--n_trials', type=int, required=True, help='Number of trials for Optuna')
    parser.add_argument('-vbs', '--val_batch_size', type=int, required=True, help='Validation batch size')
    parser.add_argument('-o', '--output_dir', type=str, 
                        default=default_output_dir, 
                        help='Output directory for results (default: %(default)s)')
    return parser.parse_args()

def main():
    args = parse_args()

    print(f"[*] Loading configuration from: {args.config_path_file}...")
    config = ConfigParser(args.config_path_file)
    print("[*] Configuration loaded successfully\n.")

    # Load, process, and scale the data
    train_dataset, valid_dataset, getBM = load_and_process_data(config, args.output_dir)

    # Optimization with Optuna using TPESampler and MedianPruner
    print("[*] Starting optimization with Optuna...")
    study_name = "DeepRBPredictor-optimization"
    optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))

    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5, n_startup_trials=5)    
    study = optuna.create_study(study_name=study_name, 
                                sampler=TPESampler(seed=config.get('seed')), 
                                direction='minimize', 
                                pruner=pruner) # Recommended budgets with this sampler (#trials: 100-1000)

    # Execute the optimization
    print(f"[*] Running {args.n_trials} trials...")
    study.optimize(lambda trial: objective(trial, config, train_dataset, valid_dataset, args.val_batch_size, getBM, args.output_dir), 
                   n_trials=args.n_trials)
    
    # Print the results of the best trial
    print(f'Best trial: {study.best_trial}, with parameters: {study.best_params} '
      f'and objective value: {study.best_value}\n')

    # Save the results to a CSV file
    df_results = study.trials_dataframe()
    print(df_results)
    results_file_path = os.path.join(args.output_dir, f'results_{study_name}.csv')
    df_results.to_csv(results_file_path, index=False)
    print(f"[*] Results saved to: {results_file_path}\n")

def load_and_process_data(config, output_dir):
    """
    Load, process, and scale the data, returning the datasets and getBM.

    Args:
        config (dict): Configuration object dictionary containing paths and parameters.
        output_dir (str): Output directory for saving processed data.

    Returns:
        tuple: train_dataset, valid_dataset, getBM
    """
    # Load the data using the data importer and getBM that relates transcript_id with gene_id info
    print("[*] Loading data...")
    data_importer = DataImporter(config.get('train_data_paths'))
    data = data_importer.load()
    getBM = pd.read_csv(config.get('getBM_path'))
    print("[*] Data loaded successfully\n.")
    # sustituir luego por esto: pipeline = DeepRBPredictorPipeline(config_path, output_dir)
    # sustituir luego por esto: data  = pipeline.import_data()

    # Filter a portion of the samples to optimize time and computational resources
    print("[*] Filtering samples to optimize time and resources...")
    subset_idx, _ = DataSplitter.split_data_class(
        data=data, 
        config=config, 
        sample_category=config.get('sample_category'), 
        test_size=config.get('sample_fraction')
    )
    data_subset = filter_data_by_sample_ids(data, subset_idx)
    print("[*] Samples filtered successfully\n.")
    # sustituir luego por esto: data_subset = pipeline.filter_samples(data)
    
    # Perform a train/validation split with the data subset
    print("[*] Splitting data into training and validation sets...")
    splitter = DataSplitter(data_subset, config)
    train_data, valid_data = splitter.split_data_sets(test_name='validation')
    print("[*] Data splitting completed\n.")
    # sustituir luego por esto: train_data, valid_data = pipeline.split_data(data_subset)
    
    # Save the training and validation samples used in the optimization to the specified output directory
    print("[*] Saving the training and validation data used ...")
    data_to_save = [
        (train_data, os.path.join(output_dir, 'data/Train'), {
            'rbp_df': 'train_RBPs_log2p_tpm.csv',
            'isoform_df': 'train_trans_log2p_tpm.csv',
            'gene_df': 'train_gn_tpm.csv',
            'metadata_df': 'train_phenotype_metadata.csv'
        }),
        (valid_data, os.path.join(output_dir, 'data/Validation'), {
            'rbp_df': 'val_RBPs_log2p_tpm.csv',
            'isoform_df': 'val_trans_log2p_tpm.csv',
            'gene_df': 'val_gn_tpm.csv',
            'metadata_df': 'val_phenotype_metadata.csv'
        })
    ]
    for data, save_path, custom_names in data_to_save: ## UNCOMMENT
         print(f"[*] Saving data to: {save_path}...")
         save_data(data, save_path, custom_names)
         print(f"[*] Data saved successfully at: {save_path}.")
    print("[*] Saved the training and validation data succesfully\n")
    # sustituir luego por esto: pipeline.save_split_data(train_data, valid_data)

    # Scale the data
    print("[*] Scaling data...")
    scaler = Scaler()
    train_data['scaled_rbp_df'] = scaler.fit_transform(train_data['rbp_df'])
    valid_data['scaled_rbp_df'] = scaler.transform(valid_data['rbp_df'])
    print("[*] Data scaled successfully\n.")
    # sustituir luego por esto:  
                        #train_data, valid_data  = pipeline.scale_data(train_data, valid_data) 
    # Create custom data sets
    print("[*] Creating custom datasets...")
    train_dataset, valid_dataset = [
        CustomTensorDataset(
            data,
            getBM,
            rbp_data_key='scaled_rbp_df', 
            gene_data_key='gene_df', 
            transcript_data_key='isoform_df',
            trans_col_name=config.get('trans_col_name'),
            gene_col_name=config.get('gene_col_name')
        ) for data in [train_data, valid_data]
    ]
    # sustituir luego por esto:  
                        #train_dataset, valid_dataset  = pipeline.create_datasets(train_data, valid_data)
    
    print("[*] Custom datasets created successfully\n.")
    return train_dataset, valid_dataset, getBM


def extract_metrics(metrics, train_metrics, val_metrics, metric_mapping=None):
    """
    Fills the `metrics` dictionary with values from `train_metrics` and `val_metrics`
    using a specific mapping for the keys.
    
    Args:
        metrics: Dictionary to fill.
        train_metrics: Dictionary containing training metrics.
        val_metrics: Dictionary containing validation metrics.
        metric_mapping: Dictionary mapping metric names to their keys in `metrics`.
    
    Returns:
        Filled `metrics` dictionary with corresponding values.
    """
    # Default metric mapping if none is provided
    if metric_mapping is None:
        metric_mapping = {
            'train_spearman_corr_general': ('spearman_corr', 'spearman_corr'),
            'train_pearson_corr': ('pearson_corr', 'pearson_corr'),
            'train_mse': ('mse', 'mse'),
            'train_r2': ('r2', 'r2'),
            'train_spearman_corr_per_gene': ('train_spearman_corr_per_gene', 'val_spearman_corr_per_gene'),
            'train_spearman_corr_per_gene_max_trans': ('train_spearman_corr_per_gene_max_trans', 'val_spearman_corr_per_gene_max_trans')
        }
    for key, (train_key, val_key) in metric_mapping.items():
        metrics[key] = train_metrics.get(train_key)
        metrics[key.replace('train_', 'val_')] = val_metrics.get(val_key)
    return metrics

def objective(trial, config, train_dataset, valid_dataset, val_batch_size, getBM, output_dir): # pipeline, y quitar: getBM, val_batch_size y config porque ya están en pipeline
    process = psutil.Process()
    mem_before = process.memory_info().rss / (1024 ** 3)  # Memory before in GB

    # Define the additional metrics to calculate
    metrics = {
        'train_spearman_corr_general': None,
        'train_pearson_corr': None,
        'train_mse': None,
        'train_r2': None,
        'val_spearman_corr_general': None,
        'val_pearson_corr': None,
        'val_mse': None,
        'val_r2': None,
        'train_spearman_corr_per_gene': None,
        'train_spearman_corr_per_gene_max_trans': None,
        'val_spearman_corr_per_gene': None,
        'val_spearman_corr_per_gene_max_trans': None,
    } 
    
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

    config.update("num_hidden_layers", num_hidden_layers)
    config.update("hidden1_nodes", hidden1_nodes)
    config.update("uniform_nodes", uniform_nodes)
    config.update("node_shrink_factor", node_shrink_factor)
    config.update("activation_func", activation_func)
    config.update("learning_rate", learning_rate)
    config.update("optimizer_name", optimizer_name)
    config.update("train_batch_size", train_batch_size)
    config.update("num_epochs", num_epochs) # JOSEBA MEJORA A TRATAR: Tal vez habria que hacer un self.config.update en model.py (en def _update_unused_variables(self):) 
                                            # después de crear el model con las variables que NO se están usando en el trial actual !!!
                                            # DESPUES DE ESO HAY QUE ACCEDER AL SISTEMA DE LOS TRIALS Y MODIFICAR LOS QUE NO SE ESTÉN USANDO!!

    # Log the hyperparameters for the current trial
    print(f"[*] Trial {trial.number}:")
    print(f"    - Number of Hidden Layers: {config.get('num_hidden_layers')}")
    print(f"    - Hidden Layer 1 Nodes: {config.get('hidden1_nodes')}")
    print(f"    - Use Uniform Nodes: {config.get('uniform_nodes')}")
    print(f"    - Node Shrink Factor: {config.get('node_shrink_factor')}")
    print(f"    - Activation Function: {config.get('activation_func')}")
    print(f"    - Learning Rate: {config.get('learning_rate')}")
    print(f"    - Optimizer: {config.get('optimizer_name')}")
    print(f"    - Training Batch Size: {config.get('batch_size')}")
    print(f"    - Number of Epochs: {config.get('num_epochs')}")
    print("\n")
    
    try:
        # Create DataLoader instances
        train_loader, val_loader = [
        DataLoader(
            dataset,
            batch_size=adjust_batch_size(dataset, (train_batch_size if idx == 0 else val_batch_size)),
            shuffle=(idx == 0),  # Solo hacer shuffle en el conjunto de entrenamiento
            drop_last=(idx == 0) # Solo drop_last en el conjunto de entrenamiento
        )
        for idx, dataset in enumerate([train_dataset, valid_dataset])   
        ]

        # sustituir luego por esto:  
                        # train_loader, val_loader = pipeline.get_loaders(train_data, valid_data)
    
        # Create model instance
        model = PredictorModel(
                input_size=next(iter(train_loader))['scaled_rbp_df'].shape[1],
                output_size=next(iter(train_loader))['isoform_df'].shape[1],
                config=config
            )
        
        # Print the model summary
        print("\n[*] Model Summary:")
        summarize_model(model, train_loader, config.get("batch_size"))
       
        # Proceed with training the model
        trainer = TrainPredictor(
                model=model,
                config=config,
                input_features=('scaled_rbp_df', 'gene_df'), 
                output_features=('isoform_df',)
        )

        train_history, val_history, _ = trainer.fit(train_loader, 
                                                    val_loader, 
                                                    epochs=config.get("num_epochs"),
                                                    path_save_results=output_dir,
                                                    optuna_trial=trial)

        # sustituir luego por esto:  
                        #train_history, val_history, trainer = pipeline.train_model(train_loader, val_loader, trial)
    
        # Create a unique plot name
        timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")   
        plot_name = f'trial_{trial.number}_{timestamp}'  
        plot_loss_curve(train_history, val_history, output_dir=output_dir, plot_name=plot_name)

        # sustituir luego por esto:  
                        #pipeline.save_model_and_history(trainer, train_history, val_history)

        # Evaluate the actual model
        preds_labels = [trainer.generate_predictions(loader) for loader in [train_loader, val_loader]]
        
        # General metrics
        train_metrics_general, val_metrics_general = [calculate_metrics(pred.flatten(), label.flatten()) for pred, label in preds_labels]  
        
        # Correlation per gene
        train_corr_per_gene, val_corr_per_gene = [calculate_spearman_corr_per_gene(train_dataset, pred, label, getBM) for pred, label in preds_labels]
        
        # luego sustituir por esto: results = self.evaluate_model(trainer, train_loader, valid_loader)

        # Add new entries from train_corr_per_gene to train_metrics_general
        train_metrics_general['train_spearman_corr_per_gene'] = train_corr_per_gene['mean_corr']
        train_metrics_general['train_spearman_corr_per_gene_max_trans'] = train_corr_per_gene['mean_corr_max_trans']
        val_metrics_general['val_spearman_corr_per_gene'] = val_corr_per_gene['mean_corr']
        val_metrics_general['val_spearman_corr_per_gene_max_trans'] = val_corr_per_gene['mean_corr_max_trans']

        # Update the metrics dictionary with training and validation metrics
        metrics = extract_metrics(metrics, train_metrics_general, val_metrics_general) # actualizar para meter la correlacion por gen en el trans y gen
        
        mem_after = process.memory_info().rss / (1024 ** 3)  # memory after in GB
        print(f"Memory used during this trial: {mem_after - mem_before:.4f} GB")

        print("\n")
        return val_history[-1]
    
    except Exception as e:
        # Catch any exception raised in the try block
        print(f"Training failed due to an error: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    
    finally:
        # Log metrics to the trial
        for key in metrics:
            trial.set_user_attr(key, metrics[key])

        # Cleanup: Free resources
        del model
        del trainer
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


 