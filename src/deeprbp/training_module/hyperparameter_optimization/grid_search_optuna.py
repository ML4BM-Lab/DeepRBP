# /src/deeprbp/training_module/hyperparameter_optimization/grid_search_optuna.py

import os
import logging
import sys
import argparse
import pandas as pd
import optuna
from optuna.samplers import TPESampler
import optuna.visualization.matplotlib as optuna_plt
import optuna.logging
import time

from ..util.logger import Logger
from torch.utils.data import DataLoader
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DataSplitter, Scaler
from .train_model import TrainPredictor
from .model import PredictorModel
from .plots import plot_loss_curve
from .evaluation import calculate_metrics, calculate_metrics_per_category

from ..util.utils import CustomTensorDataset, filter_data_by_sample_ids, save_data, adjust_batch_size

# Define the default configuration path and output directory
default_config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml'
default_output_dir = '/scratch/jsanchoz/DeepRBP/stuff/'

def main(
        config_path_file: str,
        n_trials: int,
        val_batch_size: int,
        output_dir: str
    ):
    
    print(f"[*] Loading configuration from: {config_path_file}...")
    config_parser = ConfigParser(config_path_file)
    config = config_parser.load_config()
    print("[*] Configuration loaded successfully.")

    # Load the data using the data importer and getBM that relates transcript_id with gene_id info
    print("[*] Loading data...")
    data_importer = DataImporter(config['data_paths'])
    data = data_importer.load()
    getBM = pd.read_csv(config['getBM_path'])
    print("[*] Data loaded successfully.")
    
    # Filter a portion of the samples to optimize time and computational resources
    print("[*] Filtering samples to optimize time and resources...")
    subset_idx, _ = DataSplitter.split_data_class(
        data=data, 
        config=config, 
        sample_category=config['sample_category'], 
        test_size=config['sample_fraction']
    )
    data_subset = filter_data_by_sample_ids(data, subset_idx)
    print("[*] Samples filtered successfully.")

    # Perform a train/validation split with the data subset
    print("[*] Splitting data into training and validation sets...")
    splitter = DataSplitter(data_subset, config)
    train_data, valid_data = splitter.split_data_sets(test_name='validation')
    print("[*] Data splitting completed.")

    # Save the training and validation samples used in the optimization to the specified output directory
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

    for data, save_path, custom_names in data_to_save:
        print(f"[*] Saving data to: {save_path}...")
        save_data(data, save_path, custom_names)
        print(f"[*] Data saved successfully at: {save_path}.")

    # Scale the data
    print("[*] Scaling data...")
    scaler = Scaler()
    train_data['scaled_rbp_df'] = scaler.fit_transform(train_data['rbp_df'])
    valid_data['scaled_rbp_df'] = scaler.transform(valid_data['rbp_df'])
    print("[*] Data scaled successfully.")

    # Create custom data sets
    print("[*] Creating custom datasets...")
    train_dataset, valid_dataset = [
            CustomTensorDataset(
                data,
                getBM,
                rbp_data_key='scaled_rbp_df', 
                gene_data_key='gene_df', 
                transcript_data_key='isoform_df',
                trans_col_name=config['trans_col_name'],
                gene_col_name=config['gene_col_name']
            ) for data in [train_data, valid_data]]
    print("[*] Custom datasets created successfully.")

    # Optimization with Optuna using TPESampler
    print("[*] Starting optimization with Optuna...")
    study_name = "DeepRBPredictor-optimization"
    optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
    study = optuna.create_study(study_name=study_name, 
                                sampler=TPESampler(seed=config["seed"]), 
                                direction='minimize') #pruner=pruner) 
                # Recommended budgets with this sampler (#trials: 100-1000)

    # Execute the optimization
    print(f"[*] Running {n_trials} trials...")
    study.optimize(lambda trial: objective(trial, config, train_dataset, valid_dataset, val_batch_size, output_dir), 
                   n_trials=n_trials) 
    
    # Print the results of the best trial
    print(f'Best trial: {study.best_trial}, with parameters: {study.best_params} '
      f'and objective value: {study.best_value}')

    # Save the results to a CSV file
    df_results = study.trials_dataframe()
    results_file_path = os.path.join(output_dir, f'results_{study_name}.csv')
    df_results.to_csv(results_file_path, index=False)
    print(f"[*] Results saved to: {results_file_path}")


def objective(trial, config, train_dataset, valid_dataset, val_batch_size, output_dir):
    ### Suggest Optuna: Sample hyperparameters for this Trial 
    # Select hyperparameters
    num_hidden_layers = trial.suggest_int('num_hidden_layers', 0, 4)
    hidden1_nodes = trial.suggest_categorical('hidden1_nodes', [64, 128, 256, 512, 1024, 2048, 4096])  
    uniform_nodes = trial.suggest_categorical('uniform_nodes', [True, False])
    node_shrink_factor = trial.suggest_categorical('node_shrink_factor', [2, 4, 8])
    activation_func = trial.suggest_categorical('activation_func', ["relu", "tanh", "sigmoid"])
    learning_rate = trial.suggest_categorical('learning_rate', [0.0001, 0.001, 0.01])
    optimizer_name = trial.suggest_categorical('optimizer_name', ['sgd90', 'asgd', 'adam', 'adagrad', 'adadelta', 'adamW'])
    train_batch_size = trial.suggest_categorical('batch_size', [32, 64, 128, 256, 512, 1024, 2048, 4096]) 
    num_epochs = trial.suggest_categorical('num_epochs', [50, 100, 500, 1000, 2000, 3000])  # Sugerir valores categóricos para epochs

    config["num_hidden_layers"] = num_hidden_layers
    config["hidden1_nodes"] = hidden1_nodes
    config["uniform_nodes"] = uniform_nodes
    config["node_shrink_factor"] = node_shrink_factor
    config["activation_func"] = activation_func
    config["learning_rate"] = learning_rate
    config["optimizer_name"] = optimizer_name
    config["batch_size"] = train_batch_size
    config["num_epochs"] = num_epochs

    # Create loaders
    train_loader, val_loader = [
        DataLoader(
            dataset,
            batch_size=adjust_batch_size(dataset, (train_batch_size if idx == 0 else val_batch_size)),
            shuffle=(idx == 0),  # Solo hacer shuffle en el conjunto de entrenamiento
            drop_last=(idx == 0) # Solo drop_last en el conjunto de entrenamiento
        )
        for idx, dataset in enumerate([train_dataset, valid_dataset])   
    ]

    try:
        model = PredictorModel(
                input_size=next(iter(train_loader))['scaled_rbp_df'].shape[1],
                output_size=next(iter(train_loader))['isoform_df'].shape[1],
                config=config
            )
            
        # Proceed with training the model
        trainer = TrainPredictor(
                model=model,
                config=config,
                input_features=('scaled_rbp_df', 'gene_df'), 
                output_features=('isoform_df',)
        )
        train_history, val_history, _ = trainer.fit(train_loader, 
                                                    val_loader, 
                                                    epochs=config["num_epochs"],
                                                    path_save_results=output_dir)

        # Create a unique plot name
        timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")   
        plot_name = f'trial_{trial.number}_{timestamp}'  
        plot_loss_curve(train_history, val_history, output_dir=output_dir, plot_name=plot_name)

        # Evaluate the actual model
        preds_labels = [trainer.generate_predictions(loader) for loader in [train_loader, val_loader]]
        train_metrics, val_metrics = [calculate_metrics(pred, label) for pred, label, _ in preds_labels]

        ## HERE TENDRIA QUE METER LAS NUEVAS MÉTRICAS DE: 
        ### Here we need to create the calculate of the correlation for each gen using getBM: the ranking 
        # really matters for the transcripts within each gene as opposed to all the transcripts across all genes.
        return val_history[-1]
    
    except ValueError as e:
        # Catch the ValueError raised in the PredictorModel
        print(f"Training failed due to configuration error: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Hyperparameter optimization using Optuna.')
    parser.add_argument('-c', '--config_path_file', type=str, 
                        default=default_config_path,
                        help='Path to the configuration YAML file (default: %(default)s)')   
    parser.add_argument('-n', '--n_trials', type=int, required=True, help='Number of trials for Optuna')
    parser.add_argument('-vbs', '--val_batch_size', type=int, required=True, help='Validation batch size')
    parser.add_argument('-o', '--output_dir', type=str, 
                        default=default_output_dir, 
                        help='Output directory for results (default: %(default)s)')
    args = parser.parse_args()

    main(args.config_path_file,
         args.n_trials,
         args.val_batch_size, 
         args.output_dir)
    

# The authors should provide the full table of results for the hyperparameter optimization runs
# to be able to validate the claim that more complex models (more hidden layers) are necessary.
    # Verify that the suggested trial is elegible
    # 1) si num_hidden_layers == 0 -> hidden1_nodes, uniform_nodes, node_shrink_factor y activation_func no aplican "not_used" (bien)
    # 2) si num_hidden_layers == 1 -> uniform_nodes y node_shrink_factor no aplican "not_used" (bien)
    # 3) si num_hidden_layers == 2 -> uniform_nodes no aplica "not_used" (bien)
    # 4) si num_hidden_layers == 3 -> no hay error (con uniform_nodes = False)
    # 5) si num_hidden_layers == 4 -> si hidden1_nodes es 64 y node_shrink_factor es 8, 128, 256, 512, 1024, 2048 error! (con uniform_nodes = False)
    # return model parameters (esto mejorar luego)
   
