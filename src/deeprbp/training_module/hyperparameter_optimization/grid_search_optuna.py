# /src/deeprbp/training_module/hyperparameter_optimization/grid_search_optuna.py

import os
import psutil
import argparse
import matplotlib.pyplot as plt
import torch

import lightning as L
from lightning.pytorch.callbacks import EarlyStopping
from lightning.pytorch.loggers import CSVLogger

import optuna
from optuna.samplers import TPESampler
import optuna.visualization.matplotlib as optuna_plt
from optuna.integration import PyTorchLightningPruningCallback

from ...data_loading.config_loader import ConfigParser
from ...data_preparation.data_module import DeepRBPDataModule
from ..model import PredictorModel
from ...util.utils import print_section_separator, print_gpu_memory_info, set_random_seed

def objective(trial, config, dm, args): 
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

    # Update pipeline configuration with suggested hyperparameters
    config.update("num_hidden_layers", num_hidden_layers)
    config.update("hidden1_nodes", hidden1_nodes)
    config.update("uniform_nodes", uniform_nodes)
    config.update("node_shrink_factor", node_shrink_factor)
    config.update("activation_func", activation_func)
    config.update("learning_rate", learning_rate)
    config.update("optimizer_name", optimizer_name)
    config.update("num_epochs", num_epochs) 
    config.update("batch_norm_eps", batch_norm_eps) 
    config.update("batch_norm_momentum", batch_norm_momentum) 
    
    # Update training batch size
    dm.train_batch_size = train_batch_size

    # Log the hyperparameters for the current trial
    print(f"\n[*] Trial {trial.number}:")
    print(f"    - Number of Hidden Layers: {config.get('num_hidden_layers')}")
    print(f"    - Hidden Layer 1 Nodes: {config.get('hidden1_nodes')}")
    print(f"    - Use Uniform Nodes: {config.get('uniform_nodes')}")
    print(f"    - Node Shrink Factor: {config.get('node_shrink_factor')}")
    print(f"    - Activation Function: {config.get('activation_func')}")
    print(f"    - Learning Rate: {config.get('learning_rate')}")
    print(f"    - Optimizer: {config.get('optimizer_name')}")
    print(f"    - Training Batch Size: {dm.train_batch_size}") # new
    print(f"    - Number of Epochs: {config.get('num_epochs')}")
    print(f"    - Batch Normalization Epsilon (eps): {config.get('batch_norm_eps')}")
    print(f"    - Batch Normalization Momentum: {config.get('batch_norm_momentum')}")
    print("\n")

    try:
        # Create model
        model = PredictorModel(
            config=config,
            input_size=len(dm.train_dataset.rbp_names), 
            output_size=len(dm.train_dataset.trans_names),
            gene_names=dm.train_dataset.gene_names,
            trans_names=dm.train_dataset.trans_names,
            getBM=dm.getBM)
        
        # Define lightning trainer
        early_stopping_callback = EarlyStopping(
            monitor='val_loss',  # The metric to monitor
            min_delta=args.min_delta,  # Minimum change to qualify as an improvement
            patience=args.patience,  # How many epochs to wait after the last improvement
            verbose=True,  # Print messages when stopping
            mode='min'             
        )

        pruning_callback = PyTorchLightningPruningCallback(trial, monitor='val_loss')
        logger = CSVLogger(save_dir=os.getcwd(), name=f"optuna_logs/trial_{trial.number}") # esto podría ser algo así tb.
        trainer = L.Trainer(
                accelerator="gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu",
                devices= int(os.environ.get('SLURM_NTASKS')),  
                num_nodes= int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),   
                logger=logger, 
                callbacks=[early_stopping_callback, pruning_callback],        
                max_epochs=config.get('num_epochs'),  
                deterministic=True,                
                accumulate_grad_batches=1,        
                inference_mode=True,              
                use_distributed_sampler=False,      
                detect_anomaly=False,              
                barebones=False,                  
                sync_batchnorm=True,              
        )
        # Training the model
        trainer.fit(model, dm)
        # Final validation loss
        return trainer.callback_metrics['val_loss'].item()

    except Exception as e:
        # Catch any exception raised in the try block
        print(f"Training failed due to an error: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    
    finally:
        # Update trial.params if the corresponding value is 'unused' (creating a new user attibute)
        for key in ['hidden1_nodes', 'uniform_nodes', 'node_shrink_factor', 'activation_func', 'batch_norm_eps', 'batch_norm_momentum']:
            used_config_value = model.hparams[key]
            print(f"Config value for {key}: {used_config_value}")
            print(f"Current trial parameter for {key}: {trial.params.get(key)}")

            if used_config_value == 'unused':
                trial.set_user_attr(f'real_used_{key}', used_config_value)
                print(f"Updated trial parameter for {key} to 'unused'")
            else:
                print(f"Skipping update for '{key}' because it is not marked as 'unused'")
        print("\n")
        torch.cuda.empty_cache()  # If using GPU, clear the cache
        print_section_separator()

def run_optimization(config, dm, args):
    pruner = optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=10)
    study = optuna.create_study(sampler=TPESampler(seed=42), direction='minimize', pruner=pruner) # Recommended budgets with this sampler (#trials: 100-1000)
    print(f"[*] Running {args.n_trials} trials...")
    study.optimize(lambda trial: objective(trial, config, dm, args), n_trials=args.n_trials)  
    print("Best trial:")
    trial = study.best_trial
    print(f"  Value: {trial.value}")
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")
    return study

#def test_best_model(study): # ESTO HABRIA QUE MODIFICAR CON NUESTRO CASO QUE USAMOS UN CONFIG!!
    # Getting the best hyperparameters
    # best_params = study.best_trial.params
    # # Creating the model with the best hyperparameters
    # model = MNISTClassifier(
    #     layer_1_size=best_params['layer_1_size'],
    #     layer_2_size=best_params['layer_2_size'],
    #     learning_rate=best_params['learning_rate'],
    #     dropout_rate=best_params['dropout_rate']
    # )
    # # Creating trainer instance
    # trainer = L.Trainer(max_epochs=10)
    # # Preparing the data
    # train_loader, val_loader, test_loader = prepare_data()
    # # Training the model with the best hyperparameters
    # trainer.fit(model, train_loader, val_loader)
    # # Testing the model with the test data
    # results = trainer.test(model, test_loader)
    # return results

def main():
    args = parse_args()
    # Load configuration and auxiliary file
    print('\n[grid_search_optuna] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path)

    # Load data module and prepare data for training
    print('\n[main_predictor] 🚀 Initializing DataModule...')
    dm = DeepRBPDataModule(config, args.output_dir)
    print('[main_predictor] 🚀 Preparing data for training...')
    dm.prepare_data() # Load or prepare the necessary data
    dm.setup('fit')  
    print('\n[main_predictor] ──────────────────────────────────────')




    study = run_optimization(config, dm, args)
    # Save the logged results
    results_file_path = os.path.join(args.output_dir, f'results_hyp.csv')
    study.trials_dataframe().to_csv(results_file_path, index=False)
    print(f"[*] Results saved to: {results_file_path}\n")

    # Visualize the results
    try:
      optuna.visualization.plot_optimization_history(study)
      optuna.visualization.plot_param_importances(study)
      optuna.visualization.plot_parallel_coordinate(study)
    except ImportError:
        print("Visualization requires plotly. Install with: pip install plotly")

    # Test the best model
    results = test_best_model(study)
    print(f"Test results with best hyperparameters: {results}")


def parse_args():   
    parser = argparse.ArgumentParser(description='Hyperparameter optimization using Optuna.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration YAML file')   
    parser.add_argument('--n_trials', type=int, required=True, help='Number of trials for Optuna')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for results')
    parser.add_argument('--num_workers', help='DataLoader number of workers', type=int, default=0)
    parser.add_argument('--min_delta', help='Minimum change to qualify as an improvement (for early stopping)', type=float, default=0.001)
    parser.add_argument('--patience', help='How many epochs to wait after the last improvement (for early stopping)', type=int, default=30)
    return parser.parse_args()

if __name__ == "__main__":
    if os.environ.get("LOCAL_RANK")=="0":
        print_gpu_memory_info()
    set_random_seed()
    main()

########################################################################

# Define the default configuration path and output directory
default_config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml'
default_output_dir = '/scratch/jsanchoz/DeepRBP/output/results/hyperpameter_optimization'

# en la nueva versión vamos a haber tenido ya los datos separados y escalados.
# entonces cargaremos los datos del train y val escalados y se generará el loader dentro del objective.
