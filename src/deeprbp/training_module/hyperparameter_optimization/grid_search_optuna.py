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
from ...util.utils import log_section_separator, print_gpu_memory_info, set_random_seed

def main():
    args = parse_args()

    # Load configuration and auxiliary file
    print('\n[grid_search_optuna] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path)

    # Load data module and prepare data for training
    print('\n[grid_search_optuna] 🚀 Initializing DataModule...')
    dm = DeepRBPDataModule(config, args.output_dir)
    print('[grid_search_optuna] 🚀 Preparing data for training...')
    dm.prepare_data() # Load or prepare the necessary data
    dm.setup('fit')  
    print('\n[grid_search_optuna] ──────────────────────────────────────')

    # Run Optuna optimization to find the best hyperparameters
    print('\n[grid_search_optuna] 🔍 Running hyperparameter optimization with Optuna...')
    study = run_optimization(config, dm, args)
    print('[grid_search_optuna] ✅ Hyperparameter optimization completed.')

    # Save the logged results of the optimization
    results_file_path = os.path.join(args.output_dir, f'results_hyp.csv')
    study.trials_dataframe().to_csv(results_file_path, index=False)
    print(f"\n[grid_search_optuna] 💾 Results saved to: {results_file_path}\n")

    # Generate various plots to visualize the optimization history and parameter importances
    print("\n[grid_search_optuna] 📊 Optimization history plot generated.")
    optuna.visualization.plot_optimization_history(study)
    print("\[grid_search_optuna] 📊 Parameter importances plot generated.")
    optuna.visualization.plot_param_importances(study)
    print("[grid_search_optuna] 📊 Parallel coordinates plot generated.")
    optuna.visualization.plot_parallel_coordinate(study)
   
    # Test the best model
    print("\n[grid_search_optuna] Testing the model with the best params.\n")
    dm.setup('test') 
    results = test_best_model(study, config, dm, args)
    print(f"\n [grid_search_optuna] Test results with best hyperparameters: {results}")

def run_optimization(config, dm, args):
    pruner = optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=10)
    study = optuna.create_study(sampler=TPESampler(seed=42), 
                                direction='minimize', 
                                pruner=pruner,
                                study_name='deeprbp_gridsearch_optuna') 
            # Recommended budgets with this sampler (#trials: 100-1000)
    print(f"\n[run_optimization] Running {args.n_trials} trials...")
    study.optimize(lambda trial: objective(trial, config, dm, args), n_trials=args.n_trials)  
    print("\n[run_optimization] Best trial:")
    trial = study.best_trial
    print(f"  Value: {trial.value}")
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")
    return study

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
    config.update("batch_norm_eps", batch_norm_eps) 
    config.update("batch_norm_momentum", batch_norm_momentum) 
    config.update("optimizer_name", optimizer_name)
    config.update("learning_rate", learning_rate)
    config.update("num_epochs", num_epochs) 
    
    # Update training batch size
    dm.train_batch_size = train_batch_size

    # Log the hyperparameters for the current trial
    print(f"\n[objective] Trial {trial.number}:")
    print(f"    - Number of Hidden Layers: {config.get('num_hidden_layers')}")
    print(f"    - Hidden Layer 1 Nodes: {config.get('hidden1_nodes')}")
    print(f"    - Use Uniform Nodes: {config.get('uniform_nodes')}")
    print(f"    - Node Shrink Factor: {config.get('node_shrink_factor')}")
    print(f"    - Activation Function: {config.get('activation_func')}")
    print(f"    - Learning Rate: {config.get('learning_rate')}")
    print(f"    - Optimizer: {config.get('optimizer_name')}")
    print(f"    - Training Batch Size: {dm.train_batch_size}")  
    print(f"    - Number of Epochs: {config.get('num_epochs')}")
    print(f"    - Batch Normalization Epsilon (eps): {config.get('batch_norm_eps')}")
    print(f"    - Batch Normalization Momentum: {config.get('batch_norm_momentum')}")
    print("\n")

    try:
        # Create model
        print('\n[objective] 🚀 Creating the model...')
        model = PredictorModel(
            config=config,
            input_size=len(dm.train_dataset.rbp_names), 
            output_size=len(dm.train_dataset.trans_names),
            gene_names=dm.train_dataset.gene_names,
            trans_names=dm.train_dataset.trans_names,
            getBM=dm.getBM)
        print('[objective] 🚀 Model created:', model)
        print('\n[objective] ──────────────────────────────────────')

        # Define lightning trainer
        print('\n[objective] 🚀 Defining the Lightning trainer...')
        early_stopping_callback = EarlyStopping(
            monitor='val_loss',  # The metric to monitor
            min_delta=args.min_delta,  # Minimum change to qualify as an improvement
            patience=args.patience,  # How many epochs to wait after the last improvement
            verbose=True,  # Print messages when stopping
            mode='min'             
        )
        pruning_callback = PyTorchLightningPruningCallback(trial, monitor='val_loss')
        logger=CSVLogger(f"{args.output_dir}", name=f"optuna_logs/trial_{trial.number}", version=0)
        
        trainer = L.Trainer(
                accelerator="gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu",
                devices= int(os.environ.get('SLURM_NTASKS')),  
                num_nodes= int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),   
                logger=logger, 
                callbacks=[early_stopping_callback, pruning_callback],    
                enable_checkpointing=False,    
                max_epochs=config.get('num_epochs'),  
                deterministic=True,                
                accumulate_grad_batches=1,        
                inference_mode=True,              
                use_distributed_sampler=False,      
                detect_anomaly=False,              
                barebones=False,                  
                sync_batchnorm=True,         
        )
        print('\n[objective] ──────────────────────────────────────')

        # Model training
        print('\n[objective] 🚀 Starting model training...')
        trainer.fit(model, dm)
        print('\n[objective] ──────────────────────────────────────')

        # Final validation loss
        return trainer.callback_metrics['val_loss'].item()

    except Exception as e:
        # Catch any exception raised in the try block
        print(f"\n[objective] Training failed due to an error: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    
    finally:
        # Update trial.params if the corresponding value is 'unused' (creating a new user attibute)
        print(f"\n[objective] Update trial.params if the corresponding value is 'unused'")
        for key in ['hidden1_nodes', 'uniform_nodes', 'node_shrink_factor', 'activation_func', 'batch_norm_eps', 'batch_norm_momentum']:
            used_config_value = model.hparams[key]
            print(f"\n[objective] Config value for {key}: {used_config_value}")
            print(f"[objective] Current trial parameter for {key}: {trial.params.get(key)}")

            if used_config_value == 'unused':
                trial.set_user_attr(f'real_used_{key}', used_config_value)
                print(f"\n[objective] Updated trial parameter for {key} to 'unused'")
            else:
                print(f"\n[objective] Skipping update for '{key}' because it is not marked as 'unused'")
        print('\n[objective] ──────────────────────────────────────')
        torch.cuda.empty_cache()  # If using GPU, clear the cache
        
def test_best_model(study, config, dm, args): # here we are gonna use the test data for prediction but this 
    # model didn't use all samples for training.
    log_section_separator("Testing the model with the best params")
    
    # Getting the best hyperparameters
    best_params = study.best_trial.params

    # Updating config values using the best trial parameters
    config.update("num_hidden_layers", best_params["num_hidden_layers"])
    config.update("hidden1_nodes", best_params["hidden1_nodes"])
    config.update("uniform_nodes", best_params["uniform_nodes"])
    config.update("node_shrink_factor", best_params["node_shrink_factor"])
    config.update("activation_func", best_params["activation_func"])
    config.update("batch_norm_eps", best_params["batch_norm_eps"])
    config.update("batch_norm_momentum", best_params["batch_norm_momentum"])
    config.update("optimizer_name", best_params["optimizer_name"])
    config.update("learning_rate", best_params["learning_rate"])
    config.update("num_epochs", best_params["num_epochs"])

    # Update training batch size
    dm.train_batch_size = best_params["train_batch_size"]  

    # Log the hyperparameters 
    print(f"\n[test_best_model] Test best model:")
    print(f"    - Number of Hidden Layers: {config.get('num_hidden_layers')}")
    print(f"    - Hidden Layer 1 Nodes: {config.get('hidden1_nodes')}")
    print(f"    - Use Uniform Nodes: {config.get('uniform_nodes')}")
    print(f"    - Node Shrink Factor: {config.get('node_shrink_factor')}")
    print(f"    - Activation Function: {config.get('activation_func')}")
    print(f"    - Learning Rate: {config.get('learning_rate')}")
    print(f"    - Optimizer: {config.get('optimizer_name')}")
    print(f"    - Training Batch Size: {dm.train_batch_size}")  
    print(f"    - Number of Epochs: {config.get('num_epochs')}")
    print(f"    - Batch Normalization Epsilon (eps): {config.get('batch_norm_eps')}")
    print(f"    - Batch Normalization Momentum: {config.get('batch_norm_momentum')}")
    print("\n")

    # Creating the model with the best hyperparameters
    print('\n[test_best_model] 🚀 Creating the model...')
    model = PredictorModel(
            config=config,
            input_size=len(dm.train_dataset.rbp_names), 
            output_size=len(dm.train_dataset.trans_names),
            gene_names=dm.train_dataset.gene_names,
            trans_names=dm.train_dataset.trans_names,
            getBM=dm.getBM)
    print('[test_best_model] 🚀 Model created:', model)
    print('\n[test_best_model] ──────────────────────────────────────')

    # Define lightning trainer
    print('\n[test_best_model] 🚀 Defining the Lightning trainer...')
    early_stopping_callback = EarlyStopping(
            monitor='val_loss',  # The metric to monitor
            min_delta=args.min_delta,  # Minimum change to qualify as an improvement
            patience=args.patience,  # How many epochs to wait after the last improvement
            verbose=True,  # Print messages when stopping
            mode='min'             
        )
    
    logger=CSVLogger(f"{args.output_dir}", name=f"optuna_logs/test_model_best_params", version=0)
    trainer = L.Trainer(
                accelerator="gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu",
                devices= int(os.environ.get('SLURM_NTASKS')),  
                num_nodes= int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),   
                logger=logger, 
                callbacks=[early_stopping_callback],        
                max_epochs=config.get('num_epochs'),  
                deterministic=True,                
                accumulate_grad_batches=1,        
                inference_mode=True,              
                use_distributed_sampler=False,      
                detect_anomaly=False,              
                barebones=False,                  
                sync_batchnorm=True,         
        )
    print('\n[test_best_model] ──────────────────────────────────────')
    
    # Model training
    print('\n[test_best_model] 🚀 Starting model training...')
    trainer.fit(model, dm)
    print('\n[test_best_model] ──────────────────────────────────────')
    
    # Testing the model with the test data
    print('\n[test_best_model] 🚀 Starting model testing...')
    results = trainer.test(model, dm)
    print('\n[test_best_model] ──────────────────────────────────────')
    return results

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


# run-deeprbp-predictor \
#   --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_train.yaml" \
#   --output_dir "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor" \
#   --epochs 11 \
#   --num_workers 4 \
#   --min_delta 0.001 \
#   --patience 30 \
#   --save_top_k 1

# run-hyper-optimization-optuna \
#     --config_path '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml' \
#     --n_trials 5 \
#     --output_dir '/scratch/jsanchoz/DeepRBP/output/results/hyperpameter_optimization' \
#     --num_workers 0 \
#     --min_delta 0.001 \
#     --patience 30    



