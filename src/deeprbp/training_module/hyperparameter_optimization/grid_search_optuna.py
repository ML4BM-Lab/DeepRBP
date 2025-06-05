# /src/deeprbp/training_module/hyperparameter_optimization/grid_search_optuna.py

import pickle
import os
import argparse
import matplotlib.pyplot as plt
import torch
import time

import lightning as L
from lightning.pytorch.callbacks import EarlyStopping
from lightning.pytorch.loggers import CSVLogger

import optuna
from optuna.samplers import TPESampler
import optuna.visualization.matplotlib as optuna_plt
from optuna.integration import PyTorchLightningPruningCallback
from optuna.storages import RDBStorage

from ...data_loading.config_loader import ConfigParser
from ...data_preparation.data_module import DeepRBPDataModule
from ..model import PredictorModel
from ...util.utils import print_if_main, setup_output_directory, log_section_separator, print_gpu_memory_info, set_random_seed

def parse_args():   
    parser = argparse.ArgumentParser(description='Hyperparameter optimization using Optuna.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration YAML file')   
    parser.add_argument('--n_trials', type=int, required=True, help='Number of trials for Optuna')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for results')
    parser.add_argument('--num_workers', help='DataLoader number of workers', type=int, default=0)
    parser.add_argument('--min_delta', help='Minimum change to qualify as an improvement (for early stopping)', type=float, default=0.001)
    parser.add_argument('--patience', help='How many epochs to wait after the last improvement (for early stopping)', type=int, default=30)
    parser.add_argument('--storage_path', type=str, required=True)
    parser.add_argument('--gpu_id', type=int, default=0, help='ID de la GPU usada por este proceso')
    return parser.parse_args()

def main():
    args = parse_args()
    set_random_seed()

    print_if_main(f"[grid_search_optuna] 🔧 Usando storage en: {args.storage_path}")
    config = ConfigParser(args.config_path)
    output_dir = setup_output_directory(args.output_dir)
    print_if_main('\n[grid_search_optuna] Output directory for main process: ', output_dir)

    # Load data module and prepare data for training
    print_if_main('\n[grid_search_optuna] 🚀 Initializing DataModule...')
    dm = DeepRBPDataModule(config, output_dir)
    print_if_main('[grid_search_optuna] 🚀 Preparing data for training...')
    dm.prepare_data() # Load or prepare the necessary data
    dm.setup('fit')  
    print_if_main('\n[grid_search_optuna] ──────────────────────────────────────')

    storage = RDBStorage(
        url=f"sqlite:///{args.storage_path}",
        engine_kwargs={"connect_args": {"timeout": 10}},
    )

    study = optuna.load_study(
        study_name="deeprbp_gridsearch_optuna",
        storage=storage,
    )

    print_if_main("[grid_search_optuna] 🎯 Executing optimization loop...")
    study.optimize(lambda trial: objective(trial, config, dm, args), n_trials=args.n_trials)
    print_if_main("[grid_search_optuna] ✅ Optimization completed.")

  
def objective(trial, config, dm, args): 
    # Suggest Optuna: Sample hyperparameters for this Trial. Solo el proceso principal sugiere hiperparámetros
    trial_params = {
        'num_hidden_layers': trial.suggest_int('num_hidden_layers', 0, 4),
        'hidden1_nodes': trial.suggest_categorical('hidden1_nodes', [64, 128, 256, 512, 1024, 2048, 4096]),  
        'uniform_nodes': trial.suggest_categorical('uniform_nodes', [True, False]),
        'node_shrink_factor': trial.suggest_categorical('node_shrink_factor', [2, 4, 8]),
        'activation_func': trial.suggest_categorical('activation_func', ["relu", "tanh", "sigmoid"]),
        'learning_rate': trial.suggest_categorical('learning_rate', [0.0001, 0.001, 0.01]),
        'optimizer_name': trial.suggest_categorical('optimizer_name', ['sgd90', 'asgd', 'adam', 'adagrad', 'adadelta', 'adamW']),
        'batch_size': trial.suggest_categorical('batch_size', [32, 64, 128, 256, 512, 1024, 2048]), 
        'num_epochs': trial.suggest_categorical('num_epochs', [50, 100, 500, 1000, 2000, 3000]),   
        'batch_norm_eps': trial.suggest_categorical('batch_norm_eps', [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]), 
        'batch_norm_momentum': trial.suggest_categorical('batch_norm_momentum', [0.1, 0.5, 0.9]) 
    }
    
    # Update pipeline configuration with suggested hyperparameters
    config.update("num_hidden_layers", trial_params['num_hidden_layers'])
    config.update("hidden1_nodes", trial_params['hidden1_nodes'])
    config.update("uniform_nodes", trial_params['uniform_nodes'])
    config.update("node_shrink_factor", trial_params['node_shrink_factor'])
    config.update("activation_func", trial_params['activation_func'])
    config.update("batch_norm_eps", trial_params['batch_norm_eps']) 
    config.update("batch_norm_momentum", trial_params['batch_norm_momentum']) 
    config.update("optimizer_name", trial_params['optimizer_name'])
    config.update("learning_rate", trial_params['learning_rate'])
    config.update("num_epochs", trial_params['num_epochs']) 
    dm.train_batch_size = trial_params['batch_size']

    # Log the hyperparameters for the current trial
    print_if_main(f"\n[objective] Trial {trial.number}:")
    for k, v in trial_params.items():
        print_if_main(f"    - {k}: {v}")
    print_if_main("\n")

    try:
        # Create model
        print_if_main('\n[objective] 🚀 Creating the model...')
        model = PredictorModel(
            config=config,
            input_size=len(dm.train_dataset.rbp_names), 
            output_size=len(dm.train_dataset.trans_names),
            gene_names=dm.train_dataset.gene_names,
            trans_names=dm.train_dataset.trans_names,
            getBM=dm.getBM
        )
        print_if_main('[objective] 🚀 Model created:', model)
        print_if_main('\n[objective] ──────────────────────────────────────')

        # Define lightning trainer
        print_if_main('\n[objective] 🚀 Defining the Lightning trainer...')
        logger = CSVLogger(
            save_dir=os.path.join(args.output_dir, f"optuna_logs/gpu_{args.gpu_id}"),
            name=f"trial_{trial.number}"
        )
        
        trainer = L.Trainer( # actually maybe trainer could be created outside and just change the number of epochs.
            accelerator="gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu",
            devices=1, #int(os.environ.get('SLURM_NTASKS')),  
            #num_nodes= int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),   
            logger=logger, 
            callbacks=[PyTorchLightningPruningCallback(trial, monitor='validation_loss')],
            enable_checkpointing=False,    
            max_epochs=config.get('num_epochs'),  
            deterministic=True,                
            accumulate_grad_batches=1,        
            inference_mode=True,              
            use_distributed_sampler=False,      
            detect_anomaly=False,              
            barebones=False,                  
            sync_batchnorm=True      
        )
        print_if_main('\n[objective] ──────────────────────────────────────')

        # Model training
        print_if_main('\n[objective] 🚀 Starting model training...')
        trainer.fit(model, dm)
        print_if_main('\n[objective] ──────────────────────────────────────')

        # Final validation loss
        return trainer.callback_metrics['validation_loss'].item()

    except Exception as e:
        # Catch any exception raised in the try block
        print_if_main(f"\n[objective] Training not arrived to the end due to: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    
    finally:
        # Check if the trial was pruned
        print_if_main(f"\n[objective] Checking for pruning conditions...")
        if trial.should_prune():
            trial.set_user_attr('pruned', True)
            print_if_main("⏹️ Trial was pruned!")
        else:
            trial.set_user_attr('pruned', False)

        # Update trial.params if 'unused' (else nan)
        print_if_main(f"\n[objective] Update trial.params if the corresponding value is 'unused'")
        
        for key in ['hidden1_nodes', 'uniform_nodes', 'node_shrink_factor', 
                    'activation_func', 'batch_norm_eps', 'batch_norm_momentum']:
            used_config_value = model.hparams[key]
            print_if_main(f"    - {key} → real config: {used_config_value} | original trial param: {trial.params.get(key)}")

            if used_config_value == 'unused':
                trial.set_user_attr(f'real_used_{key}', used_config_value)
                print_if_main(f"    ↳ Saved 'unused' value for {key}")
            else:
                print_if_main(f"    ↳ No update needed for {key}")
        
        print_if_main('\n[objective] ✅ Trial finalization completed.')
        torch.cuda.empty_cache()  # If using GPU, clear the cache
        
if __name__ == "__main__":
    print_gpu_memory_info()
    set_random_seed()
    main()
  



