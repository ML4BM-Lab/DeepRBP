# src/deeprbp/training_module/main_predictor.py

import argparse
import pandas as pd
import os
import torch

from ..data_loading.config_loader import ConfigParser
from ..data_preparation.data_module import DeepRBPDataModule
from .model import PredictorModel

from ..util.utils import print_gpu_memory_info, set_random_seed, find_best_checkpoint
from .evaluation import evaluate_and_visualize_metrics_by_category
from .preds2visualization import plot_all_metrics_history

import lightning as L
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from lightning.pytorch.loggers import CSVLogger

def main():
    args = parse_args()
    # Load configuration and auxiliary file
    print('\n[main_predictor] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path) 

    # Load data module and prepare data for training
    print('\n[main_predictor] 🚀 Initializing DataModule...')
    dm = DeepRBPDataModule(config, args.output_dir)
    print('[main_predictor] 🚀 Preparing data for training...')
    dm.prepare_data() # Load or prepare the necessary data
    dm.setup('fit')  
    print('\n[main_predictor] ──────────────────────────────────────')

    # Create model
    print('\n[main_predictor] 🚀 Creating the model...')
    model = PredictorModel(
            config=config,
            input_size=len(dm.train_dataset.rbp_names), 
            output_size=len(dm.train_dataset.trans_names),
            gene_names=dm.train_dataset.gene_names,
            trans_names=dm.train_dataset.trans_names,
            getBM=dm.getBM
    )
    print('[main_predictor] 🚀 Model created:', model)
    print('\n[main_predictor] ──────────────────────────────────────')

    # Define lightning trainer
    print('\n[main_predictor] 🚀 Defining the Lightning trainer...')
    callbacks = get_callbacks(args)
    trainer = L.Trainer(
            accelerator="gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu",
            devices= int(os.environ.get('SLURM_NTASKS')), # Extract GPUs per node int(os.environ.get('SLURM_NTASKS'))
            num_nodes= int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),  # Number of GPU nodes for distributed training. Default: 1. Extract number of nodes int(os.environ.get('SLURM_JOB_NUM_NODES', 1))
            logger=CSVLogger(f"{args.output_dir}/csv_logs", name="deep_rbp_predictor", version=0), 
            callbacks=callbacks,       
            max_epochs=3,    
            deterministic=True,               # Set to True for reproducibility
            accumulate_grad_batches=1,        # Accumulate gradients over multiple batches (default 1) -- In this case, if your global batch size is 20,000 and you set accumulate_grad_batches=4, each GPU will still receive 5,000 samples per mini-batch, but the optimizer will only perform an update after processing 4 mini-batches, effectively simulating a global batch size of 20,000.
            inference_mode=True,              # Whether to run in inference mode (default True) -- Whether to use torch.inference_mode() or torch.no_grad() during evaluation (validate/test/predict).
            use_distributed_sampler=False,     # Use distributed sampler (default True)
            # profiler=PyTorchProfiler(dirpath=opt.result,filename="profiler.txt"),                    # Profiler for performance tracking (default None)
            detect_anomaly=False,             # Detect anomalies in training (default False)
            barebones=False,                  # Use barebones (default False)
            # plugins=SLURMEnvironment(auto_requeue=False), # Custom plugins (default None)
            sync_batchnorm=True              # Synchronize batch normalization (default False) PROBAR CON ESTO EN FALSE
    )
    print('\n[main_predictor] ──────────────────────────────────────')

    # Model training
    print('\n[main_predictor] 🚀 Starting model training...')
    trainer.fit(model, dm)
    print('\n[main_predictor] ──────────────────────────────────────')

    # Plot training history for different metrics
    print('\n[main_predictor] 🚀 Plotting training history...')
    metrics_file_path = os.path.join(args.output_dir, 'csv_logs', 'deep_rbp_predictor', 'version_0', 'metrics.csv')
    if os.path.exists(metrics_file_path):
        metrics_df = pd.read_csv(metrics_file_path)
        plot_all_metrics_history(metrics_df, f'{args.output_dir}/metrics_training_history')
    else:
        print(f"[main_predictor] ❌ Metrics file does not exist at: {metrics_file_path}")
    print('\n[main_predictor] ──────────────────────────────────────')

    # Load model checkpoint
    print('\n[main_predictor] 🚀 Loading the best model checkpoint...')
    checkpoint_dir = f'{args.output_dir}/checkpoint_model'
    model_ckpt_path = find_best_checkpoint(checkpoint_dir)
    print(f"[main_predictor] 🚀 Using checkpoint: {model_ckpt_path}")
    model = PredictorModel.load_from_checkpoint(model_ckpt_path)
    print('\n[main_predictor] ──────────────────────────────────────')

    # Evaluation on test set (general!, se puede hacer lo mismo con el train y val?? o hay que usar para ello el .predict?) (AQUI IGUAL HAY QUE METER TB LOS GLOBALES PARA EL TRAIN Y VAL)
    print('\n[main_predictor] 🚀 Evaluating on the test set...')
    dm.setup('test') 
    trainer.test(model, dm)
    print('\n[main_predictor] ──────────────────────────────────────')

    # Evaluation of the model's performance by category
    print('\n[main_predictor] 🚀 Evaluating model performance by category...')
    datasets = [(dm.train_data, 'train'), (dm.val_data, 'val'), (dm.test_data, 'test')]

    for data, set_name in datasets:
        print(f"[main_predictor] 🚀 Evaluating on the {set_name} set...")
        evaluate_and_visualize_metrics_by_category(
            test_data=data,
            trainer=trainer,
            model=model,
            dm=dm,
            output_dir=args.output_dir,
            set_name=set_name,
            plot_results=config.get('plot_results')
        )
    print('\n[main_predictor] 🚀 Process completed.')

def get_callbacks(args):
    """Define and return the callbacks for training.
    Args:
        args: The argument parser or an object containing configuration settings.
    Returns:
        list: A list of callbacks to be used in the Trainer.
    """
    early_stopping = EarlyStopping(
        monitor='val_loss',  # The metric to monitor
        min_delta=args.min_delta,  # Minimum change to qualify as an improvement
        patience=args.patience,  # How many epochs to wait after the last improvement
        verbose=True,  # Print messages when stopping
        mode='min'             
    )
    ckpt_callback = ModelCheckpoint(
        dirpath=f'{args.output_dir}/checkpoint_model',  # Custom directory for checkpoints
        filename='deeprbp-predictor-{epoch:02d}-{val_loss:.2f}',
        monitor='val_loss',   
        verbose=True,
        save_top_k=args.save_top_k,
        mode='min',  
        enable_version_counter=False,
        save_weights_only=False
    )
    return [early_stopping, ckpt_callback]

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
    return parser.parse_args()

if __name__ == "__main__":
    if os.environ.get("LOCAL_RANK")=="0":
        print_gpu_memory_info()
    set_random_seed()
    main()