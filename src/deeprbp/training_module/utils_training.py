# src/deeprbp/training_module/utils_training.py

import pandas as pd
import os
import torch

from .model import PredictorModel

from ..util.utils import print_if_main, find_best_checkpoint
from .evaluation import evaluate_and_visualize_metrics_by_category
from .preds2visualization import plot_all_metrics_history

import lightning as L
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from lightning.pytorch.loggers import CSVLogger

def get_callbacks(args, output_dir):
    """Define and return the callbacks for training.
    Args:
        args: The argument parser or an object containing configuration settings.
    Returns:
        list: A list of callbacks to be used in the Trainer.
    """
    early_stopping = EarlyStopping(
        monitor='validation_loss',  # The metric to monitor
        min_delta=args.min_delta,  # Minimum change to qualify as an improvement
        patience=args.patience,  # How many epochs to wait after the last improvement
        verbose=True,  # Print messages when stopping
        mode='min'             
    )
    ckpt_callback = ModelCheckpoint(
        dirpath=f'{output_dir}/checkpoint_model',  # Custom directory for checkpoints
        filename='deeprbp-predictor-{epoch:02d}-{validation_loss:.2f}',
        monitor='validation_loss',   
        verbose=True,
        save_top_k=args.save_top_k,
        mode='min',  
        enable_version_counter=False,
        save_weights_only=False
    )
    return [early_stopping, ckpt_callback]

def train_and_evaluate_model(dm, args, config, output_dir):
    """
    Trains a PredictorModel using the provided data module and configuration, 
    evaluates performance on train/val/test, and plots training metrics.

    Args:
        dm (DeepRBPDataModule): The data module with train/val/test sets prepared.
        args (Namespace): Parsed command-line arguments.
        config (ConfigParser): Configuration object loaded from file.
        output_dir (str): Path where outputs and logs will be stored.

    Returns:
        model (PredictorModel): The best model loaded from checkpoint.
        trainer (L.Trainer): The PyTorch Lightning trainer used for training and evaluation.
    """
    print_if_main('[utils_training] 🚀 Preparing data for training...')
    dm.setup('fit')
    print_if_main('\n[utils_training] ──────────────────────────────────────')
    
    # Create model
    print_if_main('\n[utils_training] 🚀 Creating the model...')
    model = PredictorModel(
        input_size=len(dm.train_dataset.rbp_names),
        output_size=len(dm.train_dataset.trans_names),
        gene_names=dm.train_dataset.gene_names,
        trans_names=dm.train_dataset.trans_names,
        getBM=dm.getBM,
        verbose=args.verbose
    )
    print_if_main('[utils_training] 🚀 Model created:', model)
    print_if_main('\n[utils_training] ──────────────────────────────────────')
    
    # Define lightning trainer
    print_if_main('\n[utils_training] 🚀 Defining the Lightning trainer...')
    callbacks = get_callbacks(args, output_dir)
    trainer = L.Trainer(
        accelerator="gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu",
        devices=int(os.environ.get('SLURM_NTASKS')),
        num_nodes=int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),
        logger=CSVLogger(f"{output_dir}/csv_logs", name="deep_rbp_predictor", version=0),
        callbacks=callbacks,
        max_epochs=args.epochs,
        deterministic=True,
        accumulate_grad_batches=1,
        inference_mode=True,
        use_distributed_sampler=False,
        detect_anomaly=False,
        barebones=False,
        sync_batchnorm=True
    )
    print_if_main('\n[utils_training] ──────────────────────────────────────')
    
    # Model training
    print_if_main('\n[utils_training] 🚀 Starting model training...')
    trainer.fit(model, dm)
    print_if_main('\n[utils_training] ──────────────────────────────────────')
    
    # Evaluation on test set (only run this on rank 0 (GPU-0 or CPU))
    if trainer.global_rank == 0 or not torch.cuda.is_available():
        print('\n[utils_training] 🚀 Plotting training history...')
        metrics_file_path = os.path.join(output_dir, 'csv_logs', 'deep_rbp_predictor', 'version_0', 'metrics.csv')
        if os.path.exists(metrics_file_path):
            metrics_df = pd.read_csv(metrics_file_path)
            plot_all_metrics_history(metrics_df, f'{output_dir}/metrics_training_history')
        else:
            print(f"[utils_training] ❌ Metrics file does not exist at: {metrics_file_path}")
        print('\n[utils_training] ──────────────────────────────────────')
    
    # Load model checkpoint
    print_if_main('\n[utils_training] 🚀 Loading the best model checkpoint...')
    checkpoint_dir = callbacks[1].dirpath
    model_ckpt_path = find_best_checkpoint(checkpoint_dir)
    print_if_main(f"[utils_training] 🚀 Using checkpoint: {model_ckpt_path}")
    model = PredictorModel.load_from_checkpoint(model_ckpt_path)
    print_if_main('\n[utils_training] ──────────────────────────────────────')
    
    # Evaluation on test set
    print_if_main('\n[utils_training] 🚀 Evaluating on the test set...')
    dm.setup('test')
    trainer.test(model, dm)
    print_if_main('\n[utils_training] ──────────────────────────────────────')
    
    # Evaluation of the model's performance by category
    print_if_main('\n[utils_training] 🚀 Evaluating model performance by category...')
    datasets = [(dm.train_data, 'train'), (dm.val_data, 'val'), (dm.test_data, 'test')]
    for data, set_name in datasets:
        print_if_main(f"[main_predictor] 🚀 Evaluating on the {set_name} set...")
        evaluate_and_visualize_metrics_by_category(
            test_data=data,
            trainer=trainer,
            model=model,
            dm=dm,
            output_dir=output_dir,
            set_name=set_name,
            plot_results=config.get('plot_results')
        )
    print_if_main('\n[utils_training] 🚀 Process completed.')
    return model, trainer