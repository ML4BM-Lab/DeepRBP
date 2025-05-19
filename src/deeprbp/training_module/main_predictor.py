# pipeline_new

import argparse
import pandas as pd
import os
import torch

from ..data_loading.config_loader import ConfigParser
from ..data_preparation.prepare_data import PrepareData
from .model_new import PredictorModel
from ..util.utils import print_gpu_memory_info, set_random_seed, find_latest_checkpoint
from .evaluation import evaluate_and_visualize_metrics_by_category
from .preds2visualization import plot_all_metrics_history

import lightning as L
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from lightning.pytorch.loggers import CSVLogger

### new pipeline ##########################################################################
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_train.yaml'
output_dir = '/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor'

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP predictor training pipeline.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results.')
    parser.add_argument('--epochs', help='Training epochs for training', type=int, default=10)
    parser.add_argument('--num_workers', help='DataLoader number of workers', type=int, default=0)
    parser.add_argument('--min_delta', help='Minimum change to qualify as an improvement (for early stopping)', type=float, default=0.001)
    parser.add_argument('--patience', help='How many epochs to wait after the last improvement (for early stopping)', type=int, default=30)
    return parser.parse_args()
    
# def main(): THIS WAS OLD VERSION
args = parse_args()
set_random_seed()

if os.environ.get("LOCAL_RANK")=="0":
    # Make sure all GPUs are empty
    print_gpu_memory_info()

# Load configuration
config = ConfigParser(args.config_path)
print(config)
getBM = pd.read_csv(config.get('getBM_path'))

# Load and process data
prep_data = PrepareData(config, getBM, args.output_dir)
data = prep_data.load_data()
train_data, val_data = prep_data.split_data(data)
#self.save_split_data(train_data, valid_data)
prep_data.fit_scaler(train_data) 
scaled_train_data, scaled_val_data = prep_data.scale_train_val_data(train_data, val_data)
train_dataset, valid_dataset = prep_data.create_train_val_datasets(scaled_train_data, scaled_val_data)
train_loader, valid_loader = prep_data.create_train_val_loaders(train_dataset, valid_dataset)

# config_path y getBM lo cogemos del config y listo. 
model = PredictorModel(config=config,
                       input_size=next(iter(train_loader))['scaled_rbp_df'].shape[1], 
                       output_size=next(iter(train_loader))['isoform_df'].shape[1],
                       gene_names=train_dataset.gene_names,
                       trans_names=train_dataset.trans_names,
                       getBM=getBM)
print(model)
print(model.hparams)

# Define the EarlyStopping callback
early_stopping = EarlyStopping(
    monitor='val_loss', # The metric to monitor
    min_delta=min_delta, #args.min_delta,  # Minimum change to qualify as an improvement
    patience=patience, #args.patience, # How many epochs to wait after the last improvement
    verbose=True, # Print messages when stopping
    mode='min'             
)

# Create a checkpoint callback
# Save the model periodically by monitoring a quantity. Every metric logged with log() or log_dict()
# saves a file like: my/path/epoch=0-step=10.ckpt
ckpt_callback = ModelCheckpoint(
    dirpath=f'{output_dir}/checkpoint_model', # Custom directory for checkpoints
    filename='deeprbp-predictor-{epoch:02d}-{val_loss:.2f}',
    monitor='val_loss',   
    verbose=True,
    save_top_k=33, #1,
    mode='min',  
    enable_version_counter=False,
    save_weights_only=False
)

csv_logger = CSVLogger(f"{output_dir}/csv_logs", name="deep_rbp_predictor", version=0)  

# Train the model
accelerator = "gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu"
trainer = L.Trainer(
        accelerator=accelerator,
        devices= int(os.environ.get('SLURM_NTASKS')), # Extract GPUs per node int(os.environ.get('SLURM_NTASKS'))
        num_nodes= int(os.environ.get('SLURM_JOB_NUM_NODES', 1)),  # Number of GPU nodes for distributed training. Default: 1. Extract number of nodes int(os.environ.get('SLURM_JOB_NUM_NODES', 1))
        logger=csv_logger, # Logger for tracking experiments
        callbacks=[early_stopping, ckpt_callback],       # Add the early stopping callback
        max_epochs=10,  # Number of epochs  
        deterministic=True,               # Set to True for reproducibility
        accumulate_grad_batches=1,        # Accumulate gradients over multiple batches (default 1) -- In this case, if your global batch size is 20,000 and you set accumulate_grad_batches=4, each GPU will still receive 5,000 samples per mini-batch, but the optimizer will only perform an update after processing 4 mini-batches, effectively simulating a global batch size of 20,000.
        inference_mode=True,              # Whether to run in inference mode (default True) -- Whether to use torch.inference_mode() or torch.no_grad() during evaluation (validate/test/predict).
        use_distributed_sampler=False,     # Use distributed sampler (default True)
        # profiler=PyTorchProfiler(dirpath=opt.result,filename="profiler.txt"),                    # Profiler for performance tracking (default None)
        detect_anomaly=False,             # Detect anomalies in training (default False)
        barebones=False,                  # Use barebones (default False)
        # plugins=SLURMEnvironment(auto_requeue=False), # Custom plugins (default None)
        sync_batchnorm=True,             # Synchronize batch normalization (default False) PROBAR CON ESTO EN FALSE
    )
                
trainer.fit(model=model, train_dataloaders=train_loader, val_dataloaders=valid_loader)

# Load model checkpoint
checkpoint_dir = f'{output_dir}/checkpoint_model'
model_ckpt_path = find_latest_checkpoint(checkpoint_dir)
#model_ckpt_path = find_best_checkpoint(checkpoint_dir)
model.ckpt_path = model_ckpt_path
print(f"Using checkpoint: {model.ckpt_path}")

### AQUI JOSEBA! SE HA DEMOSTRADO QUE LA FORMA EN LA QUE KAT "CARGA EL MEJOR MODELO" NO FUNCIONA, BUSCA OTRA TIO EN LIGHTNING QUE FUNCIONE




# Plot history on Training
metrics_file_path = os.path.join(output_dir, 'csv_logs', 'deep_rbp_predictor', 'version_0', 'metrics.csv')

# Check if the metrics.csv file exists
if os.path.exists(metrics_file_path):
    # Load the metrics.csv file into a Pandas DataFrame
    metrics_df = pd.read_csv(metrics_file_path)
else:
    print(f"Metrics file does not exist at: {metrics_file_path}")
plot_all_metrics_history(metrics_df, f'{output_dir}/metrics_training_history')

# Evaluate on each tumor type
evaluate_and_visualize_metrics_by_category(
    test_data = scaled_train_data,
    sample_category = config.get('sample_category'),
    trainer = trainer,
    model = model,
    batch_size = config.get('val_batch_size'),
    output_dir = output_dir,
    set_name = 'train',
    getBM = getBM,
    plot_results = config.get('plot_results')

)

evaluate_and_visualize_metrics_by_category(
    test_data = scaled_val_data,
    sample_category = config.get('sample_category'),
    trainer = trainer,
    model = model,
    batch_size = config.get('val_batch_size'),
    output_dir = output_dir,
    set_name = 'val',
    getBM = getBM,
    plot_results = config.get('plot_results')

)

# if __name__ == "__main__":
#     main()