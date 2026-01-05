
# /src/deeprbp/explainability_module/real_knockdowns/utils_kds.py

import os
import torch
import pandas as pd
from lightning.pytorch import Trainer
from lightning.pytorch.loggers import CSVLogger

from ...data_loading.config_loader import ConfigParser
from ...data_preparation.data_module import DeepRBPDataModule

def build_datamodule(config_base, condition, processed_data_dir, base_output_dir, scaler_mode):
    """
    Creates a DeepRBPDataModule for a given experimental condition.

    Parameters:
    - config_base: ConfigParser object with base configuration.
    - condition: str, specific condition or dataset name.
    - processed_data_dir: str, directory containing processed data files.
    - base_output_dir: str, directory to save outputs specific to the condition.

    Returns:
    - DeepRBPDataModule: Configured data module for the specified condition.

    This function ensures the output directory exists and sets the test path for the condition.
    """
    output = os.path.join(base_output_dir, condition)
    os.makedirs(output, exist_ok=True)
    return DeepRBPDataModule(
        config=ConfigParser(
            **config_base.config_data,
            test_path_files=os.path.join(processed_data_dir, condition),
            scaler_mode=scaler_mode
        ),
        output_dir=output,
        verbose=1
    )

def build_trainer(config, output_dir):
    """
    Builds a PyTorch Lightning Trainer for inference/testing with a CSVLogger 
    using the same output_dir as the corresponding DataModule.

    Args:
        config (ConfigParser): Configuration object (e.g., contains 'cuda').
        output_dir (str): The output directory associated with the condition 
                          (e.g., .../GSE136366/control or .../GSE136366/knockdown).

    Returns:
        Trainer: Configured PyTorch Lightning Trainer instance.
    """
    # Create logger inside the same output_dir used by the DataModule
    logger = CSVLogger(
        save_dir=output_dir,
        name="logs",
        version=0
    )
    # Determine accelerator
    accelerator = "gpu" if config.get('cuda') and torch.cuda.is_available() else "cpu"
    # Build and return the Trainer
    return Trainer(  #Set up trainer (inference only)
        accelerator=accelerator,
        devices=int(os.environ.get("SLURM_NTASKS", 1)),
        num_nodes=int(os.environ.get("SLURM_JOB_NUM_NODES", 1)),
        logger=logger,
        max_epochs=1,
        inference_mode=True,
        deterministic=True,
        enable_progress_bar=True,
        enable_model_summary=False,
        log_every_n_steps=10,
        enable_checkpointing=False
    )