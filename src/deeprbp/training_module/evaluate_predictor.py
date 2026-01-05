# src/deeprbp/training_module/evaluate_predictor.py

import argparse
import os
import torch
import lightning as L

from ..data_loading.config_loader import ConfigParser
from ..data_preparation.data_module import DeepRBPDataModule
from .model import PredictorModel
from .evaluation import evaluate_and_visualize, evaluate_and_visualize_metrics_by_category
from ..util.utils import (
    print_if_main,
    setup_output_directory,
    print_gpu_memory_info,
    set_random_seed,
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate a trained DeepRBP predictor on a dataset (no training)."
    )
    parser.add_argument("--config_path", type=str, required=True, help="Path to the YAML config file.")
    parser.add_argument("--model_checkpoint", type=str, required=True, help="Path to a .ckpt checkpoint.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save evaluation outputs.")
    parser.add_argument("--num_workers", type=int, default=0, help="DataLoader workers (if DataModule uses it).")
    parser.add_argument("--verbose", type=int, default=1, help="Verbosity level.")
    return parser.parse_args()

def main():
    args = parse_args()
    set_random_seed()
    torch.set_float32_matmul_precision("high")

    print_if_main("\n[evaluate_predictor] 📄 Loading configuration...")
    config = ConfigParser(args.config_path)

    output_dir = setup_output_directory(args.output_dir)
    print_if_main("[evaluate_predictor] Output directory:", output_dir)

    print_if_main("\n[evaluate_predictor] 🚀 Initializing DataModule...")
    dm = DeepRBPDataModule(config, output_dir)

    print_if_main("\n[evaluate_predictor] 🚀 Loading trained model checkpoint...")
    model = PredictorModel.load_from_checkpoint(args.model_checkpoint)

    print_if_main("\n[evaluate_predictor] 🚀 Building Trainer...")
    trainer = L.Trainer(
        accelerator="gpu" if config.get("cuda") and torch.cuda.is_available() else "cpu",
        devices=1,
        logger=False,
        enable_checkpointing=False,
        deterministic=True,
    )

    print_if_main("\n[evaluate_predictor] 🧪 Setting up TEST data...")
    dm.setup("test")

    print_if_main("\n[evaluate_predictor] 🧪 Running trainer.test(...) ...")
    trainer.test(model, dm)

    print_if_main("\n[evaluate_predictor] 📊 Computing detailed metrics/plots (test split)...")
    if dm.sample_category is None:
        evaluate_and_visualize(
            trainer=trainer,
            model=model,
            dm=dm,
            output_dir=output_dir,
            set_name="test",
        )
    else:
        evaluate_and_visualize_metrics_by_category(
            test_data=dm.test_data,
            trainer=trainer,
            model=model,
            dm=dm,
            output_dir=output_dir,
            set_name="test",
            plot_results=config.get("plot_results"),
        )
    print_if_main("\n[evaluate_predictor] ✅ Done.")

if __name__ == "__main__":
    if os.environ.get("LOCAL_RANK") == "0":
        print_gpu_memory_info()
    main()
