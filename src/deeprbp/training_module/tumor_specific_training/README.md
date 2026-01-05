
# Tumor-specific versus general DeepRBP training
This module evaluates whether training a single DeepRBP predictor on all tumor types jointly 
provides better generalization than training separate, tumor-specific models.

---

## Training strategies
We compare two approaches:

- **General model**  
  A single model trained on all available tumor types using the architecture
  selected via Optuna, with a consistent training–validation split.

- **Tumor-specific models**  
  Independent models trained separately for each tumor type using the same
  optimized architecture and comparable data splits.

Batch sizes for tumor-specific models are adjusted based on the number of
available training samples in each tumor type to ensure stable optimization
and comparable convergence behavior:

- Tumor types with fewer than 100 training samples use a small batch size
  (e.g., 8).
- Tumor types with 100–300 training samples use a medium batch size
  (e.g., 32).
- Tumor types with more than 300 training samples use a larger batch size
  (e.g., 64).


---

## Evaluation
For each tumor type, we evaluate transcript abundance prediction using:

- the general model,
- the corresponding tumor-specific model.

Performance is assessed using multiple metrics, including:

- Spearman correlation
- Pearson correlation
- R²
- Mean Squared Error (MSE)

---

## Execution
The analysis is executed via a dedicated SLURM script:

```bash
sbatch slurm/training_module/tumor_specific_training/run_tcga_specific_vs_general_training.sh
```

This script:

- trains tumor-specific models with early stopping,
- loads results from the general model,
- compares both strategies across all tumor types,
- generates publication-ready plots and summary tables.

The goal is to assess whether the general model provides sufficient performance across contexts, potentially offering better generalization than models trained on individual tumor types.

### SLURM script configuration and required arguments
The tumor-specific versus general training analysis is executed via a dedicated SLURM script that 
launches a distributed PyTorch job using `torchrun`.

Before submitting the job, users must review and adapt **both the SLURM resource directives and the runtime arguments**, 
as paths and resources are cluster- and user-specific.

---
 
#### Runtime arguments passed to the training script
The SLURM script launches the following Python module:

```bash
deeprbp.training_module.tumor_specific_training.main_specific_vs_general
```

The most relevant runtime arguments are:

- `--output_base_dir`
Base directory where results for each tumor-specific model will be saved.
One subdirectory per tumor type is created automatically.

- `--all_output_dir`
Path to the output directory of the general model trained **on all tumor types together**. 
This directory must contain the
test_tumor_category_results.csv file, which is used as the reference
baseline in cross-tumor comparisons.

- `--epochs`
Maximum number of training epochs for each tumor-specific model.

- `--num_workers`
Number of DataLoader workers per process.

- `--min_delta` and `--patience`
Early stopping parameters controlling convergence behavior.

- `--save_top_k`
Number of best-performing checkpoints (based on validation loss) retained
per tumor-specific model.

**Dataset paths and batch size logic**
Dataset paths (training and test splits) and feature mappings are currently
defined inside the Python entry point and must match the user’s local
directory structure.

Batch sizes for tumor-specific models are **automatically adjusted**
according to the number of samples available for each tumor type, based on
predefined heuristics (see `BATCH_SIZE_BY_TUMOR`).

**Users do not** need to manually set batch sizes when running the script.

**Notes on distributed execution**
- The script uses **PyTorch Distributed Data Parallel (DDP)** via `torchrun`.
- All tumor-specific models are trained sequentially, while individual
training runs may internally leverage multiple GPUs.
- Result aggregation and plotting are performed only by the main process
(global rank 0).