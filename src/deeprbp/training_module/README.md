
# Training DeepRBP from scratch (optional)
This module describes how to train the **DeepRBP predictor** from scratch.

Most users will not need to train a model: if your goal is to compute attribution scores using the **pretrained model**, follow the quick start in the repository’s main `README.md`.

This documentation is intended for:
- reproducing the results reported in the manuscript,
- training DeepRBP on large cohorts (e.g., TCGA),
- adapting the feature set (RBPs/genes/transcripts) or the training configuration.

---

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [1. Downloading TCGA data](#1-downloading-tcga-data)
- [2. Feature specification](#2-feature-specification)
- [3. Gene ID annotation mapping (GENCODE v23)](#3-gene-id-annotation-mapping-gencode-v23)
- [4. Data preprocessing](#4-data-preprocessing)
- [5. Train/test split](#5-traintest-split)
- [6. Training the DeepRBP predictor](#6-training-the-deeprbp-predictor)
- [7. Hyperparameter optimization (optional)](#7-hyperparameter-optimization-optional)
- [Additional benchmarking experiments](#additional-benchmarking-experiments)
---

## Overview
A typical training workflow includes:

1. Downloading raw RNA-seq expression matrices and metadata (e.g., TCGA via UCSC Xena),
2. (Recommended) ensuring a stable and versioned gene identifier mapping for reproducibility,
3. Preprocessing and feature alignment using a **feature specification**,
4. Splitting samples into training and testing sets,
5. Training the DeepRBP predictor,
6. (Optional) hyperparameter optimization and benchmarking.

---

## Requirements
- Python ≥ 3.9
- Linux or macOS
- CUDA-enabled GPU recommended for training
- Sufficient disk space for TCGA-scale matrices

## 1. Downloading TCGA data  
To download the raw expression matrices and metadata:

```bash
bash slurm/data_preparation/download_data.sh
```

This retrieves data from the `UCSC Xena platform`, including:
- gene expression matrix (TPM, log-transformed),
- transcript (isoform) expression matrix (TPM, log-transformed),
- gene-level and transcript level expected counts (if provided),
- phenotype and sample metadata.

Files are stored under: `data/training_module/raw`

Although the downloaded files may also contain GTEx and TARGET samples, only TCGA samples are used for training in the default pipeline.

## 2. Feature specification
DeepRBP does not hard-code the set of genes, RBPs, or transcripts.
Instead, feature selection and ordering are defined by a **feature specification file**.

### Default (pretrained-compatible) feature specification
We provide:
```bash
data/feature_specs/DeepRBP_feature_spec.xlsx
```

This file defines:
- the list of RBPs,
- the list of genes,
- the corresponding multi-isoform transcripts,
- the exact feature ordering used during training.

### Pretrained vs custom feature specifications
- **Using the pretrained model:** input matrices must be aligned to `DeepRBP_feature_spec.xlsx`.
If the feature specification changes, pretrained checkpoints are not compatible.

- **Custom specification (advanced):** you may define a different feature set, but you must
preprocess data accordingly and train a new model.

## 3. Gene ID annotation mapping (GENCODE v23)
### Rationale
Several components of the project may require a consistent mapping between:

- GENCODE/Ensembl gene identifiers (release-specific),
- gene symbols and/or other identifiers used in downstream processing.

To ensure reproducibility across environments, we recommend using a frozen mapping file
that matches the reference release used throughout the project (e.g., **GENCODE v23**).

### The getBM.csv file
If your pipeline uses a mapping file such as:
```bash
data/annotations/getBM.csv
```

it should contain the fields required by your preprocessing and downstream utilities
(e.g., Ensembl gene ID ↔ gene symbol). The exact columns depend on the scripts that consume it.

### Recommended practice
If `getBM.csv` is required by multiple modules (not only training), treat it as a
versioned project asset and keep it under `data/annotations/`.

Document the expected reference release (e.g., “GENCODE v23”) to avoid mixing identifiers
across releases.

If this mapping is not required for your run (e.g., your input matrices already use the expected
identifiers consistently), this step can be skipped.

## 4. Data preprocessing
Raw matrices must be cleaned and transformed prior to training.

Preprocessing typically:
- aligns matrices to the feature specification,
- enforces a strict and reproducible feature order,
- applies required transformations (e.g., `log2(TPM + 1)` for selected inputs),
- exports standardized files used by the training and explainability pipelines.

Run preprocessing locally:

```bash
preprocess-data \
  --raw_data_dir data/training_module/raw \
  --output_dir data/training_module/processed \
  --feature_spec_file data/feature_specs/DeepRBP_feature_spec.xlsx
```
Custom feature specification:

```bash
preprocess-data \
  --raw_data_dir data/training_module/raw \
  --output_dir data/training_module/processed_custom \
  --feature_spec_file my_own_feature_spec.xlsx
```

**HPC execution (optional)**
```bash
sbatch slurm/data_preparation/preprocess_data.sh
```

## 5. Train/test split
After preprocessing, split samples into training and test sets using a reproducible strategy.
For TCGA-scale training, we recommend a **stratified split** by tumor type (or another relevant metadata field) to preserve class proportions across splits and support robust generalization assessment.

This step also supports restricting the split to a subset of tumor types via the configuration file (by default, all samples are included).

Common approaches include:
- random split with a fixed seed,
- stratified split (e.g., by tumor type),
- cohort-based split (train on one cohort, test on another).

```bash
split-and-save \
  --config_path src/deeprbp/configs/config_data_split.yaml \
  --output_dir data/training_module/splitted_datasets
```

### Configuration (config_data_split.yaml)
```yaml
path_files: "data/training_module/processed/TCGA"
sample_category: "detailed_category"   # metadata column used for stratification
select_samples: ["all"]                # or a list of specific categories
test_fraction: 0.2                     # 80/20 split by default
seed: 42                               # random seed for reproducibility
```
Key parameters:

- **path_files**: path to the processed dataset directory to be split
- **sample_category**: metadata column used for stratification/grouping (e.g., tumor type)
- **select_samples**: subset of categories to include (`["all"]` or an explicit list)
- **test_fraction**: fraction of samples assigned to the test split
- **seed**: random seed for reproducibility

**HPC execution (optional)**
```bash
sbatch slurm/data_preparation/split_and_save.sh
```

### Output structure
The command creates `Train/` and `Test/` directories, each containing:

- **gn_tpm.csv** — gene expression matrix (TPM)
- **RBPs_log2p_tpm.csv** — RBP expression matrix (log2(TPM + 1))
- **trans_log2p_tpm.csv** — transcript expression matrix (log2(TPM + 1))
- **phenotype_metadata.csv** — sample metadata used for selection/stratification

## 6. Training the DeepRBP predictor  
This section describes how to train the **DeepRBP predictor model** once the data have been preprocessed and split into training and test sets.

By default, training uses the `PredictorModel` class defined in:

```bash
src/deeprbp/training_module/model.py
```

This class implements the reference DeepRBP architecture used in the pretrained model and in the reported experiments.

**⚠️ Important**
If you wish to modify the network architecture (e.g., number of layers, layer sizes, activation functions) 
outside of hyperparameter optimization, you must edit the `__init__` method of the `PredictorModel` class directly. 
 
### 6.1 Configuration file
Training is controlled via a YAML configuration file specifying data paths, metadata columns, and runtime options.

**Example configuration (`config_model_train.yaml`)**
```yaml
train_path_files: "data/training_module/splitted_datasets/Train"
test_path_files: "data/training_module/splitted_datasets/Test"
sample_category: "detailed_category"
test_fraction: 0.2
getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
cuda: True
train_batch_size: 64
val_batch_size: 128
plot_results: True
```

*Key parameters:*
- **train_path_files / test_path_files**: directories containing the processed input matrices and metadata.
- **sample_category**: metadata column used for stratification.
- **getBM_path**: optional gene/transcript annotation mapping file.
- **cuda**: whether to use GPU acceleration.
- **train_batch_size / val_batch_size**: batch sizes for training and validation.
- **plot_results**: enable generation of diagnostic plots.

### 6.2 Running training locally
Once the configuration file is prepared, run the training command:

```bash
run-deeprbp-predictor \
  --config_path src/deeprbp/configs/config_model_train.yaml \
  --output_dir output/results/run_deeprbp_predictor \
  --epochs 2000 \
  --num_workers 4 \
  --min_delta 0.001 \
  --patience 30 \
  --save_top_k 1 \
  --verbose 1
```

*Main arguments*
- `--config_path` *(required)*: path to the YAML configuration file.
- `--output_dir` *(required)*: directory where checkpoints, logs, and plots will be saved.
- `--epochs`: maximum number of training epochs.
- `--num_workers`: number of data-loading workers.
- `--min_delta`: minimum improvement in validation MSE to be considered significant.
- `--patience`: number of epochs without improvement before early stopping.
- `--save_top_k`: number of best-performing checkpoints to retain.
- `--verbose`: verbosity level.

### 6.3 Running training on HPC (SLURM)
For large datasets (e.g., TCGA-scale cohorts), we recommend running training on an HPC cluster.
To submit a training job:

```bash
sbatch slurm/training_module/run_predictor.sh
```

## 7. Hyperparameter optimization (optional)
This section describes optional `hyperparameter optimization` of the DeepRBP predictor using `Optuna`.
This step is intended for advanced users and is `not required` to use the pretrained model or to run attribution analyses.

Hyperparameter optimization uses the `TunablePredictorModel` class, also defined in: `src/deeprbp/training_module/model.py`

Unlike `PredictorModel`, this class exposes a search space of `tunable architectural and training parameters` that can be explored automatically by Optuna.

### 7.1 Scope and rationale
Hyperparameter optimization aims to `minimize the mean squared error (MSE)` between predicted and observed transcript expression values on a validation set.

To reduce computational cost:
- only a fraction of the training data is used (default: 50%),
- samples are stratified by tumor type,
- early stopping and pruning are applied.

This procedure is designed to `run exclusively on HPC systems using SLURM`.

### 7.2 Configuration file
Optimization behavior is controlled via a dedicated YAML configuration file.

**Example configuration (`config_hyper_optimization.yaml`)**
```yaml
train_path_files: "data/training_module/splitted_datasets/Train"
sample_category: "detailed_category"
sample_fraction: 0.5
test_fraction: 0.2
getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
cuda: True
val_batch_size: 256
plot_results: False
```

Key parameters:
- **sample_fraction**: fraction of the training data used during optimization.
- **test_fraction**: fraction reserved for validation within the optimization loop.
- **val_batch_size**: batch size for validation evaluation.

### 7.3 Creating the Optuna study
Before running any hyperparameter optimization, an `Optuna study must be created explicitly`.
At this stage, the optimization strategy is fixed, including:

- the sampler (TPE),
- the pruning strategy (Median Pruner),
- and all sampler/pruner hyperparameters.

⚠️ Important
The optimization strategy is defined once, at study creation time, and is shared by all subsequent optimization runs that load this study.
If you wish to change the sampler, pruner, or their parameters, you must create a `new study`. (`src/deeprbp/training_module/hyperparameter_optimization/create_optuna_study.py`)

#### Optimization strategy (current implementation)
The study is created using:

- **Sampler**: `TPESampler`
  - `multivariate=True`
  - `group=True`
  - `n_startup_trials=10`
  - `n_ei_candidates=64`

- **Pruner**: `MedianPruner`
  - `n_startup_trials=0`
  - `n_warmup_steps=3` (epoch-based)
  - `interval_steps=1`

This configuration corresponds to a refined (phase-2–style) optimization strategy, focusing on efficient exploitation of a constrained, high-quality search space.

#### Create the study
To initialize the Optuna study and persistent storage:

```bash
create-optuna-study \
  --output_dir output/results/hyperparameter_optimization
```

This command:
- creates a local SQLite database (`optuna.db`),
- registers the study under the name `deeprbp_gridsearch_optuna`,
- fixes the sampler and pruner configuration for all subsequent runs.

### 7.4 Running hyperparameter optimization on SLURM
Hyperparameter optimization is executed on HPC systems using SLURM and relies on the
`TunablePredictorModel` class defined in: `src/deeprbp/training_module/model.py`

Only the hyperparameters explicitly exposed by this class and by the optimization script
are explored automatically.

SLURM scripts for hyperparameter optimization are located under:
```bash
slurm/training_module/hyperopt/
```

To launch the optimization job:
```bash
sbatch slurm/training_module/hyperopt/run_hyper_optimization_optuna.sh
```

#### Optimization workflow (current implementation)
- The optimization loads an **existing Optuna study created in advance**.
- A fixed number of trials (`--n_trials`) is executed.
- Each trial:
  - samples hyperparameters from a refined, constrained search space,
  - trains a `TunablePredictorModel`,
  - evaluates validation MSE,
  - may be pruned early based on **Optuna callbacks**.

The optimization objective is the validation loss (`validation_loss`) computed by the
PyTorch Lightning trainer.

#### Search space (refined)
The current implementation explores the following architecture and optimization hyperparameters, as exposed by the `TunablePredictorModel` class.

**Architecture-related hyperparameters**
- `num_hidden_layers` ∈ {2, 3, 4}
Number of fully connected hidden layers in the predictor network.

- `hidden1_nodes` ∈ {256, 512, 1024, 2048}
Number of neurons in the first hidden layer. Subsequent layers may keep the same size or shrink depending on the configuration.

- `uniform_nodes` ∈ {True, False}
If `True`, all hidden layers use the same number of neurons (`hidden1_nodes`).
If `False`, layer widths decrease progressively according to `node_shrink_factor`.

- `node_shrink_factor` ∈ {2, 4, 8}
Factor by which the number of neurons is reduced between successive layers when uniform_nodes=False.

- `activation_func` ∈ {relu, tanh}
Non-linear activation function applied after each hidden layer.

- `batch_norm_eps` ∈ {1e-6, 1e-5, 1e-4}
Numerical stability constant (`ε`) used in batch normalization layers.

- `batch_norm_momentum` ∈ {0.5, 0.9}
Momentum parameter controlling how batch normalization running statistics are updated.

**Optimization-related hyperparameters**
- `optimizer_name` ∈ {adam, adamW, adagrad}
Optimization algorithm used to update model parameters.

- `learning_rate` ∈ {1e-4, 3e-4, 7e-4, 1e-3, 3e-3, 1e-2, 3e-2}
Learning rate controlling the step size of parameter updates.
A discrete grid is used to avoid near-duplicate trials differing only by small learning rate changes.

- `batch_size` ∈ {32, 64, 128, 256}
Number of samples processed per training step.

- `num_epochs` ∈ {500, 1000, 2000}
Maximum number of training epochs for each trial.

Some hyperparameters are `conditionally inactive` depending on the sampled architecture
(e.g., `node_shrink_factor` when `uniform_nodes=True`).

Such inactive parameters are explicitly detected during training and recorded as `"unused"` in the Optuna trial metadata to ensure transparent interpretation of the results.

### 7.5 Analyzing optimization results
After optimization completes, results can be analyzed using:

```bash
analyze-optuna-results \
  --storage_path output/results/hyperparameter_optimization/optuna.db \
  --output_dir output/results/hyperparameter_optimization/analysis
```

This command:
- exports trial results to CSV,
- generates summary plots (optimization history, parameter importance, timelines),
- saves figures in both PNG and PDF formats.

## Additional benchmarking experiments
In addition to training the DeepRBP predictor, this module contains
benchmarking experiments used in the manuscript, including:

- comparisons with traditional machine learning regressors,
- tumor-specific versus general training strategies.

These experiments are documented in dedicated submodules:

- `benchmark_methods/README.md`
- `tumor_specific_training/README.md`