# DeepRBP
## Inferring RBP–gene regulatory links from bulk RNA-seq via deep learning and model attributions

**DeepRBP** is a deep learning–based framework designed to prioritize regulatory relationships between RNA-binding proteins (RBPs) and genes/transcripts from RNA-seq data.  
It combines:
- a **prediction module** that estimates transcript abundance from gene and RBP expression, and
- an **explainability module** that computes interpretable RBP–transcript and RBP–gene scores.

---

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Quick start: use the pretrained DeepRBP model](#quick-start-use-the-pretrained-deeprbp-model)
  - [1. Download the pretrained model bundle](#1-download-the-pretrained-model-bundle)
  - [2. Recommended local layout](#2-recommended-local-layout)
  - [3. Do I need data preparation or training?](#3-do-i-need-data-preparation-or-training)
- [Run the pretrained DeepRBP predictor on your data](#run-the-pretrained-deeprbp-predictor-on-your-data)
  - [Required input files](#required-input-files)
  - [Feature compatibility](#feature-compatibility)
  - [Step 1. Run the pretrained predictor](#step-1-run-the-pretrained-predictor)
  - [Step 2. Run the explainer](#step-2-run-the-explainer)
  - [Expected outputs](#expected-outputs)
- [Data preparation](#data-preparation)
  - [When this section is needed](#when-this-section-is-needed)
  - [Using DeepRBP with your own RNA-seq data](#using-deeprbp-with-your-own-rna-seq-data)
  - [Standard pipeline: kallisto → DeepRBP inputs](#standard-pipeline-kallisto--deeprbp-inputs)
- [Training a model from scratch (advanced)](#training-a-model-from-scratch-advanced)
- [Additional module-specific documentation](#additional-module-specific-documentation)
- [Citation](#citation)
- [License](#license)

---

## Overview
Alternative splicing is a key regulatory mechanism in many biological processes and diseases, particularly cancer. RNA-binding proteins (RBPs) play a central role in this regulation, but experimentally identifying RBP–target interactions (e.g. via CLIP-based assays) is costly and limited in scope.

**DeepRBP** provides a computational alternative to prioritize putative RBP–gene and RBP–transcript regulatory relationships from bulk RNA-seq data.

This repository provides:
- the full DeepRBP pipeline,
- scripts to preprocess RNA-seq data,
- a pretrained model,
- and tools to compute explainability scores.

For most users, the recommended workflow is:

1. **Download the pretrained model**
2. **Prepare your data in the expected input format**
3. **Run the predictor evaluator** to check compatibility and performance
4. **Run the explainer** to obtain TxRBP and GxRBP scores

---

## Installation

### Requirements
- Python ≥ 3.9  
- Linux or macOS  
- *(Optional but recommended)* CUDA-enabled GPU  

### Install DeepRBP
```bash
git clone https://github.com/ML4BM-Lab/DeepRBP.git
cd DeepRBP

conda create -n DeepRBP python=3.9
conda activate DeepRBP
pip install -e .
```

## Quick start: use the pretrained DeepRBP model
We provide a pretrained DeepRBP predictor hosted on Hugging Face.

🔗 **Model repository**
https://huggingface.co/ML4BM-Lab/DeepRBP

The pretrained model consists of four required files:
- `model.ckpt` — trained DeepRBP predictor checkpoint
- `scaler.joblib` — fitted input scaler
- `sigma.npy` — output scaling parameter
- `DeepRBP_feature_spec.xlsx` — feature manifest (RBPs/genes/transcripts + exact order)

⚠️ Important
All four files are required for correct inference and explainability.
Treat them as a single bundle and keep them together.

1. **Download the pretrained model bundle**

**Option A (recommended): download from Hugging Face using Git LFS**
The model checkpoint is stored using **Git Large File Storage (Git LFS)**.

**Install and initialize Git LFS**
```bash
git lfs install
```

If `git lfs` is not available, install it first:
- **macOS**: `brew install git-lfs`
- **Ubuntu/Debian**: `sudo apt install git-lfs`

**Clone the pretrained model repository**
```bash
git clone https://huggingface.co/ML4BM-Lab/DeepRBP pretrained_model
```

This downloads all required files automatically.

If the repository is private, Hugging Face may ask for credentials.
Use your Hugging Face username and an access token (not your password).

**Option B: manual download**
If you prefer not to use Git LFS:

1. Open: `https://huggingface.co/ML4BM-Lab/DeepRBP`
2. Download:
  - `model.ckpt`
  - `scaler.joblib`
  - `sigma.npy`
  - `DeepRBP_feature_spec.xlsx`
3. Place them together in a local folder.

---

2. **Recommended local layout**
```text
DeepRBP/
├── pretrained_model/
│   ├── model.ckpt
│   ├── scaler.joblib
│   ├── sigma.npy
│   └── DeepRBP_feature_spec.xlsx
├── data/
│   ├── my_dataset/
│   │   ├── RBPs_log2p_tpm.csv
│   │   ├── trans_log2p_tpm.csv
│   │   ├── gn_tpm.csv
│   │   └── phenotype_metadata.csv      # optional
│   └── annotations/
│       └── getBM.csv                   # required for some workflows
├── output/
└── src/
```

---

3. **Do I need data preparation or training?**
Use this rule of thumb:

- If you already have:
  - `RBPs_log2p_tpm.csv`
  - `trans_log2p_tpm.csv`
  - `gn_tpm.csv`
  - optional `phenotype_metadata.csv`

and they are aligned to `DeepRBP_feature_spec.xlsx`, you can run **the pretrained model directly**.
- If you do not yet have these files, go to [Data preparation](#data-preparation).
- If you want to change the feature set (RBPs/genes/transcripts) or train a new model, go to [Training a model from scratch (advanced)](#training-a-model-from-scratch-advanced).
---

## Run the pretrained DeepRBP predictor on your data
This is the recommended entry point for most users.

If your dataset is already formatted as DeepRBP expects, you do not need to retrain the model.

The recommended workflow is:

1. **Run the pretrained DeepRBP predictor first**  
   This verifies that your matrices are compatible with the pretrained model and provides standard performance metrics.

2. **Run the DeepRBP explainer second**  
   This produces transcript-level and gene-level attribution scores.

## Required input files
DeepRBP expects the following files per dataset or per experimental group:

- **RBPs_log2p_tpm.csv**
RBP expression matrix (samples × RBPs), in log2(TPM + 1)

- **trans_log2p_tpm.csv**
Transcript expression matrix (samples × transcripts), in log2(TPM + 1)

- **gn_tpm.csv**
Gene expression matrix (samples × genes), in TPM

- **phenotype_metadata.csv** (*optional*)
Sample metadata (e.g. tissue, condition, or tumor type)

If your data can be converted into this format, it can be used by DeepRBP.

## Feature compatibility
DeepRBP expects your matrices to match the feature manifest used during training:
```text
DeepRBP_feature_spec.xlsx
```

This file defines:
- the list of RBPs, genes, and transcripts,
- and their **exact order**.

⚠️ If you are using the pretrained model, your input matrices **must be aligned** to this feature specification.

If you want to use a different feature set, you must preprocess accordingly and train a new model.

---

### Annotation mapping (`getBM.csv`)

Some evaluation and downstream utilities require a frozen gene/transcript annotation mapping file:

```text
data/annotations/getBM.csv
```

For pretrained-model compatibility, this file must match the annotation reference used to build DeepRBP_feature_spec.xlsx: GENCODE v23 / Ensembl 81 / GRCh38.p3.

The recommended way to obtain it is from the pretrained model bundle:

git clone https://huggingface.co/ML4BM-Lab/DeepRBP pretrained_model

mkdir -p data/annotations
cp pretrained_model/getBM.csv data/annotations/getBM.csv

### Step 1. Run the pretrained predictor
You can run the pretrained `DeepRBP predictor` on your processed dataset. This is the recommended first step for a new dataset.

It is useful to:
- verify that your inputs are compatible with the pretrained model,
- generate transcript-level predictions,
- obtain standard metrics (Spearman/Pearson/R²/MSE).

**Minimal predictor config (YAML)**
```yaml
# src/deeprbp/configs/config_model_eval.yaml
test_path_files: "data/my_dataset"

getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
sample_category: "detailed_category"   # optional; remove if not needed
cuda: True
val_batch_size: 256
plot_results: True
seed: 42

# IMPORTANT: use the scaler bundled with the pretrained model
scaler_dir: "pretrained_model"
```

**Run locally**
```bash
run-deeprbp-evaluator \
  --config_path src/deeprbp/configs/config_model_eval.yaml \
  --model_checkpoint pretrained_model/model.ckpt \
  --output_dir output/results/eval_pretrained \
  --num_workers 4 \
  --verbose 1
```

This command runs the pretrained DeepRBP predictor on your dataset and reports standard performance metrics.

**Run on HPC / SLURM**
```bash
sbatch slurm/training_module/run_evaluate_predictor.sh
```

⚠️ **Bundled scaler (recommended)**
If you use `pretrained_model/model.ckpt`, always evaluate with the bundled `pretrained_model/scaler.joblib`.
Do not refit a new scaler when using the pretrained checkpoint.

---

### Step 2. Run the explainer
Once predictor compatibility has been confirmed, you can compute attribution scores with the DeepRBP explainer.

DeepRBP computes:
- **TxRBP scores**: transcript × RBP scores
- **GxRBP scores**: gene × RBP scores obtained by collapsing transcript scores to genes

By default, the explainer uses `DeepLIFT`.

**Minimal explainer config (single-run mode)**
Use this when:
- your dataset has no TCGA-like categories, or
- you do not provide `phenotype_metadata.csv`, or
- you simply want one run over all samples.

```yaml 
# src/deeprbp/configs/config_model_explain.yaml
test_path_files: "data/my_dataset"
getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
explanation_method: "DeepLIFT"
reference_data: "knockout_reference"
batch_reduction_method: "t-statistic"
gene_collapse_method: "max_absolute_value"
save_per_sample_scores: false
```

**Run locally**
```bash
run-deeprbp-explainer \
  --config_path src/deeprbp/configs/config_model_explain.yaml \
  --model_ckpt_path pretrained_model/model.ckpt \
  --scaler_dir pretrained_model \
  --output_dir output/results/explainer_all_samples
```

**Optional: hidden-layer attributions**
You can also compute attributions to the last hidden layer (HL × RBP) by adding:

```bash
run-deeprbp-explainer \
  --config_path src/deeprbp/configs/config_model_explain.yaml \
  --model_ckpt_path pretrained_model/model.ckpt \
  --scaler_dir pretrained_model \
  --output_dir output/results/explainer_all_samples \
  --analyze_hidden_layer
```

**Run on HPC / SLURM**
```bash
cd slurm/explainability_module

# Standard run
sbatch run_explainer.sh

# Hidden-layer attributions
ANALYZE_HL=true sbatch run_explainer.sh
```

### Expected outputs
**Predictor**
Typical outputs include:
- performance metrics,
- plots (if enabled),
 -prediction diagnostics.

**Explainer**
Single-run mode typically produces:
```bash
<output_dir>/
  df_scores_TxRBP.csv
  df_scores_GxRBP.csv
  result_table.csv
  df_scores_TxRBP_per_sample.csv   # only if save_per_sample_scores: true
  df_scores_HLxRBP.csv             # only with --analyze_hidden_layer
```

---

## Data preparation
### When this section is needed
Use this section only if you do not yet have DeepRBP-compatible input files:
- `RBPs_log2p_tpm.csv`
- `trans_log2p_tpm.csv`
- `gn_tpm.csv`
- optional `phenotype_metadata.csv`

This section explains how to convert raw RNA-seq quantifications into the format expected by the pretrained model and the explainer.

### Using DeepRBP with your own RNA-seq data
DeepRBP can be applied to any RNA-seq dataset, including:
- custom cancer cohorts,
- in-vitro experiments,
- RBP knockdown or perturbation studies.

At a high level, the required steps are:
1. **Quantify expression**
  - transcript-level TPMs,
  - gene-level TPMs (e.g. Salmon, Kallisto, or equivalent)

2. **Build DeepRBP input matrices**
  - transcript TPM matrix,
  - gene TPM matrix,
  - RBP expression matrix (subset of genes),
  - optional sample metadata

3. **Apply standard DeepRBP transformations**
  - `log2(TPM + 1)` where required

4. **Export files using the expected filenames**
  - `RBPs_log2p_tpm.csv`
  - `trans_log2p_tpm.csv`
  - `gn_tpm.csv`
  - `phenotype_metadata.csv`

### Standard pipeline: kallisto → DeepRBP inputs
If your samples were quantified with kallisto, DeepRBP provides a preprocessing pipeline that converts abundance.tsv files into DeepRBP-ready inputs using the same preprocessing logic as model training.

The pipeline expects:
- one folder per sample under kallisto_output/, each containing abundance.tsv,
- a metadata CSV with a Run column (sample IDs),
- the DeepRBP feature specification file.

#### Basic execution
```bash
preprocess-user-data \
  --kallisto_output <DATASET_DIR>/kallisto_output \
  --metadata_csv   <DATASET_DIR>/metadata.csv \
  --feature_spec   pretrained_model/DeepRBP_feature_spec.xlsx \
  --output_dir     <DATASET_DIR>/deeprbp_inputs \
  --group_col tissue_type
```

This command creates one output folder per group (e.g. per tissue or condition), each containing:
- `RBPs_log2p_tpm.csv`
- `trans_log2p_tpm.csv`
- `gn_tpm.csv`
- `phenotype_metadata.csv`

Optional arguments allow further splitting (e.g. by tumor stage or cell line).

#### Execution on HPC systems (optional)
```bash
sbatch preprocess_datauser.sh \
  --kallisto_output <DATASET_DIR>/kallisto_output \
  --metadata_csv   <DATASET_DIR>/metadata.csv \
  --feature_spec   pretrained_model/DeepRBP_feature_spec.xlsx \
  --output_dir     <DATASET_DIR>/deeprbp_inputs \
  --group_col      tissue_type
```

This is functionally equivalent to direct execution and is recommended for large-scale preprocessing.

---

## Training a model from scratch (advanced)
Most users do **not** need to train a model.

Training from scratch is only needed if you want to:
- reproduce the full training workflow,
- train on a large cohort such as TCGA,
- change the feature specification,
- or optimize a new architecture.

For that workflow, see:
```bash
src/deeprbp/training_module/README.md
```

That documentation covers:
- TCGA download,
- feature specification,
- preprocessing for training,
- train/test splits,
- predictor training,
- and optional hyperparameter optimization.

---

## Additional module-specific documentation

The repository also includes additional documentation for specialized workflows:
- `src/deeprbp/training_module/README.md`
- `src/deeprbp/explainability_module/README.md`
- `src/deeprbp/explainability_module/postar_validation/README.md`
- `src/deeprbp/explainability_module/real_knockdowns/README.md`
- `src/deeprbp/explainability_module/complex_analysis/README.md`

If you want to reproduce the **TCGA experiments reported in the manuscript**, including tumor-type/category runs and downstream validations, please follow the module-specific documentation in: `src/deeprbp/explainability_module/README.md`

## Citation
If you use DeepRBP in your work, please cite:

DeepRBP: A novel deep neural network for inferring splicing regulation
https://doi.org/10.1101/2024.04.11.589004
