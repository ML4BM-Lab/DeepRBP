# DeepRBP
## Inferring RBP–gene regulatory links from bulk RNA-seq via deep learning and model attributions

**DeepRBP** is a deep learning–based framework designed to infer regulatory relationships between RNA-binding proteins (RBPs) and genes/transcripts from RNA-seq data.  
It combines a predictor model for transcript abundance with explainability methods to produce interpretable RBP–gene and RBP–transcript scores.

## 📄 Publication
**DeepRBP: A novel deep neural network for inferring splicing regulation**  
https://doi.org/10.1101/2024.04.11.589004 

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Quick start: using a pretrained DeepRBP model](#quick-start-using-a-pretrained-deeprbp-model)
  - [Download the pretrained model](#download-the-pretrained-model)
  - [Recommended local layout](#recommended-local-layout)
- [Data preparation](#data-preparation)
  - [Required input format](#required-input-format)
  - [Feature compatibility](#feature-compatibility)
  - [Using DeepRBP with your own RNA-seq data](#using-deeprbp-with-your-own-rna-seq-data)
- [Evaluate the DeepRBP predictor](#evaluate-the-deeprbp-predictor)
- [Running the DeepRBP explainer](#running-the-deeprbp-explainer)
- [Citation](#citation)
- [License](#license)

## Overview
Alternative splicing is a key regulatory mechanism in many biological processes and diseases, particularly cancer. RNA-binding proteins (RBPs) play a central role in this regulation, but experimentally identifying RBP–target interactions (e.g. via CLIP-based assays) is costly and limited in scope.

**DeepRBP** provides a computational alternative to prioritize putative RBP–gene and RBP–transcript regulatory relationships from RNA-seq data.

DeepRBP consists of two main components:

### Prediction module
A deep neural network that predicts transcript abundance from gene and RBP expression.

### Explainability module
Attribution methods (e.g., DeepLIFT) are used to compute interpretable scores that quantify the contribution of each RBP to each gene or transcript.

This repository provides:
- the full DeepRBP pipeline,
- scripts to preprocess RNA-seq data,
- a pretrained model,
- and tools to compute explainability scores.

This README focuses on using a pretrained model to compute explainability scores, which is the recommended entry point for most users.

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

## Quick start: using a pretrained DeepRBP model
We provide a **pretrained DeepRBP predictor** hosted on Hugging Face.

🔗 **Model repository**  
https://huggingface.co/ML4BM-Lab/DeepRBP

The pretrained model consists of three required files:
- `model.ckpt` — trained DeepRBP predictor checkpoint  
- `scaler.joblib` — fitted input scaler  
- `sigma.npy` — output scaling parameter  
- `DeepRBP_feature_spec.xlsx` — feature manifest (RBPs/genes/transcripts + exact order)

⚠️ Important  
All four files are required for correct inference and explainability.  
Treat them as a single bundle and keep them together.

### Download the pretrained model

#### Option A (recommended): download from Hugging Face using Git LFS
The model checkpoint is stored using **Git Large File Storage (Git LFS)**.

##### 1. Install and initialize Git LFS (one-time setup)
```bash
git lfs install
```

If git lfs is not available, install it first:
- **macOS**: `brew install git-lfs`
- **Ubuntu/Debian**: sudo apt install git-lfs

##### 2. Clone the pretrained model repository
```bash
git clone https://huggingface.co/ML4BM-Lab/DeepRBP pretrained_model
```

This will download all required files automatically (including DeepRBP_feature_spec.xlsx).

**Note**:
If the repository is private, you may be prompted for credentials. Use your Hugging Face username and an access token (not your password).
Tokens can be created at: `https://huggingface.co/settings/tokens`

#### Option B: manual download from the Hugging Face website
If you prefer not to use Git LFS, you can download the files manually:

1. Open: `https://huggingface.co/ML4BM-Lab/DeepRBP`
2. Download:
  - `model.ckpt`
  - `scaler.joblib`
  - `sigma.npy`
  - `DeepRBP_feature_spec.xlsx`
3. Place them in a local directory (see layout below).

Suggested project structure

We recommend organizing your project as follows:

```text
DeepRBP/
├── pretrained_model/
│   ├── model.ckpt
│   ├── scaler.joblib
│   └── sigma.npy
├── data/
│   └── feature_specs/
│       └── DeepRBP_feature_spec.xlsx
├── output/      # predictions and explainability results
└── src/
```

## Data preparation
DeepRBP can be applied to TCGA datasets or to your own RNA-seq data.  
Regardless of the data source, the model always expects the same preprocessed input format.

This section explains:
- what files are required,
- how to structure your own data,
- and how to adapt existing RNA-seq datasets.

---

### 📁 Required input format
DeepRBP expects the following files per dataset:

- **`RBPs_log2p_tpm.csv`**  
  RBP expression matrix (samples × RBPs), in `log2(TPM + 1)`
- **`trans_log2p_tpm.csv`**  
  Transcript expression matrix (samples × transcripts), in `log2(TPM + 1)`
- **`gn_tpm.csv`**  
  Gene expression matrix (samples × genes), in TPM
- **`phenotype_metadata.csv`** (optional)  
  Sample metadata (e.g. tissue, condition, or tumor type)

If your data can be converted into this format, **it can be used by DeepRBP**.

---
 
**⚠️ Feature compatibility (important)**
DeepRBP expects your matrices to match the feature manifest used during training:
`DeepRBP_feature_spec.xlsx` (RBPs/genes/transcripts + exact order).

If you are using the pretrained model, you must align your data to this manifest.
If you want to use a different feature set, you’ll need to retrain the model
(see `src/deeprbp/training_module/README.md`).

### Using DeepRBP with your own RNA-seq data
DeepRBP can be applied to any RNA-seq dataset, including:

- custom cancer cohorts,
- in-vitro experiments,
- RBP knockdown datasets.

High-level steps:

1. Quantify expression:
  - transcript-level TPMs,
  - gene-level TPMs (e.g., Salmon, Kallisto, or equivalent).

2. Build the required matrices:
  - transcript TPM matrix,
  - gene TPM matrix,
  - RBP expression matrix (subset of genes),
  - optional sample metadata.

3. Apply transformations:
  - log2(TPM + 1) where required.

4. Export using the expected filenames listed above.

## Evaluate the DeepRBP predictor
You can evaluate the **pretrained DeepRBP predictor** on your processed dataset (same input format as in Data preparation). 
This is useful to sanity-check that your matrices are compatible and to obtain standard performance metrics (Spearman/Pearson/R²/MSE).

### Minimal evaluation config (YAML)

```bash
# src/deeprbp/configs/config_model_eval.yaml
test_path_files: "data/my_dataset"   # folder containing RBPs_log2p_tpm.csv, trans_log2p_tpm.csv, gn_tpm.csv, phenotype_metadata.csv (optional)

getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
sample_category: "detailed_category" # you can remove this
cuda: True
val_batch_size: 256
plot_results: True
seed: 42

# IMPORTANT: use the bundled scaler from the pretrained model folder
scaler_dir: "pretrained_model"
# scaler_mode defaults to "tcga" (recommended). Do not refit a new scaler when using the pretrained checkpoint.
```

### Run locally
```bash
run-deeprbp-evaluator \
  --config_path src/deeprbp/configs/config_model_eval.yaml \
  --model_checkpoint pretrained_model/model.ckpt \
  --output_dir output/results/eval_pretrained \
  --num_workers 4 \
  --verbose 1
```

### Run on HPC/SLURM
```bash
sbatch slurm/training_module/run_evaluate_predictor.sh
```

**⚠️ Bundled scaler (recommended)**: If you use `pretrained_model/model.ckpt`, always evaluate with the bundled `pretrained_model/scaler.joblib` (set `scaler_dir: pretrained_model`). 
Only refit a scaler if you trained your own predictor from scratch (see `src/deeprbp/training_module/README.md`).

---

## Running the DeepRBP explainer
DeepRBP computes attribution scores at the **transcript level** (TxRBP) for each sample, which can be seen as a 3D tensor
(**transcripts × RBPs × samples**). Positive values indicate that an RBP contributes to **increasing** the predicted transcript
abundance, while negative values indicate a contribution to **decreasing** it.

To obtain a single, cohort-level score per transcript–RBP pair, DeepRBP collapses per-sample attributions using a
**t-statistic** reduction (configurable via `batch_reduction_method`), yielding a 2D matrix (**transcripts × RBPs**).

Gene-level scores (GxRBP) are derived by collapsing transcript scores to genes using `getBM.csv` (Transcript_ID → Gene_ID).
With `gene_collapse_method: max_absolute_value`, DeepRBP selects, for each (gene, RBP), the transcript with the largest
absolute TxRBP score and keeps its sign. The selected representative transcript per gene–RBP is reported in `result_table.csv`.

DeepRBP computes these attributions from a trained predictor checkpoint and a compatible input dataset
(same format as in **Data preparation**: `RBPs_log2p_tpm.csv`, `trans_log2p_tpm.csv`, `gn_tpm.csv`, and optional `phenotype_metadata.csv`).

By default, DeepRBP can compute attributions with **DeepLIFT**: DeepLIFT (Shrikumar, Greenside, and Kundaje, 2017)
*Learning important features through propagating activation differences*, ICML (PMLR), pp. 3145–3153.

The repository also includes optional validation workflows (e.g., POSTAR and real knockdowns), documented here:
- `src/deeprbp/explainability_module/README.md`
- `src/deeprbp/explainability_module/postar_validation/README.md`
- `src/deeprbp/explainability_module/real_knockdowns/README.md`
- `src/deeprbp/explainability_module/complex_analysis/README.md`

If you want to reproduce the **TCGA** experiments reported in our paper (tumor-type/category runs, filters, POSTAR validation, etc.),
please follow the internal documentation in `src/deeprbp/explainability_module/README.md`.

### Single-run mode (recommended for custom datasets)
If your dataset has **no TCGA-like categories** or you do **not** provide `phenotype_metadata.csv`, DeepRBP runs in **single-run mode** (all samples, no filtering).

#### Minimal YAML (no metadata / no categories)
```yaml
# src/deeprbp/configs/config_model_explain.yaml
test_path_files: "data/my_dataset"         # folder with RBPs_log2p_tpm.csv, trans_log2p_tpm.csv, gn_tpm.csv (phenotype_metadata.csv optional)
getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
explanation_method: "DeepLIFT"
reference_data: "knockout_reference"
batch_reduction_method: "t-statistic"
gene_collapse_method: "max_absolute_value"
save_per_sample_scores: false   # optional: if True, also saves TxRBP per sample (Tx × RBP × Samples)
```

#### Option 1: Run locally (CLI)
```bash
run-deeprbp-explainer \
  --config_path src/deeprbp/configs/config_model_explain.yaml \
  --model_ckpt_path pretrained_model/model.ckpt \
  --scaler_dir pretrained_model \
  --output_dir output/results/explainer_all_samples
```

##### Optional: hidden-layer attributions (DeepLIFT only)
You can also compute attributions to the last hidden layer (HL × RBP) by adding:
```bash
  --analyze_hidden_layer
```

#### Option 2: Submit a job on HPC (SLURM)
We provide a SLURM script that can toggle hidden-layer analysis via an environment variable.

```bash
cd slurm/explainability_module

# Without hidden layer (default):
sbatch run_explainer.sh

# With hidden layer:
ANALYZE_HL=true sbatch run_explainer.sh
```

The script builds the command and adds --analyze_hidden_layer only when ANALYZE_HL=true.
Remember: hidden-layer attributions require explanation_method: "DeepLIFT".

#### Outputs
Single-run mode (all samples):
```text
<output_dir>/
  df_scores_TxRBP.csv
  df_scores_GxRBP.csv
  result_table.csv
  df_scores_TxRBP_per_sample.csv   # only if save_per_sample_scores: true
  df_scores_HLxRBP.csv             # only with --analyze_hidden_layer
```

## Citation
<!-- TODO: add BibTeX / preferred citation format -->
If you use DeepRBP in your work, please cite:

DeepRBP: A novel deep neural network for inferring splicing regulation
https://doi.org/10.1101/2024.04.11.589004

## License
<!-- TODO: add license text or link to LICENSE file -->



# QUE EL USUARIO PUEDA VER COMO PREPROCESAR CON EL 
# EJEMPLO DEL REAL KDS O LA INFO DE MARIA (ESTO HAY Q HACER) PARA VER COMO PUEDEN PRE-
# PROCESAR LOS DATOS (PERO EN TRAINIG_MODULE SOBRE TODO ES ENSEÑAR LOS PASOS BIEN)

## AQUI METERLE EL EVALUATE MODEL QUE NUNCA METÍ!!!