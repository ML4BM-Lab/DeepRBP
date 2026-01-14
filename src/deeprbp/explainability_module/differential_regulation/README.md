
# Differential regulation analysis

## Pipeline overview

```text
FASTQ
  ↓
kallisto quant
  ↓
DeepRBP preprocessing
  ↓
(1) Per-condition explainability
  ↓
(2) Differential expression (context)
  ↓
(3) Differential regulation
  ↓
(4) RBP / transcript ranking
```

This module implements **condition-aware regulatory analysis** in DeepRBP by
combining:

1. Per-condition DeepLIFT explainability scores
2. Differential expression (DE) analysis of RBPs, transcripts, and genes

The goal is to study **regulatory changes between biological conditions**
(e.g. control vs treatment, disease stages, experimental perturbations),
separating **expression effects** from **regulatory effects**.

---

## Dataset used in this module
The datasets listed below are used as **fully reproducible reference use cases**
for the differential regulation analysis presented in this work.
The same computational pipeline is applied to TCGA data (used in the main
analyses) as well as to additional public and user-provided RNA-seq datasets.

Each dataset is briefly described below. Detailed biological context and
dataset-specific documentation are available in the corresponding dataset
directories.

---

### GSE114564 — Liver disease progression
**Use case**: condition-aware regulatory rewiring across disease stages.

Paired-end RNA-seq data from human liver tissue spanning multiple stages of
liver disease progression:

- Normal liver (NL)
- Chronic hepatitis (CH)
- Liver cirrhosis (LC)
- Dysplastic nodule (DN)
- Early HCC (eHCC)
- Advanced HCC (avHCC)

This structured disease trajectory makes GSE114564 well suited for evaluating
**regulatory changes across progressive pathological states**, beyond
differential expression alone.

📄 Extended documentation:
`/data/explainability_module/differential_regulation/GSE114564/README.md`

---

### GSE101432 — Human liver transcriptome and isoform complexity
**Use case**: isoform- and splicing-aware regulatory analysis.

Paired-end RNA-seq data from a diverse collection of human liver samples,
including normal liver, benign adjacent tissue, primary and relapsed
hepatocellular carcinoma (HCC), and liver cancer cell lines.

This dataset provides a rich resource for studying isoform-level regulation
and alternative splicing, and is particularly informative for analyses
focused on post-transcriptional regulatory mechanisms mediated by RBPs.

📄 Extended documentation:
`/data/explainability_module/differential_regulation/GSE101432/README.md`

---

### RNA-seq quantification
All datasets in this module are processed using the same unified RNA-seq
quantification pipeline, ensuring full comparability across analyses.

- **kallisto**: v0.51.1
- **Reference transcriptome**: GENCODE v23
- **Mode**: paired-end
- **Bootstraps**: 100

Reference FASTA and index files are shared under `data/annotation/`.

---

#### Directory structure
```text
<DATASET_ID>/
├── SRR_Acc_List.txt
├── SraRunTable.csv
├── raw_fastq/
├── kallisto_output/
└── README.md
```

### How to run (per dataset)
The same SLURM-based pipeline is used for all datasets.
Replace `<DATASET_ID>` with the desired dataset (e.g. `GSE114564`,
`GSE101432`).

1. **Download FASTQ files**
```bash
sbatch /scratch/jsanchoz/DeepRBP/slurm/explainability_module/differential_regulation/kallisto_quant/download_fastq_array.sh <DATASET_ID>
```
2. **Run kallisto quant**
```bash
sbatch /scratch/jsanchoz/DeepRBP/slurm/explainability_module/differential_regulation/kallisto_quant/kallisto_array.sh <DATASET_ID>
```
3. **(Optional) QC summary**
```bash
bash /scratch/jsanchoz/DeepRBP/slurm/explainability_module/differential_regulation/kallisto_quant/kallisto_qc_summary.sh \
  /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/<DATASET_ID>/kallisto_output \
  /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/<DATASET_ID>/SraRunTable.csv \
  kallisto_qc
```
⚠️ QC is for exploratory assessment only; no samples are filtered at this stage.

### Data preprocessing
After RNA-seq quantification with **kallisto**, all datasets in this module are
preprocessed into **DeepRBP-compatible input matrices** following the same
standard preprocessing pipeline used during model training.

All downstream explainability and differential regulation analyses
**assume this standardized DeepRBP input format**.

> ℹ️ A full description of the required input format and preprocessing logic is
available in the main project README under Data preparation.

---

Starting from `kallisto_output/`, the pipeline generates, per **biological group or condition**, 
the following files:

- `RBPs_log2p_tpm.csv` — RBP expression matrix
- `trans_log2p_tpm.csv` — transcript expression matrix
- `gn_tpm.csv` — gene expression matrix
- `phenotype_metadata.csv` — sample metadata

All transformations (TPM handling, `log2(TPM + 1)`, feature ordering) strictly
follow the DeepRBP feature specification used during training
(`DeepRBP_feature_spec.xlsx`).

---

#### Execution (per dataset)
The same preprocessing command is used for all user-provided RNA-seq datasets
quantified with **kallisto**.  
Datasets are processed **group-wise**, where each group corresponds to a
subset of samples defined by one or more metadata columns.

By default, samples are grouped using a single metadata column specified via
`--group_col`. Optional flags allow further subdivision of samples by tumor
stage or cell line when such information is available.

#### GSE114564 — disease stage–aware preprocessing
Samples are grouped by disease stage to enable condition-aware regulatory
analysis.

##### Direct execution (interactive or local)
```bash
preprocess-user-data \
  --kallisto_output /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/GSE114564/kallisto_output \
  --metadata_csv   /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/GSE114564/kallisto_qc_with_kallisto_qc.csv \
  --feature_spec   /scratch/jsanchoz/DeepRBP/data/training_module/feature_specs/DeepRBP_feature_spec.xlsx \
  --output_dir     /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/GSE114564/processed \
  --group_col      disease_state \
  > /scratch/jsanchoz/DeepRBP/output/logs/GSE114564_preprocess.log 2>&1
```

This generates one DeepRBP input folder per stage (e.g. `NL`, `CH`, `LC`, `DN`, `eHCC`,
`avHCC`).

##### SLURM execution (recommended for HPC)
The same command can be executed on an HPC cluster using the provided SLURM
wrapper, which forwards all command-line arguments to `preprocess-user-data`.

```bash
sbatch /scratch/jsanchoz/DeepRBP/slurm/data_preparation/preprocess_datauser.sh \
  --kallisto_output /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/GSE114564/kallisto_output \
  --metadata_csv   /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/GSE114564/kallisto_qc_with_kallisto_qc.csv \
  --feature_spec   /scratch/jsanchoz/DeepRBP/data/training_module/feature_specs/DeepRBP_feature_spec.xlsx \
  --output_dir     /scratch/jsanchoz/DeepRBP/data/explainability_module/differential_regulation/GSE114564/processed \
  --group_col      disease_state
```

Job output is written to the SLURM log directory
(e.g. `/scratch/jsanchoz/DeepRBP/output/logs/preprocess_user_data_<JOBID>.out`).


#### GSE101432 — sample-type / isoform-aware preprocessing
Samples are grouped by sample type (e.g. normal, tumor, relapse, cell
line), enabling isoform- and splicing-focused regulatory analyses.

```bash
preprocess-user-data \
  --kallisto_output GSE101432/kallisto_output \
  --metadata_csv   GSE101432/metadata.csv \
  --feature_spec   pretrained_model/DeepRBP_feature_spec.xlsx \
  --output_dir     GSE101432/processed \
  --group_col      ? \
```

### Output structure
```text
<DATASET_ID>/processed/
├── <GROUP_1>/
│   ├── RBPs_log2p_tpm.csv
│   ├── trans_log2p_tpm.csv
│   ├── gn_tpm.csv
│   └── phenotype_metadata.csv
├── <GROUP_2>/
│   └── ...
```

These folders are used directly by the **DeepRBP explainability and differential regulation analysis pipelines**.

## 1. Per-condition explainability (independent)
This step:
- treats each condition independently,
- computes DeepLIFT-based regulatory scores,
- produces condition-specific explainability matrices.

For each biological condition:

- Samples belonging to that condition are loaded
- DeepRBP DeepLIFT explainability is computed
- Results are saved in a condition-specific folder

Each condition is treated as an **independent dataset**.
No cross-condition statistics are computed at this stage.

This step is implemented by: `compute_scores_per_condition.py`.

No statistical comparison between conditions is performed at this stage.

Each condition folder contains:
  - df_scores_TxRBP_per_sample.csv (if computed)
  - df_scores_TxRBP.csv
  - df_scores_GxRBP.csv
  - result_table.csv
  - (other cached/aux outputs produced by the core explainer)

These outputs are designed to be consumed by downstream differential
expression and regulatory comparison analyses.

---

### How to run  
This script can be executed locally or on an HPC cluster.
In all cases, the same arguments and configuration file are used.

#### Run locally (installed command)
If DeepRBP is installed as a package, the script is available as:

```bash
run-deeprbp-differential-regulation-scores \
  --config_path src/deeprbp/configs/config_diff_regulation.yaml \
  --model_ckpt_path /path/to/deeprbp_predictor.ckpt \
  --output_dir /path/to/output/diff_reg/Liver \
  --select_category Liver_Hepatocellular_Carcinoma \
  --conditions Primary_Tumor,Solid_Tissue_Normal
```

This is the recommended execution mode for standard and production runs.  

#### HPC execution (optional)
For batch execution on HPC systems, a reference SLURM script is provided:
```bash
sbatch slurm/explainability_module/differential_regulation/run_compute_scores_per_condition.sh
```

The SLURM script wraps the installed command and handles environment setup,
logging, and resource allocation.

### Command-line arguments
The script accepts the following arguments:

- `--config_path` (required)
Path to the explainability configuration YAML
(e.g. `config_diff_regulation.yaml`).

- `--model_ckpt_path` (required)
Path to the trained DeepRBP predictor checkpoint.

- `--output_dir` (required)
Base output directory. Results will be written under:
```text
<output_dir>/explainability_scores/<Condition>/
```

- `--select_category` (required)
Dataset category or tissue identifier.

This value must exist and be correctly defined in the dataset metadata;
otherwise the dataset cannot be loaded.

- `--conditions` (required)
Comma-separated list of condition labels to process
(e.g. `Primary_Tumor,Solid_Tissue_Normal`).

Each condition must exist in the metadata.

### Configuration file (`config_diff_regulation.yaml`)
The differential expression and differential regulation workflows rely on a
shared YAML configuration file, e.g.:
```bash
src/deeprbp/configs/config_diff_regulation.yaml
```

This file defines how the dataset is interpreted, not which conditions are
compared.

Key conventions:

- **Conditions are not defined in the YAML**
The biological conditions to compare are always provided at runtime via the
CLI argument `--conditions`.  
This allows the same configuration file to be reused across multiple contrasts.

- **test_path_files**
Points to the root dataset directory containing expression matrices and
metadata files.

- Metadata columns used for subsetting and comparison:
  - `sample_category`: column used to select a tissue or biological subset
  - `disease_condition`: column defining the conditions to compare

- The YAML defines the **analysis context shared across explainability and DE**,
including:

  - explainability parameters (e.g. DeepLIFT mode, reference definition)
  - score aggregation strategies
  - feature scaling and normalization via `scaler_dir`

 The scaler specified in `scaler_dir` must match the scaler used during
training of the predictor checkpoint, ensuring consistency between model
inference and downstream analyses.

## 2. Between-condition differential expression (limma-voom)
This section computes between-condition differential expression (DE) using a
single, canonical implementation based on **limma-voom**.

The analysis is designed to provide **expression-level context** for the
*differential regulation* results derived from explainability scores.

Differential expression is **computed independently** for three feature levels:

- **Genes**
- **Transcripts**
- **RNA-binding proteins (RBPs)** used in DeepRBP model

All levels are analyzed using the same statistical core to ensure full
comparability across results.

---

### Methodology overview
For each requested feature level:

1. Samples are subset by:
  - a user-defined **biological category** (e.g. tissue, cancer type)
  - two experimental **conditions** (e.g. Tumor vs Normal)

2. Estimated count matrices are:
  - aligned to filtered metadata
  - optionally restricted to RBPs (RBP-level only)

3. Differential expression is computed using:
  - voom mean–variance modeling
  - limma linear models
  - empirical Bayes variance moderation

4. Results are exported as:
  - tabular DE results
  - publication-ready volcano plots

```text
Python (dataset-aware orchestration)
│
│  de_expression_conditions.py
│  - loads config_diff_regulation.yaml
│  - filters metadata (category + conditions)
│  - subsets counts (genes / transcripts / RBPs)
│  - exports inputs for R
│
├── counts_<level>_features_x_samples.tsv
├── metadata_filtered.tsv
│
▼
R (statistical engine)
│
│  run_voom-limma_diff_reg.R
│   └── voom_limma_core.R
│       - design matrix
│       - voom transformation
│       - limma model + contrasts
│
│   └── plot_volcano.R
│       - consistent volcano plots
│
▼
Outputs (per level & contrast)
│
├── DE_<level>_<condB>_vs_<condA>.csv
├── volcano_<level>.png
└── volcano_<level>.pdf
```

All outputs are generated using a shared R core (`voom_limma_core.R`) to guarantee
consistent statistical behavior across datasets and analysis modes.

### 📁 Output structure
All results are written under the user-defined `output_dir`, following a
structured layout:

```php-template
output_dir/
└── de/
    ├── inputs/
    │   ├── metadata_filtered.tsv
    │   ├── counts_genes_features_x_samples.tsv
    │   ├── counts_transcripts_features_x_samples.tsv
    │   └── counts_rbps_features_x_samples.tsv
    │
    ├── genes/
    │   └── <conditionB>_vs_<conditionA>/
    │       ├── DE_genes_<conditionB>_vs_<conditionA>.csv
    │       ├── volcano_genes.png
    │       └── volcano_genes.pdf
    │
    ├── transcripts/
    │   └── <conditionB>_vs_<conditionA>/
    │       ├── DE_transcripts_<conditionB>_vs_<conditionA>.csv
    │       ├── volcano_transcripts.png
    │       └── volcano_transcripts.pdf
    │
    └── rbps/
        └── <conditionB>_vs_<conditionA>/
            ├── DE_genes_<conditionB>_vs_<conditionA>.csv
            ├── volcano_genes.png
            └── volcano_genes.pdf
```

#### 📄 DE result tables
Each `DE_*.csv` file contains the full limma output, including:
  - `logFC` – log2 fold change (conditionB − conditionA)
  - `AveExpr` – average expression
  - `t` – moderated t-statistic
  - `P.Value` – raw p-value
  - `adj.P.Val` – FDR-adjusted p-value
  - `B` – log-odds of differential expression
  - `Feature_ID` – gene / transcript / RBP identifier

### Volcano plots
For each contrast and feature level, volcano plots are generated automatically:

- points are colored by direction and magnitude of log2FC
- significance is defined by:
  - `|logFC| > logfc_thresh`  
  - `adj.P.Val or P.Value < p_cut_value`
- significant up/down counts are annotated on the plot
- feature labels default to **Feature_ID**, with optional gene/transcript names
if annotation is available

### Running locally
The analysis can be run locally via the installed CLI entrypoint:
```bash
run-deeprbp-differential-expression \
  --config_path src/deeprbp/configs/config_diff_regulation.yaml \
  --output_dir output/diff_reg/Liver \
  --select_category Liver_Hepatocellular_Carcinoma \
  --conditions Primary_Tumor,Solid_Tissue_Normal \
  --levels genes,transcripts,rbps \
  --logfc_thresh 1 \
  --p_cut_type fdr \
  --p_cut_value 0.05
```

This command will sequentially run DE for all requested feature levels.

### Running on a cluster (SLURM)
A ready-to-use SLURM script is provided. Submit the job with:
```bash
sbatch slurm/explainability_module/differential_regulation/run_de_expression_conditions.sh
```

### Relation to differential regulation analysis
These DE results are **not used to compute regulatory scores**, but instead serve
as **contextual information** to interpret them.

- In downstream analyses, they allow us to distinguish between:
- regulatory changes driven by expression shifts, and
regulatory rewiring occurring independently of expression level.

For this reason, the DE pipeline is kept separate, explicit, and fully reproducible.

# AHORA FALTARIAN LOS STEPS 4 Y 5. Que son de rankings



