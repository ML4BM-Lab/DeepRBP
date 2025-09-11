# DeepRBP: A novel deep neural network for inferring splicing regulation
#### Publication: https://doi.org/10.1101/2024.04.11.589004

<p align="center">
  <a href="https://www.python.org/downloads/release/python-390/">
    <img src="https://img.shields.io/badge/Python-3.9%2B-blue.svg" alt="Python Version">
  </a>
  <a href="#">
    <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg" alt="Platform">
  </a>
  <a href="#">
    <img src="https://img.shields.io/badge/GPU-Supported-brightgreen.svg" alt="GPU Support">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  </a>
</p>

## Description
<p align="center">
    <img src="images/methods_deepsf.png" width="700" alt="PDF Image">
</p>

Alternative splicing plays a pivotal role in various biological processes. In the context of cancer, aberrant splicing patterns can lead to disease progression and treatment resistance. Understanding the regulatory mechanisms underlying alternative splicing is crucial for elucidating disease mechanisms and identifying potential therapeutic targets.  
We present DeepRBP, a deep learning (DL)-based framework to identify potential RNA-binding protein (RBP)-Gene regulation pairs for further *in-vitro* validation. DeepRBP is composed of:  

1. **Prediction module:** A DL model that predicts transcript abundance given RBP and gene expression data.  
2. **Explainability module:** Computes informative RBP-Gene scores using DeepLIFT.

Why Use DeepRBP?
Experimental methods like CLIP are expensive, labor-intensive, and require prior knowledge of the target RBP. DeepRBP offers a computational alternative, predicting regulatory interactions without these constraints. This approach is cost-effective, scalable, and capable of uncovering complex multi-RBP interactions often missed by traditional methods.

DeepRBP has been validated on cancer datasets like TCGA, revealing potential novel regulatory relationships in AML, KICH, and HCC.
---

## Installation
To install **DeepRBP**, follow these steps:

```bash
git clone https://github.com/ML4BM-Lab/DeepRBP.git
cd DeepRBP
conda create -n DeepRBP python=3.9 # conda create --prefix  /data/jsanchoz/conda-env/DeepRBP python=3.9
conda activate DeepRBP # conda activate /data/jsanchoz/conda-env/DeepRBP
pip install -e .
```

# Prediction Module
## Datasets Information
In this project, we utilized multiple datasets, including samples from The Cancer Genome Atlas (TCGA) and the Genotype-Tissue Expression (GTEx) project. We used 80% of the TCGA samples to train the DeepRBP predictor, which learns to predict transcript abundances. The remaining TCGA samples, along with the GTEx samples, were employed to assess the model's generalization capabilities.

## Data Download
You can download the necessary datasets from the [UCSC Xena platform](https://xenabrowser.net/) (Goldman et al., 2020). To automate the process, execute (approx. 20 minutes):

```bash
sbatch slurm/download_data.sh
```
The following files will be downloaded and stored in the `/data/training_module/raw` directory:

- **Gene expression data** (TcgaTargetGtex_rsem_isoform_tpm.gz)  
  This file contains gene expression levels measured in Transcripts Per Million (TPM) as log2(tpm + 0.001) format. TPM is a normalization method that accounts for both the length of the gene and the total number of reads in a sample, allowing for comparison of gene expression levels across different samples. It includes samples across TCGA, GTEx and TARGET (not used).

- **Transcript expression data** (TcgaTargetGtex_rsem_gene_tpm.gz)  
  Similar to the gene expression file, this file provides expression levels for various transcript isoforms also in TPM format as log2(tpm + 0.001) format. Similar to the isoform data, it includes expression levels for genes from TCGA, GTEx and TARGET datasets. It enables analyses focused on specific isoforms of genes, which can have different functional roles and regulatory mechanisms.

- **Gene count expression data** (TcgaTargetGTEX_phenotype.txt)
  This file contains raw gene counts in log2(expected_count+1) from TCGA, GTEx, and TARGET, representing the number of reads mapped to each gene. Unlike TPM, raw counts do not account for gene length or total sequencing depth. They are often used in statistical methods for differential expression analysis, as they provide a direct measure of sequencing data.

- **Phenotype metadata** (TcgaTargetGTEX_gene_expected_count.gz) 
  This text file contains important clinical and biological information about the samples from the TCGA, GTEx, and TARGET datasets.

## Data Preprocessing
In this step, we will load and preprocess the raw data files to prepare input matrices for both TCGA and GTEX datasets. Specifically, we will generate and save in `output_dir`, in this case we suggest to save them in `data/training_module/processed`:

- **RBP expression matrix**: `RBPs_log2p_tpm.csv` (in log2(TPM+1)) derived from the gene expression data, with dimensions `n_patients x num_RBPs`.
- **Transcript expression matrix**: `trans_log2p_tpm.csv` (in log2(TPM+1)) with dimensions `n_patients x num_transcripts`.
- **Gene expression matrix**: `gn_tpm.csv` (in TPM) with dimensions `n_patients x num_genes`.
- **Metadata file**: `phenotype_metadata.csv`, containing phenotype information for each sample, indicating tissue or tumor type.
- **Gene count matrix**: `gn_counts.csv` (in counts) with dimensions `n_patients x num_RBPs`. For further differential expression analysis.

### Process Details
Among other tasks, this process includes:

1. Cleaning and standardizing phenotype data.
2. Cleaning gene and transcript expression and count data data by removing genome version annotations and aggregating loci.
3. Filtering out transcripts of genes with only one isoform.
4. Selecting genes and their transcripts for modeling based on either cancer-related genes or all protein-coding genes.
5. Filtering RNA-binding proteins (RBPs) for modeling and creating a subset RBP matrix from the gene matrix for expression.
6. Transforming gene expression to TPM (and counts), RBP and transcript expression to TPM. 
7. Transposing expression data so patients are rows and genes (or transcript IDs) are columns.
8. Saving processed expression and count data and phenotype metadata as CSV files in the specified output directory.

### Execution Command
To execute this, run:

```bash
preprocess-data --raw_data_dir "/scratch/jsanchoz/DeepRBP/data/training_module/raw" \
                --selected_genes_dir "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps" \
                --output_dir "/scratch/jsanchoz/DeepRBP/data/training_module/processed" \
                --transcript_expression_file "TcgaTargetGtex_rsem_isoform_tpm.gz" \
                --gene_expression_file "TcgaTargetGtex_rsem_gene_tpm.gz" \
                --gene_counts_file "TcgaTargetGTEX_gene_expected_count.gz" \
                --phenotype_data_file "TcgaTargetGTEX_phenotype.txt" \
                --chunk_size 1000 \
                --gene_selection True \
                --gene_transcript_mapping_file "getBM.csv" \
                --splicing_genes_file "Table_S5_Cancer_splicing_gene_eyras.xlsx" \
                --cancer_genes_file "Table_S6_Cancer_gene_eyras.xlsx" \
                --gene_census_file "Table_Cancer_Gene_Census.tsv" \
                --rbp_genes_file "Table_S2_list_RBPs_eyras.xlsx"
```

### Command Arguments
- **raw_data_dir (str)**: Directory containing raw data files.
- **selected_genes_dir (str)**: Directory with lists of RNA-binding proteins (RBPs) and selected genes for modeling.
- **output_dir (str)**: Directory for saving processed files.
- **transcript_expression_file (str)**: Filename for transcript expression data (transcripts x n_patients) in log2(tpm+0.001).
- **gene_expression_file (str)**: Filename for gene expression data (genes x n_patients) in log2(tpm+0.001).
- **gene_counts_file (str)**: Filename for gene-level expected counts data ('genes x n_patients' matrix) in log2(expected_count+1).
- **phenotype_data_file (str)**: Filename for phenotype data (patients x phenotype features).
- **chunk_size (int)**: Rows to process per chunk for memory efficiency.
- **gene_selection (bool)**: Flag to indicate gene selection; True uses cancer and alternative splicing-related genes. False uses all protein-coding genes with more than one isoform.
- **gene_transcript_mapping_file (str)**: Output file mapping transcript IDs/names to gene IDs/names and biotypes.
- **splicing_genes_file (str)**: Excel file with genes implicated in alternative splicing in cancer.
- **cancer_genes_file (str)**: Excel file listing ~900 predicted cancer-driver genes based on mutations or copy number alterations.
- **gene_census_file (str)**: TSV file with Cancer Gene Census data.
- **rbp_genes_file (str)**: Excel file listing RNA-binding proteins (RBPs).

### HPC Execution
Alternatively, you can submit this command on an HPC system with Slurm:

```bash
sbatch slurm/preprocess_data.sh
```

## Model Training from Scratch  
### Selecting Tumor Samples and Stratifying Processed Data into Training and Testing Sets
To ensure that different tumor types are equally represented in both the training and testing sets, we will perform a stratified split. This method maintains the proportion of each class in the splits, providing a more reliable evaluation of the model's generalization capabilities.

We will use the processed data from The Cancer Genome Atlas (TCGA) for this task. The training set will consist of 80% of the data, while the remaining 20% will be reserved for testing the model's generalization and explainability module. The training set will later be utilized in hyperparameter optimization with Optuna.

In this step before splitting data we can select specific tumor types (defined in config file). By default we use all the samples.

To execute this, run:

```bash
split-and-save --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_data_split.yaml" \
               --output_dir "/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets" 
                
```
#### Command Arguments
- **config_path (str)**: Path to the config file with the processed data files, sample selection, tumor types (categories), train-test fraction and source name.
- **output_dir (str)**: Directory to save the splitted datasets.

The configuration file should look like this:

#### **Example Configuration File (`config_data_split.yaml`)**

```yaml
path_files: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA"
sample_category: "detailed_category"  # The column in metadata to stratify on 
select_samples: ["all"]  # Can be ['all'] or a list of specific sample types: ['Lung_Squamous_Cell_Carcinoma', 'Rectum_Adenocarcinoma']
test_fraction: 0.2
seed: 42
```

#### HPC Execution
Alternatively, you can submit this command on an HPC system with Slurm:

```bash
sbatch slurm/split_and_save.sh
```

As a result, in the `--output_dir`, you will find two folders: **Train** and **Test**. Each folder contains the following `.csv` files:

1. **`gn_tpm.csv`**: This dataframe contains the expression levels of the genes involved in this study, measured in TPM (Transcripts Per Million).
2. **`RBPs_log2p_tpm.csv`**: This dataframe presents the expression levels of RNA-binding proteins (RBPs) in log2(TPM + 1) for the selected samples.
3. **`trans_log2p_tpm.csv`**: This file contains the expression levels of transcripts in log2(TPM + 1).
4. **`phenotype_metadata.csv`**: This file includes metadata for the samples used in the study.


### Hyperparameter Optimization for DeepRBPredictor with Optuna 
In this section, we will implement hyperparameter optimization for the DeepRBP predictor using Optuna, a hyperparameter optimization framework designed for machine learning. 
This process aims to find the best set of hyperparameters that maximize model performance. To do so we are using a subset of the training data. Exactly the 50\% with a stratified split based on tumor 
type to enhance efficiency and optimize computational resources. The training data is sourced from the directory located at `/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Train`. 

You can use the following `config` file:

#### **Example Configuration File (`config_hyper_optimization.yaml`)**
```yaml
train_path_files: '/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Train'
sample_category: "detailed_category"   
sample_fraction: 0.5
test_fraction: 0.2
getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
cuda: True
val_batch_size: 256
plot_results: False
```
where,
- **`train_path_files`**: Paths to input folder containing the data to model:
  - **`rbp_path`**: File containing RBP expression data in log2(TPM+1).  
  - **`isoform_expr_path`**: File containing isoform expression data in log2(TPM+1).  
  - **`metadata_path`**: Metadata file (sample-level information).  
  - **`gene_expr_path`**: File with gene expression data per isoform in TPM.  
  
this data it is just used to test the model using the best parameters selected by the MSE between the log2p TPM values 
between predicted and real values on a validation set.
- **`sample_category`**: Metadata column used to stratify samples.  
- **`sample_fraction`**: Portion of the training data used for the optimization process to save time and computational resources.
- **`test_fraction`**: Portion of the remained training data used for validation.
- **`getBM_path`**: File for selected gene-RBP mappings. 

- **`gene_col_name`**: The name of the column in getBM representing Gene_ID.
- **`trans_col_name`**: The name of the column in getBM representing Transcript_ID.
- **`cuda`**: Whether to use cude or not.
- **`val_batch_size`**: Batch size used for the validation data loader.
- **`plot_results`**: Enable or disable visualization of results. 


#### Step 1: Create the Optuna Study  
Before running any hyperparameter optimization, you must first initialize an Optuna study where the results will be stored. This step sets up the study with a `TPESampler` and a `MedianPruner`, and creates a local `SQLite` database to store the optimization history.

To create the study, simply run the following command from the terminal:

```bash
create-optuna-study --output_dir /scratch/jsanchoz/DeepRBP/output/results/hyperparameter_optimization_SLURM
```
where,
- **`--output_dir`: Path where the Optuna storage (optuna.db) will be saved.

#### Step 2: Submit the Hyperparameter Optimization Job via SLURM  
This step is intended to run exclusively on an HPC cluster (via SLURM) due to the computational requirements; we do not provide a command-line entry point.

**Two-stage optimization.** We split the search into two phases. **Phase 1** explored a broad space with random sampling, deliberately testing extreme configurations (e.g., zero hidden layers, alternative optimizers, large batch sizes) and discarding clearly suboptimal extremes based on validation MSE. **Phase 2** refined the search with Optuna’s `TPESampler`. In total, we ran **280 trials** (80 in Phase 1 and 200 in Phase 2), which falls within Optuna’s recommended range for effective exploration.

**Hierarchical search space.** The space is conditional: when `num_hidden_layers = 0`—and in some cases when it is 1–2—parameters such as `uniform_nodes`, `node_shrink_factor`, hidden-layer `activation_function`, and `batch_normalization` are inactive. We report all Optuna-suggested values per trial and **explicitly mark unused hyperparameters** given each configuration.

**Phase 2 sampler settings.** We enabled `multivariate=True` and `group=True` in the `TPESampler`. Briefly, `multivariate=True` samples parameters jointly to capture inter-parameter dependencies, and `group=True` decomposes the space into conditionally active groups and samples from their joint distribution—well suited to the hierarchical design. The **first 10 trials** in Phase 2 were randomly sampled as a warm-up for the sampler. Additionally, we **switched the learning rate to a categorical grid** in Phase 2 to reduce its search width and focus exploitation on the remaining hyperparameters.

**Early stopping with pruning.** We used Optuna’s `MedianPruner` to accelerate the search and reduce computation on underperforming configurations. Pruning was **disabled until at least 25 trials** had completed. Within each trial, pruning was further delayed until **470 training steps (epochs)** had been reached, after which **checks were performed every 47 steps**. In **Phase 2**, we **reduced this to 3 warm-up steps with checks once per epoch**, enabling evaluation of more configurations in less time, in line with Optuna best practices.

**SLURM execution.** We submit Phase 1 and Phase 2 as separate SLURM array jobs. 

To launch the optimization job, simply run the following command (adjust the SLURM script as needed for your HPC setup):
```bash
sbatch slurm/run_hyper_optimization_optuna.sh
```

Below is a description of the main arguments used in the hyperparameter optimization step:

* `--storage_path`: Full path to the Optuna SQLite database file where all trial information will be stored.
* `--n_trials`: Number of parameter combinations to try using Optuna. 💡 We run 250 trials per GPU (total 1000 trials), as recommended by Optuna's TPE sampler (100–1000 trials for best results).
* `--config_path`: Path to the YAML config file defining model architecture and training parameters.
* `--output_dir`: Directory to save all results, logs, and intermediate f  iles during the optimization.
* `--num_workers`: Number of worker threads for data loading. `0` is safe for most environments; increase if your system allows.
* `--min_delta`: Minimum performance improvement threshold to continue training. If the monitored metric improves less than this value, it may trigger early stopping.
* `--patience`: Number of epochs to wait without improvement before stopping training.
* `--gpu_id`: The GPU index to use for the trial batch. 🎯 This allows running one optimization job per GPU in parallel.

The hyperparameter that are going to be optimised are:
- **`num_hidden_layers`**: Number of hidden layers.
- **`hidden1_nodes`**: Number of nodes in the first hidden layer.
- **`uniform_nodes`**:  Whether to use uniform nodes across layers.
- **`node_shrink_factor`**: Factor to reduce nodes in layers.
- **`activation_func`**: Activation function to use (e.g., 'relu', 'tanh').
- **`learning_rate`**: Learning rate for the optimizer.
- **`optimizer_name`**:  Name of the optimizer (e.g., 'adamW').
- **`batch_norm_eps`**: Epsilon value for batch normalization.
- **`batch_norm_momentum`** Momentum value for batch normalization.
- **`batch_size`** Number of samples processed before the model parameters are updated.
- **`num_epochs`** Number of full passes over the training dataset during training.
 

#### Step 3: Analyze the Hyperparameter Optimization results  
Use the following command to analyze the results of your Optuna hyperparameter search and generate summary plots:

```bash
analyze-optuna-results --storage_path /scratch/jsanchoz/DeepRBP/final_results/hyperparameter_optimization/phase2_refined/hyperparameter_optimization_SLURM_gpu0_20250819_220917_definitive_edition/optuna_backups/optuna.db \
                       --output_dir /scratch/jsanchoz/DeepRBP/final_results/hyperparameter_optimization/phase2_refined/analyze_results
```

This will:
* Load the Optuna study from the SQLite database.
* Export trial results to a CSV file.
* Generate and save informative plots (e.g., optimization history, parameter importance, timeline, etc.) in both PNG and PDF formats.

### Trying alternative Machine Learning benchmark methods  
In this section, we benchmark our **DeepRBP predictor**—a deep learning-based model—against a series of traditional machine learning regressors.
We evaluate the following algorithms using a `MultiOutputRegressor` setup to predict isoform abundances: `svr`, `decision_tree`, `elastic_net`, `ridge`.

Each model predicts isoform-level abundances, which are then scaled by their corresponding gene expression (TPM) values to produce transcript TPMs. As in our Deep Learning model, we apply a log2(TPM + 1) transformation before computing metrics.

We calculate the following evaluation metrics:

- **R²**
- **Mean Squared Error (MSE)**
- **Pearson correlation**
- **Spearman correlation (Computed for both **all transcripts** and **aggregated per gene**.)**

To ensure a fair comparison with our deep learning models (which were optimized via Optuna), we reuse the **same training and validation splits** from the previous pipeline step. Final results are saved in a CSV table for all benchmark algorithms.

---

#### Configuration
To run the benchmark, define a configuration YAML file like the following:

```yaml
# src/deeprbp/configs/config_benchmark_methods.yaml
pre_split_train_path: '/scratch/jsanchoz/DeepRBP/final_results/hyp_optimization_optuna/job1/data/Train'  
pre_split_val_path: '/scratch/jsanchoz/DeepRBP/final_results/hyp_optimization_optuna/job1/data/Validation'
#sample_category: "detailed_category"  #esto en principio no necesito
getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
plot_results: False

```

#### Run via SLURM
To launch the benchmark experiments on your HPC cluster:

```bash
sbatch slurm/run_benchmark_ridge.sh
```

## Executing DeepRBP Predictor (using optimized hyperparameters)
There are three options:
* Running the Python script 
* Submitting a job to a HPC queue
* Running with Docker
---

### **Option 1: Running the Python Script**  
To execute DeepRBP on the **TCGA** dataset, use a `.yaml` configuration file. Below is an example configuration file:

#### **Example Configuration File (`config_model_train.yaml`)**

```yaml
train_path_files: '/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Train'
test_path_files: '/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Test'
sample_category: "detailed_category"  # The column in metadata to stratify on 
test_fraction: 0.2
getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
cuda: True
train_batch_size: 64
val_batch_size: 128
plot_results: True
```

Once the `config` file is ready, execute the script as follows, specifying the path `output_dir` where you want to save the results:

```bash
run-deeprbp-predictor \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_model_train.yaml" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor" \
  --epochs 2000 \
  --num_workers 4 \
  --min_delta 0.001 \
  --patience 30 \
  --save_top_k 1 \
  --verbose 1
```

### Explanation of Arguments
- **`--config_path`** (`type: str`, `required: True`): 
This argument specifies the path to the configuration file (config.yaml) that contains the settings and parameters for the DeepRBP training process. Ensure the path is correctly set to avoid errors during execution.

- **`--output_dir`** (`type: str`, `required: True`): 
This argument defines the directory where you want to save the results of the training process.  

- **`--epochs`** (`type: int`, `default: 10`):  
This argument sets the number of training epochs for the model. 

- **`--num_workers`** (`type: int`, `default: 0`):
This argument specifies the number of worker threads to use for loading the data in the DataLoader. Increasing the number of workers can speed up data loading, especially with larger datasets. By default, it is set to 0, meaning that data loading will occur in the main thread.

- **`--min_delta`** (`type: float`, `default: 0.001`): 
This argument sets the minimum change in the monitored quantity (like validation loss) that qualifies as an improvement for the purpose of early stopping. If the change is less than this value, the training will not be considered improved.

- **`--patience`** (`type: int`, `default: 30`): 
This argument determines how many epochs to wait after the last improvement before stopping the training process early. If no improvement is observed for the specified number of epochs, training will be halted to prevent overfitting.

- **`--save_top_k`** (`type: int`, `default: 1`):
This argument controls how many of the best models, according to the MSE validation, will be saved. If save_top_k is set to 0, no models will be saved. If it is set to -1, all models will be saved. Adjust this parameter based on your needs for model retention and evaluation.

### Option 2: Submit a Job in a HPC
If the number of training datasets or the total number of samples is high, we recommend submitting the job using the provided `run_predictor.sh` script from the cluster directory. 
This script is adapted to Slurm, but can be easily modified to work on SGE. 
The specific parameters should be adapted depending on the specifications of the HPC.

```bash
cd slurm
sbatch run_predictor.sh
```
### Option 3: Running with Docker
## (work to do here)
---

### Evaluate DeepRBP Predictor (a día de hoy no he hecho un notebook para esto)
### **Option 1: Running the Python Script**  
# para esto haz un jupyter notebook para que el usuario pueda usar el modelo sobre su propio data si quiere.
 



### Tumor-Specific vs General Training Benchmark
In this section, we evaluate whether training `DeepRBPredictor` on all tumor types combined improves transcript expression prediction performance compared to using a model trained on a single tumor type.

We compare two training strategies:

- *General model*: Trained on all available tumor types using the best architecture selected through Optuna, employing a consistent training-validation split for all tumor types (via `run_predictor.sh`).
- *Tumor-specific models*: Individually trained models for each tumor type using the same optimized architecture. The training-validation splits are consistent with those used for the general model.

Batch sizes for tumor-specific models are determined by the number of samples available:
- For tumor types with fewer than 100 training samples: A small batch size, such as 8, is used.
- For tumor types with 100-300 training samples: A medium batch size, such as 32, is employed.
- For tumor types with more than 300 training samples: A larger batch size, such as 64, is utilized.
<!-- # programa aun que yo diga cual es el batch size que hay que utilizar para cada uno que creo que eso no está puesto -->

For each tumor type, we evaluate isoform expression prediction accuracy using:
- The general model trained across all tumors.
- The tumor-specific model trained only on that tumor type.

This allows us to quantify the generalization capabilities of the model trained on all tumor types versus the specialized training approach.

To run this analysis:

```bash
sh /scratch/jsanchoz/DeepRBP/slurm/run_tcga_specific_vs_general_training.sh
```

This script performs the following steps:

- Trains one model per tumor type using a previously optimized architecture, with early stopping applied to ensure model convergence.
- Loads performance results from the best general model, trained on all tumor types together (via `run_predictor.sh`).
- Evaluates and compares both strategies (tumor-specific vs general) based on transcript expression prediction across all cancer types using multiple performance metrics.
- Generates publication-ready plots and summary tables:
  - The main manuscript plot includes only selected tumor types: Liver, Kidney, and AML.
  - The supplementary figure includes all tumor types.
- All plots are generated for multiple metrics (e.g., Spearman, Pearson, R², MSE...) to allow flexibility and completeness in the analysis.

The goal is to assess whether the general model provides sufficient performance across contexts, potentially offering better generalization than models trained on individual tumor types.




 




### aqui joseba (cuando todo esto esté ejecutado bien puedes borrar lo de aquí)
#TODO: 

-3)	Ángel me ha dicho una idea sobre: DeepLIFT de los RBPs que están el mismo complejo debería estar más correlado que los que no. Correlacion complejo > Correlacion no complejo, para ver si detectamos complejos y familias de rbp que se autoregulan. Lo saco de está página web:

La adición de los hidden layers puede hacer que la red aprendar interacciones y oposum ver familias.
https://mips.helmholtz-muenchen.de/corum/
descarga de aquí : https://mips.helmholtz-muenchen.de/corum/?query=(fcg_id=1)
















---

# Explainability Module
This module uses the already trained DeepRBP Predictor to compute TxRBP (transcript-by-RBP) and GxRBP (gene-by-RBP) scores using DeepLIFT (Shrikumar, Greenside, and Kundaje, 2017) [Learning important features through propagating activation differences, International Conference on Machine Learning, PMLR, pages 3145–3153].

With DeepLIFT, the contribution of each RBP-Transcript pair is determined for every sample in the input data, resulting in a three-dimensional score matrix with the following dimensions:
- Number of transcripts 
- Number of RBPs
- Number of samples.

Positive scores indicate activation of the transcript, while negative scores indicate transcript inhibition. To obtain a single score indicative of the general behavior of each RBP-Transcript pair, we collapse the scores across samples by computing the t-statistic (labeled as “t-stat”), calculated by the formula:
\[ \text{t-stat} = \frac{\text{mean}}{\left(\frac{\sigma}{\sqrt{n}}\right)} \]
where \( n \) represents the number of samples. This results in a score matrix of size: number of transcripts by number of RBPs.

Specific TCGA samples (not presented in the training process) are used to calculate the scores. This module is validated using a binary matrix indicating experimental evidence of regulation in POSTAR3 (Zhao et al., 2022) [POSTAR3: an updated platform for exploring post-transcriptional regulation coordinated by RNA-binding proteins, Nucleic Acids Research, volume 50, D1, pages D287–D294]. POSTAR3 is a comprehensive Post-Transcriptional Regulation database that provides protein binding sites on RNA obtained from CLIP experiments.

Additionally, we have applied our model in in-silico knockdown experiments.

This module aims to provide insights into how RBPs regulate gene expression. Below is an overview of the validation process and instructions to access the required data.

---

## Data Access (EJECUTA ESTO OTRA VEZ POR SI ACASO)
The necessary raw data for running this module is available through the provided Zenodo link (#TODO: AFTER YOU GET THE COMMUNITY PERMISSION UPLOAD THE LINK: https://zenodo.org/uploads/15337302) 

. Below is a description of the files:  

- **Events_Regions_gc23_400nt.RData**: Contains detailed information about the genomic regions associated with the events.  
- **EventsFound_gencode23.txt**: Provides metadata about the events, including genomic positions, event types, names, and IDs.  
- **human.txt**: Contains POSTAR3 data, including RNA-binding protein (RBP) binding sites on RNA. This file is tissue-specific and includes data from CLIP experiments as POSTAR peaks.  

To generate a tissue-specific POSTAR matrix, you can use the script `create_gene_rbp_postar_matrix.R`. This script processes the provided files and outputs a matrix customized for a specific tissue.  

**Command example:**  

```bash
module load R/4.3.2
Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R \
    --input_path /data/jsanchoz/DeepRBP/data/explainability_module/postar3 \
    --output_path /scratch/jsanchoz/DeepRBP/data/explainability_module/postar3/processed \
    --output_file_name human_liver \
    --postar_file human.txt \
    --events_regions_file Events_Regions_gc23_400nt.RData \
    --events_gencode_file EventsFound_gencode23.txt \
    --selected_tissue_cell_line HepG2,Huh7 \
    --getBM_file getBM.csv
```
*Argument details*:
- **`--input_path`**: The directory containing the necessary input files for processing. This should include the POSTAR file, events regions file, and any other relevant data files.
- **`--output_path`**: The directory where the processed output files will be saved. Ensure that this path exists or will be created by the script.
- **`--output_file_name`**: The base name of the output file that will be generated. The resulting file will be named `<output_file_name>_GxRBP.csv`.
- **`--postar_file`**: The POSTAR file containing information about RNA-binding protein (RBP) binding. This file is critical for identifying RBP interactions in the data.
- **`--events_regions_file`**: The file specifying the genomic regions associated with events of interest. This file should be in RData format and contain the necessary genomic range data.
- **`--events_gencode_file`**: A text file containing metadata about events, including identifiers (IDs) and their genomic positions. This information is used to correlate events with RBPs.
- **`--selected_tissue_cell_line`**: A comma-separated list of cell lines from the POSTAR experiments to include in the analysis. This allows for filtering the data based on specific tissues or cell types.
- **`--getBM_file`**: The name of the CSV file containing mapping information to relate gene names with their corresponding gene IDs. This file is essential for converting RBP names to gene IDs in the output matrix.

By following these steps, you can generate a POSTAR matrix tailored to your specific tissue and experimental needs. For acute myeloid leukemia (AML) use the K562 cell-line; for kidney chromophobe (KICH) use HEK293 cell-line, and for 
liver hepatocellular carcinoma (HCC) use HepG2 and Huh7 cell-lines.

```bash
Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R \
    --input_path /data/jsanchoz/DeepRBP/data/explainability_module/postar3 \
    --output_path /scratch/jsanchoz/DeepRBP/data/explainability_module/postar3/processed \
    --output_file_name human_aml \
    --postar_file human.txt \
    --events_regions_file Events_Regions_gc23_400nt.RData \
    --events_gencode_file EventsFound_gencode23.txt \
    --selected_tissue_cell_line K562 \
    --getBM_file getBM.csv
```
```bash
Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R \
    --input_path /data/jsanchoz/DeepRBP/data/explainability_module/postar3 \
    --output_path /scratch/jsanchoz/DeepRBP/data/explainability_module/postar3/processed \
    --output_file_name human_kidney \
    --postar_file human.txt \
    --events_regions_file Events_Regions_gc23_400nt.RData \
    --events_gencode_file EventsFound_gencode23.txt \
    --selected_tissue_cell_line HEK293 \
    --getBM_file getBM.csv -->
```

## Run Executing DeepRBP Explainer
To execute the explainability pipeline on our model trained with DeepLIFT, you have the following options:
* Running the Python script 
* Submitting a job to an HPC queue
* Running with Docker

---

### **Option 1: Running the Scripts**  
To execute DeepRBP Explainer with a specific tumor type on the **TCGA** test dataset, you need: An explainability config (e.g., config_model_explain_deeplift_knock_t_stat.yaml)

New (recommended): pass --select_category from the command line for one or multiple tumor types.
If --select_category is provided, it overrides any select_category in the YAML.
If it’s not provided, the explainer falls back to select_category in the YAML (if present).


#### **Example Configuration File (`config_model_explain_deeplift_knock_t_stat.yaml`)**

```yaml
test_path_files: '/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets/Test'
sample_category: "detailed_category"  # The column in metadata to stratify on
disease_condition: "sample_type" # The column in metadata to stratify on
# select_category: "Liver_Hepatocellular_Carcinoma"  # (optional if you pass it via CLI)
select_condition: 
  - "Primary_Tumor"
  #- "Solid_Tissue_Normal"
explanation_method: "DeepLIFT"
reference_data: "knockout_reference"
batch_reduction_method: "t-statistic"
gene_collapse_method: "max_absolute_value"
getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
```
where,
- **`sample_category`**: Column name in metada that refers to samples' tumor type to help filtering samples.
- **`disease_condition`**: The column in the metadata to stratify on.
- **`select_category`**: (optional in YAML): tumor type(s) to select. If you pass --select_category via CLI, this YAML key is ignored.
- **`select_condition`**: Conditions for filtering the samples, such as "Primary_Tumor" or potentially including "Solid_Tissue_Normal".
by tumor type. You could select only primary tumors or also include normal samples. You could calculate scores for 
normal samples and then run for tumor samples in a new execution to study the differences.
- **`explanation_method`**: The method you will use to compute the TxRBP explainability scores. Alternatively, you can use the `Pseudoknockdown` method, which simulates a knockdown or knockup.
- **`reference_data`**: Required only by DeepLIFT to perform calculations.
- **`batch_reduction_method`**: Specifies how we collapse the dimension of the samples once we have calculated the scores for each sample. An alternative option could be `sum_scores`.
- **`gene_collapse_method`**: Set to `"max_absolute_value"`, which is the method used to collapse the TxRBP scores matrix to GxRBP. This method takes the highest absolute value score among the transcripts of a gene and retains its sign.
- **`scaler_dir`**: Directory where the scaler is saved (in a joblib format) along with the sigma value.

##### **Run examples**
*A) Single tumor type (from CLI)*
```bash
run-deeprbp-explainer \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml" \
  --model_ckpt_path "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt" \
  --scaler_dir "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/explainer_dl_kout_t_stat_EXAMPLE" \
  --select_category "Liver_Hepatocellular_Carcinoma"
```

*B)Multiple tumor types (comma-separated list)*
```bash
run-deeprbp-explainer \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml" \
  --model_ckpt_path "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt" \
  --scaler_dir "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/explainer_dl_kout_t_stat_MULTI" \
  --select_category "Liver_Hepatocellular_Carcinoma,Acute_Myeloid_Leukemia,Kidney_Chromophobe"
```

Output layout (per category)
<output_dir>/<select_category>/df_scores_TxRBP.csv
<output_dir>/<select_category>/df_scores_GxRBP.csv
<output_dir>/<select_category>/result_table.csv


Tip: You can leave select_category commented out in the YAML (recommended), and drive categories entirely from CLI for reproducible multi-run sweeps. If you keep it in the YAML, the CLI flag still takes precedence.

### **Option 2: Submit a Job in a HPC**   
```bash
cd slurm
sbatch run_explainer.sh
```

### Option 3: Running with Docker
## (work to do here)
--- -->

As a result, you will receive three CSV files containing the scores of the RBPs at the transcript level (size: number of transcripts x RBPs) and at the gene level (size: number of genes x RBPs). Additionally, there will be a results table for each RBP-Gene (transcript) interaction with the following fields: RBP ID (attached RBP), RBP name, Gene ID (selected Gene), Gene name, Transcript ID (selected transcript of that Gene), Transcript name, Transcript biotype, Score, and the number of transcripts per gene (indicating how many transcripts the gene has).

## Evaluation of Explainability Scores Using POSTAR
For evaluating explainability scores against POSTAR experimental data, the `run-postar-validator`command facilitates this process by integrating explainability scores with POSTAR data, allowing for validation and analysis of the results.

To run the POSTAR validator, use the following command:

```bash  
run-postar-validator \
  --postar_matrix_dir "/scratch/jsanchoz/DeepRBP/data/explainability_module/postar3/processed" \
  --postar_file "human_liver_GxRBP.csv" \
  --scores_result_dir "/scratch/jsanchoz/DeepRBP/output/results/explainability_deeplift_knock_t_stat/results" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/explainability_deeplift_knock_t_stat/results"
```

### Command-Line Arguments
The following command-line arguments are required to execute the script:

- `--postar_matrix_dir` (str): Directory path to the POSTAR matrix. This directory should contain the POSTAR data files needed for validation.
- `--postar_file` (str): Filename of the POSTAR binary matrix. This file should be in the format of genes x RNA Binding Proteins (RBPs).
- `--scores_result_dir` (str): Directory path for the scores matrices and results table. This directory should contain the explainability scores that you want to evaluate.
- `--output_dir` (str): Directory path where validation results will be saved. The results of the validation will be stored in this directory.
- `--verbose` (int, default=1): Verbosity level for logging messages. Adjust this parameter to control the amount of information logged during the execution (default is 1 for minimal logging).

### Results Visualization (locally)
After running the DeepRBP Explainer, you can visualize the results to compare scores derived from explainability techniques, such as DeepLIFT, against POSTAR experimental data.

Among other things, this script performs a Wilcoxon test at the gene level and a test at the RBP level to see if the difference between the POSTAR values of 0 ("Not Binding") and 1 ("Binding") is statistically significant. Then it plots the results.

Use the following command to generate the visualization:

```bash
module load R/4.3.2
Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/postar_validation/run_postar_plot_generation.R \
  --input_path "/scratch/jsanchoz/DeepRBP/output/results/explainability_deeplift_knock_t_stat/results/postar_validation" \
  --results_filename "result_table_completed.csv" \
  --count_genes_per_rbp_file "count_genes_per_rbp.csv" \
  --count_rbps_per_gen_file "count_rbps_per_gen.csv" \
  --getBM_path "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv" \
  --getBM_filename "getBM.csv" \
  --output_path "/scratch/jsanchoz/DeepRBP/output/results/explainability_deeplift_knock_t_stat/results/postar_validation" \
  --output_filename "plot_score_results.pdf" \
  --save_plot TRUE \
  --index_start 1 \
  --index_end 4 \
  --max_iterations 7
```

This command generates box plots illustrating the distribution of scores across different RNA Binding Proteins (RBPs) and genes, utilizing POSTAR labels for classification. Additionally, the plot displays the results of conducting a Wilcoxon test across RBPs or genes between the groups 0 and 1 of POSTAR, allowing you to assess whether the difference in medians is statistically significant. You can see the results obtained in the Wilcoxon test in `stat_test_rbps.csv`and `stat_test_genes.csv`.

- **Binding**: Indicates that the gene is regulated by the RBP (Postar_Score = 1).
- **Not Binding**: Indicates that the gene is not regulated by the RBP (Postar_Score = 0).
- **Unknown**: Indicates that there is insufficient data to determine the binding status (Postar_Score = N/A).

## Parameters
The parameters, in order, are as follows:
- **`input_path`**:  
  A character string specifying the path to the input files.
- **`results_filename`**:  
  A character string indicating the name of the results file (CSV) that contains the calculated scores.
- **`count_genes_per_rbp_file`**:  
  A character string indicating the name of the file (CSV) that contains the number of genes per RBP in Postar ordered by positive rbp 
  gene interactions in decreasing order
- **`count_rbps_per_gen_file`**:  
  A character string indicating the name of the file (CSV) that contains the number of rbps per Gene in Postar ordered by positive rbp
  gene interactions in decreasing order
- **`getBM_path`**:  
  A character string indicating the path of the getBM file (CSV) that contains gene ID mappings.
- **`getBM_filename`**:  
  A character string indicating the name of the getBM file (CSV) that contains gene ID mappings.
- **`output_path`**:  
  A character string specifying the path where the output files will be saved.
- **`output_filename`**:  
  A character string indicating the name of the output file (including the extension).
- **`save_plot`**:  
  A logical value indicating whether to save the plot as a PDF. Default is `TRUE`.
- **`index_start`**:  
  An integer specifying the starting index for slicing the RBP and gene lists for plotting. Default is `1`.
- **`index_end`**:  
  An integer specifying the ending index for slicing the RBP and gene lists for plotting. Default is `4`.
- **`max_iterations`**:  
  An integer specifying the maximum number of iterations for processing RBPs and genes for plotting. Default is `5`.


## Complex-level Correlation Analysis
In this section we perform a exploratory analysis of RBP explainability scores through the lens of protein complexes to investigate whether RBPs that form part of the same protein complex (e.g., CORUM) exhibit higher correlation in their explainability score patterns across genes compared to unrelated RBPs.

A correlation matrix of shape n_RBPs × n_RBPs is computed from the DeepRBP score matrix (n_RBPs × n_genes), measuring similarity between RBP profiles. Then we:

- reorder the correlation matrix to group complex RBPs together.
- generate a minimal pheatmap-style heatmap, focusing on the upper triangle to simplify large complexes.
- perform a Wilcoxon rank-sum test to compare within-complex vs outside-complex correlation values.
- Create a boxplot contrasting correlation distributions (within vs outside).

Finally we collect all p-values and summarize them using Stouffer’s method.

To run this module:

```bash
sh /scratch/jsanchoz/DeepRBP/slurm/run_corum_complex_analysis.sh
```
### NMF-based Complex Detection Analysis
In this section, we apply Non-negative Matrix Factorization (NMF) to discover latent protein complexes from gene × RBP explainability score matrices. Starting from the number of known CORUM complexes, we vary the number of components (complexes) and evaluate reconstruction error to find a meaningful decomposition.

For each tested number of components, the script:
- factorizes the gene × RBP score matrix into component loadings (genes) and feature loadings (RBPs),
- saves heatmaps of the feature loading matrix (components vs RBPs),
- plots the reconstruction error trend across tested component numbers.

To run:

```bash
sh /scratch/jsanchoz/DeepRBP/slurm/run_nmf_complex_analysis.sh
```


## Evaluation of Explainability Scores Using real RBP knockdown data
This section describes how to validate explainability scores using real RNA-binding protein (RBP) knockdown RNA-seq data.

*Environment Setup*
First, create a conda environment with Kallisto installed:

```bash
conda create -n kallisto_env -c bioconda kallisto
conda activate kallisto_env
```

### Dataset Download and Quantification
To evaluate explainability using real knockdown data, we need RNA-seq quantification data at the transcript level. For this, we use publicly available RBP knockdown datasets, which we process using Kallisto, a fast and accurate tool for pseudoalignment-based quantification.

For each dataset, we provide a dedicated script that:

1. Downloads paired-end RNA-seq FASTQ files from the corresponding GEO accession.
2. Downloads the GENCODE transcriptome (e.g., v23) if not already present.
3. Builds a Kallisto index (only once per transcriptome version).
4. Runs Kallisto quantification with 100 bootstraps per sample.
5. Saves the results in a structured output folder, where each sample has its own subdirectory with:
  - **abundance.tsv**: transcript-level quantification
  - **abundance.h5**: binary HDF5 version of the same data

**Example: GSE136366**
To process the dataset GSE136366 (TDP-43 knockdown), run the following script:

```bash
sbatch /scratch/jsanchoz/DeepRBP/slurm/download_GSE136366_kd_data.sh
```
This will download and quantify all samples, and save the output in:

`/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE136366/kallisto_output/`
 
Each sample will have its own folder inside this directory, containing both `abundance.tsv` and `abundance.h5`.

### Differential Expression with voom-limma  
After quantifying all samples with Kallisto, use the following command to run differential expression analysis between 
the two experimental conditions: control and knockdown:

```bash
module load R/4.3.2
Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/real_knockdowns/run_voom-limma.R \
  --path=/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE136366 \
  --output_dir=/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns/GSE136366 \
  --path_getBM=/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv \
  --condition_control=Rescued_tdp43 \
  --condition_kd=tdp43_ko
```

*Input arguments*
- **--path**: Path to the directory containing the Kallisto quantification output and the info_samples.txt file.
- **--output_dir**: Directory where the vool-limma results will be saved.
- **--path_getBM**: Path to the CSV file with transcript-to-gene mapping. Used to map transcripts to genes during aggregation.
- **--condition_control**: Label that identifies the control samples in info_samples.txt. Must match the second column (V2) of that file.
- **--condition_kd**: Label that identifies the knockdown samples in info_samples.txt.
Must also match what's in column V2.

The script will generate a new folder inside the dataset directory: `voom-limma_output/`

This folder contains:

- **DE_genes_<condition_kd>_vs_<condition_control>.csv**: Differential expression results at the gene level
- **DE_transcripts_<condition_kd>_vs_<condition_control>.csv**: Differential expression results at the transcript level
Each file includes statistics such as log-fold change, adjusted p-values, and B-statistics from the voom-limma analysis.

In the dataset folder, make sure you have the following file:
- **info_samples.txt**: a tab-separated file containing metadata for each sample. In the 'run_accession' column, you can read the sample name (for example SRR10045016), and if you go to the folder with the same name in path_data, you will find the data containing abundance and ``expression data in TPM``.

### Data Preprocessing: Read and process raw samples to generate the input data for DeepRBP 
Once the Kallisto quantification step is completed and `abundance.tsv` files are available for all samples, the next step is to preprocess the expression data to generate the model input matrices required by DeepRBP.
For this step we need to go back to the `/data/jsanchoz/conda-env/DeepRBP` environment.

This step includes:

1. Merging transcript and gene-level expression matrices (TPM) across samples.
2. Filtering and selecting relevant genes and transcripts.
3. Extracting RNA-binding protein (RBP) expression matrices.
4. Performing a log2(TPM + 1) transformation.
5. Transposing matrices so that samples are rows and features are columns.
6. Saving the process ed matrices to disk.

*Running the preprocessing script*
To run the preprocessing pipeline on the knockdown dataset `GSE136366`, execute:

```bash
preprocess-realkd-data \
    --path_dataset '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE136366' \
    --output_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE136366/processed' \
    --condition_control 'Rescued_tdp43' \
    --condition_knockdown 'tdp43_ko' 
```
<!-- preprocess-realkd-data \
    --path_dataset '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE75491' \
    --output_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE75491/processed' \
    --condition_control 'Control' \
    --condition_knockdown 'RBM47 KD' -->

<!-- preprocess-realkd-data \
    --path_dataset '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/PRJEB39343' \
    --output_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/PRJEB39343/processed_MBNL1' \
    --condition_control 'HFE145siNeg_1' \
    --condition_knockdown 'HFE145siMBNL1' -->

<!-- FUS, TAF15, TARDBP -->

<!-- preprocess-realkd-data \ # kd1
    --path_dataset '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE77702' \
    --output_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE77702/processed_FUS' \
    --condition_control 'Control' \
    --condition_knockdown 'FUS KD' -->

<!-- preprocess-realkd-data \ # kd2
    --path_dataset '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE77702' \
    --output_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE77702/processed_TAF15' \
    --condition_control 'Control' \
    --condition_knockdown 'TAF15 KD' -->

<!-- preprocess-realkd-data \ # kd3
    --path_dataset '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE77702' \
    --output_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE77702/processed_TARDBP' \
    --condition_control 'Control' \
    --condition_knockdown 'TARDBP KD' -->

This will:
- Read the sample metadata from `info_samples.txt` inside the dataset folder.
- Load the transcript-to-gene mapping file `getBM.csv` from the selected_genes_dir.
- Identify the sample IDs for each condition.
- Process and save the following matrices for both `control` and `knockdown` conditions:
  - **RBPs_log2p_tpm.csv**: RBP expression matrix
  - **trans_log2p_tpm.csv**: transcript expression matrix
  - **gn_tpm.csv**: gene expression matrix


In this case all outputs are saved in : `/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/raw/GSE136366/processed/<condition>/`

where <condition> is either `control` or `knockdown`.
 
*Input Arguments*
- `--path_dataset`: Path to the dataset root (must contain `kallisto_output/` and `info_samples.txt`).
- `--output_dir`: Directory where the processed matrices will be saved.
- `--condition_control`: Label used in `info_samples.txt` to identify `control` samples.
- `--condition_knockdown`: Label used to identify `knockdown` samples.

**📝 Notes**
- The `info_samples.txt` file must contain a column with sample IDs (`run_accession`) and a column with condition labels (e.g., `tdp43_ko`, `Rescued_tdp43`).
- This script uses **TPM-normalized transcript and gene expression values**.
- The processing pipeline uses **selected gene sets related to cancer and splicing** if `--gene_selection` is set to `True` (default).
- This step is required **before running explainability analyses or training a model** on real knockdown data.


### Executing DeepRBP Evaluation and Explainability on Knockdown Data
This section outlines the process to execute the main script for evaluating the DeepRBP model and calculating explainability scores using real RBP knockdown data.

The script leverages a configuration file to define key parameters such as the path to the trained scaler, the model checkpoint, the explainability method, and other necessary settings.
It processes input data and labels for both control and knockdown conditions.

**Configuration File**
The configuration file (`config_real_knockdowns.yaml`) includes:
- **scaler_dir**: Path where the trained scaler is stored. This scaler is used to normalize input data to the same scale as during the model's training phase.
- **model_checkpoint_path**: Specifies the location of the trained model's checkpoint file.
- **explanation_method**: Defines the method used for explainability, such as "DeepLIFT".
- **reference_data**: Reference data type for the explainability analysis.
- **batch_reduction_method**: Method used to reduce batch-level explainability results.
- **gene_collapse_method**: Approach for aggregating transcript-level data to gene-level.
- **getBM_path**: Path to the CSV file for transcript-to-gene mapping, essential for data aggregation.
- **gene_col_name** and **trans_col_name**: Columns in the mapping file that identify gene and transcript IDs, respectively.

*Example Configuration File*
Below is an example of how the `config_real_knockdowns.yaml` might be structured:

```YAML
scaler_dir: "/scratch/jsanchoz/DeepRBP/output/results/stuff"
model_checkpoint_path: "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=09-validation_loss=4.40.ckpt"
explanation_method: "DeepLIFT"
reference_data: "knockout_reference"
batch_reduction_method: "t-statistic"
gene_collapse_method: "max_absolute_value"
getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
```

#### Running the Script*
To execute the script on the dataset `GSE136366`, use the following command:
```bash
run-deeprbp-realkd \
    --config_path '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_real_knockdowns.yaml' \
    --processed_data_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE136366/processed' \
    --output_dir '/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns/GSE136366'
```

*Running via SLURM*
If you prefer to execute the script using SLURM for job scheduling, follow these steps:

```bash
cd slurm
sbatch run_realkds.sh
```

Steps Performed by the Script:
- **Load and Scale Data**: The script loads processed data from `processed_data_dir` and scales it using the pre-trained scaler defined in the configuration file.
- **Model Evaluation**: It evaluates the DeepRBP trained predictor on both `control` and `knockdown` datasets, making predictions and calculating performance metrics in `log2(TPM + 1)` format.
- **Explainability Analysis**: The script computes explainability scores specifically for the `control` dataset, utilizing the method specified in the configuration file (e.g., "DeepLIFT").
- **Save Results**: Finally, it saves all evaluation metrics and explainability scores to the specified `output_dir`.


### Comparing RBP Explainability Scores Between DE and Non-DE Genes/Transcripts
The script `generate_realkd_plot.R`  generates visualizations of explainability scores for RNA-binding proteins (RBPs) across multiple real knockdown experiments. It compares scores between differentially expressed (DE) and non-differentially expressed genes/transcripts, using statistical testing.

*Pipeline Overview*
1. **Load scores from real knockdown results:** Transcript-level and gene-level explainability scores and DE results from voom-limma.
2. **Filter significant features:** Keep only those with `adj.P.Val < 0.05` from DE analysis
3. **Prepare melted dataframes:** Annotate DE vs non-DE and Extract scores for each RBP using getBM mapping.
4. **Statistical testing:** Mann-Whitney-Wilcoxon test between DE and non-DE with bonferroni adj.
5. **Generate plots:** Boxplots with jittered points, Annotated with significance levels, Legends and titles formatted
6. **Arrange and export:** All transcript plots combined into a single PDF, All gene plots combined into another PDF

*How to Run*
```bash
cd /scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/real_knockdowns # change as needed
Rscript generate_realkd_plot.R \
  --output_dir "/Users/joseba/Desktop/real_knockdowns/figures" \
  --path_getBM "/Users/joseba/Desktop/data:real_kds/getBM.csv" \
  --path_results "/Users/joseba/Desktop/real_knockdowns" \
  --experiment_list "PRJEB39343,GSE75491,GSE77702-FUS,GSE77702-TAF15,GSE77702-TARDBP,GSE136366" \
  --rbp_interest_list "MBNL1,RBM47,FUS,TAF15,TARDBP,TARDBP"
```
*Input Arguments*
- **output_dir**: Directory to save output plots (PDFs)
- **path_getBM**:	Path to CSV file mapping Gene_ID to Gene_name (used to match RBPs)
- **path_results**:	Path to folder containing real knockdown results (with subfolders for each experiment)
- **experiment_list**: Comma-separated list of experiment names
- **rbp_interest_list**: Comma-separated list of RBP names (must match getBM$Gene_name) corresponding to experiments
⚠️ The number of experiments and RBPs must match one-to-one.

*Output*
Two arranged PDF figures saved to `output_dir`: 
- **transcripts_arranged.pdf**: Boxplots of scores for DE vs non-DE transcripts
- **genes_arranged.pdf**: Boxplots of scores for genes with DE transcripts vs others

Each figure contains multiple subplots (one per experiment), with statistical comparisons (Wilcoxon test) annotated.





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# asi se veía mi modelo viejo: 
<!-- DeepRBP(
  (linear0): Linear(in_features=1282, out_features=1024, bias=True)
  (bn0): BatchNorm1d(1024, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
  (a_function0): ReLU()
  (linear1): Linear(in_features=1024, out_features=128, bias=True)
  (bn1): BatchNorm1d(128, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
  (a_function1): ReLU()
  (final_layer): Linear(in_features=128, out_features=11462, bias=True)
) -->



# JOSEBA CUANDO ACABES TODAS LAS EJECUCIONES Y EL MANUSCRITO!!!!! 
# CUANDO ACABE TODO EL PROYECTO VAMOS A CREAR UN README GENERAL QUE EXPLIQUE TODO EL PROYECTO.
# UN README PARA EL TRAINING (CON APARTADOS EXTRA DE ENTRENAR FROM SCRATCH POSIBLEMENTE EN OTRO README)
# UN README PARA EL EXPLAINABILITY (CON DIVISION POR SUBSECCIONES)


<!-- # nuevo organigrama !!! (puede estar aun sujeto a muchos cambios) ACTUALIZA ESTO BROTHER!!!
/DeepRBP
├── data (esto hay que actualizar)
│   ├── training_module                       
│   │   ├── raw                                # Datos crudos descargados de TCGA y GTEx
│   │   │   ├── TcgaTargetGtex_rsem_isoform_tpm.gz   # Datos de transcritos en log2(tpm+0.001) (TCGA y GTEx)
│   │   │   ├── TcgaTargetGtex_rsem_gene_tpm.gz      # Datos de genes en log2(tpm+0.001) (TCGA y GTEx)
│   │   │   ├── TcgaTargetGTEX_phenotype.txt         # Datos de fenotipo de TCGA, GTEx y TARGET
│   │   │   ├── TcgaTargetGTEX_gene_expected_count.gz  # Datos de genes cuentas raw de TCGA, GTEx y TARGET
│   │   │
│   │   ├── selected_genes_rbps                # Listas seleccionadas de genes y RBPs relevantes
│   │   │   ├── Table_Cancer_Gene_Census.tsv        # Tabla con el censo de genes de cáncer
│   │   │   ├── Table_S2_list_RBPs_eyras.xlsx      # Tabla S2 con lista de RBPs (Eyras)
│   │   │   ├── Table_S5_Cancer_splicing_gene_eyras.xlsx   # Tabla S5 con genes de splicing en cáncer (Eyras)
│   │   │   ├── Table_S6_Cancer_gene_eyras.xlsx            # Tabla S6 con genes de cáncer (Eyras)
│   │   │   ├── getBM.csv                          # Relaciona genes id con su trans id correspondiente
│   │   │
│   │   ├── processed   # Datos procesados del raw data listos para ser usados por el modelo (después de ser transformados)
│   │   │   ├── TCGA / GTEx                      # Datos de TCGA o GTEX(sin normalizar, escalar ni dividir en Train/Test)
│   │   │   │   ├── RBPs_tpm.csv           # Expresión de RBPs en datos de TCGA
│   │   │   │   ├── RBPs_log2p_counts.csv  # Counts de RBPs en datos de TCGA
│   │   │   │   ├── gn_tpm         # Expresión de genes en datos de TCGA
│   │   │   │   ├── trans_tpm.csv          # Expresión de transcritos en datos de TCGA
│   │   │   │   ├── phenotype_metadata.csv       # Datos de fenotipo de TCGA o GTEX
│   │   │   │
│
│   ├── explainability_module                   # Módulo dedicado a la validación y explicación del modelo
│   │   ├── postar3                             # Datos de POSTAR3 para validación
│   │   │   ├── raw
│   │   │   │   ├── Events_Regions_gc23_400nt.RData # las regiones de los eventos
│   │   │   │   ├── EventsFound_gencode23.txt # info de los eventos: posición, tipo de evento, nombre, id, etc...
│   │   │   │   ├── human.txt: Postar3 of the selected tissue: information of the RBPs attatch in genome.
                        (/data/jsanchoz/DeepRBP/data/explainability_module/postar3 - los primeros 3 docs)
│   │   │   │   ├── human_postar3_cell_line_info.csv

                # meter aquí los distintos human.txt que se creen
│   │   │   │ 
│   │   │   ├── processed
│   │   │   │   ├── 
│   │   └── real_kds                            # Datos de experimentos de knockdown (KD)
│ 
├── output/                           # Carpeta para almacenar resultados y modelos entrenados
│   ├── checkpoints/                  # Almacena checkpoints del modelo durante el entrenamiento
│   ├── results/                      # Almacena los resultados de las evaluaciones y análisis
│   │   ├── analysis/                 # Sección dedicada a análisis y resultados
│   │   │   ├── TCGA-Lung-Breast-2024-10-09/   # Identificador único para esta corrida
│   │   │   │   ├── train_prediction_model/    # Resultados del módulo de entrenamiento

│   │   │   │   │ data/  
│   │   │   │   │   ├── scaler_trained/  
│   │   │   │   │   │   │   ├── scaler.joblib 
│   │   │   │   │   │   │   ├── sigma.npy # Valor de desviación estándar utilizado en el escalado

│   │   │   │   │   ├── train_data/     
│   │   │   │   │   │   │   ├── rbp_expr_df.csv         # Expresión de RBPs en log2(tpm+1)
│   │   │   │   │   │   │   ├── scaled_rbp_expr_df.csv  # Expresión de RBPs scaled 0-1
│   │   │   │   │   │   │   ├── gene_expr_df.csv        # Expresión de genes en tpm
│   │   │   │   │   │   │   ├── trans_expr_df.csv       # Expresión de transcritos en log2(tpm+1)
│   │   │   │   │   │   │   └── metadata_df.csv         # Metadata de los samples, indicando el set
│   │   │   │   │   │   │
│   │   │   │   │   ├── valid_data/   
│   │   │   │   │   ├── test_data/   
│   └── logs/                # Carpeta general para logs de todo el proyecto

<!-- ├── notebooks  # Notebooks de análisis y pruebas
│   ├── Tutorial_predict_transcript_expression.ipynb  # Tutorial para predecir expresión de transcriptos
│   ├── Tutorial_replicate_postar3.ipynb  # Tutorial para replicar los resultados en POSTAR3
│   └── Tutorial_replicate_real_kds.ipynb  # Tutorial para replicar knockdown experiments -->

<!-- src/  # Main code for the DeepRBP package
│
├── deeprbp/
│   ├── __init__.py                      # Initialize the DeepRBP package
│   │
│   ├── training_module/ 
│   │   ├── main_predictor.py              # Main function to execute the predictor training pipeline 
│   │   ├── pipeline.py                    # Class DeepRBPredictorPipeline
│   │   ├── train_model.py                 # Class TrainPredictor for training the model
│   │   ├── model.py                       # Defines the prediction model class (PredictorModel)
│   │   ├── plots.py                       # Functions for plotting results
│   │   ├── evaluation.py            # Functions to evaluate the performance of the predictive model
│   │
│   │── explainability_module/ 
│   │   ├── main_explainer.py              # Main function to execute the explainer postar pipeline 
│   │   ├── postar_pipeline.py   # Class DeepRBPostarExplainabilityPipeline
│   │   ├── postar_validator.py  # Class PostarValidator to validate results against POSTAR experimental data.
│   │   ├── deeplift_handler.py            # Class DeepLiftHandler for DeepLift calculations at transcript and gene levels.
│   │   ├── pseudokd_handler.py            # (In Progress)
│   │   ├── model.py             # Defines the explanatory model class (ExplainerModel)
│   │   ├── results_visualization/  # Visualizacion de los resultados de postar/real kds
│   │   │   │   └── generate_postar_plots.R  # function
│   │   │   │   └── run_postar_plot_generation.R  # main
│   │   │   │   └── plots.py                       # Functions for plotting results
│   │ 
│   │── pretrained_model/  # Contiene el modelo preentrenado y sus archivos asociados 
│   │   ├── config.json  # Configuración del modelo preentrenado
│   │   ├── model.pt  # Modelo preentrenado
│   │   ├── scaler_sfs.joblib  # Escaladores usados en el preprocesamiento
│   │       └── sigma_sfs.txt  # Parámetros adicionales del modelo

│   ├── configs/  # Configuraciones
│   │   └── config_tcga_train.yaml 
│   │   └── config_gtex.yaml 
│   │   └── config_tcga_explain.yaml 
│   │   └── config_deg.yaml 
│   │ 
│   ├── data_loading/  
│   │   └── config_loader.py                  # Config class and load_config to load configurations from YAML
│   │   └── data_loader.py  # Classes for data loading, filtering, splitting, etc.
│   │
│   ├── util/  
│   │   └── utils.py               # Utility functions
│   │   └── logger.py              # Class for logging print, warning, and errors

│   ├── data_preprocessing/  # Raw data preprocessing
│   │   └── prep_model_inputs.py  # Prepares TCGA/GTEx data for model input
│   │   └── create_gene_rbp_postar_matrix.R # Creates the GxRBP matrix for specific tissues
│   │
│   ├── differential_expression_analysis/   
│   │   └── deg_analysis.py    
│   │ 
│   └── tests/  # Tests unitarios para el paquete DeepRBP
│       └── test_data_loader.py  # Test unitario para la clase DataLoader (por ejemplo)

├── slurm  # Scripts ejecutables
│   ├── download_data.sh  # Script para descargar datos (TCGA, GTEx)
│   ├── generate_model_inputs.sh  # Script para procesar los datos descargados y generar matrices de input
│   ├── run_predictor_pipeline.sh  # Script para entrenar y evaluar el predictor DeepRBP
│   └── run_explainer_postar_pipeline.sh  # Script para ejecutar el módulo de explainability

├── images  # Imágenes para visualización (por ejemplo, diagramas o ejemplos de resultados.
├── README.md  # Instrucciones y documentación del proyecto
├── .gitignore  # Archivos y carpetas a ignorar en el control de versiones
└── setup.py  # Script de instalación para el paquete DeepRBP


# ME QUEDA LUEGO HACER UN GET_POTENTIAL_CANDIDATES A PARTIR DE UN DF_SUMMARY. Y RESULTADOS DE DEG.
# ME QUEDA PODER USAR DEEPRBPEXPLAINER PARA REAL KDS DATA (EL GET INPUT DATA HAY Q HACERLO TB PARA ESTOS)
# ME QUEDA QUE EL PREP_MODEL_INPUTS COJA LOS COUNTS DE LOS TRANS Y TODOS LOS GENES.
# hacer el pseucode.py y que el deeplift_handler pueda trabajar con todos los casos.
# hacer merge de este branch en git y publicar la versión!
# hacer notebooks! --> -->
