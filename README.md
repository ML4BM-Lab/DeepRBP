
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
In this project, we have used several datasets, including TCGA and GTEx samples. The TCGA samples are used to train the DeepRBP predictor that learns transcript abundances, while GTEx samples are used to evaluate generalization.

## Data Download
You can download the necessary datasets from the [UCSC Xena platform](https://xenabrowser.net/) (Goldman et al., 2020). To automate the process, execute (approx. 20 minutes):

```bash
sbatch slurm/download_data.sh
```
The datasets will be automatically saved to the following directory: `/data/training_module/raw`

## Data Preprocessing
In this step, we will load and preprocess the raw data files: gene expression (`TcgaTargetGtex_rsem_gene_tpm.gz`), transcript expression (`TcgaTargetGtex_rsem_isoform_tpm.gz`), and phenotype metadata (`TcgaTargetGTEX_phenotype.txt`). These will be used to prepare input matrices for both TCGA and GTEX datasets. Specifically, we will generate:

- **RBP expression matrix**: `RBPs_log2p_tpm.csv` (log2(tpm+1)), derived from the gene expression data, with dimensions `n_patients x num_RBPs`.
- **Transcript expression matrix**: `trans_log2p_tpm.csv` (log2(tpm+1)), with dimensions `n_patients x num_transcripts`.
- **Gene expression matrix**: `gn_expr_each_iso_tpm.csv` (in TPM) with dimensions `n_patients x num_transcripts`.
- **Metadata file**: `phenotype_metadata.csv`, containing phenotype information for each sample, indicating tissue or tumor type.

### Process Details
Among other tasks, this process includes:

1. Cleaning and standardizing phenotype data.
2. Cleaning gene and transcript expression data by removing genome version annotations and aggregating loci.
3. Filtering out transcripts of genes with only one isoform.
4. Selecting genes and their transcripts for modeling based on either cancer-related genes or all protein-coding genes.
5. Filtering RNA-binding proteins (RBPs) for modeling and creating a subset RBP matrix from the gene matrix.
6. Transforming gene expression to TPM, and RBP and transcript expression to log2(tpm+1).
7. Transposing expression data so patients are rows and genes (or transcript IDs) are columns.
8. Saving processed expression data and phenotype metadata as CSV files in the specified output directory.

## Execution Command
To execute this, run:

```bash
prepare-model-inputs --raw_data_dir "/scratch/jsanchoz/DeepRBP/data/training_module/raw" \
                     --selected_genes_dir "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps" \
                     --output_dir "/scratch/jsanchoz/DeepRBP/data/training_module/processed" \
                     --transcript_expression_file "TcgaTargetGtex_rsem_isoform_tpm.gz" \
                     --gene_expression_file "TcgaTargetGtex_rsem_gene_tpm.gz" \
                     --phenotype_data_file "TcgaTargetGTEX_phenotype.txt" \
                     --chunk_size 1000 \
                     --gene_selection True \
                     --gene_transcript_mapping_file "getBM.csv" \
                     --splicing_genes_file "Table_S5_Cancer_splicing_gene_eyras.xlsx" \
                     --cancer_genes_file "Table_S6_Cancer_gene_eyras.xlsx" \
                     --gene_census_file "Table_Cancer_Gene_Census.tsv" \
                     --rbp_genes_file "Table_S2_list_RBPs_eyras.xlsx"
```

# Command Arguments
- **raw_data_dir (str)**: Directory containing raw data files.
- **selected_genes_dir (str)**: Directory with lists of RNA-binding proteins (RBPs) and selected genes for modeling.
- **output_dir (str)**: Directory for saving processed files.
- **transcript_expression_file (str)**: Filename for transcript expression data (transcripts x n_patients) in log2(tpm+0.001).
- **gene_expression_file (str)**: Filename for gene expression data (genes x n_patients) in log2(tpm+0.001).
- **phenotype_data_file (str)**: Filename for phenotype data (patients x phenotype features).
- **chunk_size (int)**: Rows to process per chunk for memory efficiency.
- **gene_selection (bool)**: Flag to indicate gene selection; True uses cancer and alternative splicing-related genes.
- **gene_transcript_mapping_file (str)**: Output file mapping transcript IDs/names to gene IDs/names and biotypes.
- **splicing_genes_file (str)**: Excel file with genes implicated in alternative splicing in cancer.
- **cancer_genes_file (str)**: Excel file listing ~900 predicted cancer-driver genes based on mutations or copy number alterations.
- **gene_census_file (str)**: TSV file with Cancer Gene Census data.
- **rbp_genes_file (str)**: Excel file listing RNA-binding proteins (RBPs).

## HPC Execution
Alternatively, you can submit this command on an HPC system with Slurm:

```bash
sbatch slurm/generate_model_inputs.sh
```

# Darle una vuelta a estos dos comentarios:
# (explicar más en detalle), que entra que sale, como , porque  decir como se consiguen los RBPs, que procesos hacemos etc

# filter TCGA and GTEx samples and save the gene and transcript expression matrices for each dataset. From the gene expression matrix, we create a subset containing the expression of RNA-binding proteins (RBPs) genes, which will serve as the primary input for our model.

## Executing DeepRBP Predictor
There are three options:
* Running the Python script 
* Submitting a job to a HPC queue
* Running with Docker

---

### **Option 1: Running the Python Script**  
To execute DeepRBP on the **TCGA** dataset, use a `.yaml` configuration file. Below is an example configuration file:

#### **Example Configuration File (`config.yaml`)**

```yaml
source_name: "TCGA"

# Paths for the data files
data_paths:
  rbp_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/RBPs_log2p_tpm.csv"
  isoform_expr_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/trans_log2p_tpm.csv"
  metadata_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/phenotype_metadata.csv"
  gene_expr_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/gn_expr_each_iso_tpm.csv"
  getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"

# Model predictor configuration
model:
  input_size: 1348               # Number of input features
  output_size: 11459             # Number of output isoforms
  num_hidden_layers: 2           # Number of hidden layers
  max_node: 1024                 # Maximum nodes per layer
  uniform_nodes: true            # Uniform node distribution across layers
  node_shrink_factor: 2          # Factor to reduce nodes per layer
  activation_func: "relu"        # Activation function
  learning_rate: 0.001           # Optimizer learning rate
  optimizer_name: "adamW"        # Optimizer to use
  cuda: true                     # Set to true for GPU acceleration

# Training configuration
training:
  epochs: 1000                   # Total number of training epochs
  batch_size: 128                # Batch size for training
  print_every: 10                # Print progress every N epochs
  train_test_split: true         # Enable train-test split
  train_val_split: true          # Enable train-validation split
  test_fraction: 0.2             # Fraction of data reserved for testing
  val_fraction: 0.15             # Fraction of training data for validation
  seed: 0                        # Random seed for reproducibility

# Sample selection
sample_category: "detailed_category"   # Metadata column used for stratification
select_samples:                       # Specify sample categories to include
  - "Lung_Adenocarcinoma"
  - "Breast_Invasive_Carcinoma"

output_dir: ""                        # Output directory (generated automatically if not specified)
plot_results: True                    # Enable visualization of results
```

where,  
- **`source_name`**: Name of the dataset being used.  
- **`data_paths`**: Paths to input data files:  
  - **`rbp_path`**: File containing RBP expression data.  
  - **`isoform_expr_path`**: File containing isoform expression data.  
  - **`metadata_path`**: Metadata file (sample-level information).  
  - **`gene_expr_path`**: File with gene expression data per isoform.  
  - **`getBM_path`**: File for selected gene-RBP mappings.  
- **`model`**: Configuration for the neural network model, including the number of layers, activation function, and optimizer.  
- **`training`**: Training hyperparameters, including epochs, batch size, and data splits.  
- **`sample_category`**: Metadata column used to stratify samples.  
- **`select_samples`**: Specify sample categories to include during training. Use `"all"` to include all samples.  
- **`output_dir`**: Directory where results will be saved. If not specified, it is generated automatically.  
- **`plot_results`**: Enable or disable visualization of results.  

**Note**: You can also include a configuration for an external dataset where evaluations will be performed using the already trained model. This is in addition to the test fraction defined from the source dataset. The configuration for this external dataset would look as follows:  

#### **Example: Configuration for External Dataset Evaluation (`config_gtex.yaml`)** 

```yaml
source_name: "GTEX"

# Paths for the data files
data_paths:
  rbp_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/RBPs_log2p_tpm.csv"
  isoform_expr_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/trans_log2p_tpm.csv"
  metadata_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/phenotype_metadata.csv"
  gene_expr_path: "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/gn_expr_each_iso_tpm.csv"
  getBM_path: "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv" 

# Training configuration
training:
  train_test_split: false
  train_val_split: false
  test_fraction: 0
  val_fraction: 0

# Sample selection
sample_category: "detailed_category"  # The column in metadata to stratify on
select_samples: ["all"]  # Use 'all' to include all samples
output_dir: ""  # Automatically generated if not specified
seed: 0
plot_results: True
```

Once the config.yaml file is ready, execute the script as follows:

```bash
run-deeprbp-predictor \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml" \
  --external_config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_gtex.yaml"
```

### Option 2: Submit a Job in a HPC
If the number of training datasets or the total number of samples is high, we recommend submitting the job using the provided `run_predictor_pipeline.sh` script from the cluster directory. 
This script is adapted to Slurm, but can be easily modified to work on SGE. 
The specific parameters should be adapted depending on the specifications of the HPC.

```bash
cd slurm
sbatch run_predictor_pipeline.sh
```
### Option 3: Running with Docker

---


### Explainability Module
This module uses the already trained DeepRBP Predictor to compute TxRBP (transcript-by-RBP) and GxRBP (gene-by-RBP) scores using DeepLIFT (Shrikumar, Greenside, and Kundaje, 2017) [Learning important features through propagating activation differences, International Conference on Machine Learning, PMLR, pages 3145–3153].

With DeepLIFT, the contribution of each RBP-Transcript pair is determined for every sample in the input data, resulting in a three-dimensional score matrix with the following dimensions:
- Number of transcripts 
- Number of RBPs
- Number of samples.

Positive scores indicate activation of the transcript, while negative scores indicate transcript inhibition. To obtain a single score indicative of the general behavior of each RBP-Transcript pair, we collapse the scores across samples by computing the t-statistic (labeled as “t-stat”), calculated by the formula:
\[ \text{t-stat} = \frac{\text{mean}}{\left(\frac{\sigma}{\sqrt{n}}\right)} \]
where \( n \) represents the number of samples. This results in a score matrix of size: number of transcripts by number of RBPs.

Specific TCGA samples (not presented in the training process) are used to calculate the scores. This module is validated primarily using a binary matrix indicating experimental evidence of regulation in POSTAR3 (Zhao et al., 2022) [POSTAR3: an updated platform for exploring post-transcriptional regulation coordinated by RNA-binding proteins, Nucleic Acids Research, volume 50, D1, pages D287–D294]. POSTAR3 is a comprehensive Post-Transcriptional Regulation database that provides protein binding sites on RNA obtained from CLIP experiments.

Additionally, we have applied our model in in-vitro knockdown experiments.

This module aims to provide insights into how RBPs regulate gene expression. Below is an overview of the validation process and instructions to access the required data.

---

#### Data Access  
The necessary data for running this module is available through the provided Zenodo link. Below is a description of the files:  

- **Events_Regions_gc23_400nt.RData**: Contains detailed information about the genomic regions associated with the events.  
- **EventsFound_gencode23.txt**: Provides metadata about the events, including genomic positions, event types, names, and IDs.  
- **human.txt**: Contains POSTAR3 data, including RNA-binding protein (RBP) binding sites on RNA. This file is tissue-specific and includes data from CLIP experiments as POSTAR peaks.  

To generate a tissue-specific POSTAR matrix, you can use the script `create_gene_rbp_postar_matrix.R`. This script processes the provided files and outputs a matrix customized for a specific tissue.  

**Command example:**  

```bash
Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R \
    /data/jsanchoz/DeepRBP/data/explainability_module/postar3 \
    /scratch/jsanchoz/DeepRBP/data/explainability_module/postar3/processed \
    human_liver \
    human.txt \
    Events_Regions_gc23_400nt.RData \
    EventsFound_gencode23.txt \
    HepG2,Huh7
```
*Argument details*:

- **/path/to/input**: Directory containing the input files.
- **/path/to/output**: Directory where the processed output will be saved.
- **human_liver**: The name of the output file.
- **human.txt**: The POSTAR file containing RBP binding information.
- **Events_Regions_gc23_400nt.RData**: File specifying the genomic regions of the events.
- **EventsFound_gencode23.txt**: File with metadata on events, including IDs and positions.
- **HepG2,Huh7**: Specifies the cell lines from POSTAR experiments to include in the matrix.

By following these steps, you can generate a POSTAR matrix tailored to your specific tissue and experimental needs.
For AML use the K562 cell-line.

## Executing DeepRBP Explainer
There are three options:
* Running the Python script 
* Submitting a job to a HPC queue
* Running with Docker

---

### AQUIIII !!!!!!!!!!!!!!!!!!!!!!!!!!!!

### **Option 1: Running the Python Script**  
To execute DeepRBP on the **TCGA** dataset, use a `.yaml` configuration file. Below is an example configuration file:

#### **Example Configuration File (`config.yaml`)**

```yaml
source_name: "TCGA"




#Lo que antes era Liver_GxRBP.csv ahora se llama human_liver_GxRBP.csv












# nuevo organigrama !!! (puede estar aun sujeto a muchos cambios) ACTUALIZA ESTO BROTHER!!!
/DeepRBP
├── data
│   ├── training_module                       
│   │   ├── raw                                # Datos crudos descargados de TCGA y GTEx
│   │   │   ├── TcgaTargetGtex_rsem_isoform_tpm.gz   # Datos de transcritos en log2(tpm+0.001) (TCGA y GTEx)
│   │   │   ├── TcgaTargetGtex_rsem_gene_tpm.gz      # Datos de genes en log2(tpm+0.001) (TCGA y GTEx)
│   │   │   ├── TcgaTargetGTEX_phenotype.txt         # Datos de fenotipo de TCGA, GTEx y TARGET
│   │   │
│   │   ├── selected_genes_rbps                # Listas seleccionadas de genes y RBPs relevantes
│   │   │   ├── Table_Cancer_Gene_Census.tsv        # Tabla con el censo de genes de cáncer
│   │   │   ├── Table_S2_list_RBPs_eyras.xlsx      # Tabla S2 con lista de RBPs (Eyras)
│   │   │   ├── Table_S5_Cancer_splicing_gene_eyras.xlsx   # Tabla S5 con genes de splicing en cáncer (Eyras)
│   │   │   ├── Table_S6_Cancer_gene_eyras.xlsx            # Tabla S6 con genes de cáncer (Eyras)
│   │   │   ├── getBM.csv                          # Relaciona genes id con su trans id correspondiente
│   │   │
│   │   ├── processed   # Datos procesados y listos para ser usados en el modelo
│   │   │   ├── TCGA / GTEx                      # Datos de TCGA o GTEX(sin normalizar, escalar ni dividir en Train/Test)
│   │   │   │   ├── RBPs_log2p_tpm.csv           # Expresión de RBPs en datos de TCGA
│   │   │   │   ├── gn_expr_each_iso_tpm         # Expresión de genes en datos de TCGA
│   │   │   │   ├── trans_log2p_tpm.csv          # Expresión de transcritos en datos de TCGA
│   │   │   │   ├── phenotype_metadata.csv       # Datos de fenotipo de TCGA o GTEX
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

├── notebooks  # Notebooks de análisis y pruebas
│   ├── Tutorial_predict_transcript_expression.ipynb  # Tutorial para predecir expresión de transcriptos
│   ├── Tutorial_replicate_postar3.ipynb  # Tutorial para replicar los resultados en POSTAR3
│   └── Tutorial_replicate_real_kds.ipynb  # Tutorial para replicar knockdown experiments

├── src  # Código principal del paquete DeepRBP
│   ├── deeprbp
│   │   ├── __init__.py  # Inicialización del paquete DeepRBP
│   │   ├── config_loader.py  # Clase Config y load_config para cargar configuraciones desde un YAML
│   │   ├── main_predictor.py  # (PENDIENTE DE TRASLOCO) función main para ejecutar la pipeline de training del predictor.
│   │   ├── model.py  # Definición de la clase modelo de predicción y explainer (PredictorModel y ExplainerModel)
│   │   ├── processing.py  # clases de procesamiento de datos, responsable de cargar, filtrar, dividir, transformar y escalar los datos
│   │   ├── train_predictor.py  # la clase para entrenar el modelo predictor
│   │   ├── train_explainer.py  # la clase para entrenar el modelo explainer (EN OBRAS)
│   │   ├── utils.py  # Funciones auxiliares
│   │   └── pretrained_model/  # Contiene el modelo preentrenado y sus archivos asociados 
│   │       ├── config.json  # Configuración del modelo preentrenado
│   │       ├── model.pt  # Modelo preentrenado
│   │       ├── scaler_sfs.joblib  # Escaladores usados en el preprocesamiento
│   │       └── sigma_sfs.txt  # Parámetros adicionales del modelo

│   ├── data_preprocessing/  # Preprocesamiento de datos crudos
│   │   └── prep_model_inputs.py  # Preprocesa los datos TCGA/GTEx para generar matrices de input
│   │   └── create_gxrbp.R  # Creates the GxRBP matrix for specific tissues





│   └── tests/  # Tests unitarios para el paquete DeepRBP
│       └── test_data_loader.py  # Test unitario para la clase DataLoader (por ejemplo)

├── slurm  # Scripts ejecutables
│   ├── download_data.sh  # Script para descargar datos (TCGA, GTEx)
│   ├── generate_model_inputs.sh  # Script para procesar los datos descargados y generar matrices de input
│   ├── run_DeepRBP_predictor.sh  # Script para entrenar y evaluar el predictor DeepRBP
│   └── run_explainability.sh  # Script para ejecutar el módulo de explainability

├── images  # Imágenes para visualización (por ejemplo, diagramas o ejemplos de resultados)

├── README.md  # Instrucciones y documentación del proyecto

├── .gitignore  # Archivos y carpetas a ignorar en el control de versiones

└── setup.py  # Script de instalación para el paquete DeepRBP
