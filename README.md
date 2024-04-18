# DeepRBP
## A novel deep neural network for inferring splicing regulation
#### Publication: https://doi.org/10.1101/2024.04.11.589004

<!-- ABOUT THE PROJECT -->
## Description
href https://github.com/ML4BM-Lab/DeepRBP/blob/main
 <p align="center"><img src="images/methods_deepsf.pdf" width="700" alt="Logo"></p>

 <!-- <p align="center"><a href=https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(23)00333-X/fulltext>Discovering the mechanism of action of drugs with a sparse explainable network<a></p> -->

Alternative splicing plays a pivotal role in various biological processes. In the context of cancer, aberrant splicing patterns can lead to disease progression and treatment resistance. 
Understanding the regulatory mechanisms underlying alternative splicing is crucial for elucidating disease mechanisms and identifying potential therapeutic targets.
We present DeepRBP, a deep learning (DL) based framework to identify potential RNA-binding proteins (RBP)-Gene regulation pairs for further in-vitro validation. DeepRBP is composed of a DL model 
that predicts transcript abundance given RBP and gene expression data coupled with an explainability module that computes informative RBP-
Gene scores using DeepLIFT. We show that the proposed framework is able to identify known RBP-Gene regulations, demonstrating its applicability to identify new ones.

This project includes instructions for:
* training the DeepRBP predictor
* making transcript abundance predictions and calculate RBP-Gene pair scores using the trained DeepRBP DL model and DeepLIFT 

#### Publication: 

### Built With

*   <a href="https://www.python.org/">
      <img src="images/python.png" width="110" alt="python" >
    </a>
*   <a href="https://pytorch.org/">
      <img src="images/pytorch.png" width="105" alt="pytorch" >
    </a>

## Clone repository & create the Conda environment
To crate a Conda environment using the provided environment.yml, follow these steps:

```bash
git clone https://github.com/ML4BM-Lab/DeepRBP.git
cd DeepRBP
conda env create -f environment.yml
conda activate DeepRBP
```

## Datasets information
In this project, we have used several databases. On one hand, we have used a cohort that contains, among others, samples from TCGA and GTEx. The samples from the former have been used to train the DeepRBP predictor that learns transcript abundances, and the samples from GTEx have been used to evaluate how well the predictive model generalizes.

On the other hand, the DeepRBP explainer has been validated using TCGA samples from a tumor type to calculate the GxRBP scores and a binary matrix with shape GxRBP indicating whether the evidence in POSTAR3 experiments indicates regulation or not. Finally, the DeepRBP explainer has been tested in real knockdown experiments. Below, you are instructed on how to download these data.

### Downloading the datasets
To download TCGA and GTeX data you have to run the *download_script.sh*.

```bash
cd DeepRBP
bash data/TCGA_GTeX/raw/download_script.sh
```
And then to process these data you have to run *create_data.sh*
```bash
cd DeepRBP
bash create_data.sh
```
After running the shell script, a folder named *splitted_datasets* will be generated in the directory *./data/TCGA_GTeX*. This folder contains processed samples from both TCGA and GTeX datasets, divided by tumor type or tissue. 

Within the TCGA dataset, samples are divided for each tumor type, allocating 80% for training and 20% for testing. Each of these tumor type folders contains the two main inputs of the model and output variable:

- **_.gn_expr_each_iso_tpm.csv_**:  Gene expression in TPM for each transcript with shape n_samples x n_transcripts.
- **_.RBPs_log2p_tpm.csv_**: RBP expression in log2(TPM+1) with shape n_samples x n_rbps
- **_.trans_log2p_tpm.csv_**: Transcripts expression in log2(TPM+1) with shape n_samples x n_transcripts

During the data generation process, the tables located in *./data/selected_genes_rbps* are utilized to choose the RBPs and genes under consideration for this study. An already trained model with TCGA is given in *./tutorials/model* with the .pt file and a .json with the model configuration and the scaler used in training.

Furthermore, within *./data/extra*, you'll find a getBM DataFrame that contains information linking transcript names to each gene, along with a reduced version that includes only the genes ultimately considered in this study.

To download the Postar3 human database containing information on CLIP experiments you have to run the *download_script.sh*
```bash
cd DeepRBP
bash data/postar3/input_data/download_script.sh
```
To download real knockdown experiments raw data, *./tutorials/tutorial_0* contains a .Rmd file that explains how to download the fastq files, run kallisto and voom-limma for a given experiment.

### Tutorials
In ./tutorials, several tutorials are presented to test the different functionalities of DeepRBP:
- Tutorial 1 : DeepRBP on real knockdown data using a trained model with TCGA.
- Tutorial_2: DeepRBP on POSTAR3 using a trained model with TCGA.
- Tutorial_3: How to train the DeepRBP predictor from scratch.
