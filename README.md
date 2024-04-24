# DeepRBP
## A novel deep neural network for inferring splicing regulation
#### Publication: https://doi.org/10.1101/2024.04.11.589004

<!-- ABOUT THE PROJECT -->
## Description
<p align="center">
    <img src="https://github.com/ML4BM-Lab/DeepRBP/blob/main/images/methods_deepsf.pdf" width="700" alt="PDF Image">
</p>

 <!-- <p align="center"><a href=https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(23)00333-X/fulltext>Discovering the mechanism of action of drugs with a sparse explainable network<a></p> -->

Alternative splicing plays a pivotal role in various biological processes. In the context of cancer, aberrant splicing patterns can lead to disease progression and treatment resistance. 
Understanding the regulatory mechanisms underlying alternative splicing is crucial for elucidating disease mechanisms and identifying potential therapeutic targets.
We present DeepRBP, a deep learning (DL) based framework to identify potential RNA-binding proteins (RBP)-Gene regulation pairs for further in-vitro validation. DeepRBP is composed of a DL model 
that predicts transcript abundance given RBP and gene expression data coupled with an explainability module that computes informative RBP-
Gene scores using DeepLIFT. We show that the proposed framework is able to identify known RBP-Gene regulations, demonstrating its applicability to identify new ones.

This project includes instructions for:
 1) Apply already trained DeepRBP to calculate RBP-Gene pair scores (explainability module) on your data using a trained model
 2) Train your own DeepRBP predictor using TCGA or another data.
 3) Replicate paper results using 3 POSTAR3 datasets and 6 real knockdown experiments 

#### Publication: 

### Built With
*   <a href="https://www.python.org/">
      <img src="https://www.python.org/static/community_logos/python-logo.png" width="110" alt="python" >
    </a>
*   <a href="https://pytorch.org/">
      <img src="https://pytorch.org/assets/images/pytorch-logo.png" width="105" alt="pytorch" >
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

On the other hand, the DeepRBP explainer has been validated using TCGA samples from a tumor type to calculate the GxRBP scores and a binary matrix with shape GxRBP indicating whether the evidence in POSTAR3 experiments indicates regulation or not. Also, the DeepRBP explainer has been tested in real knockdown experiments. Below, you are instructed on how to download these data.

## Data
To download and process the TCGA and GTeX data used in this project you need to execute the following shell scripts:

```bash
data/input_create_model/raw/download_script.sh
data/input_create_model/processed/create_data.sh
```
With this data, you will be able to create a model from scratch and/or check the performance of the prediction of the transcripts, such as performing explainability.
On the other hand, if you want to identify the RBPs that regulate a gene for your experiment, as a previous step in the next tutorial, *./data/data_real_kds/tutorial_download_kd_experiments*, we explain how to download the fastq files and run kallisto. We also explain how to use voom-limma for differential expression analysis between two conditions as further validation of the scores obtained with DeepRBP."

After running the shell script, a folder named *splitted_datasets* will be generated in the directory *./data/input_create_model/processed*. This folder contains processed samples from both TCGA and GTeX datasets, divided by tumor type or tissue. 

Within the TCGA dataset, samples are divided for each tumor type, allocating 80% for training and 20% for testing. Each of these tumor type folders contains the two main inputs of the model and output variable:

- **_.gn_expr_each_iso_tpm.csv_**:  Gene expression in TPM for each transcript with shape *n_samples x n_transcripts*.
- **_.RBPs_log2p_tpm.csv_**: RBP expression in log2(TPM+1) with shape *n_samples x n_rbps*.
- **_.trans_log2p_tpm.csv_**: Transcripts expression in log2(TPM+1) with shape *n_samples x n_transcripts*.

To download the POSTAR3 data with information on RBP binding sites in the genome from CLIP experiments, you need to go to the POSTAR3 website and request access to the data.

## Apply the model
In the `./apply_the_model` folder, you will find tutorials on how to use a DeeRBP predictor model trained for transcript prediction and explainability.
In particular, there's a tutorial to verify the results obtained in the DeepRBP explain module on Postar3 data using TCGA data and in a real knockout experiment. Alternatively, you can once again run the .sh scripts to directly apply a trained model for calculating DeepLIFT scores in experiments of your interest.

## Create the model
In the following folder `./model`, you will be able to create your own predictive model using the DeepRBP architecture, using the samples generated from TCGA (or others that the user wants). A Jupyter notebook tutorial is provided for you to familiarize yourself with the method. Alternatively, you can call the shell script main_predictor.sh to train the DL model according to your preferences.





