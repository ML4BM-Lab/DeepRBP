
In this tutorial, we will show you an example of the application of DeepRBP in the GSE136366 experiment where the RBP TARDBP (TDP43) was knocked out. In Tutorial1_real_kds.ipynb, the code is provided to perform the following steps:

1) Process raw samples to generate the input data for DeepRBP (JUAN is doing an Rmarkdown for this). From this code, a folder named datasets is created that contains the necessary inputs to run the DeepRBP model for both control and knockout samples (RBP and gene expression matrices, as well as those of the transcripts).
2) Pass the control data through the predictive model trained on TCGA and make predictions.
3) Calculate Pearson correlation, Spearman, and MSE.
4) Pass the data through the explainability module.
5) Calculate the RBPxTranscripts and RBPxGenes matrices.
6) Some plots.

## Data information
In the *experiments* folder, you will find:
The 6 samples from this experiment containing abundance and expression data in TPM.
A txt file with the description of each sample (whether is "control" or "knockout".
In the limma folder, there is a table with results from the differential expression analysis ("Transcript_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Gene_ID").
This information is used to plot the results once the scores have been calculated in step 5. You also have the R script to reproduce this analysis in limma_analysis_real_kds.R.

In the *tutorials/model* folder, you will find a model already trained with TCGA, along with information on the training configuration and scaler used.

You also have a getBM that links transcript identifiers to their genes.

This notebook uses functions from *DeepRBP/utils* scripts for processing and plotting, as well as calculates the explainability scores using *DeepRBP/calculate_deeplift_values.py*.

As a result, a results folder will be created containing the explainability scores at the transcript and gene levels, as well as plots comparing the scores of DE transcripts and genes against non-DE using Limma results.