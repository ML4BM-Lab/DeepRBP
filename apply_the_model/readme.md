# Model Directory (./apply_the_model)

The `./apply_the_model` directory contains all the necessary code and utilities to apply your custom model using preprocessed data.

## Directory Structure
  
- `calculate_deeplift_values.py`: Module with functions for calculating deeplift scores.
  
- `main_evaluate_predictions.py`: Script to make predictions using a trained DeepRBP model and processed data.
  
- `main_explain_postar3.py`: Script for explainability analysis with DeepRBP using custom samples and comparison with a Postar3 matrix.
  
- `main_explain_real_kds.py`: Script for explainability analysis with DeepRBP using control samples from a knockdown experiment and comparing results with DE transcripts or genes from a limma experiment.

## Tutorials
There are three tutorials available to help you get familiar with the method:

- `apply_the_model/Tutorial_predict_transcript_expression.ipynb`
  
- `apply_the_model/Tutorial_replicate_postar3.ipynb`
  
- `apply_the_model/Tutorial_replicate_real_kds.ipynb`

## Model Output
In the model output, we will find 3 folders with the results obtained from the prediction of transcripts, in postar3, or in real knockdown experiments.

#### Scripts
This directory contains 6 example shell scripts with real knockdowns of an RBP, 3 shell scripts with 3 Postar matrices, and one script to evaluate the prediction of transcripts from your model.

- `apply_the_model/scripts/execute_explain_real_kds_lines_GSE136366_TARDBP_kd.sh`
  
- `apply_the_model/scripts/execute_explain_real_kds_lines_GSE77702_kd1_FUS_kd.sh`
  
- `apply_the_model/scripts/execute_explain_real_kds_lines_GSE77702_kd2_TAF15_kd.sh`
  
- `apply_the_model/scripts/execute_explain_real_kds_lines_GSE77702_kd3_TARDBP_kd.sh`
  
- `apply_the_model/scripts/execute_explain_real_kds_lines_HFE_MBNL1_kd.sh`
  
- `apply_the_model/scripts/execute_explain_real_kds_lines_SRR296_RBM47_kd.sh`
  
- `apply_the_model/scripts/run_DeepRBP_evaluate_predictions.sh`
  
- `apply_the_model/scripts/run_DeepRBP_explain_postar3_AML.sh`
  
- `apply_the_model/scripts/run_DeepRBP_explain_postar3_Kidney.sh`
  
- `apply_the_model/scripts/run_DeepRBP_explain_postar3_Liver.sh`
