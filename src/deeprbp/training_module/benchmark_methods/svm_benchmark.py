
# src/deeprbp/training_module/benchmark_methods/svm_benchmark.py

import argparse
import pandas as pd
import os
import numpy as np
 
from ...data_loading.config_loader import ConfigParser
from ...data_preparation.data_module import DeepRBPDataModule

from .utils_benchmark import generate_X_y_data, create_multi_output_regressor
from ...util.utils import print_if_main, setup_output_directory
from ..evaluation import calculate_general_metrics, spearmanr_per_gene

# maybe we could just use one alternative_model_benchmark y en el decirle que modelo alternativo queremos probar
# y nos acortamos de usar dos scripts para casi lo mismo.
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/config_alternative/config_model_train_svm_benchmark.yaml'
output_dir = '/scratch/jsanchoz/DeepRBP/output/results/run_svm_benchmark'
selected_algorithm = 'svr' # 'decision_trees'
calculate_abundance = True

# vamos a probar dos approaches: due possibili metodi
        # - o calcular directamente la expression del transcrito
        # - o calcular el isoforma abundance como nuestro modelo
# en ambos casos al final calculamos las métricas que calculamos tb para nuestro modelo 
# en log2tpm+1

args = parse_args()
    
# Load configuration and auxiliary file
print_if_main('\n[main_predictor] 🚀 Loading configuration...')
config = ConfigParser(args.config_path) 

# Determine the output directory based on gpu rank or cpu device
output_dir = setup_output_directory(args.output_dir)
print_if_main('\n[main_predictor] Output directory for main process: ', output_dir)

# Load data module and prepare data for training
print_if_main('\n[main_predictor] 🚀 Initializing DataModule...')
dm = DeepRBPDataModule(config, output_dir)
print_if_main('[main_predictor] 🚀 Preparing data for training...')
dm.setup('fit')  
dm.setup('test') 
print_if_main('\n[main_predictor] ──────────────────────────────────────')

# Get data
X_train, y_train = generate_X_y_data(
        dm.train_dataset.to_numpy(), 
        calculate_abundance=calculate_abundance) #args.calculate_abundance)

X_val, y_val = generate_X_y_data(
        dm.val_dataset.to_numpy(), 
        calculate_abundance=calculate_abundance) #args.calculate_abundance)

X_test, y_test = generate_X_y_data(
        dm.test_dataset.to_numpy(), 
        calculate_abundance=calculate_abundance) #args.calculate_abundance)

# Definir el modelo como un regressor de múltiples salidas
model = create_multi_output_regressor(args.selected_algorithm)

# Ajustar el modelo

# print(X_train.shape, X_val.shape, X_test.shape)
# print(y_train.shape, y_val.shape, y_test.shape)
# input_columns = [0, 1, 2]  # Ejemplo de columnas para X
# output_columns = [0, 1]

# # Crear subconjuntos de las matrices X e y
# X_train = X_train[:, input_columns]
# X_val = X_val[:, input_columns]
# X_test = X_test[:, input_columns]
# y_train = y_train[:, output_columns]
# y_val = y_val[:, output_columns]
# y_test = y_test[:, output_columns]
# #############

multioutput_model.fit(X_train, y_train)

# Realizar predicciones
y_pred_val = multioutput_model.predict(X_val) 
y_pred_test = multioutput_model.predict(X_test)

# Aplicar la transformación logarítmica a las predicciones para obtener la y_pred_final
if args.calculate_abundance:
    print("Performing the transformation to obtain final predictions in log2(TPM + 1).")
    y_pred_val = np.log2((y_pred_val * dm.val_dataset.to_numpy()['gene_df']) + 1)
    y_pred_test = np.log2((y_pred_test * dm.test_dataset.to_numpy()['gene_df']) + 1)

# Calcular métricas de evaluación


# Ejemplo de uso
# Supongamos que dataset, getBM, gene_names, trans_names, outputs, y labels están definidos
calculate_metrics(outputs, labels, getBM, gene_names, trans_names):

# VAL 
val_metrics_general = calculate_general_metrics(dm.val_dataset.to_numpy()['isoform_df'].flatten(), y_pred_val.flatten())
val_metrics_per_gene = spearmanr_per_gene(
                        dm.val_dataset.gene_names, 
                        dm.getBM,
                        dm.val_dataset.trans_names,
                        y_pred_val, dm.val_dataset.to_numpy()['isoform_df'])

# TEST
test_metrics_general = calculate_general_metrics(dm.test_dataset.to_numpy()['isoform_df'].flatten(), y_pred_test.flatten())
test_metrics_per_gene = spearmanr_per_gene(
                        dm.test_dataset.gene_names, 
                        dm.getBM,
                        dm.test_dataset.trans_names,
                        y_pred_test, dm.test_dataset.to_numpy()['isoform_df'])
    
  

def parse_args():
    parser = argparse.ArgumentParser(description='Run the predictor using other traditional Machine Learning algorithms as a benchmark.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results.')
    parser.add_argument('--selected_algorithm', type=str, choices=['svr', 'decision_trees'], default='svr',
                        help='Algorithm to be used for prediction. Choices are: "svr" or "decision_trees".')
    parser.add_argument('--calculate_abundance', type=bool, default=True,
                        help='If True, fit the model with transcript abundance; if False, use transcript log2p(TPM) expression.')
    return parser.parse_args()


if __name__ == "__main__":
    if os.environ.get("LOCAL_RANK")=="0":
        print_gpu_memory_info()
    set_random_seed()
    main()




#-2)	Idoia: comparar el Predictor con un decisión tree o SVM como otro baseline. Mira multi output regressor.


# esto habria que amoldarlo para que acepte models que no sean nuestro Deep Learning model

# Evaluation of the model's performance by category
        # print('\n[main_predictor] 🚀 Evaluating model performance by category...')
        # datasets = [(dm.train_data, 'train'), (dm.val_data, 'val'), (dm.test_data, 'test')]

        # for data, set_name in datasets:
        #     print(f"[main_predictor] 🚀 Evaluating on the {set_name} set...")
        #     evaluate_and_visualize_metrics_by_category(
        #         test_data=data,
        #         trainer=trainer,
        #         model=model,
        #         dm=dm,
        #         output_dir=output_dir,
        #         set_name=set_name,
        #         plot_results=config.get('plot_results')
        #     )
        # print('\n[main_predictor] 🚀 Process completed.')
