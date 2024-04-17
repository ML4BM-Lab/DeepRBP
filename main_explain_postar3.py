import argparse
import os
import json

ruta_utils = '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/utils'
os.environ['PYTHONPATH'] = ruta_utils + ':' + os.environ.get('PYTHONPATH', '')
print('utils')

from utils import Utils as utils 
from utils.Config import Config
from calculate_deeplift_values import main as calculate_deeplift_values

def main_explain(args, config_obj, path_train_files, path_data, tumor2tissue, path_save):
    ### 1) Prepare INPUTS to perform the in-silico validation of our DL model.
    data_inputs = utils.prepare_inputs_explainability(args, config_obj, path_data, path_train_files)
    print('[explainability] 1) Prepare INPUTS to perform the in-silico validation of our DL model ... -> DONE')

    ### 2) Load POSTAR experimenental data with GxRBP relationships
    df_val_GxRBP = utils.read_postar3_data(PATH_DATA, args.tumor_type, tumor2tissue, args.cell_line)
    print('[explainability] 2) Load Postar data ... -> DONE')
    
    ### 3) Perform DEEPLIFT method
    print('[explainability] 3) Do DeepLIFT ... -> Initializing')
    result = calculate_deeplift_values(df_scaled_test=data_inputs.df_scaled_test, 
                                      test_labels=data_inputs.test_labels, 
                                      test_gn=data_inputs.test_gn, 
                                      model=data_inputs.model, 
                                      path_save=path_save, 
                                      df_val_GxRBP=df_val_GxRBP,
                                      path_data=path_data)

    df_tstat_scores_GxRBP_list = result.df_tstat_scores_GxRBP_list
    df_sum_scores_GxRBP_list = result.df_sum_scores_GxRBP_list
    deep_studies_list = result.deep_studies_list
    print('[explainability] 3) DeepLIFT ... -> DONE')

 
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Explainability parameters')
  parser.add_argument('--tumor_type', type=str, default='Liver_Hepatocellular_Carcinoma', help='The test data is obtained for this tumor to perform explainability')
  parser.add_argument('--cell_line', type=str, default='all', help='The postar GxRBP data is obtained from this cell-line', required=False)
  parser.add_argument('--source_explain', default='TCGA', choices=['TCGA', 'GTEX', 'all'], help='Data source')
  parser.add_argument('--model_selected', type=str, default='1024N_2HL_8f', help='Trained model selected')
  parser.add_argument('--path_result', type=str, default='/scratch/jsanchoz/ML4BM-Lab/DeepRBP/results', help='Path for the results directory')
  parser.add_argument('--path_data', type=str, default='/scratch/jsanchoz/data', help='Path for the data directory')
  args = parser.parse_args()

  PATH_SAVE = args.path_result
  PATH_DATA = args.path_data

  tumor2tissue = { 
            "Kidney_Chromophobe": "Kidney_embryo_GxRBP",
            "Liver_Hepatocellular_Carcinoma": ["Liver_huh7_GxRBP", "Liver_hepg2_GxRBP", "Liver_all_GxRBP"], 
            "Acute_Myeloid_Leukemia": "Myeloid_GxRBP"
            }

  print("Welcome again! You're about to perform the explainability of the DeepRBP model")
  print(f"The experiments will be shown in {PATH_SAVE}")
  print(f"The path_data used for this explainability is: {PATH_DATA}")
  print(f"The selected model to perform explainability is {args.model_selected}")

  # if args.cell_line is not None:
    # path_save = os.path.join(PATH_SAVE, 'explainability', f'model_{args.model_selected}', args.tumor_type, args.cell_line)
  # else:
  path_save = os.path.join(PATH_SAVE, 'explainability', f'model_{args.model_selected}', args.tumor_type)
  utils.check_create_new_directory(path_save)
  
  path_train_files = PATH_DATA+f'/trained_models/model_{args.model_selected}'
  print(f'[run_experiment] Path train files is: {path_train_files}')
  with open(f'{path_train_files}/config.json', 'r') as file:
    config_dict = json.load(file)

  # Select the tumor type to be tested in the explainability and load the final config_obj
  config_dict['test_tumor_types'] = args.tumor_type
  config_obj = Config(**config_dict)
  print(config_obj.get_config())

  # Call the main_explain function with all required arguments
  main_explain(args, config_obj, path_train_files, PATH_DATA, tumor2tissue, path_save)
