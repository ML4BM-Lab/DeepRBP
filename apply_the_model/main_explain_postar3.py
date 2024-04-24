import argparse
import os
import json
import torch

from utils.Utils import prepare_inputs_explainability, read_postar3_data
from utils.Config import Config
from DLModelClass.modelsNN import DeepRBP
from calculate_deeplift_values import main as calculate_deeplift_values

def main_explain(args, config_obj, path_train_files, PATH_DATA, tumor2tissue, path_save):
    ### 1) Prepare INPUTS to perform the in-silico validation of our DL model.
    data_inputs = prepare_inputs_explainability(args.tumor_type, args.source_explain, config_obj, PATH_DATA, path_train_files)
    print('[explainability] 1) Prepare INPUTS to perform the in-silico validation of our DL model ... -> DONE')

    ### 2) Load the trained model
    model = DeepRBP(n_inputs=data_inputs.df_scaled_test.shape[1], n_outputs=data_inputs.test_labels.shape[1], config=config_obj, device=data_inputs.device)
    model.load_state_dict(torch.load(f'{path_train_files}/model.pt', map_location=data_inputs.device))      
    model.eval()
    
    ### 3) Load POSTAR experimental data with GxRBP relationships   
    df_val_GxRBP = read_postar3_data(PATH_DATA+'/../../data_postar3', args.tumor_type, tumor2tissue)
    print('[explainability] 2) Load Postar data ... -> DONE')
    
    ### 4) Perform DEEPLIFT method
    print('[explainability] 3) Do DeepLIFT ... -> Initializing')
    result = calculate_deeplift_values(df_scaled_test=data_inputs.df_scaled_test, 
                                      test_labels=data_inputs.test_labels, 
                                      test_gn=data_inputs.test_gn, 
                                      model=model, 
                                      path_save=path_save, 
                                      df_val_GxRBP=df_val_GxRBP,
                                      path_data=PATH_DATA)
    
    df_tstat_scores_GxRBP_list = result.df_tstat_scores_GxRBP_list
    df_sum_scores_GxRBP_list = result.df_sum_scores_GxRBP_list
    deep_studies_list = result.deep_studies_list
    print('[explainability] 3) DeepLIFT ... -> DONE')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Explainability parameters')
    parser.add_argument('--path_save', type=str, help='Path to save results', required=True)
    parser.add_argument('--path_data', type=str, help='Path to the data directory', required=True)
    parser.add_argument('--tumor_type', type=str, help='Tumor type for explainability', required=True)
    parser.add_argument('--source_explain', type=str, help='Source for explainability', required=True)
    parser.add_argument('--path_train_files', type=str, help='Path to trained files', required=True)
    
    args = parser.parse_args()

    PATH_SAVE = args.path_save
    PATH_DATA = args.path_data

    tumor2tissue = { 
            "Kidney_Chromophobe": "Kidney_embryo_GxRBP",
            "Liver_Hepatocellular_Carcinoma": "Liver_GxRBP", 
            "Acute_Myeloid_Leukemia": "Myeloid_GxRBP"
            }
    
    print("Welcome again! You're about to perform the explainability of the DeepRBP model")
    print(f"The experiments will be shown in {PATH_SAVE}")
    print(f"The path_data used for this explainability is: {PATH_DATA}")
    
    with open(f'{args.path_train_files}/config.json', 'r') as file:
        config_dict = json.load(file)
    print(config_dict)
    
    # Select the tumor type to be tested in the explainability and load the final config_obj
    config_dict['test_tumor_types'] = args.tumor_type
    config_obj = Config(**config_dict)
    print(config_obj.get_config())
    print(f'Samples selected to perform the explainability with POSTAR3: {config_obj.test_tumor_types}')

    main_explain(args, config_obj, args.path_train_files, PATH_DATA, tumor2tissue, PATH_SAVE)
