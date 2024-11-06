# /scratch/jsanchoz/DeepRBP/src/deeprbp/config_loader.py

import yaml
from datetime import datetime

def generate_unique_id(config):
    timestamp = datetime.now().strftime("%Y-%m-%d")
    source_train = config["source_name"]
    tumor_types = '-'.join([t.split('_')[0] for t in config["select_samples"]])
    return f"{source_train}_{tumor_types}_{timestamp}"

def load_config(yaml_path):
    with open(yaml_path, 'r') as file:
        config = yaml.safe_load(file)
    if "output_dir" not in config:
        unique_id = generate_unique_id(config)
        config["output_dir"] = f'/scratch/jsanchoz/DeepRBP/output/results/analysis/{unique_id}/train_prediction_model'
    return config



### version vieja igual borrar: 
# import datetime
# def generate_unique_id_1(config):
#     timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
#     source_train = config["train_prediction_model"]["source_name"]
#     tumor_types = '-'.join([t.split('_')[0] for t in config["train_prediction_model"]["train_tumor_types"]])
#     return f"{source_train}_{tumor_types}_{timestamp}"
# def get_default_config(): 
#     config = {
#         "train_prediction_model": { 
#             "source_name": 'TCGA', # source dataset used for training (ESTA CREO QUE VA A SOBRAR, E IGUAL SE PUEDE METER UN TOY=FALSE PARA OPTIMIZATION)
#             "train_test_split": True, # whether we wanna split data in train and test
#             "test_frac": 0.2, # by tumor-type-stratified test data fraction
#             "train_val_split": True, # from remaining train data or data if we wanna split in train and val
#             "val_frac": 0.15,
#             "seed": 0, 
#             "batch_size": 128,
#             "cuda": True,
#             "epochs": 1000,
#             "train_tumor_types": ['Lung_Adenocarcinoma', 'Breast_Invasive_Carcinoma'],  # list with tumor type names or 'all' (string) used for modelling
#             "sample_category": 'detailed_category' # column in the metadata DataFrame that contains the classification of each patient or sample according to the specific type of tumor they have (OPTIONAL)
#         },
#         "val_prediction_model": {
#             "source_name": 'GTEX', # source dataset used for training (ESTA CREO QUE VA A SOBRAR, E IGUAL SE PUEDE METER UN TOY=FALSE PARA OPTIMIZATION)
#             "train_test_split": False, # whether we wanna split data in train and test
#             "train_val_split": False, # from remaining train data or data if we wanna split in train and val
#             "seed": 0, 
#             "batch_size": 128,
#             "cuda": True
#         # "explain_model" luego explain_module con parametros como deeplift method, sample batch dimension reduction
#         #"data": {
#         #    "training_module": {
#         #        "TCGA": {
#         #            "rbp_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/RBPs_log2p_tpm.csv",
#         #            "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/trans_log2p_tpm.csv",
#         #            "metadata_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/phenotype_metadata.csv",
#         #            "gene_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/gn_expr_each_iso_tpm.csv"
#         #        }
#         #        "GTEX": {
#         #            "rbp_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/RBPs_log2p_tpm.csv",
#         #            "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/trans_log2p_tpm.csv",
#         #            "metadata_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/phenotype_metadata.csv",
#         #            "gene_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/GTEX/gn_expr_each_iso_tpm.csv"
#         #        }
#         #    }
#         #}
#     }
#     unique_id = generate_unique_id_1(config)
#     # path to save results in this section
#     config["train_prediction_model"]["output_dir"] = f'/scratch/jsanchoz/DeepRBP/output/results/analysis/{unique_id}/train_prediction_model'
#     return config
# def get_config(user_conf={}):
#     config = get_default_config()
#     for key in user_conf.keys():
#         if key in config.keys():
#             config[key] = {**config[key], **user_conf[key]} if isinstance(config[key], dict) else user_conf[key]
#         else:
#             config[key] = user_conf[key]
#     return config
