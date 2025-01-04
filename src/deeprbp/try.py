### TRYING THE NEW CODE!!

import os
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')
from captum.attr import DeepLift
from config_loader import ConfigParser
from processing import DataImporter, DatasetLoader, Scaler
from model import PredictorModel
from utils import *
from logger import Logger
from .deeplift_handler import DeepLiftHandler

# from .utils import *
# from .config_loader import ConfigParser
# self.config_parser = ConfigParser(config_path)
# self.base_config = self.config_parser.get_base_config()
# self.explain_config = self.config_parser.get_explainability_config()

### 1) Prepare INPUTS to perform the in-silico validation of our DL model (old title list).
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_explain.yaml'
config_parser = ConfigParser(config_path)
base_config = config_parser.get_base_config()
explain_config = config_parser.get_explainability_config()

# Define paths for saving data and results
#self.path_save_data = os.path.join(self.base_config['output_dir'], 'data')
#self.path_save_results = os.path.join(self.base_config['output_dir'], 'results')
path_save_data = os.path.join(base_config['output_dir'], 'data')
path_save_results = os.path.join(
    base_config['output_dir'], 
    'results',
    f"{explain_config['explanation_method']}_{explain_config['reference_data']}_{explain_config['batch_reduction_method']}_{explain_config['gene_collapse_method']}"
)
ensure_directory_exists(path_save_data)
ensure_directory_exists(path_save_results)

# Initialize DataImporter for primary data
#self.data_importer = DataImporter(self.base_config['data_paths'])
#self.scaler = Scaler()
data_importer = DataImporter(base_config['data_paths'])
data_loader = DatasetLoader(data_importer, base_config)
scaler = Scaler.load(explain_config['scaler_path'])

#def load_and_process_data(self, loader):
data = data_loader.load_data()
#return data
data['scaled_rbp_expr_df'] = scaler.transform(data['rbp_expr_df'])

### 2) Load the trained model (old title list).
config_path_train = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml'
config_train_parser = ConfigParser(config_path_train)
training_config = config_train_parser.get_model_training_config()

model = PredictorModel.load_model(
        path_to_weights = os.path.join(explain_config['trained_model_path'], explain_config['model_file']),
        config=training_config
        )

### 3) Load POSTAR experimental data with GxRBP relationships   
df_val_GxRBP = pd.read_csv(
        os.path.join(explain_config['postar_matrix_path'], explain_config['postar_file']), 
        index_col=0
)

### 4) Perform DEEPLIFT method
# Create an instance of DeepLiftHandler
deeplift_handler = DeepLiftHandler(model, data, explain_config)

# Prepare RBP tensors
scaled_rbp_tensor, gn_tensor, reference_rbp_tensor = deeplift_handler.prepare_rbp_tensors()

# Compute attribution scores
list_batch_scores = deeplift_handler.compute_attribution_scores(scaled_rbp_tensor, reference_rbp_tensor, gn_tensor)

# Reduce batch dimension (RBP x T)
df_deeplift_scores_TxRBP = deeplift_handler.reduce_batch_dimension(list_batch_scores)

# Set low-expressed genes' scores (mean < 1TPM) to 0.
df_deeplift_scores_TxRBP = deeplift_handler.filter_scores_for_low_expressed_genes(
    deeplift_scores=df_deeplift_scores_TxRBP,
    gene_expr_df=data['gene_expr_df'],
    threshold=1
)

print(df_deeplift_scores_TxRBP)

####
####

# 4) Collapse scores to genes (RBP x G) ME HE QUEDADO AQUI BROTHER !!!
#### AQUII !!! explain_config['gene_collapse_method'] getBM_path = base_config['data_paths'].get('getBM_path', None)
getBM = pd.read_csv(getBM_path, index_col=0)
#path_save_results
df_deeplift_scores_GxRBP = collapse_transcript_scores_to_genes(
            df_deeplift_scores_TxRBP, 
            getBM)


# Sorting the getBM DataFrame to ensure proper alignment
getBM_sorted = getBM.sort_values(by='Transcript_ID').reset_index(drop=True)
# Ensure that df_score_TxRBP's index matches Transcript_ID from getBM
df_deeplift_scores_TxRBP = df_deeplift_scores_TxRBP.loc[getBM_sorted['Transcript_ID']]
  
# COMPROBACIONES QUE HAY QUE HACER:
#•	En Postar3 NO puede haber NaN en los genes, verificar al construir la matriz de Postar si esos genes porque no se matchean que tenemos NaNs! A ÁNGEL le sorprendia igualmenter que a at gene level los NaN sean mas bajos tb.
#•	Referencias DeepLIFT usar todo 0.5 como referencia? 