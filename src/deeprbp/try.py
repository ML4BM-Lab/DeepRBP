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
getBM_path = base_config['data_paths'].get('getBM_path', None)
getBM = pd.read_csv(getBM_path, index_col=0)
#path_save_results

# steps:
explain_config['explanation_method']
explain_config['reference_data']
explain_config['batch_reduction_method']
explain_config['gene_collapse_method']

### calculate_deeplift_values.py (DeepLiftHandler Class)
# from here down it could go inside a new script "calculate_deeplift_values.py" or a better name

##
rbps_id = list(data['scaled_rbp_expr_df'].columns)
trans_id = list(data['trans_expr_df'].columns)

def prepare_rbp_tensors(data, explain_config, logger): # the logger we need to remove it later on when we introduce this function inside a bigger class.
    """
    Prepare scaled RBP tensor and reference RBP tensor based on the configuration.

    Parameters:
        data (dict): Dictionary containing the scaled RBP expression DataFrame under 'scaled_rbp_expr_df'
                     and gene expression DataFrame under 'gene_expr_df'.
        explain_config (dict): Dictionary containing the configuration. 
                               Expected key: 'reference_data' with options:
                               - 'median_reference': Median of the scaled RBP expression.
                               - 'knockdown_reference': Zeroed tensor (knockdown).
                               - 'half_reference': Half-zeroed tensor.
        logger (Logger): Instance of the Logger class for logging messages.

    Returns:
        tuple: A tuple containing:
            - scaled_rbp_tensor (Tensor): Tensor of scaled RBP expression values.
            - gn_tensor (Tensor): Tensor of gene expression values.
            - reference_rbp_tensor (Tensor): Tensor used as a reference for RBP expression.
    """
    # Convert scaled RBP and gene expression data to tensors in float32
    scaled_rbp_tensor = torch.tensor(data['scaled_rbp_expr_df'].values, dtype=torch.float32)
    gn_tensor = torch.tensor(data['gene_expr_df'].values, dtype=torch.float32)
    logger.log("[prepare_rbp_tensors] Converted scaled RBP expression and gene expression to tensors (float32).")
    
    # Determine reference tensor based on configuration
    reference_type = explain_config.get('reference_data')
    if reference_type == 'median_reference':
        logger.log("[prepare_rbp_tensors] Using median of RBP expression as reference.")
        reference_rbp_tensor = torch.tensor(np.median(data['scaled_rbp_expr_df'], axis=0))
        reference_rbp_tensor = torch.reshape(reference_rbp_tensor, (1, reference_rbp_tensor.size(0)))
        logger.log("[prepare_rbp_tensors] Median reference tensor created.")

    elif reference_type == 'knockdown_reference':
        logger.log("[prepare_rbp_tensors] Using knockdown (zeroed) RBP expression as reference.")
        reference_rbp_tensor = torch.zeros(1, data['scaled_rbp_expr_df'].shape[1])
        logger.log("[prepare_rbp_tensors] Knockdown reference tensor created.")

    elif reference_type == 'half_reference':
        logger.log("[prepare_rbp_tensors] Using half-zeroed RBP expression as reference.")
        reference_rbp_tensor = torch.ones(1, data['scaled_rbp_expr_df'].shape[1]) * 0.5
        logger.log("[prepare_rbp_tensors] Half-zeroed reference tensor created.")

    else:
        logger.error(f"[prepare_rbp_tensors] Unknown reference type: '{reference_type}'", ValueError)
    return scaled_rbp_tensor, gn_tensor, reference_rbp_tensor

def compute_attribution_scores_for_batch(target_node, deeplift_explainer, scaled_rbp_inputs, reference_rbp_inputs, gene_expression_inputs):
    """
    Compute attribution scores for a batch of input data using the DeepLift explainer.
    
    Parameters:
        target_node (int): Index of the target node or output for which scores are calculated.
        deeplift_explainer (DeepLift): Instance of the DeepLift class used to attribute importance.
        scaled_rbp_inputs (Tensor): Tensor of scaled RBP input data.
        reference_rbp_inputs (Tensor): Baseline tensor (reference) for comparison during attribution.
        gene_expression_inputs (Tensor): Additional input arguments (e.g., gene expression data).

    Returns:
        scores (Tensor): Importance scores attributed to input features.
    """
    scores = deeplift_explainer.attribute(
        inputs=scaled_rbp_inputs,
        baselines=reference_rbp_inputs,
        target=target_node,
        additional_forward_args=gene_expression_inputs,
    )
    return scores

def reduce_batch_dimension(
    list_batch_scores: list, 
    explain_config: dict, 
    logger, 
    trans_id: list, 
    rbps_id: list) -> pd.DataFrame:
    """
    Reduce the batch dimension of DeepLIFT scores to generate final TxRBP scores using a specified method.

    Parameters:
        list_batch_scores (list): List of computed attribution scores for a batch using DeepLift.
        explain_config (dict): Dictionary containing the configuration.
                               Expected key: 'batch_reduction_method' with options:
                               - 't-statistic': Computes the t-statistic across samples.
                               - 'sum_scores': Sums the scores across samples.
        logger (Logger): Instance of the Logger class for logging messages.
        trans_id (list): List of transcript IDs.
        rbps_id (list): List of RBP IDs.

    Returns:
        pd.DataFrame: DataFrame containing the reduced DeepLIFT scores.
    """
    batch_reduction_type = explain_config.get('batch_reduction_method')
    logger.log(f"[reduce_batch_dimension] The reduce method selected is {batch_reduction_type}.")
   
    # Validate the reduction method
    if batch_reduction_type not in ['t-statistic', 'sum_scores']:
        logger.error(f"[reduce_batch_dimension] Unknown batch reduction method: '{batch_reduction_type}'", ValueError)
        raise ValueError(f"Invalid batch reduction method: {batch_reduction_type}")
    
    # Convert all tensors in the list to NumPy arrays upfront
    batch_scores_np = [tensor.detach().numpy() for tensor in list_batch_scores]
    if batch_reduction_type == 't-statistic':
        logger.log("[reduce_batch_dimension] Calculating t-statistic.")
        
        # Stack tensors for efficient calculations
        scores_stack = np.stack(batch_scores_np, axis=0)  # Shape: (num_batches, num_samples, num_features)
        mean = np.mean(scores_stack, axis=1)  # Mean along the sample dimension
        std = np.std(scores_stack, axis=1)   # Std along the sample dimension
        num_samples = scores_stack.shape[1]
        with np.errstate(divide='ignore', invalid='ignore'):
            t_stat_scores = np.where(std != 0, mean / (std / np.sqrt(num_samples)), 0)

        # Detect where std == 0
        zero_std_indices = np.argwhere(std == 0)  # Get indices where std is zero
        if zero_std_indices.size > 0:
            logger.log(f"[reduce_batch_dimension] Found {len(zero_std_indices)} positions with std=0.")
            for idx in zero_std_indices:
                trans = trans_id[idx[0]]  # Get the transcript ID
                rbp = rbps_id[idx[1]]    # Get the RBP ID
                logger.log(f"[reduce_batch_dimension] Zero std at Transcript: {trans}, RBP: {rbp}.")
        
        df_deeplift_TxRBP = pd.DataFrame(t_stat_scores, index=trans_id, columns=rbps_id)
        logger.log("[reduce_batch_dimension] T-statistic reduction completed.")

    elif batch_reduction_type == 'sum_scores':
        logger.log("[reduce_batch_dimension] Calculating sum of scores.")
        scores_stack = np.stack(batch_scores_np, axis=0)  # Shape: (num_batches, num_samples, num_features)
        sum_scores = np.sum(scores_stack, axis=1)  # Sum along the sample dimension
        df_deeplift_TxRBP = pd.DataFrame(sum_scores, index=trans_id, columns=rbps_id)
        logger.log("[reduce_batch_dimension] Sum reduction completed.")
    return df_deeplift_TxRBP

def filter_scores_for_low_expressed_genes(
        deeplift_scores: pd.DataFrame, 
        gene_expr_df: pd.DataFrame, 
        logger,
        threshold: float = 1
        ) -> pd.DataFrame:
    """
    Set low-expressed genes (mean expression < threshold) to 0 in the TxRBP scores DataFrame.

    Parameters:
        deeplift_scores (pd.DataFrame): DataFrame with DeepLIFT scores (TxRBP).
        gene_expr_df (pd.DataFrame): DataFrame with gene expression values.
        threshold (float): Expression threshold (default=1 TPM).

    Returns:
        pd.DataFrame: Updated DeepLIFT scores with low-expressed genes set to 0.
    """
    # Identify genes with mean expression below the threshold
    low_expr_genes = gene_expr_df.mean() < threshold
    # Set their corresponding rows in DeepLIFT scores to 0
    deeplift_scores.loc[low_expr_genes, :] = 0
    logger.log(f"Low-expressed genes (mean expression < {threshold} TPM) have been excluded from the TxRBP scores.")
    return deeplift_scores

### main
logger = Logger(verbose=1)
# this is the old '### Calculate the reference of the RBP expression data'
scaled_rbp_tensor, gn_tensor, reference_rbp_tensor = prepare_rbp_tensors(data, explain_config, logger)

# Step 1: Create an instance of the DeepLift explainer
deeplift_explainer = DeepLift(model)
len_out_features = len(trans_id)

# Suppress specific user warnings
warnings.filterwarnings("ignore", message="Input Tensor .* did not already require gradients")
warnings.filterwarnings("ignore", message="Setting forward, backward hooks and attributes on non-linear activations")

# Step 2: Compute attribution scores batch-wise for RBP x S x T
list_batch_scores = []
for target_node in tqdm(range(len_out_features), desc="Calculating scores", unit="node"):
    scores_batch = compute_attribution_scores_for_batch(
        target_node, deeplift_explainer, scaled_rbp_tensor, reference_rbp_tensor, gn_tensor
    )
    list_batch_scores.append(scores_batch)
    print('[compute_attribution_scores_for_batch] Calculate deeplift scores for all nodes ... -> DONE')

# 2) Reduce batch dimension (RBP x T)
df_deeplift_scores_TxRBP = reduce_batch_dimension(list_batch_scores, explain_config, logger, trans_id[:119], rbps_id)

# 3) Set low-expressed genes's scores (mean < 1TPM) to 0.
df_deeplift_scores_TxRBP = filter_scores_for_low_expressed_genes(
        deeplift_scores = df_deeplift_scores_TxRBP, 
        gene_expr_df = data["gene_expr_df"],
        threshold = 1,
        logger = logger)




# 4) Collapse scores to genes (RBP x G) ME HE QUEDADO AQUI BROTHER !!!
#### AQUII !!!

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