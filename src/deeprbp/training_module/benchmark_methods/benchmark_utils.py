# src/deeprbp/training_module/benchmark_methods/benchmark_utils.py

import numpy as np
import random
import os
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, ElasticNet, Ridge
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.multioutput import MultiOutputRegressor

def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)

def generate_X_y_data(data_np, calculate_abundance=True):
    """
    Extract features and target data from the input dataset.

    Parameters:
    - data_np (dict): A dictionary containing the dataset with the following keys:
        - 'scaled_rbp_df' (numpy.ndarray): A 2D array of shape (n_samples, n_features)
          containing the log2p RBP expression data scaled and transformed in a range [0,1].
        - 'isoform_df' (numpy.ndarray): A 2D array of shape (n_samples, n_isoforms)
          containing the target isoform expression data in log2tpm+1.
        - 'gene_df' (numpy.ndarray): A 2D array of shape (n_samples, n_genes)
          containing gene expression data in TPM.
          
    - calculate_abundance (bool): If True, calculate the target `y` as abundances using
      the formula (2 ** transcript_data - 1) / gene_data, setting y to 0 where gene_data is 0.
      If False, set `y` directly to transcript_data.

    Returns:
    - X_data (numpy.ndarray): The feature data as a 2D array.
    - y (numpy.ndarray): The target data as a 2D array.
    """
    # Extract the RBP, isoform, and gene data
    X_data = data_np['scaled_rbp_df']  # RBP Features
    transcript_data = data_np['isoform_df']  # Target
    gene_data = data_np['gene_df']  # Gene data
    # Calculate y based on the calculate_abundance flag
    if calculate_abundance:
        # Calculate y, setting y to 0 where gene_data is 0
        print("Calculating y as abundances using the formula: (2 ** transcript_data - 1) / gene_data")
        tpm_transcript = 2 ** transcript_data - 1
        y = np.zeros_like(tpm_transcript, dtype=np.float32)
        np.divide(tpm_transcript, gene_data, out=y, where=gene_data != 0)
        #y = np.where(gene_data == 0, 0, (2 ** transcript_data - 1) / gene_data)
    else:
        print("Setting y directly to transcript_data.")
        y = transcript_data
    return X_data, y

# podriamos probar de aquí la que no hemos probado: 
# vamos a probar dos approaches: due possibili metodi
        # - o calcular directamente la expression del transcrito (sin usar la expresion de los genes para entrenar)
        # - o calcular el isoforma abundance como nuestro modelo y luego posteriormente ya multiplicamos por la expresion del gen
# en ambos casos al final calculamos las métricas que calculamos tb para nuestro modelo 
# en log2tpm+1

def create_multi_output_regressor(selected_algorithm: str, **kwargs):
    """
    Create a MultiOutputRegressor with the specified regression algorithm.

    Parameters:
    selected_algorithm (str): The name of the regression algorithm to use. 
                              Available options are:
                              - 'linear_regression'
                              - 'svr'
                              - 'decision_tree'
                              - 'random_forest'
                              - 'gradient_boosting'
                              - 'xgboost'
                              - 'lightgbm'
                              - 'knn'
                              - 'elastic_net'
                              - 'ridge'
    
    **kwargs: Additional keyword arguments to pass to the regression model. 
               This can include parameters specific to the chosen model class.
    
    Returns:
    MultiOutputRegressor: An instance of MultiOutputRegressor wrapping the specified regression model.
    
    Raises:
    ValueError: If an unsupported algorithm name is provided.
    
    Example:
    >>> model = create_multi_output_regressor('random_forest', n_estimators=100)
    """
    # Select the base model based on the specified algorithm
    if selected_algorithm == 'svr':
        model = SVR(kernel='rbf', **kwargs)   
    elif selected_algorithm == 'decision_tree':
        model = DecisionTreeRegressor(**kwargs)
    elif selected_algorithm == 'random_forest':
        model = RandomForestRegressor(**kwargs)
    elif selected_algorithm == 'gradient_boosting':
        model = GradientBoostingRegressor(**kwargs)
    elif selected_algorithm == 'xgboost':
        model = XGBRegressor(objective='reg:squarederror', **kwargs)
    elif selected_algorithm == 'lightgbm':
        model = LGBMRegressor(**kwargs)
    elif selected_algorithm == 'knn':
        model = KNeighborsRegressor(**kwargs)
    elif selected_algorithm == 'elastic_net':
        model = ElasticNet(**kwargs)
    elif selected_algorithm == 'ridge':
        model = Ridge(**kwargs)
    else:
        raise ValueError(f"Unsupported algorithm: {selected_algorithm}")
    multioutput_model = MultiOutputRegressor(model)
    return multioutput_model

