# src/deeprbp/training_module/benchmark_methods/utils_benchmark.py

from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, ElasticNet, Ridge
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.multioutput import MultiOutputRegressor
import numpy as np

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
        y = np.where(gene_data == 0, 0, (2 ** transcript_data - 1) / gene_data)
        print("Calculating y as abundances using the formula: (2 ** transcript_data - 1) / gene_data")
    else:
        y = transcript_data
        print("Setting y directly to transcript_data.")
    return X_data, y

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
    if selected_algorithm == 'linear_regression':
        model = LinearRegression(**kwargs)
    elif selected_algorithm == 'svr':
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
