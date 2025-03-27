# /src/deeprbp/training_module/hyperparameter_optimization/grid_search_optuna.py

import os
import argparse
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DataSplitter

# Define the default configuration path
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml'

#args = parse_args()

# Load configuration settings from the specified YAML file
print("Loading configuration from:", config_path)
config_parser = ConfigParser(config_path)
config = config_parser.load_config()
print("Configuration loaded successfully.")

# Load the data using the data importer
print("Loading data...")
data_importer = DataImporter(config['data_paths'])
data = data_importer.load()
print("Data loaded successfully")

# Filter half the samples to optimize time and computational resources
# train_idx, _ = DataSplitter.split_data_class(data=data, 
#                                              config=config, 
#                                              sample_category='detailed_category', 
#                                              test_size=0.5)
# # Filter the data by the selected samples
# train_data = filter_data_by_sample_ids(data, train_idx)

#   # For performing a train/validation split
#     train_data, valid_data = splitter.split_data_sets(test_name='validation')