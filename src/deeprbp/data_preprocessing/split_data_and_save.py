# /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/split_data_and_save.py

import os
import argparse
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DataSplitter
from ..util.utils import select_sample_ids_by_type, filter_data_by_sample_ids, save_data

# Define the default configuration path
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_data_split.yaml'

def main():
    args = parse_args()

    # Load configuration settings from the specified YAML file
    print("Loading configuration from:", args.config_path)
    config_parser = ConfigParser(args.config_path)
    config = config_parser.load_config()
    print("Configuration loaded successfully.")

    # Load the data using the data importer
    print("Loading data...")
    data_importer = DataImporter(config['data_paths'])
    data = data_importer.load()
    print("Data loaded successfully")

    # Filter data by selecting the samples based on the provided configuration
    print("Filter data by selecting the samples based on the provided configuration ...")
    sel_sample_ids = select_sample_ids_by_type(
            metadata_df=data['metadata_df'],
            sample_category_col=config['sample_category'],
            sample_types=config['select_samples']
        )
    # Filter the data by the selected samples
    data = filter_data_by_sample_ids(data, sel_sample_ids)
    print('Data filtered successfully')

    # Split the data into training and test sets
    print("Splitting data into training and test sets...")
    splitter = DataSplitter(data, config)
    train_data, _, test_data = splitter.split_data_sets()
    print("Data split successfully.")

    # Save the training and test datasets to the specified output directory
    print(f"Saving training data to: {os.path.join(args.output_dir, 'Train')}")
    save_data(train_data, os.path.join(args.output_dir, 'Train'))
    print("Training data saved successfully.")

    print(f"Saving test data to: {os.path.join(args.output_dir, 'Test')}")
    save_data(test_data, os.path.join(args.output_dir, 'Test'))
    print("Test data saved successfully.")

def parse_args():
    parser = argparse.ArgumentParser(description=(
        'This script loads configuration settings from a specified YAML file, '
        'imports the necessary data, select desired tumor types and splits it into training and test datasets. '
        'The resulting datasets are then saved to the specified output directory. '
        'The configuration file should include paths to the data files, sample selection criteria, '
        'tumor types (categories), train-test split fraction, and other necessary parameters.'
    ))

    parser.add_argument('--config_path', type=str, default=config_path, 
        help='Path to the config file with the processed data files, sample selection, tumor types (categories), train-test fraction and source name.')
    parser.add_argument('--output_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets', 
                        help='Directory to save the splitted datasets.')
    return parser.parse_args()

if __name__ == '__main__':
    main()