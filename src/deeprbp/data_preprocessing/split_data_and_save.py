# /src/deeprbp/data_preprocessing/split_data_and_save.py

import os
import argparse
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DataSplitter
from ..util.utils import select_sample_ids_by_type, filter_data_by_sample_ids, save_data

def run_split_data_and_save(config_parser, output_dir: str):
    """
    Executes the data loading, filtering, splitting, and saving steps.

    Parameters:
    - config_parser (ConfigParser): A ConfigParser object with loaded config.
    - output_dir (str): Directory where the Train/Test datasets will be saved.
    """
    # Load the data using the data importer
    print("Loading data...")
    data_importer = DataImporter(config_parser.get('path_files'))
    data = data_importer.load()
    print("Data loaded successfully")
    # Filter data by selecting the samples based on the provided configuration
    print("Filter data by selecting the samples based on the provided configuration ...")
    sel_sample_ids = select_sample_ids_by_type(
            metadata_df=data['metadata_df'],
            sample_category_col=config_parser.get('sample_category'),
            sample_types=config_parser.get('select_samples')
        )
    # Filter the data by the selected samples
    data = filter_data_by_sample_ids(data, sel_sample_ids)
    print('Data filtered successfully')
    # Split the data into training and test sets
    print("Splitting data into training and test sets...")
    splitter = DataSplitter(data, config_parser.get('sample_category'))
    train_data, test_data = splitter.split_data_sets(config_parser.get('test_fraction'))
    print("Data splitted successfully.")
    # Save the training and test datasets to the specified output directory
    print(f"Saving training data to: {os.path.join(output_dir, 'Train')}")
    save_data(train_data, 
              os.path.join(output_dir, 'Train'), 
              custom_names = {
                            'rbp_df': 'RBPs_log2p_tpm.csv',
                            'isoform_df': 'trans_log2p_tpm.csv',
                            'gene_df': 'gn_tpm.csv',
                            'metadata_df': 'phenotype_metadata.csv'}
                            )
    print("Training data saved successfully.")
    print(f"Saving test data to: {os.path.join(output_dir, 'Test')}")
    save_data(test_data, 
              os.path.join(output_dir, 'Test'), 
              custom_names = {
                            'rbp_df': 'RBPs_log2p_tpm.csv',
                            'isoform_df': 'trans_log2p_tpm.csv',
                            'gene_df': 'gn_tpm.csv',
                            'metadata_df': 'phenotype_metadata.csv'}
                            )
    print("Test data saved successfully.")

def main():
    args = parse_args()
    print("Loading configuration from:", args.config_path)
    config_parser = ConfigParser(args.config_path)
    print("Configuration loaded successfully.")
    run_split_data_and_save(config_parser=config_parser, output_dir=args.output_dir)

def parse_args():
    parser = argparse.ArgumentParser(description=(
        'This script loads configuration settings from a specified YAML file, '
        'imports the necessary data, select desired tumor types and splits it into training and test datasets. '
        'The resulting datasets are then saved to the specified output directory. '
        'The configuration file should include paths to the data files, sample selection criteria, '
        'tumor types (categories), train-test split fraction, and other necessary parameters.'
    ))

    parser.add_argument('--config_path', type=str, default='/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_data_split.yaml', 
        help='Path to the config file with the processed data files, sample selection, tumor types (categories), train-test fraction and source name.')
    parser.add_argument('--output_dir', type=str, default='/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets', 
                        help='Directory to save the splitted datasets.')
    return parser.parse_args()

if __name__ == '__main__':
    main()

