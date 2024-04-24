# create_data.py
import timeit 
import datetime
import argparse
import os

# print the current working directory
#print('Original working directory:', os.getcwd())
# change the working directory
#os.chdir('/Users/joseba/Downloads/ML4BM-Lab2/DeepRBP/model/utils')
# print the updated working directory
#print('Updated working directory:', os.getcwd())

from utils.Utils import check_create_new_directory
import utils.Generate_Datasets as gd


def create_data(args, path_data):
    """
    Function to create training and test datasets from raw data.
    
    Args:
    - args (argparse.Namespace): Arguments parsed from command line.
    - path_data (str): Path to the data directory.
    
    Returns:
    None
    """
    check_create_new_directory(path_data+'processed/splitted_datasets')
    print('[create_data] Creating the Training and Test folders and dividing the samples per tumor_type')
    gd.create_datasets(path_data=path_data, chunksize=args.chunksize, select_genes = args.select_genes)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Datasets Tool')
    parser.add_argument('--chunksize', 
                    type=int, 
                    help='Number of samples to process per chunk. '
                         'This determines the size of each subset of patient data that will be processed at a time. '
                         'The processing steps include the generation of the RBP expresion, gene expression and '
                         'transcript expression matrices : '
                    
                        '  1) Combine different versions of a gene or loci. '
                        '  2) Convert gene expression matrix to TPMs and apply log2(tpm+1) transformation for transcript expression matrix.'
                        '  3) Categorize data based on tumor type or tissue. '
                        '  4) Eliminate genes that have only one transcript. '
                        '  5) Choose genes according to specific criteria, such as those related to cancer. '
                        '  6) Create the RBP matrix using a list of identified RBPs and transform the data to log2 scale. '
                        '  7) Create a reduced getBM data relating transcript ids to its gene for the selected genes for modelling '
                        '  8) Split the data into categories, reserving 20*%* for testing. '

                         'Larger chunk sizes may require more memory but can speed up processing.')
    parser.add_argument('--select_genes', 
                    type=str, 
                    help='Genes selected for modelling. '
                         'Provide a list of specific genes to be used in the model. '
                         'If not specified, all available genes will be used.')
    parser.add_argument('--path_data', 
                    type=str, 
                    help='Path to the data directory. '
                         'This should point to the directory where the raw and processed data files are stored. '
                         'The tool will read from and write to this directory.')
    args = parser.parse_args()
    
    print(f'path_data is: {args.path_data}')
    start_time = timeit.default_timer()
    
    create_data(args, args.path_data)
    duration = datetime.timedelta(seconds=timeit.default_timer() - start_time)
    print('[create_data] Create the Data ... -> DONE')
    print('[create_data] The time of execution is:', str(duration))