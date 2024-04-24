import os
import re
import json
import joblib
from glob import glob
import pandas as pd
import numpy as np
from collections import namedtuple
import torch
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import timeit
import datetime
from scipy import stats
from scipy.stats import linregress, spearmanr, pearsonr
from collections import namedtuple
from sklearn.metrics import roc_curve, auc

#ruta_utils = '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/utils'
ruta_utils = './../utils'
os.environ['PYTHONPATH'] = ruta_utils + ':' + os.environ.get('PYTHONPATH', '')
from utils import Plots as plots
#from utils.Config import Config
#ruta_modelsNN = '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/modelsNN'
#os.environ['PYTHONPATH'] = ruta_modelsNN + ':' + os.environ.get('PYTHONPATH', '')
#from modelsNN.modelsNN import DeepRBP

#ruta_utils = '/Users/joseba/Downloads/ML4BM-Lab2/DeepRBP/utils'
#ruta_modelsNN = '/Users/joseba/Downloads/ML4BM-Lab2/DeepRBP/modelsNN'
#os.environ['PYTHONPATH'] = ':'.join([ruta_utils, ruta_modelsNN]) + ':' + os.environ.get('PYTHONPATH', '')
#from utils import Plots as plots
#from utils.Config import Config
#from modelsNN.modelsNN import DeepRBP

def check_create_new_directory(path):
    """ 
    Function that checks if a directory is created and if not create the new directory
    Parameters
    ----------
    path: string with the path
    """
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print(f'[Utils] The new directory {path} is created!')

def check_if_files_exist(folder_path):
    for root, dirs, files in os.walk(folder_path):
        if files:
            return True  
    return False 

#####################################  DeepRBP predictor modelling ################################

def check_existing_trained_model(use_trained_model, path_save_files, config):
    """
    Checks if there is an existing trained model with the same settings.
    
    Args:
    - use_trained_model (bool): Indicates whether to use a pre-trained model.
    - path_save_files (str): Path to the directory where model files are saved.
    
    Raises:
    - ValueError: If there is no model saved with the current settings and use_trained_model is True.
    - FileNotFoundError: If the directory specified by path_save_files does not exist.
    """
    trained_model_available = False
    if use_trained_model:
        if os.path.exists(path_save_files):
            print(f"[Utils] The directory {path_save_files} already exists.")
            print("[Utils] Checking if there is a model saved with the actual settings")
            config_file_path = os.path.join(path_save_files, 'config.json')
            model_file_path = os.path.join(path_save_files, 'model.pt')

            # Check if both config.json and model.pt exist
            if os.path.exists(config_file_path) and os.path.exists(model_file_path):

                # Load config from the existing model
                with open(config_file_path, 'r') as config_file:
                    config_data_trained_model = json.load(config_file)
                
                keys_to_compare = ['max_node', 'num_hidden_layers', 'learning_rate', 'batch_size', 'epochs',
                                    'same_num_nodes', 'node_divition_factor', 'activation_layer', 'optimizer',
                                    'train_tumor_types', 'test_tumor_types']
                config_data_trained_model_subset = {key: config_data_trained_model[key] for key in keys_to_compare}
                config_subset = {key: config.get_config()[key] for key in keys_to_compare}
                
                if config_data_trained_model_subset == config_subset:
                    print("[Utils] We found an already trained model with the actual settings")
                    trained_model_available = True
            else:
                raise ValueError("[Utils] Either config.json or model.pt is missing. You need to set use_trained_model to False")
        else:
            raise FileNotFoundError(f"[Utils] The directory {path_save_files} does not exist.")
    else:
        # We are going to create the .json for the model we're gonna train with the actual settings
        check_create_new_directory(path_save_files)
        with open(os.path.join(path_save_files, 'config.json'), 'w') as file:
            json.dump(config.get_config(), file)
    return trained_model_available

def split_training_test_sets(df_sfs, df_trans, df_gns_each_trans, test_size=0.2):
    """
    Split the dataset into training and test sets.

    Args:
        df_sfs (DataFrame): DataFrame containing RNA-binding proteins (RBPs) data.
        df_trans (DataFrame): DataFrame containing transcript (labels) data.
        df_gns_each_trans (DataFrame): DataFrame containing gene expression data for each isoform.
        test_size (float, optional): Proportion of the dataset to include in the test split. Defaults to 0.2.

    Returns:
        DataSplit: Namedtuple containing the split datasets.
            - df_train (DataFrame): Training dataset.
            - train_labels (DataFrame): Labels for the training dataset.
            - train_gn (DataFrame): Gene expression data for the training dataset.
            - df_test (DataFrame): Test dataset.
            - test_labels (DataFrame): Labels for the test dataset.
            - test_gn (DataFrame): Gene expression data for the test dataset.
    """
    df_train, df_test = train_test_split(df_sfs, test_size=test_size, random_state=0)  
    train_labels = df_trans.loc[df_train.index]
    test_labels = df_trans.loc[df_test.index]

    train_gn = df_gns_each_trans.loc[df_train.index]
    test_gn = df_gns_each_trans.loc[df_test.index]
    DataSplit = namedtuple('DataSplit', ['df_train', 'train_labels', 'train_gn',
                                         'df_test', 'test_labels', 'test_gn'])
    data_split = DataSplit(df_train, train_labels, train_gn, df_test, test_labels, test_gn)          
    return data_split   

def get_data(config, path, set_mode='training', toy_set=False, frac = None, perform_train_test_split=True):  
    """
    This function retrieves and preprocesses data. It selects data based on tumor types and sources, reads the relevant files, 
    and organizes them into data frames. It can also perform train-test splits if needed. Finally, it returns the processed data.

    Args:
        config (object): Object containing tumor types and source configurations.
        path (str): Path to the directory containing the data.
        set_mode (str, optional): Mode of the data ('training' or 'test'). Defaults to 'training'.
        toy_set (bool, optional): Whether to use a subset of the samples for training. Defaults to False.
        frac (float, optional): Proportion of samples to use if toy_set is True. Defaults to None.
        perform_train_test_split (bool, optional): Whether to perform train-test split. Defaults to True.

    Returns:
        Namedtuple: Namedtuple containing the processed data.
            - df_train (DataFrame): Training dataset.
            - df_val (DataFrame): Validation dataset if perform_train_test_split is True.
            - train_labels (DataFrame): Labels for the training dataset.
            - valid_labels (DataFrame): Labels for the validation dataset if perform_train_test_split is True.
            - train_gn (DataFrame): Gene expression data for the training dataset.
            - valid_gn (DataFrame): Gene expression data for the validation dataset if perform_train_test_split is True.
            - df_raw (DataFrame): Test dataset if set_mode is 'test'.
            - raw_labels (DataFrame): Labels for the test dataset if set_mode is 'test'.
            - raw_gn (DataFrame): Gene expression data for the test dataset if set_mode is 'test'.
    """

    ### 0) Set the set type: training data or test data
    if set_mode == 'training':
        select_tumor_types = config.train_tumor_types
        select_source = config.source_train
        # Initialize result dataframe
        if perform_train_test_split:
            df_train, df_val, train_labels, valid_labels, train_gn, valid_gn = [pd.DataFrame() for _ in range(6)]
        else:
            df_train, train_labels, train_gn = [pd.DataFrame() for _ in range(3)]

    else: # set = test
        select_tumor_types = config.test_tumor_types
        select_source = 'all' #config.source_pred # We are testing in both TCGA and GTEX (aunque obviamente separadamente)
        # Initialize result dataframe
        df_raw, raw_labels, raw_gn = [pd.DataFrame() for _ in range(3)]
    print(f'[Utils][get_data] set_mode is {set_mode}')

    ### 1) Select tumor types from config
    if select_tumor_types == 'all': # We ARE selecting all tumor types
        if select_source == 'TCGA' or select_source == 'GTEX':
            path_raw = glob(f'{path}/splitted_datasets/{select_source}/{set_mode}/*')
        else: # else = config.source == 'all':
            path_raw_tcga = glob(f'{path}/splitted_datasets/TCGA/{set_mode}/*')
            path_raw_gtex = glob(f'{path}/splitted_datasets/GTEX/{set_mode}/*')
            path_raw = path_raw_tcga+path_raw_gtex
    
    else: # We are not selecting all tumor types
        print('[Utils][get_data] We are not selecting all the tumor types')
        tumor_types = select_tumor_types.split(',')
        print(tumor_types)
        if select_source == 'TCGA' or select_source == 'GTEX':
            path_all = glob(os.path.join(path, "splitted_datasets", select_source, set_mode, "*"))

            if select_source == 'TCGA':
                print('The source is TCGA')
                path_raw = [p for p in path_all if os.path.basename(p) in tumor_types]
            if select_source == 'GTEX':
                print('The source is GTEX')
                tumor_types = list(set([val.split("_")[0] for val in tumor_types]))
                path_raw = [p for p in path_all if os.path.basename(p).split("_")[0] in tumor_types]
            
        else: # else = config.source == 'all':
            path_all_tcga = glob(f'{path}/splitted_datasets/TCGA/{set_mode}/*')
            path_all_gtex = glob(f'{path}/splitted_datasets/GTEX/{set_mode}/*')
            path_raw_tcga = [p for p in path_all_tcga if os.path.basename(p) in tumor_types]
            path_raw_gtex = [p for p in path_all_gtex if os.path.basename(p).split("_")[0] in list(set([val.split("_")[0] for val in tumor_types]))]
            path_raw = path_raw_tcga+path_raw_gtex

    ### 2) Get actually the data from the path_raw list:
    for i in path_raw:
        print('[Utils][get_data] path_raw: ', i)
        file_list = os.listdir(i)
        temp_rbps_raw = pd.read_csv(f'{i}/'+ next((f for f in file_list if re.compile(r'RBP').search(f)), None), index_col=0)
        temp_tpm_raw = pd.read_csv(f'{i}/'+ next((f for f in file_list if re.compile(r'trans').search(f)), None), index_col=0)
        temp_gn_raw = pd.read_csv(f'{i}/'+ next((f for f in file_list if re.compile(r'gn_expr').search(f)), None), index_col=0)

        # Add column with the name of the tumor type and the source of the sample.               
        temp_rbps_raw['tumor_type'] = i.split('/')[-1]
        temp_rbps_raw['source'] = temp_rbps_raw.index.str.split("-").str[0]
        print('[Utils][get_data] tumor_type & source columns added to this set ... -> DONE')

        if set_mode == 'training':
            if toy_set: # If we are using a toy set with a subset of the samples per tumor type
                print(f'[Utils][get_data] You are taking {int(frac*100)}% of the samples for {i.split("/")[-1]} tumor type')
                temp_rbps_raw = temp_rbps_raw.sample(frac=frac, random_state=33)
                temp_tpm_raw = temp_tpm_raw.sample(frac=frac, random_state=33)
                temp_gn_raw = temp_gn_raw.sample(frac=frac, random_state=33)

            if perform_train_test_split:
                print(f'[Utils][get_data] You are performing a train-test split')
                # In Training mode Divide Raw Data in Training/Validation for each tumor type
                data_split = split_training_test_sets(df_sfs=temp_rbps_raw, df_trans=temp_tpm_raw, df_gns_each_trans=temp_gn_raw)
                print(f'[Utils][get_data] Divide data in {i.split("/")[-1]} in Training/Validation ... -> DONE')
                                                                       
                df_train = pd.concat([df_train, data_split.df_train], axis = 0)
                train_labels = pd.concat([train_labels, data_split.train_labels], axis = 0)
                train_gn = pd.concat([train_gn, data_split.train_gn], axis = 0)
                df_val = pd.concat([df_val, data_split.df_test], axis = 0)
                valid_labels = pd.concat([valid_labels, data_split.test_labels], axis = 0)
                valid_gn = pd.concat([valid_gn, data_split.test_gn], axis = 0)

            else:
                print(f'[Utils][get_data] You are using all raw data to train')
                df_train = pd.concat([df_train, temp_rbps_raw], axis = 0)
                train_labels = pd.concat([train_labels, temp_tpm_raw], axis = 0)
                train_gn = pd.concat([train_gn, temp_gn_raw], axis = 0)

        else: # set_mode == 'test':
            df_raw = pd.concat([df_raw, temp_rbps_raw], axis = 0)
            raw_labels = pd.concat([raw_labels, temp_tpm_raw], axis = 0)
            raw_gn = pd.concat([raw_gn, temp_gn_raw], axis = 0)

    # 3) Return result ordered
    if set_mode == 'training':
        df_train = df_train.sort_index(axis=1)
        train_labels = train_labels.sort_index(axis=1)
        train_gn = train_gn.sort_index(axis=1)

        if perform_train_test_split == False:
            DataRaw = namedtuple('DataRaw', ['df_train', 'train_labels', 'train_gn'])
            data_raw = DataRaw(df_train, train_labels, train_gn)
        else:
            df_val = df_val.sort_index(axis=1)
            valid_labels = valid_labels.sort_index(axis=1)
            valid_gn = valid_gn.sort_index(axis=1)
            DataRaw = namedtuple('DataRaw', ['df_train', 'df_val', 'train_labels', 'valid_labels', 'train_gn', 'valid_gn'])
            data_raw = DataRaw(df_train, df_val, train_labels, valid_labels, train_gn, valid_gn)

    else:
        df_raw = df_raw.sort_index(axis=1)
        raw_labels = raw_labels.sort_index(axis=1)
        raw_gn = raw_gn.sort_index(axis=1)
    
        DataRaw = namedtuple('DataRaw', ['df_raw', 'raw_labels', 'raw_gn'])
        data_raw = DataRaw(df_raw, raw_labels, raw_gn)
    return data_raw

def scale_rbp_input_data(df_train, df_validation):
    """ 
    Function that scales the training and validation set expression applying a standard scale, then clips the max/min 
    values in +-2sigma and then the data values between 0 and 1.
    Parameters
    ----------
    df_train: df, train set
    df_validation: df, validation set

    Returns: data_scaled, named tuple.
    -------
    - scaledTrain_df: df, scaled training set
    - scaledValidation_df: df, scaled validation set
    - scaler: the scale object
    - sigma: the standard deviation used in the scaling process
    """
    scaler_sfs = StandardScaler()  # Initialize
    # We put the content inside the scaler. For each feature mean and std.
    scaler_sfs.fit(df_train)
    scaledTrain_df = pd.DataFrame(scaler_sfs.transform(df_train), index=df_train.index, columns=df_train.columns)
    scaledValidation_df = pd.DataFrame(scaler_sfs.transform(df_validation), index=df_validation.index, columns=df_validation.columns)
    sigma_sfs = np.std(scaledTrain_df.values)
    # Fijamos los valores menores de -2sigma_sfs en -2sigma_sfs y los mayores de +2sigma_sfs en +2sigma_sfs
    scaledTrain_df = np.clip(scaledTrain_df, -2*sigma_sfs, 2*sigma_sfs, axis=1)
    scaledValidation_df = np.clip(scaledValidation_df, -2*sigma_sfs, 2*sigma_sfs, axis=1)
    # Después, desplazamos todos los datos 2sigma_sfs sumándoles 2
    scaledTrain_df += 2*sigma_sfs
    scaledValidation_df += 2*sigma_sfs
    # Finalmente, dividimos todos los datos entre 4sigma_sfs
    scaledTrain_df /= 4*sigma_sfs
    scaledValidation_df /= 4*sigma_sfs
    DataScaled = namedtuple('DataScale', ['scaledTrain_df', 'scaledValidation_df', 'scaler_sfs', 'sigma_sfs'])
    data_scaled = DataScaled(scaledTrain_df, scaledValidation_df, scaler_sfs, sigma_sfs)
    return data_scaled 

def create_data_loader(scaled_rbps_df, labels_df, gn_df, config, set_mode='training'):
    """
    Function to create data loaders for training or testing.

    Args:
        scaled_rbps_df (DataFrame): Scaled RNA-binding proteins (RBPs) data.
        config (Config): Configuration object containing model settings.
        labels_df (DataFrame): Labels data.
        gn_df (DataFrame): Gene data.
        set_mode (str): Mode for data loader, 'training' or 'test'.

    Returns:
        DataLoader: PyTorch data loader.
    """
    # Convert to PyTorch dataset
    ds = TensorDataset(torch.tensor(scaled_rbps_df.values, dtype=torch.float32),
                        torch.tensor(labels_df.values, dtype=torch.float32),
                        torch.tensor(gn_df.values, dtype=torch.float32))
    if set_mode == 'training':
        print('[Utils][create_data_loader] Creating a training data_loader')
        data_loader = DataLoader(ds, config.batch_size, shuffle=True, drop_last=True, pin_memory=True)
    else:
        print('[Utils][create_data_loader] Creating a test data_loader')
        data_loader = DataLoader(ds, config.batch_size * 2, drop_last=True, pin_memory=True)
    print('[Utils][create_data_loader] Convert to PyTorch dataset ... -> DONE ')
    return data_loader

def get_training_validation_loaders(path_data, config, frac=None, toy_set=False, return_data_scale=False):
    """
    Function to prepare training and validation data loaders.

    Args:
        path_data (str): Path to the data.
        config (Config): Configuration object containing model settings.
        frac (float): Fraction of data to use for training.
        toy_set (bool): Flag to indicate whether to use a toy dataset.
        return_data_scale (bool): Flag to indicate whether to return data scaling information.

    Returns:
        Data_Loader: Namedtuple containing training and validation data loaders.
    """

    ### 1) Read Raw Data and Divide Raw Data in Training/Validation
    data_raw = get_data(config=config, path=path_data, toy_set=toy_set, frac=frac)
    print('[Utils][get_data_loaders] Read Raw Data and Divide Raw Data in Training/Validation ... -> DONE')
    df_train = data_raw.df_train
    df_val = data_raw.df_val
    df_train = df_train.drop(['source', 'tumor_type'], axis=1)
    df_val = df_val.drop(['source','tumor_type'], axis=1)
    print(f'[Utils][get_data_loaders] Shape of Training Data: {df_train.shape}')
    print(f'[Utils][get_data_loaders] Shape of Validation Data: {df_val.shape}')
    train_labels = data_raw.train_labels
    valid_labels = data_raw.valid_labels
    train_gn = data_raw.train_gn
    valid_gn = data_raw.valid_gn

    ### 2) Scale the RBP input data
    data_scale = scale_rbp_input_data(df_train, df_val)
    print('[Utils][get_data_loaders] Scale the SF input data ... -> DONE')
    print(data_scale.scaledTrain_df) 
    print('rbps:', list(data_scale.scaledTrain_df.columns)[0:33])
    print('trans:', list(train_labels.columns)[0:33])
    print('genes:', list(train_gn.columns)[0:33])
    
    ### 3) Create Training/Validation data loaders
    train_loader = create_data_loader(scaled_rbps_df=data_scale.scaledTrain_df, labels_df=train_labels, gn_df=train_gn, config=config)
    val_loader = create_data_loader(scaled_rbps_df=data_scale.scaledValidation_df, labels_df=valid_labels, gn_df=valid_gn, config=config, set_mode='test')
    
    if return_data_scale:
        print('[Utils][get_data_loaders] Returning the data scale') 
        Data_Loader = namedtuple('Data_Loader', ['train_loader', 'val_loader', 'data_scale'])
        data_loader = Data_Loader(train_loader, val_loader, data_scale)
    else:
        Data_Loader = namedtuple('Data_Loader', ['train_loader', 'val_loader'])
        data_loader = Data_Loader(train_loader, val_loader)  
    return data_loader

def build_optimizer(model, optimizer, learning_rate):
    """ 
    Function that creates the optimizer object
    
    Parameters
    ----------
    model
    optimizer: string, name of the selected optimizer
    learning_rate
    
    Returns: optimizer object
    -------
    """ 
    if optimizer == 'sgd90':
        optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate, momentum=0.9)
    elif optimizer == 'sgd70':
        optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate, momentum=0.7)
    elif optimizer == 'sgd50':
        optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate, momentum=0.5)
    elif optimizer == 'asgd':
        optimizer = torch.optim.ASGD(model.parameters(), lr=learning_rate, lambd=0.0001, alpha=0.75)
    elif optimizer == 'adam':
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    elif optimizer == 'adagrad':
        optimizer = torch.optim.Adagrad(model.parameters(), lr=learning_rate)
    elif optimizer == 'adadelta':
        optimizer = torch.optim.Adadelta(model.parameters(), lr=learning_rate)
    elif optimizer == 'adamW':
        optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate)
    return optimizer

def load_or_train_model(trained_model_available, model, path_save_files, device, path_data, config, optimizer=None, save_model=None):
    """
    Load or train a model based on availability of pre-trained model.

    Args:
        model: PyTorch model object.
        path_save_files (str): Path to directory where model files are saved.
        device: Torch device (CPU or GPU).
        path_data (str): Path to the data directory.
        config: Configuration object containing model settings.
        save_model (bool): Whether to save the trained model. Defaults to True.

    Returns:
        model: PyTorch model object (loaded or trained).
        data_scale: Namedtuple containing scaler and sigma used in model training.
    """
    if trained_model_available:
        print(f'[Utils] Warning: trained_model_available = {trained_model_available}')
        print('[Utils] Loading the weights and scaler of the already trained models')
        model.load_state_dict(torch.load(path_save_files+'/model.pt', map_location=device))      
        model.eval()
        print('[Utils] Loading the weights of the already trained models ... -> DONE\n')
        # Load the scaler and sigma used in the model training
        scaler_sfs = joblib.load(path_save_files+'/scaler_sfs.joblib')
        with open(path_save_files+'/sigma_sfs.txt', 'r') as f:
            sigma_sfs = f.readline().strip()
        sigma_sfs = np.float128(sigma_sfs)
        Data_Scale = namedtuple('Data_Scale', ['scaler_sfs', 'sigma_sfs'])
        data_scale = Data_Scale(scaler_sfs, sigma_sfs)
        print('[Utils] Load the scaler and sigma used in the model training ... -> DONE')
        print('[Utils] Warning: We are using the original scaling to scale the test data')
    else:
        print(f'[Utils] Warning: trained_model_available = {trained_model_available}')
        print('[Utils] There is not an already trained model')
        print('[Utils] Get the Training & Validation data loaders and the scaler to train the model ...')
        data_loader = get_training_validation_loaders(path_data, config, return_data_scale=True)
        train_loader = data_loader.train_loader
        val_loader = data_loader.val_loader
        data_scale = data_loader.data_scale
        if save_model:
            print(f'[Utils] The model weights and training configuration is saved in: {path_save_files}')
        optimizer = build_optimizer(model, config.optimizer, config.learning_rate)
        #### Train the model
        start_time = timeit.default_timer()
        model.fit(
                epochs=config.epochs,     
                train_loader=train_loader, 
                val_loader=val_loader, 
                optimizer=optimizer,      
                save_model=save_model, 
                path=path_save_files, 
                model_name='model.pt'
                )
        duration = datetime.timedelta(seconds=timeit.default_timer() - start_time)
        print('[Utils] Train the model ... -> DONE')
        print('[Utils] The time of execution is:', str(duration))
        ### Save the scaler and sigma used in Training to a file
        filename_scaler = os.path.join(path_save_files, 'scaler_sfs.joblib')
        filename_sigma = os.path.join(path_save_files, 'sigma_sfs.txt')
        joblib.dump(data_scale.scaler_sfs, filename_scaler)
        np.savetxt(filename_sigma, [data_scale.sigma_sfs])
        print(f'[Utils] scaler and sigma used in Training saved in {path_save_files}') 
        print('[Utils] Save the scaler and sigma used in Training to a file ... -> DONE\n')
    return model, data_scale

#####################################  DeepRBP predictor model's performance evaluation ################################

def do_predictions(config, data_test, data_scale, model, device, pred_source, path_getBM, path_save=None, plot_results=False,):
    """
    Perform predictions on test data and evaluate the model performance.

    Args:
        config: Configuration object containing model settings.
        data_test: Test data object.
        data_scale: Namedtuple containing scaler and sigma used in model training.
        model: PyTorch model object.
        device: Torch device (CPU or GPU).
        pred_source (str): Source for predictions (e.g., TCGA or GTEX).
        path_getBM (str): Path to the file containing gene expression data.
        path_save (str, optional): Path to directory where results will be saved. Defaults to None.
        plot_results (bool, optional): Whether to plot the results. Defaults to False.

    Returns:
        ResultsPrediction: Namedtuple containing Spearman correlation, Pearson correlation, and mean squared error for each tumor type.
    """

    # 1) Initialize results frames:
    df_spear_cor = pd.DataFrame()
    df_pears_cor = pd.DataFrame()
    df_error = pd.DataFrame()
    
    # 2) Get Test Data to do the predictions
    df_test = data_test.df_raw
    test_labels = data_test.raw_labels
    test_gn = data_test.raw_gn
    print('[utils][do_predictions] Get Test Data to do the predictions ... -> DONE')

    # 3) Select just the samples of the desired source and remove column source
    df_test = df_test[df_test.source == pred_source]
    print(f'[utils][do_predictions] Selected source to make the predictions: {df_test.source.unique()}')
    test_labels = test_labels.loc[df_test.index]
    test_gn = test_gn.loc[df_test.index]
    print('[utils][do_predictions] Select just the samples of the desired source and remove column source... -> DONE')
    df_test = df_test.drop('source', axis=1)
    tumor_types = list(df_test.tumor_type.unique())

    for i in tumor_types:
        print(f'[utils][do_predictions] Calculating predictions for {i} tumor type')
        # 4) Select just tumor type i samples
        df_rbps = df_test[df_test.tumor_type == i]
        df_labels = test_labels.loc[df_rbps.index]
        df_gn = test_gn.loc[df_rbps.index]
        # 4.1) Remove tumor_type column
        df_rbps = df_rbps.drop('tumor_type', axis=1)
        print(f'[utils][do_predictions] Shape of {i} tumor type test data')
        print('rbp test data:', df_rbps.shape)
        print('test label data:', df_labels.shape)
        print('test_gn data:', df_gn.shape)

        # 5) Scale Data 
        df_scaled_rbps = get_scaled_rbp_test_data(df_test=df_rbps, data_scale=data_scale)
        print(f'[utils][do_predictions] Scale Data of {i} tumor type test data... -> DONE')

        # 6) Make Predictions and Calculate the total Spearman Correlation and the Pearson Correlation
        trans_names = list(df_labels.columns)
        index_names = list(df_labels.index)

        model.eval()
        with torch.no_grad():
            # Move input tensors to the same device as the model
            input_rbps = torch.Tensor(df_scaled_rbps.values.astype(np.float64)).to(device)
            input_gn = torch.Tensor(df_gn.values.astype(np.float64)).to(device)
            pred = model(input_rbps, input_gn).detach().cpu().numpy()     

        pred_df = pd.DataFrame(pred, columns=trans_names, index=index_names) # Plot Pred vs Label expression ratio histogram

        pred = pred.flatten()
        labels = df_labels.values.flatten()
        
        spear_cor = stats.spearmanr(pred, labels)[0] # total correlation between real and pred
        mse = mean_squared_error(pred, labels)
        pear_cor = stats.pearsonr(pred, labels)[0]
        print(f'[utils][do_predictions] Calculate the total Spearman Correlation and the Pearson Correlation of {i} tumor type test data... -> DONE')
        print(f"[utils][do_predictions] Results -> spear cor: {round(spear_cor, 3)}, pear_cor: {round(pear_cor, 3)} & mse: {round(mse, 3)}")

        # 7) Plot results
        if plot_results:
            if path_save is None:
                raise ValueError("If plot_results is True, path_save must be provided.")
            path_save_scat = path_save + f'/{pred_source}/scat_plot_real_vs_pred_value'
            path_save_hist = path_save + f'/{pred_source}/pred_label_expression_ratio_histogram'
            check_create_new_directory(path_save_scat)
            check_create_new_directory(path_save_hist)
    
            # Plot real vs pred
            plots.plot_real_vs_pred(config, i, pred, labels, spear_cor, pear_cor, path_save_scat, config.source_train, pred_source) 
            # Checking that the sum of predicted transcripts approaches theoretical gene expression: ratio histogram.
            plots.plot_expression_ratio_histogram(path_getBM, df_trans_pred=pred_df, df_trans_label=df_labels, df_gn=df_gn, 
                                                  tumor_name=i, path_save=path_save_hist, source=config.source_train, source_pred=pred_source)
            # plots.pred_label_expression_ratio_histogram(path_getBM, pred_df, df_gn, i, path_save_hist, config.source_train, pred_source) 
        print(f'[utils][do_predictions] Plot results for {i} tumor type test data... -> DONE')

        # 5) Save the results  
        df_spear_cor.loc[0, i] = spear_cor
        df_pears_cor.loc[0, i] = pear_cor
        df_error.loc[0, i] = mse
        
    ResultsPrediction = namedtuple('ResultsPrediction', ['df_spear_cor', 'df_pears_cor', 'df_error'])
    results_prediction = ResultsPrediction(df_spear_cor, df_pears_cor, df_error)
    return results_prediction

def get_scaled_rbp_test_data(df_test, data_scale):
    """
    Scale test data with the scaler used in Training
    Args:
        df_test (DataFrame): DataFrame containing the test data
        data_scale (object): Object containing the scaler used for scaling during training
        
    Returns:
        scaledTest_df (DataFrame): Scaled test data
    """
    # Scale test data with the scaler used in Training
    scaledTest_df = pd.DataFrame(data_scale.scaler_sfs.transform(df_test),  # Standard Scale
                                index=df_test.index, 
                                columns=df_test.columns)
    # Apply clip 2sigma
    scaledTest_df = np.clip(scaledTest_df, -2*data_scale.sigma_sfs, 2*data_scale.sigma_sfs, axis=1)
    # Move the data 2sigma
    scaledTest_df += 2*data_scale.sigma_sfs
    # Divide the data by 4sigma
    scaledTest_df /= 4*data_scale.sigma_sfs 
    print('[utils][get_scaled_rbp_test_data] Scale test data with the scaler used in Training ... -> DONE ')
    return scaledTest_df

def perform_predictions_on_test_data(config, data_test, data_scale, model, device, path_data, path_save, plot_results=False):
    """
    Perform predictions on test data.

    Args:
        config: Configuration object containing model settings.
        data_test: Test data.
        data_scale: Namedtuple containing scaler and sigma used in model training.
        model: PyTorch model object.
        device: Torch device (CPU or GPU).
        plot_results (bool): Whether to plot and save prediction results. Defaults to False.
        
    Returns:
        results_tcga: Results of predictions on TCGA data.
        results_gtex: Results of predictions on GTEX data.
    """

    check_create_new_directory(path_save)
    path_getBM = path_data+'../extra/getBM_reduced.csv'

    results_tcga = do_predictions(
        config=config,
        data_test=data_test,
        data_scale=data_scale,
        model=model,
        device=device,
        path_save=path_save,
        plot_results=plot_results,
        pred_source='TCGA',
        path_getBM=path_getBM
    )

    results_gtex = do_predictions(
        config=config,
        data_test=data_test,
        data_scale=data_scale,
        model=model,
        device=device,
        path_save=path_save,
        plot_results=plot_results,
        pred_source='GTEX',
        path_getBM=path_getBM
    )
    return results_tcga, results_gtex

########################################## in-silico validation of DeepRBP using Postar 3 ##############################

def prepare_inputs_explainability(tumor_type, source_explain, config, path_data, path_train_files):
    """
    Prepare inputs for POSTAR3 explainability analysis.

    Args:
        tumor_type: String specifying the tumor type.
        source_explain: String specifying the source for explainability.
        config: Configuration object, possibly containing hyperparameters and other settings.
        path_data: Path to the data folder.
        path_train_files: Path to the directory containing training files, including the trained model and scalers.

    Returns:
        namedtuple: Inputs including scaled test data, test labels, gene names, and device.
    """
    ### Get the Test data of the desired tumor type
    print(f"[prepare_inputs_explainability] Explainability is going to be performed for: {tumor_type} of source {source_explain}")
    data_test = get_data(config, path=path_data, set_mode='test')
    df_test = data_test.df_raw
    test_labels = data_test.raw_labels
    test_gn = data_test.raw_gn
    print("[prepare_inputs_explainability] Read Test data ... -> DONE")
    ## In the Test data select just the samples of the desired source and remove column source
    if source_explain != 'all':
        df_test = df_test[df_test.source == source_explain]
        test_labels = test_labels.loc[df_test.index]
        test_gn = test_gn.loc[df_test.index]
    ## Drop the 'source' and 'tumor_type' columns
    df_test = df_test.drop(['tumor_type', 'source'], axis=1)  
    print(f"[prepare_inputs_explainability] Select just the samples of the {source_explain} source and remove column 'source' & 'tumor_type'... -> DONE")
    ### Load the scaler and sigma used in the SF data in the model training and scale the Test data
    scaler_sfs = joblib.load(path_train_files+'/scaler_sfs.joblib')
    with open(path_train_files+'/sigma_sfs.txt', 'r') as f:
        sigma_sfs = f.readline().strip()
    sigma_sfs = np.float128(sigma_sfs)
    Data_Scale = namedtuple('Data_Scale', ['scaler_sfs', 'sigma_sfs'])
    data_scale = Data_Scale(scaler_sfs, sigma_sfs)
    df_scaled_test = get_scaled_rbp_test_data(df_test=df_test, data_scale=data_scale)
    print(f"[prepare_inputs_explainability] Load the scaler and sigma vale used for SF data in model training and scale the Test data ... -> DONE")
    ### Load the trained model
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    print('[prepare_inputs_explainability] The device selected is:', device)
    #Inputs = namedtuple('Inputs', ['df_scaled_test', 'test_labels', 'test_gn', 'model', 'device'])
    Inputs = namedtuple('Inputs', ['df_scaled_test', 'test_labels', 'test_gn', 'device'])
    return Inputs(df_scaled_test, test_labels, test_gn, device)

def read_postar3_data(path_data, select_tumor_type, tumor2tissue):
    """
    Read the Postar3 data for the selected tissue.

    Args:
        path_data (str): Path to the directory containing the data.
        select_tumor_type (str): The selected tumor type for explainability.
        tumor2tissue (dict): Dictionary mapping tumor types to tissues.

    Returns:
        pandas.DataFrame: DataFrame containing the Postar3 data.
    """
    ## Read the GxRBP Validation data for the selected tissue
    matching_postar_file = tumor2tissue.get(select_tumor_type)
    print('[explainability][read_validation_data] The selected tissue test used to do the explainability is:', select_tumor_type)
    print('[explainability][read_validation_data] For this tissue you have just this Postar data available:', matching_postar_file)
    print('[explainability][read_validation_data] Loading Postar data ...')
    tissue_path = f"{path_data}/create_gxRBP/results/{matching_postar_file}.csv"
    df_val_GxRBP = pd.read_csv(tissue_path, index_col=0)
    return df_val_GxRBP

def analyze_postar_matrix(df_val_GxRBP, path_save): 
    """
    Analyze the POSTAR matrix and save the results.

    Parameters:
        df_val_GxRBP (DataFrame): The POSTAR matrix DataFrame containing gene-RBP scores.
        path_save (str): The path to save the results.
    """
    df_val_GxRBP_filtered = df_val_GxRBP.copy()
    
    # Number of RBPs per Gene
    df_count_rbps_per_gen = pd.DataFrame()
    df_count_rbps_per_gen['Class 0'] = df_val_GxRBP_filtered.apply(lambda x: (x == 0).sum(), axis=1)
    df_count_rbps_per_gen['Class 1'] = df_val_GxRBP_filtered.apply(lambda x: (x == 1).sum(), axis=1)
    df_count_rbps_per_gen['Class NaN'] = df_val_GxRBP_filtered.apply(lambda x: x.isna().sum(), axis=1)
    df_count_rbps_per_gen['Genes'] = df_count_rbps_per_gen.index
    df_count_rbps_per_gen = df_count_rbps_per_gen.reset_index(drop=True)

    # Number of genes per RBP
    df_count_genes_per_rbp = pd.DataFrame()
    df_count_genes_per_rbp['Class 0'] = df_val_GxRBP_filtered.apply(lambda x: (x == 0).sum(), axis=0)
    df_count_genes_per_rbp['Class 1'] = df_val_GxRBP_filtered.apply(lambda x: (x == 1).sum(), axis=0)
    df_count_genes_per_rbp['Class NaN'] = df_val_GxRBP_filtered.apply(lambda x: x.isna().sum(), axis=0)
    df_count_genes_per_rbp['RBPs'] = df_count_genes_per_rbp.index
    df_count_genes_per_rbp = df_count_genes_per_rbp.reset_index(drop=True)
    
    path_save = path_save+'/Postar'
    check_create_new_directory(path_save)
    
    # Write the results of the analysis in csv files & plot the analysis
    name_save1 = 'df_count_rbps_per_gen_nan_included'
    name_save2 = 'df_count_genes_per_rbp_nan_included'
    
    df_count_rbps_per_gen.to_csv(f'{path_save}/{name_save1}_all_Postar.csv')
    df_count_genes_per_rbp.to_csv(f'{path_save}/{name_save2}_all_Postar.csv')
        
    df_count_rbps_per_gen = df_count_rbps_per_gen.drop(columns=['Class NaN'])
    df_count_genes_per_rbp = df_count_genes_per_rbp.drop(columns=['Class NaN'])
    plots.plot_analyze_postar_matrix(df_count_rbps_per_gen, df_count_genes_per_rbp, path_save)

def analyze_results_per_rbp_or_gene(df_score_GxRBP, df_val_nan_included_GxRBP, path_save, path_data, getBM): 
    """
    Analyze the results per RNA binding protein (RBP) or gene, including plotting distributions, ROC curves, and identifying NaN candidates.

    Parameters:
        df_score_GxRBP (DataFrame): DataFrame containing scores data for GxRBP relationships.
        df_val_nan_included_GxRBP (DataFrame): DataFrame containing NaN-included values for GxRBP relationships.
        path_save (str): Path to save the results.
        path_data (str): Path to data directory.
        getBM (DataFrame): DataFrame containing gene information.

    Returns:
        None
    """
    
    path = '/'.join(path_save.rsplit('/')[0:9])
    df_count_genes_per_rbp = pd.read_csv(f'{path}/Postar/df_count_genes_per_rbp_nan_included_all_Postar.csv', index_col=0)
    list_rbps_postar = df_count_genes_per_rbp.sort_values(by='Class 1', ascending=False)['RBPs'].values.tolist()
    num_trans_per_gene = getBM.Gene_ID.value_counts()
    
    # Reorder the score dataframes 
    df_score_top_rbps_GxRBP = df_score_GxRBP.copy()
    df_val_top_rbps_GxRBP = df_val_nan_included_GxRBP.copy()
    df_score_top_rbps_GxRBP = df_score_top_rbps_GxRBP.loc[:,list_rbps_postar]
    df_val_top_rbps_GxRBP = df_val_top_rbps_GxRBP.loc[:,list_rbps_postar]
    
    # Melt the information of scores and postar classes
    df_score_melted = df_score_top_rbps_GxRBP.reset_index().melt(id_vars=['index'], var_name='RBP_name', value_name='Scores')
    df_val_melted = df_val_top_rbps_GxRBP.reset_index().melt(id_vars=['index'], var_name='RBP_name', value_name='Postar')
    df_score_melted.rename(columns={'index': 'Gene_ID'}, inplace=True)
    df_val_melted.rename(columns={'index': 'Gene_ID'}, inplace=True)
    df_combined1 = pd.merge(df_score_melted, df_val_melted, on=['Gene_ID','RBP_name'], how='inner')
    df_combined1['Num_Trans_Per_Gene'] = df_combined1['Gene_ID'].map(num_trans_per_gene)
    
    ####### 2) Genes with more RBPs class 1s relationships in Postar
    df_count_rbps_per_gen = pd.read_csv(f'{path}/Postar/df_count_rbps_per_gen_nan_included_all_Postar.csv', index_col=0)
    list_genes_postar = df_count_rbps_per_gen.sort_values(by='Class 1', ascending=False)['Genes'].values.tolist()

    # Filter the score dataframes 
    df_score_top_genes_GxRBP = df_score_GxRBP.copy()
    df_val_top_genes_GxRBP = df_val_nan_included_GxRBP.copy()
    df_score_top_genes_GxRBP = df_score_top_genes_GxRBP.loc[list_genes_postar,:]
    df_val_top_genes_GxRBP = df_val_top_genes_GxRBP.loc[list_genes_postar,:]
    
    # Melt the information of scores and postar classes
    df_score_melted = df_score_top_genes_GxRBP.reset_index().melt(id_vars=['index'], var_name='RBP_name', value_name='Scores')
    df_val_melted = df_val_top_genes_GxRBP.reset_index().melt(id_vars=['index'], var_name='RBP_name', value_name='Postar')
    df_score_melted.rename(columns={'index': 'Gene_ID'}, inplace=True)
    df_val_melted.rename(columns={'index': 'Gene_ID'}, inplace=True)
    df_combined2 = pd.merge(df_score_melted, df_val_melted, on=['Gene_ID','RBP_name'], how='inner')
    
    # Save Results
    check_create_new_directory(path_save)
    df_combined1.sort_values(by='Scores', ascending=False).to_csv(f'{path_save}/ordered_GxS_scores_of_k_rbps_with_more_genes_in_all_Postar.csv')
    df_combined2.sort_values(by='Scores', ascending=False).to_csv(f'{path_save}/ordered_GxS_scores_of_k_genes_with_more_rbps_in_all_Postar.csv')
    
    print('[analyze_results_per_rbp_or_gene] Plotting results ...')
    plots.plot_score_results(list_rbps_postar, list_genes_postar, df_combined1, df_combined2, path_save)
    print('[analyze_results_per_rbp_or_gene] Plotting results ... -> DONE\n')

    print('[analyze_results_per_rbp_or_gene] Calculating RBP thresholds ...')
    df_combined1_without_nan = df_combined1[df_combined1.Postar.isin([0,1])] 
    optimal_thresholds_df = calculate_optimal_thresholds_per_rbp(df_combined1_without_nan, path_save)
    plots.plot_one_ktop_rbp_with_thresholds(df_scores=df_combined1, rbp_name=df_combined1.iloc[0].RBP_name, optimal_thresholds_df=optimal_thresholds_df, path_save=path_save, getBM=getBM)
    print('[analyze_results_per_rbp_or_gene] Calculating RBP thresholds ... -> DONE\n')

    print('[analyze_results_per_rbp_or_gene] Giving some candidades from NaN Postar ...')
    df_combined1_nan_candidates = df_combined1[df_combined1['Postar'].isna()].copy()
    max_calculated_thr = optimal_thresholds_df.Optimal_Score_Threshold.max()

    ### 1) NaNs in Known RBPs in Postar
    df_combined1_well_nan_candidates = pd.merge(df_combined1_nan_candidates, optimal_thresholds_df, left_on='RBP_name', right_on='RBP_name', how='right')
    
    # Create a new column indicating whether the sample exceeds the threshold and add info about its well-characterization
    df_combined1_well_nan_candidates['Threshold_exceeded'] = df_combined1_well_nan_candidates['Scores'] > df_combined1_well_nan_candidates['Optimal_Score_Threshold']
    df_combined1_well_nan_candidates['Known_RBP_in_Postar'] = True
    
    ### 2) NaNs in UnKnown RBPs in Postar
    df_combined1_new_nan_candidates = df_combined1_nan_candidates[~df_combined1_nan_candidates.RBP_name.isin(optimal_thresholds_df.index.tolist())].copy()
    df_combined1_new_nan_candidates['Optimal_Score_Threshold'] = max_calculated_thr
    df_combined1_new_nan_candidates['Threshold_exceeded'] = df_combined1_new_nan_candidates['Scores'] > df_combined1_new_nan_candidates['Optimal_Score_Threshold']
    df_combined1_new_nan_candidates['Known_RBP_in_Postar'] = False

    df_nan_candidates_results = pd.concat([df_combined1_new_nan_candidates, df_combined1_well_nan_candidates], ignore_index=True)
    df_nan_candidates_results = df_nan_candidates_results.sort_values(by=['Known_RBP_in_Postar', 'Threshold_exceeded', 'Scores'], ascending=[False, False, False])
    df_nan_candidates_results.to_csv(f'{path_save}/results_table_rbp_nan_candidates.csv')
    print('[analyze_results_per_rbp_or_gene] Giving some candidades from NaN Postar ... -> DONE')

def calculate_optimal_thresholds_per_rbp(df, path_save):    
    """
    Calculate optimal thresholds for well-characterized RBPs in Postar and plot distributions and ROC curves.

    Parameters:
        df (DataFrame): DataFrame containing scores and Postar classes.
        path_save (str): Path to save the results.

    Returns:
        DataFrame: DataFrame containing optimal thresholds for each RBP.
    """
    print("[calculate_optimal_thresholds_per_rbp] Calculating optimal thresholds per RBP...")
    # Create a list to store optimal thresholds for each RBP
    optimal_thresholds = []
    # Get the unique list of RBP_names in the DataFrame
    unique_rbps = df['RBP_name'].unique()
    # Iterate over each RBP and calculate the optimal threshold
    for rbp in unique_rbps:
        print(rbp)
        # Filter the DataFrame for the current RBP
        filtered_df = df[(df['RBP_name'] == rbp)]
        # Skip if there are no rows for the current RBP
        if filtered_df.empty:
            continue
        # Calculate the ROC curve
        fpr, tpr, thresholds = roc_curve(filtered_df['Postar'], filtered_df['Scores'])
        # Find the index of the optimal threshold (maximizes the sum of sensitivity and specificity)
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        # Store the optimal threshold in the list
        optimal_thresholds.append({'RBP_name': rbp, 'Optimal_Score_Threshold': optimal_threshold})
        # Plot distributions and ROC curve for the first 33 RBPs
        if len(optimal_thresholds) <= 37: # ERRROR!!
            plots.plot_distributions_and_roc(filtered_df, rbp, optimal_threshold, fpr, tpr, optimal_idx, path_save)
    # Create a DataFrame with the results
    optimal_thresholds_df = pd.DataFrame(optimal_thresholds)
    csv_path = f'{path_save}/optimal_thresholds.csv'
    optimal_thresholds_df.to_csv(csv_path, index=False)
    print(f'Results saved to: {csv_path}')
    return optimal_thresholds_df


########################################## in-silico validation of DeepRBP using real kds ################################
def check_data_exists(path_data, experiment):
    """
    Checks if data related to a specific experiment already exists in the given directory.

    Args:
        path_data (str): The path where the data is located.
        experiment (str): The name of the experiment to check for.

    Returns:
        bool: True if the data does not exist, False if the data exists.

    """
    folder_path = f'{path_data}/datasets'
    if any(
        os.path.isdir(os.path.join(folder_path, subfolder)) and experiment in subfolder and any(
            os.path.isfile(os.path.join(folder_path, subfolder, filename))
            for filename in os.listdir(os.path.join(folder_path, subfolder))
        )
        for subfolder in os.listdir(folder_path)
    ):
        print('[utils] The files are already created! \n')
        print('[utils] The already created data is: \n')
        for subfolder in os.listdir(folder_path):
            subfolder_path = os.path.join(folder_path, subfolder)
            if os.path.isdir(subfolder_path) and experiment in subfolder:
                print(f'Subdirectory: {subfolder}')
                files = [filename for filename in os.listdir(subfolder_path) if os.path.isfile(os.path.join(subfolder_path, filename))]
                print(f'Files: {files}\n')
        return False
    else:
        print('[utils] Warning: The data for the experiment you are asking is not created')
        return True

def create_vec_exp_matrix_from_abundance(path_data, id_sample):
    """
    Create vectors for transcript expression and gene expression from abundance data.

    Parameters:
    - path_data (str): Path to the data directory.
    - id_sample (str): ID of the sample.

    Returns:
    - df_trans_tpm (DataFrame): DataFrame with transcript TPM values.
    - df_trans_counts (DataFrame): DataFrame with transcript count values.
    - df_genes_tpm (DataFrame): DataFrame with gene TPM values.
    """
    ### Read the abundance.tsv file
    df = pd.read_csv(f'{path_data}/{id_sample}/abundance.tsv', sep='\t')
    print(f'[create_vec_exp_matrix_from_abundance] Read the abundance.tsv file of sample: {id_sample} ... -> DONE')
    ### Expand detailed information from target_id column and reorder columns
    df[['Transcript_ID', 'Gene_ID', 'OTTHUMG', 'OTTHUMT', 'Gene_Symbol', 'Gene_Name', 'Length', 'Biotype']] = df['target_id'].str.rstrip('|').str.split('|', expand=True)
    df.drop(columns=['target_id', 'Length'], inplace=True)
    df = df[['Transcript_ID', 'Gene_ID', 'OTTHUMT', 'OTTHUMG', 'Gene_Symbol', 'Gene_Name', 'Biotype', 'length', 'eff_length', 'est_counts', 'tpm']]
    print('[create_vec_exp_matrix_from_abundance] Expand detailed information from target_id column and reorder columns ... -> DONE')
    ### Remove version name from Transcript_ID and Gene_ID
    df['Transcript_ID'] = df['Transcript_ID'].str.split('.').str[0]
    df['Gene_ID'] = df['Gene_ID'].str.split('.').str[0]
    print('[create_vec_exp_matrix_from_abundance] Remove version name from Transcript_ID and Gene_ID ... -> DONE')
    ### Save the final abundance file (with counts)
    df.to_csv(f'{path_data}/{id_sample}/processed_abundance.csv')
    ### Create the final vectors of transcript expression and gene expression for this sample
    # trans (tpm & counts)
    df_trans_tpm = pd.DataFrame(df['tpm'].values, index=df['Transcript_ID'].values, columns=[f'{id_sample}']).reset_index()
    df_trans_counts = pd.DataFrame(df['est_counts'].values, index=df['Transcript_ID'].values, columns=[f'{id_sample}']).reset_index() 
    # genes
    sum_tpm_df = df.groupby('Gene_ID')['tpm'].sum().reset_index() # Create the grouping by 'Gene_ID' and calculate the sum of the 'tpm' column
    df_genes_tpm = pd.DataFrame(sum_tpm_df['tpm'].tolist(), index=sum_tpm_df['Gene_ID'].tolist(), columns=[f'{id_sample}']).reset_index()
    df_trans_tpm.rename(columns={'index': 'sample'}, inplace=True)
    df_trans_counts.rename(columns={'index': 'sample'}, inplace=True)
    df_genes_tpm.rename(columns={'index': 'sample'}, inplace=True)
    print('[create_vec_exp_matrix_from_abundance] Create the final vectors of transcript expression and gene expression for this sample ... -> DONE')
    print('\n')
    return df_trans_tpm, df_trans_counts, df_genes_tpm

def process_and_save_expression_data(path_exp, id_samples, condition, path_data):
    """
    Process expression data from multiple samples, aggregate, and save the results.

    Parameters:
    - path_exp (str): Path to the experiment data directory.
    - id_samples (list): List of sample IDs.
    - condition (str): Condition or tissue type.
    - path_data (str): Path to the data directory.

    Saves:
    - CSV files containing processed expression data for RBPs, transcripts, counts, and genes.

    Returns:
    None
    """

    # Initialize cumulative DataFrames
    df_trans_tpm_cumulative = pd.DataFrame()
    df_trans_counts_cumulative = pd.DataFrame()
    df_genes_tpm_cumulative = pd.DataFrame()

    # 1. Merge each sample
    print(f'[process_and_save_expression_data] 1. Merging the samples ... ')
    for sample in id_samples:
        print(f'[process_and_save_expression_data] We are processing the sample: {sample} from tissue: {condition}\n')
        # Obtain data for the current sample
        df_trans_tpm, df_trans_counts, df_genes_tpm = create_vec_exp_matrix_from_abundance(path_exp, sample)
        
        # Merge in the cumulative DataFrames on the 'sample' column
        if df_trans_tpm_cumulative.empty:
            # If the cumulative DataFrame is empty, simply copy the data
            df_trans_tpm_cumulative = df_trans_tpm.copy()
            df_trans_counts_cumulative = df_trans_counts.copy()
            df_genes_tpm_cumulative = df_genes_tpm.copy()
        else:
            # Merge in the cumulative DataFrames on the 'sample' column
            df_trans_tpm_cumulative = pd.merge(df_trans_tpm_cumulative, df_trans_tpm, on='sample', how='outer')
            df_trans_counts_cumulative = pd.merge(df_trans_counts_cumulative, df_trans_counts, on='sample', how='outer')
            df_genes_tpm_cumulative = pd.merge(df_genes_tpm_cumulative, df_genes_tpm, on='sample', how='outer')
    print(f'[process_and_save_expression_data] 1. Merging the samples ... -> DONE\n')
    
    # 2. Data Cleaning and aggregation of different loci
    # Give to different loci the same gen_id and trans_id 
    df_trans_tpm_cumulative['sample'] = df_trans_tpm_cumulative['sample'].str.replace('R', '0')
    df_trans_counts_cumulative['sample'] = df_trans_counts_cumulative['sample'].str.replace('R', '0')
    df_genes_tpm_cumulative['sample'] = df_genes_tpm_cumulative['sample'].str.replace('R', '0')
    df_trans_tpm_cumulative = df_trans_tpm_cumulative.iloc[:, 1:].groupby(df_trans_tpm_cumulative['sample']).sum()
    df_trans_counts_cumulative = df_trans_counts_cumulative.iloc[:, 1:].groupby(df_trans_counts_cumulative['sample']).sum()
    df_genes_tpm_cumulative = df_genes_tpm_cumulative.iloc[:, 1:].groupby(df_genes_tpm_cumulative['sample']).sum()      
    print('[process_and_save_expression_data] 2. Data Cleaning and aggregation of different loci ... -> DONE\n')
    
    # 3. Log2 transformation of Transcript data
    df_trans_log2p_tpm_cumulative = np.log2(df_trans_tpm_cumulative + 1)
    print('[process_and_save_expression_data] 3. Log2 transformation of Transcript data ... -> DONE\n')
    
    # 4. Obtaining RBP expression dataframe in log2(tpm+1)
    getBM = pd.read_csv(f'{path_data}/extra/getBM_total.csv')
    df_RBPs = pd.read_excel(f'{path_data}/selected_genes_rbps/Table_S2_list_RBPs_eyras.xlsx', skiprows=2) # List of splicing genes
    ids_RBPs = getBM[getBM['Gene_name'].isin(df_RBPs['HGNC symbol'].tolist())]['Gene_ID'].unique().tolist()
    df_rbp_tpm_cumulative = df_genes_tpm_cumulative.loc[ids_RBPs,:]
    df_rbp_log2p_tpm_cumulative = np.log2(1+df_rbp_tpm_cumulative)
    gene_id_name_map = {row['Gene_ID']: row['Gene_name'] for index, row in getBM[getBM.Gene_ID.isin(ids_RBPs)].drop_duplicates(subset='Gene_name').iterrows()} # Change the colnames of rbp set by gene_name
    df_rbp_log2p_tpm_cumulative = df_rbp_log2p_tpm_cumulative.loc[gene_id_name_map.keys(), :] #the are 4 rbp_ids that have the same rbp_name
    df_rbp_log2p_tpm_cumulative = df_rbp_log2p_tpm_cumulative.rename(index=gene_id_name_map)
    print('[process_and_save_expression_data] 4. Obtaining RBP expression dataframe in log2(tpm+1) ... -> DONE\n')
    
    # 5. Gene selection depending on model used (cancer genes or protein coding genes) and get the dataset with the gene expression for each transcript
    getBM_mini = pd.read_csv(f'{path_data}/extra/getBM_reduced.csv', index_col=0)
    list_transcripts = list(getBM_mini.Transcript_ID)
    list_genes = list(getBM_mini.Gene_ID)
    df_trans_log2p_tpm_cumulative = df_trans_log2p_tpm_cumulative.loc[list_transcripts, :] # Get the dataset with the expression of the transcripts
    df_gn_each_iso_tpm_cumulative = df_genes_tpm_cumulative.loc[list_genes,:] # Get the dataset with the gene expression for each transcript
    print('[process_and_save_expression_data] 5. Gene selection depending on model used and get the dataset with the gene expression for each transcript ... -> DONE\n')
    
    # 6. Transpose DataFrames
    df_rbp_log2p_tpm_cumulative = df_rbp_log2p_tpm_cumulative.T
    df_trans_log2p_tpm_cumulative = df_trans_log2p_tpm_cumulative.T
    df_gn_each_iso_tpm_cumulative = df_gn_each_iso_tpm_cumulative.T
    df_gn_each_iso_tpm_cumulative.columns = list_transcripts # We put to the gene_expr_each_iso dataframe the colnames of the transcripts they are related to
    print('[process_and_save_expression_data] 6. Transpose DataFrames ... -> DONE\n')
    
    # 7. Save created files in folder
    path = f'{path_exp}/datasets/{condition}'
    check_create_new_directory(path)
    path_rbp = f'{path}/{condition}_RBPs_log2p_tpm.csv'
    path_trans =  f'{path}/{condition}_trans_log2p_tpm.csv'
    path_trans_counts = f'{path}/{condition}_trans_est_counts.csv'
    path_gn = f'{path}/{condition}_gn_expr_each_iso_tpm.csv'
    df_rbp_log2p_tpm_cumulative.to_csv(path_rbp) #rbps
    df_trans_log2p_tpm_cumulative.to_csv(path_trans) #trans
    df_trans_counts_cumulative.to_csv(path_trans_counts)
    df_gn_each_iso_tpm_cumulative.to_csv(path_gn) #genes for each iso
    print('[process_and_save_expression_data] 7. Save created files in folder ... -> DONE \n')
    print(f'[process_and_save_expression_data] {condition} files were saved in {path}')
    

def identify_samples_from_tissue(df_filereport, condition):
    """
    Identify sample IDs based on the given condition from the file report DataFrame.

    Parameters:
    - df_filereport (DataFrame): DataFrame containing file report data.
    - condition (str): Condition or tissue type.

    Returns:
    - id_samples (list): List of sample IDs related to the given condition.
    """
    
    print(f'[identify_samples_from_tissue] The condition selected is {condition}')
    if condition == 'HFE_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('HFE145siNeg')]['run_accession'].tolist()
    elif condition == 'HFE_MBNL1_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('HFE145siMBNL1')]['run_accession'].tolist()
    elif condition == 'AGS_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('AGSsiNeg')]['run_accession'].tolist()
    elif condition == 'AGS_PTB1_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('AGSsiPTBP1')]['run_accession'].tolist()
    elif condition == 'AGS_ESRP2_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('AGSsiESRP2')]['run_accession'].tolist()
    elif condition == 'SRR296_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('Control')]['run_accession'].tolist()
    elif condition == 'SRR296_RBM47_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('RBM47 KD')]['run_accession'].tolist()
    elif condition == 'GSE77702-kd1_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('Control')]['run_accession'].tolist()
    elif condition == 'GSE77702-kd1_FUS_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('FUS KD')]['run_accession'].tolist()
    elif condition == 'GSE77702-kd2_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('Control')]['run_accession'].tolist()
    elif condition == 'GSE77702-kd2_TAF15_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('TAF15 KD')]['run_accession'].tolist()
    elif condition == 'GSE77702-kd3_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('Control')]['run_accession'].tolist()
    elif condition == 'GSE77702-kd3_TARDBP_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('TARDBP KD')]['run_accession'].tolist()
    elif condition == 'GSE136366_control':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('Rescued_tdp43')]['run_accession'].tolist()
    elif condition == 'GSE136366_TARDBP_kd':
        id_samples = df_filereport[df_filereport['sample_title'].str.contains('tdp43_ko')]['run_accession'].tolist()
    else:
        print(f"Condition not known: {condition}")
        id_samples = []
    print(f'[identify_samples_from_tissue] The samples related to this condition are: {id_samples}\n')
    return id_samples

def generate_data(df_filereport, condition1, condition2, path_exp, path_data):
    """
    Generates data for the given conditions

    Args:
        df_filereport (DataFrame): DataFrame containing file report data.
        condition1 (str): Name of the first condition.
        condition2 (str): Name of the second condition.
        path_data (str): Path to the general data directory.
        path_exp (str): Path to the experiment data folder
        
    Returns:
        None
    """
    print('[utils] Generating the data ...')
    id_samples_normal = identify_samples_from_tissue(df_filereport, condition=condition1)
    id_samples_kd = identify_samples_from_tissue(df_filereport, condition=condition2)
    print('[utils] Identify the samples for this tissue ... -> DONE')
    
    process_and_save_expression_data(path_exp, id_samples=id_samples_normal, condition=condition1, path_data=path_data)
    process_and_save_expression_data(path_exp, id_samples=id_samples_kd, condition=condition2, path_data=path_data)
    print('[utils] Using the id_samples create the inputs,targets for DeepSF model and save in path ... -> DONE')
    print('[utils] Generating the data ... -> DONE \n')

def scale_inputs_real_knockout_explainability(path_model, folder_path, condition): 
    """ 
    Loads the data related to the specified condition, sorts the columns,
    scales the data using the scaler and sigma used in the SF data in the model training,
    and returns the scaled data along with the labels and gene expression data.

    Args:
        path_model (str): The path to the directory containing the model files.
        folder_path (str): The path to the folder containing the data related to the condition.
        condition (str): The condition for which the data needs to be prepared.

    Returns:
        Inputs: A namedtuple containing the scaled RBP data, labels, and gene expression data.
    """
    ### Load Data & Sort columns
    df_rbps = pd.read_csv(f'{folder_path}/{condition}/{condition}_RBPs_log2p_tpm.csv', index_col=0)
    df_labels = pd.read_csv(f'{folder_path}/{condition}/{condition}_trans_log2p_tpm.csv', index_col=0)
    df_gns = pd.read_csv(f'{folder_path}/{condition}/{condition}_gn_expr_each_iso_tpm.csv', index_col=0)
    df_rbps = df_rbps.sort_index(axis=1)
    df_labels = df_labels.sort_index(axis=1)
    df_gns = df_gns.sort_index(axis=1)
    print('[utils] Prepare inputs real knockout explainability ...')
    print('[utils] Load Data & Sort columns ... -> DONE \n')
    ### Scale Data: Load the scaler and sigma used in the SF data in the model training and scale
    scaler_sfs = joblib.load(path_model+'/scaler_sfs.joblib')
    with open(path_model+'/sigma_sfs.txt', 'r') as f:
        sigma_sfs = f.readline().strip()
    sigma_sfs = np.float128(sigma_sfs)
    Data_Scale = namedtuple('Data_Scale', ['scaler_sfs', 'sigma_sfs'])
    data_scale = Data_Scale(scaler_sfs, sigma_sfs)
    df_scaled_rbps = get_scaled_rbp_test_data(df_test=df_rbps, data_scale=data_scale)
    print('[utils] Scale Data: Load the scaler and sigma used in the SF data in the model training and scale \n')
    Inputs = namedtuple('Inputs', ['df_scaled_rbps', 'df_labels', 'df_gns'])
    return Inputs(df_scaled_rbps, df_labels, df_gns)

#def load_model(path_model, df_rbps_control, df_labels_control, device):
#    """
#    Loads a trained model from the specified path and returns it.
#
#    Args:
#        path_model (str): Path to the directory containing the model files.
#        df_rbps_control (pd.DataFrame): DataFrame containing the control RBP data.
#        df_labels_control (pd.DataFrame): DataFrame containing the control labels data.
#        device (torch.device): The device on which to load the model.
#
#    Returns:
#        torch.nn.Module: The loaded model.
#    """
#    with open(f'{path_model}/config.json', 'r') as file:
#        config_dict = json.load(file)
#    config_obj = Config(**config_dict)
#    model = DeepRBP(n_inputs=df_rbps_control.shape[1], n_outputs=df_labels_control.shape[1], config=config_obj, device=device)
#    model.load_state_dict(torch.load(path_model+'/model.pt', map_location=device))      
#    model.eval()
#    return model

def evaluate_model(model, df_rbps_control, df_gns_control, df_labels_control, df_rbps_kd, df_gns_kd, df_labels_kd, device, condition1, condition2):
    """
    Evaluates a trained model's performance on control and knockdown data.

    Args:
        model (torch.nn.Module): The trained model to evaluate.
        df_rbps_control (DataFrame): DataFrame containing scaled RBP features for control data.
        df_gns_control (DataFrame): DataFrame containing gene expression features for control data.
        df_labels_control (DataFrame): DataFrame containing labels for control data.
        df_rbps_kd (DataFrame): DataFrame containing scaled RBP features for knockdown data.
        df_gns_kd (DataFrame): DataFrame containing gene expression features for knockdown data.
        df_labels_kd (DataFrame): DataFrame containing labels for knockdown data.
        device (torch.device): The device on which to perform evaluations.
        condition1 (str): Name of the first condition (e.g., control condition).
        condition2 (str): Name of the second condition (e.g., knockdown condition).

    Returns:
        None
    """
    pred_control = model(torch.Tensor(df_rbps_control.values.astype('float64')).to(device), torch.Tensor(df_gns_control.values.astype('float64')).to(device)).detach().cpu().numpy()
    pred_kd = model(torch.Tensor(df_rbps_kd.values.astype('float64')).to(device), torch.Tensor(df_gns_kd.values.astype('float64')).to(device)).detach().cpu().numpy()
    label_control = df_labels_control.values.flatten()
    label_kd = df_labels_kd.values.flatten()
    
    spear_control = stats.spearmanr(pred_control.flatten(), label_control)[0] 
    mse_control = mean_squared_error(pred_control.flatten(), label_control)
    pear_control = stats.pearsonr(pred_control.flatten(), label_control)[0]
    
    spear_kd = stats.spearmanr(pred_kd.flatten(), label_kd)[0]
    mse_kd = mean_squared_error(pred_kd.flatten(), label_kd)
    pear_kd = stats.pearsonr(pred_kd.flatten(), label_kd)[0]
    print(f"[utils] Results in all {condition1} data -> spear cor: {round(spear_control, 3)}, pear_cor: {round(pear_control, 3)} & mse: {round(mse_control, 3)}")
    print(f"[utils] Results in all {condition2} data -> spear cor: {round(spear_kd, 3)}, pear_cor: {round(pear_kd, 3)} & mse: {round(mse_kd, 3)}")
    
def get_significant_samples(rbp_interest, experiment, path_data, df_deeplift_scores_TxRBP):
    """
    Retrieves significant samples in transcripts and genes from a Limma result file
    if it exists.

    Parameters:
    - rbp_interest (str): The RNA-binding protein of interest.
    - experiment (str): The name of the experiment.
    - path_data (str): The path to the data directory.
    - df_deeplift_scores_TxRBP (DataFrame): DataFrame containing DeepLIFT scores for transcripts.

    Returns:
    - significant_samples (list): List of significant samples in transcripts.
    - significant_samples_gns (list): List of significant samples in genes.

    Raises:
    - FileNotFoundError: If the specified Limma result file does not exist.
    """
    file_path = f"{path_data}/limma/limma_{experiment}_{rbp_interest}_results.csv"
    if os.path.isfile(file_path):
        print("There is a Limma result dataframe saved in your directory")
        print(file_path)
        df_limma = pd.read_csv(file_path, index_col=0)
        df_limma.set_index('Transcript_ID', inplace=True)
        df_significant = df_limma.loc[df_deeplift_scores_TxRBP.index].copy()
        df_significant = df_significant[df_significant['adj.P.Val'] < 0.05]
        df_significant['Transcript_ID'] = df_significant.index
        df_significant.reset_index(inplace=True, drop=True)
        significant_samples = df_significant.Transcript_ID.tolist()
        significant_samples_gns = df_significant.Gene_ID.unique().tolist()
        return significant_samples, significant_samples_gns
    else:
        # If the file does not exist, print a message indicating it was not found and raise an exception
        print("There is NOT a Limma result dataframe saved in your directory")
        raise FileNotFoundError("The specified file does not exist")