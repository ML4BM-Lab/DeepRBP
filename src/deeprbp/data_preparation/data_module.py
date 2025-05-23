
import pandas as pd
import lightning as L

from ..data_loading.config_loader import ConfigParser
from .prepare_data import PrepareData
from ..util.utils import adjust_batch_size

class DeepRBPDataModule(L.LightningDataModule):
    """
    A PyTorch Lightning DataModule for preparing and loading data related to RNA-binding proteins (RBPs), genes, and transcripts.

    This DataModule handles the loading, processing, and splitting of datasets for training, validation, and testing. 
    It leverages configuration settings provided through a ConfigParser object and manages data preparation steps 
    such as scaling and creating tensor datasets.

    Attributes:
        config (ConfigParser): Configuration parser containing all relevant settings for data loading and processing.
        output_dir (str): Directory where processed data and outputs will be saved.
        num_workers (int): Number of subprocesses to use for data loading. Default is 0 (no multiprocessing).
        verbose (int): Verbosity level for logging. Higher values increase the amount of output generated. 
                       Default is 1 (basic logging).
        getBM (pd.DataFrame): DataFrame containing mapping information for transcripts and genes, loaded from a CSV file.
        trans_col_name (str): Name of the column in the DataFrame that contains transcript identifiers. 
                              Default is "Transcript_ID".
        gene_col_name (str): Name of the column in the DataFrame that contains gene identifiers. 
                             Default is "Gene_ID".
        train_path_files (str): Path(s) to the training data files.
        val_fraction (float): Fraction of the training data to be used for validation.
        test_path_files (str): Path(s) to the test data files, if available.
        train_batch_size (int): Batch size for the training dataset.
        val_batch_size (int): Batch size for the validation dataset.
        sample_category (str): Colname in metadata referred to sample category used for stratied splitting of data.
        toy_sample_fraction (float or None): Fraction of samples to be used for a toy dataset.
        prep_data (PrepareData): Instance of PrepareData for handling the loading and preparation of datasets.
    """
    def __init__(self, config: ConfigParser, output_dir: str, num_workers: int = 0, verbose: int = 1):
        super().__init__()
        self.config = config
        self.output_dir = output_dir
        self.num_workers = num_workers
        self.verbose = verbose
        
        # Load the mapping DataFrame from the specified path in the configuration
        self.getBM = pd.read_csv(self.config.get('getBM_path'))
        self.trans_col_name=self.config.get('trans_col_name', default="Transcript_ID") # just used on create_tensor_dataset
        self.gene_col_name=self.config.get('gene_col_name', default="Gene_ID")
        
        # Load paths for training and testing data
        self.train_path_files = self.config.get('train_path_files')
        self.val_fraction = self.config.get('test_fraction')
        self.test_path_files = self.config.get('test_path_files', None)
        
        # Load batch sizes for training and validation
        self.train_batch_size = self.config.get('train_batch_size', None) 
        self.val_batch_size = self.config.get('val_batch_size', 256)
        
        # Load sampling information
        self.sample_category = self.config.get('sample_category') 
        self.toy_sample_fraction = self.config.get('sample_fraction', None) # Used for filtering samples
        
        # Initialize the PrepareData instance for data preparation
        self.prep_data = PrepareData(self.getBM, self.trans_col_name, self.gene_col_name, self.output_dir, self.verbose) 

        # Flags to check if data has already been prepared or set up
        self.is_data_prepared = False
        self.is_train_setup = False
        self.is_test_setup = False
 
    def prepare_data(self):
        """Load data"""
        if not self.is_data_prepared:
            self.training_all_data = self.prep_data.load_data(self.train_path_files)
            if self.test_path_files:
                self.test_data = self.prep_data.load_data(self.test_path_files)
            self.is_data_prepared = True
        else:
            print("[DeepRBPDataModule] 🔍 Data already prepared; skipping loading data.")

    def setup(self, stage=None):
        """Split, transform data and create tensor datasets. Called on each GPU separately - stage defines if we are at fit or test step.
        Setup is called from every process across all the nodes. Setting state here is recommended.
        """
        print('ESTOY LLAMANDO A SETUP BRO')
        # we set up only relevant datasets when stage is specified (automatically set by Pytorch-Lightning)
        if stage == 'fit' or stage is None:
            print('ESTOY INTENTANDO ENTRAR MÁS')
            if not self.is_train_setup:
                if self.toy_sample_fraction:
                    print(f"[DeepRBPDataModule] 🧩 Using toy sample fraction: {self.toy_sample_fraction}")
                    self.training_all_data = self.prep_data.filter_samples(self.training_all_data, self.sample_category, self.toy_sample_fraction)
                                                                           
                print("[DeepRBPDataModule] 🛠 Setting up training and validation data...")
                # Split data on train and validation
                self.train_data, self.val_data = self.prep_data.split_data(self.training_all_data, self.sample_category, self.val_fraction)
                self.prep_data.save_split_data(self.train_data, self.val_data)
                
                # Transform data
                self.prep_data.fit_scaler(self.train_data) 
                self.train_data, self.val_data = self.prep_data.scale_train_val_data(self.train_data, self.val_data)
                
                # Create the tensor dataset
                self.train_dataset = self.prep_data.create_tensor_dataset(self.train_data)
                self.val_dataset = self.prep_data.create_tensor_dataset(self.val_data)
                
                # Print the shapes of the resulting tensors and layers details only on rank 0  
                if self.trainer is not None and hasattr(self.trainer, 'is_global_zero') and self.trainer.is_global_zero:
                    self.train_dataset.print_tensor_shapes() 
                    self.val_dataset.print_tensor_shapes() 
                self.is_train_setup = True 
            else:
                print("[DeepRBPDataModule] 🔍 Training and validation data already set up; skipping setup.")
   
        if stage == 'test' or stage is None:
            if not self.is_test_setup:
                print("[DeepRBPDataModule] 🛠 Setting up test data...")
                if self.prep_data.scaler is None:
                    self.prep_data.load_scaler(self.output_dir+'/data')

                self.test_data = self.prep_data.scale_data(self.test_data)
                self.test_dataset = self.prep_data.create_tensor_dataset(self.test_data)
                
                if self.trainer is not None and hasattr(self.trainer, 'is_global_zero') and self.trainer.is_global_zero:
                    self.test_dataset.print_tensor_shapes() 

                self.is_test_setup = True
            else:
                print("[DeepRBPDataModule] 🔍 Test data already set up; skipping setup.")

    def train_dataloader(self):
        """Returns loader for training set"""
        print('ESTOY LLAMANDO A train_dataloader BRO')
        train_data_loader = self.prep_data.create_data_loader(
                    self.train_dataset, batch_size=adjust_batch_size(self.train_dataset, self.train_batch_size), 
                    shuffle = True, drop_last = True, num_workers = self.num_workers
        )
        return train_data_loader
 
    def val_dataloader(self):
        """Returns loader for validation set"""
        print('ESTOY LLAMANDO A val_dataloader BRO')
        val_data_loader = self.prep_data.create_data_loader(
                    self.val_dataset, batch_size=adjust_batch_size(self.val_dataset, self.val_batch_size), 
                    shuffle = False, drop_last = False, num_workers = self.num_workers
        )
        return val_data_loader
 
    def test_dataloader(self):
        """Returns loader for test set"""
        print('ESTOY LLAMANDO A test_dataloader BRO')
        test_data_loader = self.prep_data.create_data_loader(
                    self.test_dataset, batch_size=len(self.test_dataset), shuffle = False, drop_last = False, num_workers = self.num_workers
        )
        # Use a large batch size to ensure all samples are included for predictions.
        return test_data_loader
   
    def predict_dataloader(self, mode='test', custom_loader=None):
        """Returns loader for prediction set based on the specified mode.
        
        Args:
            mode (str): The mode can be 'train', 'val', 'test', or 'predict'. Default is 'test'.
            custom_loader (DataLoader, optional): A custom DataLoader to use when mode is 'predict'.

        Returns:
            DataLoader: The appropriate DataLoader for predictions.
        
        Raises:
            ValueError: If an invalid mode is provided. The mode must be one of 'train', 'val', 'test', or 'predict'.
     
        Note:
            This method is designed to be used with the `predict` method of the PyTorch Lightning `Trainer`. 

            ```python
            # Assuming 'trainer' is an instance of the Trainer class and 'model' is your PyTorch Lightning model:
            
            # To get predictions on the training dataset:
            data_module = DeepRBPDataModule(config, output_dir)
            train_predictions = trainer.predict(model, datamodule=data_module.predict_dataloader(mode='train'))

            # To get predictions on the validation dataset:
            val_predictions = trainer.predict(model, datamodule=data_module.predict_dataloader(mode='val'))

            # To get predictions on the test dataset:
            test_predictions = trainer.predict(model, datamodule=data_module.predict_dataloader(mode='test'))
        ```
        """
        print('ESTOY LLAMANDO A predict_dataloader BRO!!!!!')
        if mode == 'train':
            return self.train_dataloader()
        elif mode == 'val':
            return self.val_dataloader()
        elif mode == 'test':
            return self.test_dataloader()
        elif mode == 'predict':
            if custom_loader is None:
                raise ValueError("A custom DataLoader must be provided when mode is 'predict'.")
            return custom_loader
        else:
            raise ValueError("Invalid mode. Choose 'train', 'val', 'test', or 'predict'.")
    


    
   
