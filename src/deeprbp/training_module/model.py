# src/deeprbp/training_module/model.py

from typing import List, Optional
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score
# Pytorch-Lightning
import lightning as L
from torchmetrics import MeanMetric

from ..util.logger import Logger
from ..data_loading.config_loader import ConfigParser
from .evaluation import spearmanr_per_gene

import warnings
# Suppress FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)

class BaseLightningModule(L.LightningModule):  
    """
    Base Lightning module for predictor model training and evaluation
    Args:
        gene_names (List[str]): A list of gene IDs used in the model.
        trans_names (List[str]): A list of transcript IDs used in the model.
        getBM (pd.DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.
        verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during the pipeline execution.
                                 - 0: No logging.
                                 - 1: Basic logging (show progress and essential logs).
                                 - 2: Detailed logging (show additional information).
                                 Default is 1.
    """
    def __init__(self, gene_names: List[str], trans_names: List[str], getBM: pd.DataFrame,
                 input_features: Optional[str] = None, output_features: Optional[str] = None, 
                 verbose: int = 0):
        super().__init__()
        self.input_features = input_features if input_features is not None else ('scaled_rbp_df', 'gene_df')
        self.output_features = output_features if output_features is not None else ('isoform_df',)  
        self.criterion = nn.MSELoss()
        self.gene_names = gene_names
        self.trans_names = trans_names
        self.getBM = getBM
        self.verbose = verbose
        self.debugging = Logger(verbose=self.verbose)
        self.initialize_metrics()
    def initialize_metrics(self):
        """Initializes metrics for training, validation, and testing. 
        This method sets up the following metrics using `MeanMetric` from `torchmetrics`:
        - **Mean Squared Error (MSE)**: It's the loss functionm, measures the average squared difference between log2(tpm+1) predictions vs actual values.
        - **Pearson Correlation**: Assesses linear correlation between predictions and actual values, ranging from -1 to 1.
        - **Spearman Correlation**: Evaluates rank correlation, indicating the strength of monotonic relationships, also ranging from -1 to 1.
        - **R² (Coefficient of Determination)**: Indicates the proportion of variance explained by the model, with values from 0 to 1.
        """ 
        # Initialize train metrics
        self.train_loss = MeanMetric()
        self.train_corr_spearman = MeanMetric()
        self.train_corr_pearson = MeanMetric()
        self.train_r2 = MeanMetric()
        # Initialize validation metrics
        self.validation_loss = MeanMetric()
        self.validation_corr_spearman = MeanMetric()
        self.validation_corr_pearson = MeanMetric()
        self.validation_r2 = MeanMetric()
        # Initialize test metrics
        self.test_loss = MeanMetric()
        self.test_corr_spearman = MeanMetric()
        self.test_corr_pearson = MeanMetric()
        self.test_r2 = MeanMetric()
        self.test_corr_spearman_per_gene = MeanMetric()
        self.test_corr_spearman_per_gene_max = MeanMetric()
    ###
    def _prepare_batch(self, batch):
        """Prepares the inputs and targets from the batch.
        
        Args:
            batch (dict): Batch of data containing input and target features.
        
        Returns:
            tuple: A tuple containing (inputs, targets).
        """
        inputs = [batch[feature].float() for feature in self.input_features]
        targets = [batch[feature].float() for feature in self.output_features]
        targets = torch.stack(targets).squeeze(0)
        # Print the devices of inputs and targets
        for i, input_tensor in enumerate(inputs):
            self.debugging.log(f"[_prepare_batch] Input feature {i} device: {input_tensor.device}", level=2)
        self.debugging.log(f"[_prepare_batch] Targets device: {targets.device}", level=2)
        return inputs, targets
    ###
    def training_step(self, train_batch, batch_idx):  
        """Performs a training step."""
        inputs, labels = self._prepare_batch(train_batch)
        rbp_expr, gen_expr = inputs
        outputs = self(rbp_expr, gen_expr)  # Model predictions
        loss = self.criterion(labels, outputs)
        self._update_metrics("train", loss, outputs, labels)
        return loss
    ###
    def validation_step(self, val_batch, batch_idx):
        inputs, labels = self._prepare_batch(val_batch)
        rbp_expr, gen_expr = inputs
        outputs = self(rbp_expr, gen_expr)  
        loss = self.criterion(labels, outputs)
        self._update_metrics("validation", loss, outputs, labels)
    ###
    def test_step(self, test_batch, batch_idx):
        """Performs a test step."""
        inputs, labels = self._prepare_batch(test_batch)
        rbp_expr, gen_expr = inputs
        outputs = self(rbp_expr, gen_expr)   
        loss = self.criterion(labels, outputs)
        self._update_metrics("test", loss, outputs, labels)
    ###
    def _update_metrics(self, data_type, loss, outputs, labels):
        """Updates the metrics dynamically based on data_type: 'train', 'validation', 'test'."""
        self.debugging.log(f"[_update_metrics] Outputs device: {outputs.device}", level=2)
        self.debugging.log(f"[_update_metrics] Labels device: {labels.device}", level=2)
        self.debugging.log(f"[_update_metrics] Loss device: {loss.device}", level=2)
        # Update loss
        getattr(self, f"{data_type}_loss")(loss.cpu())
        # Convert to numpy for sklearn/scipy metrics
        labels_np = labels.flatten().cpu().numpy()
        outputs_np = outputs.flatten().cpu().detach().numpy()
        # Update correlation and R² metrics
        getattr(self, f"{data_type}_corr_pearson")(pearsonr(labels_np, outputs_np)[0])
        getattr(self, f"{data_type}_corr_spearman")(spearmanr(labels_np, outputs_np)[0])
        getattr(self, f"{data_type}_r2")(r2_score(labels_np, outputs_np))
        # Only for test: per-gene Spearman
        if data_type == "test":
            results = spearmanr_per_gene(
                    gene_names=self.gene_names, 
                    getBM=self.getBM, 
                    trans_names=self.trans_names, 
                    outputs=outputs.cpu().detach().numpy(), 
                    labels=labels.cpu().numpy()
            )
            self.test_corr_spearman_per_gene(results["mean_corr"])
            self.test_corr_spearman_per_gene_max(results["mean_corr_max"])
    ###
    def _log_metrics(self, data_type):
        """Logs the metrics for the specified dataset type and returns the computed values."""
        base_metrics = ["loss", "corr_pearson", "corr_spearman", "r2"]
        # additional test metrics
        if data_type == "test":
            base_metrics += ["corr_spearman_per_gene", "corr_spearman_per_gene_max"]
        metrics_dict = {
            f"{data_type}_{metric_name}": getattr(self, f"{data_type}_{metric_name}").compute()
            for metric_name in base_metrics
        }
        self.log_dict(metrics_dict, on_epoch=True, logger=True)
        return metrics_dict
    ###
    def _reset_metrics(self, data_type):
        """Resets metrics for the specified dataset type."""
        base_metrics = ["loss", "corr_pearson", "corr_spearman", "r2"]
        if data_type == "test":
            base_metrics += ["corr_spearman_per_gene", "corr_spearman_per_gene_max"]
        for metric_name in base_metrics:
            getattr(self, f"{data_type}_{metric_name}").reset()
    ###
    def on_validation_epoch_end(self):
        """Called at the end of the validation epoch."""
        # Skip logging and printing during validation sanity check
        if self.trainer.sanity_checking:
            return
        # Log metrics and retrieve computed values
        train_metrics = self._log_metrics("train")
        val_metrics = self._log_metrics("validation")
        device = self.device  # Get the current device
        epoch = self.current_epoch # Get the current epoch number
        # Print metrics and device information in a single line
        if self.trainer.local_rank==0:
            print(f"Epoch end {epoch} | "
            f"Device: {device} | "
            f"Training Loss: {train_metrics['train_loss']:.4f} | "
            f"Training Pearson Corr.: {train_metrics['train_corr_pearson']:.4f} | "
            f"Training Spearman Corr.: {train_metrics['train_corr_spearman']:.4f} | "
            f"Training R²: {train_metrics['train_r2']:.4f} | "
            f"Validation Loss: {val_metrics['validation_loss']:.4f} | "
            f"Validation Pearson Corr.: {val_metrics['validation_corr_pearson']:.4f} | "
            f"Validation Spearman Corr.: {val_metrics['validation_corr_spearman']:.4f} | "
            f"Validation R²: {val_metrics['validation_r2']:.4f} | ")
        # Reset metrics for the next epoch
        self._reset_metrics("train")
        self._reset_metrics("validation")
    ###
    def on_test_epoch_end(self):
        """Called at the end of the test epoch."""
        device = self.device  # Get the current device
        # Log and reset test metrics
        test_metrics = self._log_metrics("test")
        # Print metrics and device information in a single line
        if self.trainer.local_rank==0:
            print(f"Device: {device} | "
                f"Test Loss: {test_metrics['test_loss']:.4f} | "
                f"Test Pearson Corr.: {test_metrics['test_corr_pearson']:.4f} | "
                f"Test Spearman Corr.: {test_metrics['test_corr_spearman']:.4f} | "
                f"Test R²: {test_metrics['test_r2']:.4f} | "
                f"Test Spearman Corr. per Gene: {test_metrics['test_corr_spearman_per_gene']:.4f} | "
                f"Test Spearman Corr. per Gene (max trans): {test_metrics['test_corr_spearman_per_gene_max']:.4f} | ")
        self._reset_metrics("test")
    ###
    def predict_step(self, batch, batch_idx):
        """Performs a prediction step and returns both predictions and targets.
        
        Args:
            batch (dict): Batch of data containing input features and targets.
            batch_idx (int): Index of the batch.
        
        Returns:
            tuple: A tuple containing (predictions, targets).
        """
        inputs, targets = self._prepare_batch(batch)  # Get both inputs and targets
        rbp_expr, gen_expr = inputs
        predictions = self(rbp_expr, gen_expr)
        return predictions, targets  # Return both predictions and targets

# 
class TunablePredictorModel(BaseLightningModule):
    """
    A PyTorch Lightning module designed for hyperparameter optimization of the isoform predictor.
    Args:
        input_size (int): Size of the input features.
        output_size (int): Size of the output features.
        config (object): Class object with the model hyperparameter settings.
        gene_names (List[str]): List of gene names.
        trans_names (List[str]): List of transcript names.
        getBM (pd.DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.
        input_features (Optional[str]): Feature colum in loader to use as input.
        output_features (Optional[str]): Feature column in loader to use as predict.
        verbose (int): Verbosity level.
    
    Hyperparameters:
        - num_hidden_layers (int): Number of hidden layers.
        - hidden1_nodes (int): Number of nodes in the first hidden layer.
        - uniform_nodes (bool): Whether to use uniform nodes across layers.
        - node_shrink_factor (float): Factor to reduce nodes in layers.
        - activation_func (str): Activation function to use (e.g., 'relu', 'tanh').
        - batch_norm_eps (float): Epsilon value for batch normalization.
        - batch_norm_momentum (float): Momentum value for batch normalization.
        - optimizer_name (str): Name of the optimizer (e.g., 'adamW').
        - learning_rate (float): Learning rate for the optimizer. 
    """
    def __init__(self, input_size: int, output_size: int, config: ConfigParser,
                 gene_names: List[str], trans_names: List[str], getBM: pd.DataFrame,
                 input_features: Optional[str] = None, output_features: Optional[str] = None, 
                 verbose: int = 0):
        super().__init__(gene_names, trans_names, getBM, input_features, output_features, verbose)
        self.save_hyperparameters(ignore=['verbose']) # save all the variables passed to init simply by calling 
        self.input_size = input_size
        self.output_size = output_size
        # Hyperparameters
        self.num_hidden_layers = config.get('num_hidden_layers')
        self.hidden1_nodes = config.get('hidden1_nodes')
        self.uniform_nodes = config.get('uniform_nodes')
        self.node_shrink_factor = config.get('node_shrink_factor') 
        self.activation_name = config.get('activation_func')
        self.batch_norm_eps = config.get('batch_norm_eps', 1e-5)   
        self.batch_norm_momentum = config.get('batch_norm_momentum', 0.1) 
        # Training parameters
        self.optimizer_name = config.get('optimizer_name')
        self.learning_rate = config.get('learning_rate')
        # Initialize the variable usage tracking
        self.variable_usage = {
            'hidden1_nodes': False,
            'uniform_nodes': False,
            'node_shrink_factor': False,
            'activation_func': False,
            'batch_norm_eps': False,
            'batch_norm_momentum': False
        }
        # Configure model layers
        self._configure_layers()
        # Update unused variables before saving hyperparameters
        self._update_unused_variables()
        # Save specific hyperparameters into the save_hyperparameters
        self.hparams.num_hidden_layers = self.num_hidden_layers
        self.hparams.hidden1_nodes = self.hidden1_nodes
        self.hparams.uniform_nodes = self.uniform_nodes
        self.hparams.node_shrink_factor = self.node_shrink_factor
        self.hparams.activation_func = self.activation_name
        self.hparams.batch_norm_eps = self.batch_norm_eps
        self.hparams.batch_norm_momentum = self.batch_norm_momentum
        self.hparams.optimizer_name = self.optimizer_name
        self.hparams.learning_rate = self.learning_rate
    def _configure_layers(self):
        """Configures the model layers based on configuration."""
        if self.num_hidden_layers > 0:
            node_count = self.hidden1_nodes  # Number of nodes for the first hidden layer
            self._mark_used_variables()
            # First hidden layer
            self.add_module('hidden_linear_0', nn.Linear(self.input_size, node_count)) # Input size to first layer
            self.add_module('batch_norm_0', nn.BatchNorm1d(node_count, eps=self.batch_norm_eps, momentum=self.batch_norm_momentum))
            self.add_module('activation_0', self._get_activation_module())
            # Subsequent hidden layers
            for i in range(1, self.num_hidden_layers):
                input_size = node_count  # Use the output size of the previous layer
                if self.uniform_nodes:
                    if i == self.num_hidden_layers - 1: # Last layer: apply shrink factor
                        node_count = round(node_count / self.node_shrink_factor)
                else:
                    # Reduce the number of nodes in each layer if `uniform_nodes` is False
                    node_count = round(node_count / self.node_shrink_factor)
                # Add the linear layer --> batch normalization --> activation
                layer = nn.Linear(input_size, node_count)
                self.add_module(f'hidden_linear_{i}', layer)
                self.add_module(f'batch_norm_{i}', nn.BatchNorm1d(layer.out_features, eps=self.batch_norm_eps, momentum=self.batch_norm_momentum))
                self.add_module(f'activation_{i}', self._get_activation_module())
        else:
            node_count = self.input_size  # No hidden layers, use input size directly
        # Configure the output layer
        self.linear_output = nn.Linear(node_count, self.output_size)
        self.add_module('linear_output', self.linear_output)
        # Add the activation layer for the output
        self.output_activation = nn.Sigmoid()
        self.add_module('output_activation', self.output_activation)
    ###
    def _mark_used_variables(self):
        """Marks the variables as used based on the current configuration if
        number of hidden layers is greater to zero."""
        self.variable_usage['hidden1_nodes'] = True   
        self.variable_usage['activation_func'] = True 
        self.variable_usage['batch_norm_eps'] = True 
        self.variable_usage['batch_norm_momentum'] = True 
        if self.num_hidden_layers > 1:
            self.variable_usage['node_shrink_factor'] = True   
        if self.num_hidden_layers >= 3:
            self.variable_usage['uniform_nodes'] = True 
    ###
    def _update_unused_variables(self):
        """Updates unused variables with the string 'unused'."""
        for var in self.variable_usage:
            if not self.variable_usage[var]:
                setattr(self, var, 'unused')
    ###                    
    def _get_activation_module(self):
        """Returns the activation layer based on the given name."""
        if self.activation_name == "relu":
            return nn.ReLU()
        elif self.activation_name == "tanh":
            return nn.Tanh()
        elif self.activation_name == "sigmoid":
            return nn.Sigmoid()
        else:
            self.log("Invalid activation_layer. Supported options are 'relu', 'tanh', and 'sigmoid'.")
    ###
    def configure_optimizers(self):
        """Configures the optimizer based on the provided name and learning rate.
        Returns:
            torch.optim.Optimizer: Configured optimizer instance.
        """
        if self.optimizer_name == 'sgd90':
            return torch.optim.SGD(self.parameters(), lr=self.learning_rate, momentum=0.9)
        elif self.optimizer_name == 'asgd':
            return torch.optim.ASGD(self.parameters(), lr=self.learning_rate, lambd=0.0001, alpha=0.75)
        elif self.optimizer_name == 'adam':
            return torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        elif self.optimizer_name == 'adagrad':
            return torch.optim.Adagrad(self.parameters(), lr=self.learning_rate)
        elif self.optimizer_name == 'adadelta':
            return torch.optim.Adadelta(self.parameters(), lr=self.learning_rate)
        elif self.optimizer_name == 'adamW':
            return torch.optim.AdamW(self.parameters(), lr=self.learning_rate)
        else:
            self.log(f"Unsupported optimizer '{self.optimizer_name}'. Valid options: ['sgd90', 'asgd', 'adam', 'adagrad', 'adadelta', 'adamW']")
    ###
    def forward(self, rbp_expr, gen_expr): # this was updated to work with the new modules
        """Defines the forward pass of the model.
        
        Args:
             bp_expr (torch.Tensor): Input features (e.g., RBP).
            gen_expr (torch.Tensor): Additional gene input for the final output scaling.
        
        Returns:
            torch.Tensor: Predicted transcript expression in log2(TPM+1).
        """
        x = rbp_expr  
        # Pass through all hidden layers
        for i in range(self.num_hidden_layers):
            x = self._modules[f'hidden_linear_{i}'](x) # Linear layer
            x = self._modules[f'batch_norm_{i}'](x) # Batch normalization
            x = self._modules[f'activation_{i}'](x) # Activation
        # Final output layer
        out = self._modules['linear_output'](x)
        out = self.output_activation(out)
        # Log-transform and scale
        out = torch.log2((out * gen_expr) + 1)
        return out



class PredictorModel(BaseLightningModule): # Despues de los resultados de la optimizacion estoy hay que cambiar.
    """A PyTorch Lightning module used to train the isoform predictor and to serve as the 
    reference model for downstream DeepRBP explainability analysis.
    
    Args:
        input_size (int): Size of the input features.
        output_size (int): Size of the output features.
        gene_names (List[str]): List of gene names.
        trans_names (List[str]): List of transcript names.
        getBM (pd.DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.
        input_features (Optional[str]): Feature colum in loader to use as input.
        output_features (Optional[str]): Feature column in loader to use as predict.
        verbose (int): Verbosity level.
    
    Hyperparameters:
        - learning_rate (float): Learning rate for the optimizer. 
    """
    def __init__(self, input_size: int, output_size: int, 
                 gene_names: List[str], trans_names: List[str], getBM: pd.DataFrame,
                 input_features: Optional[str] = None, output_features: Optional[str] = None, 
                 verbose: int = 0):
        
        super().__init__(gene_names, trans_names, getBM, input_features, output_features, verbose)
        self.save_hyperparameters(ignore=['verbose']) # save all the variables passed to init simply by calling 
        
        self.input_size = input_size
        self.output_size = output_size
        
        ### (CHANGE THIS TO WRITE THE FINAL MODEL)
        self.learning_rate = 0.0001

        # Define the actual neural network
        self.abundance_estimator = nn.Sequential(
            nn.Linear(input_size, 128),
            nn.BatchNorm1d(128, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.ReLU(),
            nn.Linear(64, output_size),
            nn.Sigmoid()
        )

    def configure_optimizers(self):
        """Configures the optimizer"
        Returns:
            torch.optim.Optimizer: Configured optimizer instance.
        """
        return torch.optim.AdamW(self.parameters(), lr=self.learning_rate)

    def forward(self, rbp_expr, gen_expr):
        """Defines the forward pass of the model.
            Args:
            rbp_expr (torch.Tensor): Input features (e.g., RBP).
            gen_expr (torch.Tensor): Additional gene input for the final output scaling.
            
            Returns:
                torch.Tensor: Predicted transcript expression in log2(TPM+1).
        """
        x = self.abundance_estimator(rbp_expr)
        out = torch.log2((x * gen_expr) + 1)
        return out


