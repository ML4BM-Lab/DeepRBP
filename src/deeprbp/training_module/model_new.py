
# voy a intentar aplicar el model con lighting

import pandas as pd
# Pytorch modules
import torch
import torch.nn as nn
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score

# Pytorch-Lightning
import lightning as L
from torchmetrics import MeanMetric
#from torchmetrics.regression import SpearmanCorrCoef, PearsonCorrCoef, R2Score

from .evaluation import spearmanr_per_gene

from typing import List, Any, Optional
import warnings
# Suppress FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)

class PredictorModel(L.LightningModule):
    """
    Args:
        output_folder (str):  
        config (object): Class object with the model hyperparameter settings.
        input_size (int): Size of the input features.
        output_size (int): Size of the output features.
        gene_names (List[str]): A list of gene IDs used in the model.
        trans_names (List[str]): A list of transcript IDs used in the model.
        getBM (pd.DataFrame): DataFrame containing mapping of Gene_ID to Transcript_ID.

    Hyperparameters:
        - activation_func (str): Activation function to use (e.g., 'relu', 'tanh').
        - batch_norm_eps (float): Epsilon value for batch normalization.
        - batch_norm_momentum (float): Momentum value for batch normalization.
        - hidden1_nodes (int): Number of nodes in the first hidden layer.
        - learning_rate (float): Learning rate for the optimizer.
        - node_shrink_factor (float): Factor to reduce nodes in layers.
        - num_hidden_layers (int): Number of hidden layers.
        - optimizer_name (str): Name of the optimizer (e.g., 'adamW').
        - uniform_nodes (bool): Whether to use uniform nodes across layers.
    """
    def __init__(self, config: Any, input_size: int, output_size: int, 
                 gene_names: List[str], trans_names: List[str], getBM: pd.DataFrame, 
                 input_features: Optional[str] = None, output_features: Optional[str] = None):  
        super().__init__()
        # Set the output folder path
        #self.output_folder = output_folder + self.fold + "/"
        #self.ckpt_path = ckpt_path
        # Model configuration
        self.input_size = input_size 
        self.output_size = output_size
        self.input_features = input_features if input_features is not None else ('scaled_rbp_df', 'gene_df')
        self.output_features = output_features if output_features is not None else ('isoform_df',)   
        self.gene_names = gene_names
        self.trans_names = trans_names
        self.getBM = getBM
        # Hyperparameters
        self.num_hidden_layers = config.get('num_hidden_layers')
        self.hidden1_nodes = config.get('hidden1_nodes')
        self.uniform_nodes = config.get('uniform_nodes')
        self.node_shrink_factor = config.get('node_shrink_factor') 
        self.activation_name = config.get('activation_func')
        self.batch_norm_eps = config.get('batch_norm_eps', 1e-5)   
        self.batch_norm_momentum = config.get('batch_norm_momentum', 0.1) 
        self.criterion = nn.MSELoss()
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
        # Save specific hyperparameters
        self.hparams.num_hidden_layers = self.num_hidden_layers
        self.hparams.hidden1_nodes = self.hidden1_nodes
        self.hparams.uniform_nodes = self.uniform_nodes
        self.hparams.node_shrink_factor = self.node_shrink_factor
        self.hparams.activation_func = self.activation_name
        self.hparams.batch_norm_eps = self.batch_norm_eps
        self.hparams.batch_norm_momentum = self.batch_norm_momentum
        self.hparams.optimizer_name = self.optimizer_name
        self.hparams.learning_rate = self.learning_rate
        # Intialize metrics
        self.metrics = self._initialize_metrics()
    ###
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
    def _initialize_metrics(self):
        """Initializes metrics for training, validation, and testing. 
        This method sets up the following metrics using `MeanMetric` from `torchmetrics`:
        - **Mean Squared Error (MSE)**: It's the loss functionm, measures the average squared difference between log2(tpm+1) predictions vs actual values.
        - **Pearson Correlation**: Assesses linear correlation between predictions and actual values, ranging from -1 to 1.
        - **Spearman Correlation**: Evaluates rank correlation, indicating the strength of monotonic relationships, also ranging from -1 to 1.
        - **R² (Coefficient of Determination)**: Indicates the proportion of variance explained by the model, with values from 0 to 1.
        """ 
        return {
            "train": {
                "loss": MeanMetric(),
                "corr_pearson": MeanMetric(),
                "corr_spearman": MeanMetric(),
                "r2": MeanMetric(),
            },
            "val": { # igual mejor quitarlo del val que ocupa mucho el cabron de computacion
                "loss": MeanMetric(),
                "corr_pearson": MeanMetric(),
                "corr_spearman": MeanMetric(),
                "r2": MeanMetric(),
            },
            "test": {
                "loss": MeanMetric(),
                "corr_pearson": MeanMetric(),
                "corr_spearman": MeanMetric(),
                "r2": MeanMetric(),
                "corr_spearman_per_gene": MeanMetric(),
                "corr_spearman_per_gene_max": MeanMetric(),
            }
        }
    ###
    def forward(self, rbp_expr, gen_expr): # this was updated to work with the new modules
        """Defines the forward pass of the model.
        
        Args:
             bp_expr (torch.Tensor): Input features (e.g., RBP).
            gen_expr (torch.Tensor): Additional gene input for the final output scaling.
        
        Returns:
            torch.Tensor: Predicted transcript abundance.
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
        return inputs, targets
    ###
    def _flatten_and_to_numpy(self, tensor):
        """Convert a tensor to a flattened numpy array."""
        return tensor.flatten().cpu().detach().numpy()
    ###
    def training_step(self, train_batch, batch_idx):  
        """Performs a training step."""
        # Prepare the input and labels
        inputs, labels = self._prepare_batch(train_batch)
        rbp_expr, gen_expr = inputs
        # Forward pass
        outputs = self(rbp_expr, gen_expr)  # Model predictions
        loss = self.criterion(outputs, labels)
        # Update training metrics
        self._update_metrics("train", loss, outputs, labels)
        return loss
    ###
    def validation_step(self, val_batch, batch_idx):
        # Prepare the input and labels
        inputs, labels = self._prepare_batch(val_batch)
        rbp_expr, gen_expr = inputs
        # Forward pass
        outputs = self(rbp_expr, gen_expr)  # Model predictions
        loss = self.criterion(outputs, labels)
        # Update validation metrics
        self._update_metrics("val", loss, outputs, labels)
    ###
    def test_step(self, test_batch, batch_idx):
        """Performs a test step."""
        inputs, labels = self._prepare_batch(test_batch)
        rbp_expr, gen_expr = inputs
        outputs = self(rbp_expr, gen_expr)  # Model predictions
        loss = self.criterion(outputs, labels)
        # Update test metrics
        self._update_metrics("test", loss, outputs, labels)
    ###
    def _update_metrics(self, data_type, loss, outputs, labels):
        """Updates the specified metrics based on the type."""
        self.metrics[data_type]["loss"](loss)
        self.metrics[data_type]["corr_pearson"](pearsonr(self._flatten_and_to_numpy(outputs), self._flatten_and_to_numpy(labels))[0])
        self.metrics[data_type]["corr_spearman"](spearmanr(self._flatten_and_to_numpy(outputs), self._flatten_and_to_numpy(labels))[0])
        self.metrics[data_type]["r2"](r2_score(self._flatten_and_to_numpy(outputs), self._flatten_and_to_numpy(labels)))
        if data_type == "test":
            results = spearmanr_per_gene(gene_names=self.gene_names, getBM=self.getBM, trans_names=self.trans_names, 
                                          outputs=outputs.cpu().detach().numpy(), labels=labels.cpu().numpy())
            self.metrics[data_type]["corr_spearman_per_gene"](results["mean_corr"])
            self.metrics[data_type]["corr_spearman_per_gene_max"](results["mean_corr_max"])
    ###
    def _log_metrics(self, data_type):
        """Logs the metrics for the specified dataset type and returns the computed values."""
        metrics_dict = {
            f"{data_type}_loss": self.metrics[data_type]["loss"].compute(),
            f"{data_type}_corr_pearson": self.metrics[data_type]["corr_pearson"].compute(),
            f"{data_type}_corr_spearman": self.metrics[data_type]["corr_spearman"].compute(),
            f"{data_type}_r2": self.metrics[data_type]["r2"].compute(),
        }
        # Check if we're logging test metrics
        if data_type == "test":
            metrics_dict[f"{data_type}_corr_spearman_per_gene"] = self.metrics[data_type]["corr_spearman_per_gene"].compute()
            metrics_dict[f"{data_type}_corr_spearman_per_gene_max"] = self.metrics[data_type]["corr_spearman_per_gene_max"].compute()
        self.log_dict(metrics_dict, on_epoch=True, logger=True)
        return metrics_dict
    ###
    def _reset_metrics(self, data_type):
        """Resets metrics for the specified dataset type."""
        for metric in self.metrics[data_type].values():
            metric.reset()
    ###
    def on_validation_epoch_end(self):
        """Called at the end of the validation epoch."""
        # Skip logging and printing during validation sanity check
        if self.trainer.sanity_checking:
            return
        # Log metrics and retrieve computed values
        train_metrics = self._log_metrics("train")
        val_metrics = self._log_metrics("val")
        # Get the current device
        device = self.device  # Get the current device
        # Get the current epoch number
        epoch = self.current_epoch 
        # Print metrics and device information in a single line
        print(f"Epoch end {epoch} | "
          f"Device: {device} | "
          f"Training Loss: {train_metrics['train_loss']:.4f} | "
          f"Training Pearson Corr.: {train_metrics['train_corr_pearson']:.4f} | "
          f"Training Spearman Corr.: {train_metrics['train_corr_spearman']:.4f} | "
          f"Training R²: {train_metrics['train_r2']:.4f} | "
          f"Validation Loss: {val_metrics['val_loss']:.4f} | "
          f"Validation Pearson Corr.: {val_metrics['val_corr_pearson']:.4f} | "
          f"Validation Spearman Corr.: {val_metrics['val_corr_spearman']:.4f} | "
          f"Validation R²: {val_metrics['val_r2']:.4f} | ")
        # Reset metrics for the next epoch
        self._reset_metrics("train")
        self._reset_metrics("val")
    ###
    def on_test_epoch_end(self):
        """Called at the end of the test epoch."""
        device = self.device  # Get the current device
        # Log and reset test metrics
        test_metrics = self._log_metrics("test")
        # Print metrics and device information in a single line
        print(f"Device: {device} | "
              f"Test Loss: {test_metrics['train_loss']:.4f} | "
              f"Test Pearson Corr.: {test_metrics['train_corr_pearson']:.4f} | "
              f"Test Spearman Corr.: {test_metrics['train_corr_spearman']:.4f} | "
              f"Test R²: {test_metrics['train_r2']:.4f} | "
              f"Test Spearman Corr. per Gene: {test_metrics['corr_spearman_per_gene']:.4f} | "
              f"Test Spearman Corr. per Gene (max trans): {test_metrics['corr_spearman_per_gene_max']:.4f} | ")
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


   


