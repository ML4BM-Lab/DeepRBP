# src/deeprbp/training_module/model.py

import os
import torch
import torch.nn as nn
import torch.nn.functional as F
#from typing import Any, Dict
import yaml
from colorama import Fore, init

from ..util.logger import Logger

class PredictorModel(nn.Module):
    """Create a neural network with multiple hidden layers allowing for flexible configuration of each layer's node count. 
    This network architecture is defined by specifying the number of hidden layers and the number of nodes in the 
    first hidden layer 'hidden1_nodes'. Each subsequent hidden layer's node count decreases by a specified division factor 'node_shrink_factor'. 

    Additionally, the user can opt for all hidden layers to have the same number of nodes by setting `uniform_nodes` to True. 
    If False, each hidden layer will have half the number of nodes of the previous layer by default, with the option to adjust 
    this division factor. 

    The activation function for the hidden layers can be customized.

    This neural network, DeepRBP, aims to predict relative transcript abundances (as a percentage) using the expression of the RBPs.
    To achieve this, the RBP expression data is passed through fully-connected layers followed by batch normalization layers, as 
    described in (Ioffe and Szegedy, 2015), and activation layers. At the output layer, the sigmoid activation function is applied 
    to obtain the percentage of the gene expression that contributes to each transcript.

    The training of the model is performed by minimizing the Mean Squared Error (MSE) between the predicted and real values. 
    A caveat to this formulation is the inaccurate representation of low-expressed genes, which can harm the training of the model. 
    For this reason, we also introduce the expression of the genes as input, and compute instead transcript abundances, i.e., the transcripts’ 
    expression. In particular, the predicted percentages are multiplied by the corresponding gene expression (in TPMs), yielding 
    transcript abundances also in TPMs. 

    Args:
        config (object class): Class object with the model configuration.
        input_size (int): Size of the input features.
        output_size (int): Size of the output features.

    Attributes:
        input_size (int): Number of input features.
        output_size (int): Number of output features.
        num_hidden_layers (int): Number of hidden layers.
        hidden1_nodes (int): Number of nodes in the first hidden layer.
        uniform_nodes (bool): Whether all hidden layers (except the last one) have the same number of nodes.
        node_shrink_factor (float): Factor by which the number of nodes decreases per layer (if `uniform_nodes` is False).
        activation_func (str): Activation function to use in hidden layers (e.g., 'relu', 'tanh').
        optimizer_name (str): Name of the optimizer to use (e.g., 'adam', 'sgd90').
        learning_rate (float): Learning rate for the optimizer.
        batch_norm_eps (float): Epsilon value for numerical stability in batch normalization, preventing division by zero.
        batch_norm_momentum (float): Momentum for the running mean and variance in batch normalization, controlling how quickly these values are updated.
    """
    def __init__(self, config, input_size: int = None, output_size: int = None, verbose: int = 1):
        super(PredictorModel, self).__init__()
        self.logger = Logger(verbose=verbose)
        init(autoreset=True)
        # Model configuration
        self.config = config
        self.input_size = input_size 
        self.output_size = output_size
        self.num_hidden_layers = self.config.get('num_hidden_layers')
        self.hidden1_nodes = self.config.get('hidden1_nodes')
        self.uniform_nodes = self.config.get('uniform_nodes')
        self.node_shrink_factor = self.config.get('node_shrink_factor') 
        self.activation_func = self.config.get('activation_func')
        self.learning_rate = self.config.get('learning_rate')
        self.optimizer_name = self.config.get('optimizer_name')
        self.batch_norm_eps = self.config.get('batch_norm_eps', 1e-5)   
        self.batch_norm_momentum = self.config.get('batch_norm_momentum', 0.1) 
        # Initialize the variable usage tracking
        self.variable_usage = {
            'hidden1_nodes': False,
            'uniform_nodes': False,
            'node_shrink_factor': False,
            'activation_func': False,
            'batch_norm_eps': False,
            'batch_norm_momentum': False
        }
        # Configure layers and optimizer
        self.logger.log("Initializing layers and optimizer...")
        self._configure_layers()
        #self.optimizer = self._configure_optimizer(self.optimizer_name, self.learning_rate)
        self.optimizer = None 
        # Check the model layers to ensure there are no layers with out_features equal to 0
        self._check_layer_outputs()
        # Print and update unused variables and print used variables
        self._print_unused_variables()
        self._print_used_variables()
        self._update_unused_variables()
        self._print_model_architecture()
    def _configure_layers(self):
        """Configures the hidden and output layers based on model configuration."""
        if self.num_hidden_layers > 0:
            node_count = self.hidden1_nodes  # Number of nodes for the first hidden layer
            self._mark_used_variables()
            # First hidden layer
            self.add_module('hidden_linear_0', nn.Linear(self.input_size, node_count)) # Input size to first layer
            self.add_module('batch_norm_0', nn.BatchNorm1d(node_count, eps=self.batch_norm_eps, momentum=self.batch_norm_momentum))
            self.add_module('activation_0', self._get_activation_module(self.activation_func))
            # Subsequent hidden layers
            for i in range(1, self.num_hidden_layers):
                input_size = node_count  # Use the output size of the previous layer
                if self.uniform_nodes:
                    if i == self.num_hidden_layers - 1: # Last layer: apply shrink factor
                        node_count = round(node_count / self.node_shrink_factor)
                else:
                    # Reduce the number of nodes in each layer if `uniform_nodes` is False
                    node_count = round(node_count / self.node_shrink_factor)
                # Add the layer, batch normalization, and activation
                layer = nn.Linear(input_size, node_count)
                self.add_module(f'hidden_linear_{i}', layer)
                self.add_module(f'batch_norm_{i}', nn.BatchNorm1d(layer.out_features, eps=self.batch_norm_eps, momentum=self.batch_norm_momentum))
                self.add_module(f'activation_{i}', self._get_activation_module(self.activation_func))
        else:
            node_count = self.input_size  # No hidden layers, use input size directly
        # Configure the output layer
        self.linear_output = nn.Linear(node_count, self.output_size)
        self.add_module('linear_output', self.linear_output)
        # Add the activation layer for the output
        self.output_activation = nn.Sigmoid()
        self.add_module('output_activation', self.output_activation)
    def _print_model_architecture(self):
        """Prints a high-level overview of the model architecture."""
        print("Model architecture:")
        print(self)
    def _check_layer_outputs(self):
        """
        Checks that no layer in the model has out_features equal to 0.

        Raises:
            ValueError: If any layer has out_features equal to 0.
        """
        for name, layer in self.named_children():
            if isinstance(layer, nn.Linear) and layer.out_features == 0:
                raise ValueError(
                    f"Error: Layer '{name}' has out_features equal to 0, which is invalid for training. "
                    f"Please try increasing 'hidden1_nodes' or decreasing 'node_shrink_factor' "
                    f"based on the number of hidden layers ('num_hidden_layers') you are using."
                )    
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
    def _print_unused_variables(self):
        """Prints the variables that are not being used."""
        init(autoreset=True)
        not_used = [var for var, used in self.variable_usage.items() if not used]  
        title = "📊 Unused Variables Report 📊"
        separator = "-----------------------------------"
        # Print title and separator
        print(Fore.CYAN + title)
        print(Fore.CYAN + separator)
        if not_used:
            print(Fore.RED + "⚠️ The following variables are not being used:")
            for var in not_used:
                print(Fore.LIGHTYELLOW_EX + f"🔹 {var}")
        else:
            print(Fore.GREEN + "✅ All relevant variables are in use.")
        print(Fore.CYAN + separator)
    def _print_used_variables(self):
        """Prints the variables that are being used."""
        used_vars = {  # Only include variables that are marked True in variable_usage
            'hidden1_nodes': self.hidden1_nodes if self.variable_usage['hidden1_nodes'] else None,
            'uniform_nodes': self.uniform_nodes if self.variable_usage['uniform_nodes'] else None,
            'node_shrink_factor': self.node_shrink_factor if self.variable_usage['node_shrink_factor'] else None,
            'activation_func': self.activation_func if self.variable_usage['activation_func'] else None,
            'batch_norm_eps': self.batch_norm_eps if self.variable_usage['batch_norm_eps'] else None,
            'batch_norm_momentum': self.batch_norm_momentum if self.variable_usage['batch_norm_momentum'] else None,
            'num_hidden_layers': self.num_hidden_layers,  # Included without tracking
            'optimizer_name': self.optimizer_name,  # Included without tracking
            'learning_rate': self.learning_rate  # Included without tracking
        }
        title = "✅ Used Variables Report ✅"
        separator = "-----------------------------------"
        # Print title and separator
        print(Fore.CYAN + title)
        print(Fore.CYAN + separator)
        for var, value in used_vars.items():
            if value is not None:  # Only print values that are not None
                print(Fore.LIGHTGREEN_EX + f"🔹 {var}: {value}")
        print(Fore.CYAN + separator)
    def _update_unused_variables(self):
        """Updates unused variables with a placeholder character."""
        # Map variable names to their corresponding attributes
        variable_map = {
            'hidden1_nodes': 'hidden1_nodes',
            'uniform_nodes': 'uniform_nodes',
            'node_shrink_factor': 'node_shrink_factor',
            'activation_func': 'activation_func',
            'batch_norm_eps': 'batch_norm_eps',
            'batch_norm_momentum': 'batch_norm_momentum'
        }
        for var in self.variable_usage:
            if not self.variable_usage[var]:  # If the variable is marked as unused ###
                if var in self.config.config_data:  # Use the ConfigParser's config_data
                    self.config.update(var, 'unused')  # Update using the ConfigParser update method
                else:
                    print(f"Warning: {var} not found in configuration.")
                # Update the corresponding attribute to None
                if var in variable_map:
                    setattr(self, variable_map[var], None)
    def configure_optimizer(self):
        """Configures the optimizer after the model has been moved to the appropriate device."""
        if self.optimizer is None:  # Only configure if not already done
            self.optimizer = self._configure_optimizer(self.optimizer_name, self.learning_rate)
    def _configure_optimizer(self, optimizer_name, learning_rate):
        """Configures the optimizer based on the provided name and learning rate.
        Args:
            optimizer_name (str): Name of the optimizer to use ('sgd90', 'sgd70', etc.).
            learning_rate (float): Learning rate for the optimizer.
        Returns:
            torch.optim.Optimizer: Configured optimizer instance.
        """
        if optimizer_name == 'sgd90':
            return torch.optim.SGD(self.parameters(), lr=learning_rate, momentum=0.9)
        elif optimizer_name == 'asgd':
            return torch.optim.ASGD(self.parameters(), lr=learning_rate, lambd=0.0001, alpha=0.75)
        elif optimizer_name == 'adam':
            return torch.optim.Adam(self.parameters(), lr=learning_rate)
        elif optimizer_name == 'adagrad':
            return torch.optim.Adagrad(self.parameters(), lr=learning_rate)
        elif optimizer_name == 'adadelta':
            return torch.optim.Adadelta(self.parameters(), lr=learning_rate)
        elif optimizer_name == 'adamW':
            return torch.optim.AdamW(self.parameters(), lr=learning_rate)
        else:
            self.logger.error(f"Unsupported optimizer '{optimizer_name}'. Valid options: ['sgd90', 'asgd', 'adam', 'adagrad', 'adadelta', 'adamW']")
    def _get_activation_module(self, activation_name):
        """Returns the activation layer based on the given name."""
        if activation_name == "relu":
            return nn.ReLU()
        elif activation_name == "tanh":
            return nn.Tanh()
        elif activation_name == "sigmoid":
            return nn.Sigmoid()
        else:
            self.logger.error("Invalid activation_layer. Supported options are 'relu', 'tanh', and 'sigmoid'.")
    def forward(self, xb, gb):
        """Defines the forward pass of the model.
        
        Args:
            xb (torch.Tensor): Input features (e.g., RBP and gene expression).
            gb (torch.Tensor): Additional gene input for the final output scaling.
        
        Returns:
            torch.Tensor: Predicted transcript abundance.
        """
        x = xb
        # Pass through all modules except the last activation
        for _, module in list(self.named_children())[:-1]:  # Exclude the last module (output_activation)
            x = module(x)
        # Apply the output activation function and the specified transformation
        out = torch.log2((self.output_activation(x) * gb) + 1)
        return out
    def train_step(self, inputs, targets, gen_expr):
        """Performs a single training step (forward + backward pass).
        Args:
            inputs (torch.Tensor): Input features.
            targets (torch.Tensor): Ground truth labels.
            gen_expr (torch.Tensor): Gene expression data.
        Returns:
            torch.Tensor: Training loss (detached).
        """
        out = self(inputs, gen_expr) # Forward pass
        loss = F.mse_loss(out, targets) # Calculate loss
        loss.backward() # Backpropagation
        self.optimizer.step() # Update weights
        self.optimizer.zero_grad() # Reset gradients
        return loss.detach()
    def validate_step(self, inputs, targets, gen_expr):
        """Performs a single validation step (forward pass only).
        Args:
            inputs (torch.Tensor): Input features.
            targets (torch.Tensor): Ground truth labels.
            gen_expr (torch.Tensor): Gene expression data.
        Returns:
            torch.Tensor: Validation loss (detached).
        """
        with torch.no_grad():
            out = self(inputs, gen_expr) # Forward pass
            val_loss = F.mse_loss(out, targets)  # Calculate loss
        return val_loss.detach()
    def save_model(self, output_dir, model_name):
        """Saves the model to the specified directory."""
        self.logger.log(f"Saving model to {output_dir} with name {model_name}...")
        os.makedirs(output_dir, exist_ok=True)
        # Save model weights
        torch.save(self.state_dict(), os.path.join(output_dir, model_name))
        self.logger.log(f"Model saved successfully in {output_dir}")
    def save_config(self, output_dir):
        """Saves the model configuration to a YAML file in the specified directory."""
        self.logger.log(f"Saving configuration to {output_dir}...")
        os.makedirs(output_dir, exist_ok=True)
        # Save configuration to a YAML file
        config_path = os.path.join(output_dir, 'config.yaml')
        with open(config_path, 'w') as config_file:
            yaml.dump(self.config.config_data, config_file)  # Use config_data from ConfigParser
        self.logger.log(f"Configuration saved successfully in {config_path}")
    @classmethod
    def load_model(cls, path_to_weights, config, input_size, output_size, device=None):
        """
        Class method to initialize the model and load pre-trained weights.

        Args:
            path_to_weights (str): Path to the file containing the pre-trained weights.
            config Configuration class for the model.
            input_size (int): Input size of the model.
            output_size (int): Output size of the model.
            device (torch.device, optional): The device to load the model onto (CPU or GPU). If None, it will be set automatically.

        Returns:
            PredictorModel: An instance of PredictorModel with pre-trained weights loaded.
        """
        # Determine the device if not provided
        if device is None:
            #device = torch.device('cuda:0' if config.get('cuda') and torch.cuda.is_available() else 'cpu')
            device = torch.device('cuda' if config.get('cuda') and torch.cuda.is_available() else 'cpu')
            print(f"⚠️ No device specified. Using default device: {device}.")  # Inform the user about the default device
        # Initialize the model with the provided config
        logger = Logger(config.get('verbose', 1))  # Initialize a logger
        logger.log(f"Loading model from weights at {path_to_weights}...")
        # Initialize the model with the provided config
        model = cls(config, input_size, output_size)
        # Load pre-trained weights
        try:
            state_dict = torch.load(path_to_weights, map_location=device, weights_only=True)   
            model.load_state_dict(state_dict)
            logger.log("Pre-trained weights loaded successfully.")
        except Exception as e:
            logger.error(f"Failed to load weights: {e}", exception_type=RuntimeError)
        return model
    