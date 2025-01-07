#class DeepRBPModel
import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Any

from logger import Logger
from config_loader import ConfigParser
from processing import DataImporter, DatasetLoader, Scaler
from deeplift_handler import DeepLiftHandler
from utils import ensure_directory_exists

class PredictorModel_new(nn.Module): # version new, esta es con la que nos vamos a quedar
    """Create a neural network with multiple hidden layers allowing for flexible configuration of each layer's node count. 
    This network architecture is defined by specifying the number of hidden layers and the maximum number of nodes in the 
    first hidden layer. Each subsequent hidden layer's node count decreases by a specified division factor. 

    Additionally, the user can opt for all hidden layers to have the same number of nodes by setting `same_num_nodes` to True. 
    If False, each hidden layer will have half the number of nodes of the previous layer by default, with the option to adjust 
    this division factor. 

    The activation function for the hidden layers can be customized, with ReLU being the default.

    This neural network, DeepRBP, is designed for predicting transcript abundance given RNA-binding protein (RBP) and gene expression data.

    Args:
        config (Any): Dictionary-like object with the model configuration.
        input_size (int, optional): Size of the input features. Overrides `config['input_size']` if provided.
        output_size (int, optional): Size of the output features. Overrides `config['output_size']` if provided.

    Attributes:
        input_size (int): Number of input features.
        output_size (int): Number of output features.
        num_hidden_layers (int): Number of hidden layers.
        max_nodes (int): Maximum number of nodes in the first hidden layer.
        uniform_nodes (bool): Whether all hidden layers have the same number of nodes.
        node_shrink_factor (float): Factor by which the number of nodes decreases per layer (if `uniform_nodes` is False).
        activation_func (str): Activation function to use in hidden layers (e.g., 'relu', 'tanh').
        optimizer_name (str): Name of the optimizer to use (e.g., 'adam', 'sgd90').
        learning_rate (float): Learning rate for the optimizer.
    """
  
    def __init__(self, config: Any, input_size: int = None, output_size: int = None, verbose: int = 1):
        super(PredictorModel_new, self).__init__()
        self.logger = Logger(verbose=verbose)

        # Model configuration
        self.input_size = input_size if input_size is not None else config['input_size']
        self.output_size = output_size if output_size is not None else config['output_size']
        self.num_hidden_layers = config.get('num_hidden_layers', 1)
        self.max_nodes = config['max_node']
        self.uniform_nodes = config['uniform_nodes']
        self.node_shrink_factor = config['node_shrink_factor']
        self.activation_func = config['activation_func']
        self.learning_rate = config['learning_rate']
        self.optimizer_name = config['optimizer_name']

        # Configure layers and optimizer
        self.logger.log("Initializing layers and optimizer...")
        self._configure_layers()
        self.optimizer = self._configure_optimizer(self.optimizer_name, self.learning_rate)

    def _configure_layers(self):
        """Configures the hidden and output layers based on model configuration."""
        if self.num_hidden_layers > 0:
            self.logger.log(f"Using a model with {self.num_hidden_layers} hidden layers")

            # Initialize the number of nodes for the first hidden layer
            node_count = self.max_nodes
            
            for i in range(self.num_hidden_layers):
                # Determine the input size for the current layer
                input_size = self.input_size if i == 0 else node_count

                # Configure hidden layers
                if self.uniform_nodes:
                    if i == self.num_hidden_layers - 1: # Last layer: apply shrink factor
                        node_count = round(node_count / self.node_shrink_factor)
                    layer = nn.Linear(input_size, node_count)
                else:
                    # Reduce the number of nodes in each layer if `uniform_nodes` is False
                    layer = nn.Linear(input_size, node_count)
                    node_count = round(node_count / self.node_shrink_factor)
                
                # Add the layer, batch normalization, and activation
                self.add_module(f'hidden_linear_{i}', layer)
                bn_layer = nn.BatchNorm1d(layer.out_features)
                self.add_module(f'batch_norm_{i}', bn_layer)

                # Add activation function (e.g., ReLU, Tanh, etc.)
                activation_name = self._get_activation_module(self.activation_func)
                self.add_module(f'activation_{i}', activation_name)

            # The input size for the final layer is the number of nodes in the last hidden layer
            final_hidden_layer_input_size = node_count
        else:
            self.logger.log("Using a model with zero hidden layers")
            final_hidden_layer_input_size = self.input_size  # No hidden layers, use input size directly

        # Configure the output layer
        self.linear_output = nn.Linear(final_hidden_layer_input_size, self.output_size)
        self.add_module('linear_output', self.linear_output)

        # Add the activation layer for the output
        self.output_activation = nn.Sigmoid()
        self.add_module('output_activation', self.output_activation)

    def _configure_optimizer(self, optimizer_name, learning_rate):
        """Configures the optimizer based on the provided name and learning rate.
        Args:
            optimizer_name (str): Name of the optimizer to use ('sgd90', 'sgd70', etc.).
            learning_rate (float): Learning rate for the optimizer.
        Returns:
            torch.optim.Optimizer: Configured optimizer instance.
        """
        self.logger.log(f"Configuring optimizer: {optimizer_name}")
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
        self.logger.log(f"Configuring activation function: {activation_name}")
        if activation_name == "relu":
            return nn.ReLU()
        elif activation_name == "tanh":
            return nn.Tanh()
        elif activation_name == "sigmoid":
            return nn.Sigmoid()
        else:
            raise ValueError("Invalid activation_layer. Supported options are 'relu', 'tanh', and 'sigmoid'.")
        
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
        torch.save(self.state_dict(), os.path.join(output_dir, model_name))
        self.logger.log(f"Model saved successfully in {output_dir}")

    @classmethod
    def load_model(cls, path_to_weights, config: Any, input_size: int = None, output_size: int = None):
        """
        Class method to initialize the model and load pre-trained weights.

        Args:
            path_to_weights (str): Path to the file containing the pre-trained weights.
            config (Any): Configuration for the model.
            input_size (int, optional): Input size of the model.
            output_size (int, optional): Output size of the model.

        Returns:
            PredictorModel: An instance of PredictorModel with pre-trained weights loaded.
        """
        # Initialize the model with the provided config
        logger = Logger(config.get('verbose', 1))  # Initialize a logger
        logger.log(f"Loading model from weights at {path_to_weights}...")
    
        # Initialize the model with the provided config
        model = cls(config, input_size, output_size)

        # Load pre-trained weights
        try:
            state_dict = torch.load(path_to_weights, map_location=torch.device('cpu'), weights_only=True)  # Asegura compatibilidad con CPU
            model.load_state_dict(state_dict)
            logger.log("Pre-trained weights loaded successfully.")
        except Exception as e:
            logger.error(f"Failed to load weights: {e}", exception_type=RuntimeError)
        return model

class ExplainerModel(): 
    def __init__(self, config_path_explain, config_path_train):
        self.logger = Logger(verbose=1)
        
        # Load configurations
        self.config_parser = ConfigParser(config_path_explain)
        self.base_config = self.config_parser.get_base_config()
        self.explain_config = self.config_parser.get_explainability_config()
        self.config_train_parser = ConfigParser(config_path_train)
        self.training_config = self.config_train_parser.get_model_training_config()
        
        # Define paths for saving data and results
        self.path_save_results = os.path.join(
            self.base_config['output_dir'], 
            'results',
            f"{self.explain_config['explanation_method']}_{self.explain_config['reference_data']}_{self.explain_config['batch_reduction_method']}_{self.explain_config['gene_collapse_method']}"
        )
        ensure_directory_exists(self.path_save_results)
        
        # Initialize data processing components
        self.data_importer = DataImporter(self.base_config['data_paths'])
        self.data_loader = DatasetLoader(self.data_importer, self.base_config)
        self.scaler = Scaler.load(self.explain_config['scaler_path'])
        
        # Initialize attributes for data and model (to be loaded later)
        self.data = None
        self.model = None
        self.explainer_handler = None
    
    def load_and_process_data(self):
        """Load and preprocess the data."""
        self.logger.log("Loading and processing data...")
        self.data = self.data_loader.load_data()
        self.data['scaled_rbp_expr_df'] = self.scaler.transform(self.data['rbp_expr_df'])
        self.logger.log("Data loaded and processed successfully.")
        return self.data
    
    def load_model(self):
        """Load the trained model."""
        self.logger.log("Loading the trained model...")
        self.model = PredictorModel.load_model(
            path_to_weights=os.path.join(self.explain_config['trained_model_path'], self.explain_config['model_file']),
            config=self.training_config
        )
        self.logger.log("Model loaded successfully.")
        return self.model
    
    def _initialize_explainer_handler(self):
        """Initialize the ExplainerHandler with the loaded model and data."""
        if self.model is None or self.data is None:
            self.logger.error("Model and data must be loaded before initializing ExplainerHandler.")
        
        # Determine explanation method based on configuration
        reference_type = self.explain_config.get('explanation_method')
        if reference_type == "DeepLIFT":
            self.explainer_handler = DeepLiftHandler(self.model, self.data, self.base_config, self.explain_config)
            self.logger.log("DeepLiftHandler initialized successfully.")
        elif reference_type == "Pseudoknocking":  # yet to develop (work in progress)
            pass  # Handle other explanation methods as needed
    
    def perform_explainer(self):
        """Perform the explanation method based on the configured explainer handler and return the results."""
        self._initialize_explainer_handler()
        # Prepare RBP tensors
        scaled_rbp_tensor, gn_tensor, reference_rbp_tensor = self.explainer_handler.prepare_rbp_tensors()
        
        # Compute attribution scores
        list_batch_scores = self.explainer_handler.compute_attribution_scores(scaled_rbp_tensor, reference_rbp_tensor, gn_tensor)
        
        # Reduce batch dimension (RBP x T)
        df_scores_TxRBP = self.explainer_handler.reduce_batch_dimension(list_batch_scores)
        
        # Filter scores for low-expressed genes
        df_scores_TxRBP = self.explainer_handler.filter_scores_for_low_expressed_genes(
            deeplift_scores=df_scores_TxRBP,
            gene_expr_df=self.data['gene_expr_df'],
            threshold=1
        )
        
        # Collapse scores to genes (RBP x G)
        result_table, df_scores_GxRBP = self.explainer_handler.collapse_transcript_scores_to_genes(df_scores_TxRBP)
        
        # Save results as CSV files
        self.save_results(df_scores_TxRBP, df_scores_GxRBP, result_table)
        return {
            'df_scores_TxRBP': df_scores_TxRBP,
            'df_scores_GxRBP': df_scores_GxRBP,
            'result_table': result_table
        }
    
    def save_results(self, df_scores_TxRBP, df_scores_GxRBP, result_table):
        """Save the results as CSV files."""
        try:
            df_scores_TxRBP.to_csv(os.path.join(self.path_save_results, 'df_scores_TxRBP.csv'), index=False)
            df_scores_GxRBP.to_csv(os.path.join(self.path_save_results, 'df_scores_GxRBP.csv'), index=False)
            result_table.to_csv(os.path.join(self.path_save_results, 'result_table.csv'), index=False)
            self.logger.log("Results saved successfully.")
        except Exception as e:
            self.logger.error(f"Error saving results: {e}")









#################################
class PredictorModel(nn.Module): # version old (borrar)
    """Create a neural network with multiple hidden layers allowing for flexible configuration of each layer's node count. 
    This network architecture is defined by specifying the number of hidden layers and the maximum number of nodes in the 
    first hidden layer. Each subsequent hidden layer's node count decreases by a specified division factor. 

    Additionally, the user can opt for all hidden layers to have the same number of nodes by setting `same_num_nodes` to True. 
    If False, each hidden layer will have half the number of nodes of the previous layer by default, with the option to adjust 
    this division factor. 

    The activation function for the hidden layers can be customized, with ReLU being the default.

    This neural network, DeepRBP, is designed for predicting transcript abundance given RNA-binding protein (RBP) and gene expression data.

    Args:
        config (Any): Dictionary-like object with the model configuration.
        input_size (int, optional): Size of the input features. Overrides `config['input_size']` if provided.
        output_size (int, optional): Size of the output features. Overrides `config['output_size']` if provided.

    Attributes:
        input_size (int): Number of input features.
        output_size (int): Number of output features.
        num_hidden_layers (int): Number of hidden layers.
        max_nodes (int): Maximum number of nodes in the first hidden layer.
        uniform_nodes (bool): Whether all hidden layers have the same number of nodes.
        node_shrink_factor (float): Factor by which the number of nodes decreases per layer (if `uniform_nodes` is False).
        activation_func (str): Activation function to use in hidden layers (e.g., 'relu', 'tanh').
        optimizer_name (str): Name of the optimizer to use (e.g., 'adam', 'sgd90').
        learning_rate (float): Learning rate for the optimizer.
    """
  
    def __init__(self, config: Any, input_size: int = None, output_size: int = None, verbose: int = 1):
        super(PredictorModel, self).__init__()
        self.logger = Logger(verbose=verbose)

        # Model configuration
        self.input_size = input_size if input_size is not None else config['input_size']
        self.output_size = output_size if output_size is not None else config['output_size']
        self.num_hidden_layers = config.get('num_hidden_layers', 1)
        self.max_nodes = config['max_node']
        self.uniform_nodes = config['uniform_nodes']
        self.node_shrink_factor = config['node_shrink_factor']
        self.activation_func = config['activation_func']
        self.learning_rate = config['learning_rate']
        self.optimizer_name = config['optimizer_name']

        # Configure layers and optimizer
        self.logger.log("Initializing layers and optimizer...")
        self._configure_layers()
        self.optimizer = self._configure_optimizer(self.optimizer_name, self.learning_rate)

    def _configure_layers(self):
        """Configures the hidden and output layers based on model configuration."""
        if self.num_hidden_layers > 0:
            self.logger.log(f"Using a model with {self.num_hidden_layers} hidden layers")

            # Initialize the number of nodes for the first hidden layer
            node_count = self.max_nodes
            
            for i in range(self.num_hidden_layers):
                # Determine the input size for the current layer
                input_size = self.input_size if i == 0 else node_count

                # Configure hidden layers
                if self.uniform_nodes:
                    if i == self.num_hidden_layers - 1: # Last layer: apply shrink factor
                        node_count = round(node_count / self.node_shrink_factor)
                    layer = nn.Linear(input_size, node_count)
                else:
                    # Reduce the number of nodes in each layer if `uniform_nodes` is False
                    layer = nn.Linear(input_size, node_count)
                    node_count = round(node_count / self.node_shrink_factor)
                
                # Add the layer, batch normalization, and activation
                self.add_module(f'hidden_linear_{i}', layer)
                bn_layer = nn.BatchNorm1d(layer.out_features)
                self.add_module(f'batch_norm_{i}', bn_layer)

                # Add activation function (e.g., ReLU, Tanh, etc.)
                activation_name = self._get_activation_module(self.activation_func)
                self.add_module(f'activation_{i}', activation_name)

            # The input size for the final layer is the number of nodes in the last hidden layer
            final_hidden_layer_input_size = node_count
        else:
            self.logger.log("Using a model with zero hidden layers")
            final_hidden_layer_input_size = self.input_size  # No hidden layers, use input size directly

        # Configure the output layer
        self.output_layer = nn.Linear(final_hidden_layer_input_size, self.output_size)
        self.add_module('output_layer', self.output_layer)

    def _configure_optimizer(self, optimizer_name, learning_rate):
        """Configures the optimizer based on the provided name and learning rate.
        Args:
            optimizer_name (str): Name of the optimizer to use ('sgd90', 'sgd70', etc.).
            learning_rate (float): Learning rate for the optimizer.
        Returns:
            torch.optim.Optimizer: Configured optimizer instance.
        """
        self.logger.log(f"Configuring optimizer: {optimizer_name}")
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
        self.logger.log(f"Configuring activation function: {activation_name}")
        if activation_name == "relu":
            return nn.ReLU()
        elif activation_name == "tanh":
            return nn.Tanh()
        elif activation_name == "sigmoid":
            return nn.Sigmoid()
        else:
            raise ValueError("Invalid activation_layer. Supported options are 'relu', 'tanh', and 'sigmoid'.")
        
    def forward(self, xb, gb):
        """Defines the forward pass of the model.
        Args:
            xb (torch.Tensor): Input features (e.g., RBP and gene expression).
            gb (torch.Tensor): Additional gene input for the sigmoid scaling.
        Returns:
            torch.Tensor: Predicted transcript abundance.
        """
        x=xb
        for _, module in self.named_children():
            x = module(x)
        out = torch.log2((torch.sigmoid(x) * gb) + 1)
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
        torch.save(self.state_dict(), os.path.join(output_dir, model_name))
        self.logger.log(f"Model saved successfully in {output_dir}")

    @classmethod
    def load_model(cls, path_to_weights, config: Any, input_size: int = None, output_size: int = None):
        """
        Class method to initialize the model and load pre-trained weights.

        Args:
            path_to_weights (str): Path to the file containing the pre-trained weights.
            config (Any): Configuration for the model.
            input_size (int, optional): Input size of the model.
            output_size (int, optional): Output size of the model.

        Returns:
            PredictorModel: An instance of PredictorModel with pre-trained weights loaded.
        """
        # Initialize the model with the provided config
        logger = Logger(config.get('verbose', 1))  # Initialize a logger
        logger.log(f"Loading model from weights at {path_to_weights}...")
    
        # Initialize the model with the provided config
        model = cls(config, input_size, output_size)

        # Load pre-trained weights
        try:
            state_dict = torch.load(path_to_weights, map_location=torch.device('cpu'), weights_only=True)  # Asegura compatibilidad con CPU
            model.load_state_dict(state_dict)
            logger.log("Pre-trained weights loaded successfully.")
        except Exception as e:
            logger.error(f"Failed to load weights: {e}", exception_type=RuntimeError)
        return model

#class ExplainerModel():