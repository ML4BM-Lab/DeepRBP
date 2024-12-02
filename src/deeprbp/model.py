#class DeepRBPModel
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, Any

class PredictorModel(nn.Module):
    """Create a neural network with multiple hidden layers allowing for flexible configuration of each layer's node count. 
    This network architecture is defined by specifying the number of hidden layers and the maximum number of nodes in the 
    first hidden layer. Each subsequent hidden layer's node count decreases by a specified division factor. 

    Additionally, the user can opt for all hidden layers to have the same number of nodes by setting `same_num_nodes` to True. 
    If False, each hidden layer will have half the number of nodes of the previous layer by default, with the option to adjust 
    this division factor. 

    The activation function for the hidden layers can be customized, with ReLU being the default.

    This neural network, DeepRBP, is designed for predicting transcript abundance given RNA-binding protein (RBP) and gene expression data.

    Args:
      
    Attributes:
  
    """
    def __init__(self, config: Any, input_size: int = None, output_size: int = None):
        super(PredictorModel, self).__init__()
        # Model configuration
        config = config['model']
        self.input_size = input_size if input_size is not None else model_config['input_size']
        self.output_size = output_size if output_size is not None else model_config['output_size']
        self.num_hidden_layers = config.get('num_hidden_layers', 1)
        self.max_nodes = config['max_node']
        self.uniform_nodes = config['uniform_nodes']
        self.node_shrink_factor = config['node_shrink_factor']
        self.activation_func = config['activation_func']
        self.learning_rate = config['learning_rate']
        self.optimizer_name = config['optimizer_name']
        self._configure_layers()
        self.optimizer = self._configure_optimizer(self.optimizer_name, self.learning_rate)

    def _configure_layers(self):
        """Configures the hidden and output layers based on model configuration."""
        if self.num_hidden_layers > 0:
            print(f'Using a model with {self.num_hidden_layers} hidden layers')
            # Inicializa el número de nodos
            node_count = self.max_nodes
            for i in range(self.num_hidden_layers):
                # Determina el tamaño de entrada para cada capa
                input_size = self.input_size if i == 0 else node_count
                # Configura las capas ocultas
                if self.uniform_nodes:
                    if i == self.num_hidden_layers - 1:  # Última capa: reducir nodos con el shrink factor
                        node_count = round(node_count / self.node_shrink_factor)
                    layer = nn.Linear(input_size, node_count)
                else:
                    # Si uniform_nodes es False, reducimos el número de nodos en cada capa
                    layer = nn.Linear(input_size, node_count)
                    node_count = round(node_count / self.node_shrink_factor)
                # Añade la capa, normalización y activación
                self.add_module(f'hidden_linear_{i}', layer)
                bn_layer = nn.BatchNorm1d(layer.out_features)
                self.add_module(f'batch_norm_{i}', bn_layer)
                activation_name = self._get_activation_module(self.activation_func)
                self.add_module(f'activation_{i}', activation_name)
            # El tamaño de entrada para la capa final es el número de nodos en la última capa oculta
            final_hidden_layer_input_size = node_count
        else:
            print('Using a model with zero hidden layers')
            final_hidden_layer_input_size = self.input_size  # Si no hay capas ocultas, el tamaño final es igual a input_size
        # Configura la capa de salida
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
            raise ValueError(f"Unsupported optimizer: {optimizer_name}") 
        
    def _get_activation_module(self, activation_name):
        """Returns the activation layer based on the given name."""
        if activation_name == "relu":
            return nn.ReLU()
        elif activation_name == "tanh":
            return nn.Tanh()
        elif activation_name == "sigmoid":
            return nn.Sigmoid()
        else:
            raise ValueError("Invalid activation_layer. Supported options are 'relu', 'tanh', and 'sigmoid'.")
        
    def forward(self, xb, gb):
        x=xb
        for _, module in self.named_children():
            x = module(x)
        out = torch.log2((torch.sigmoid(x) * gb) + 1)
        return out
    
    def train_step(self, inputs, targets, gen_expr):
        out = self(inputs, gen_expr) # Forward pass
        loss = F.mse_loss(out, targets) # Calculate loss
        loss.backward() # Backpropagation
        self.optimizer.step() # Update weights
        self.optimizer.zero_grad() # Reset gradients
        return loss.detach()

    def validate_step(self, inputs, targets, gen_expr):
        with torch.no_grad():
            out = self(inputs, gen_expr) # Forward pass
            val_loss = F.mse_loss(out, targets)  # Calculate loss
        return val_loss.detach()

#class ExplainerModel(nn.Module):