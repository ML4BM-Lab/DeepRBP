import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from tqdm import tqdm
from collections import OrderedDict
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import seaborn as sns
import matplotlib.pyplot as plt
import torch.nn.functional as F
from datetime import datetime
import re
import os
import optuna
import pandas as pd

class DeepRBP(nn.Module):
    """Create a neural network with multiple hidden layers allowing for flexible configuration of each layer's node count. 
    This network architecture is defined by specifying the number of hidden layers and the maximum number of nodes in the 
    first hidden layer. Each subsequent hidden layer's node count decreases by a specified division factor. 

    Additionally, the user can opt for all hidden layers to have the same number of nodes by setting `same_num_nodes` to True. 
    If False, each hidden layer will have half the number of nodes of the previous layer by default, with the option to adjust 
    this division factor. 

    The activation function for the hidden layers can be customized, with ReLU being the default.

    This neural network, DeepRBP, is designed for predicting transcript abundance given RNA-binding protein (RBP) and gene expression data.

    Args:
        config: Configuration object containing neural network settings.
        device: Device to run the computation (e.g., 'cpu' or 'cuda').
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Attributes:
        n_inputs (int): Number of input features. Defaults to 1282.
        n_outputs (int): Number of output features. Defaults to 11462.
        device (str): Device to run the computation (e.g., 'cpu' or 'cuda').
        num_hidden_layers (int): Number of hidden layers.
        max_node (int): Maximum number of nodes in the first hidden layer.
        same_num_nodes (bool): Indicates whether all hidden layers have the same number of nodes.
        node_division_factor (int): Factor by which the number of nodes decreases from one hidden layer to the next.
        activation_layer (str): Name of the activation function used in the hidden layers.

    Methods:
        forward: Performs a forward pass through the network.
        training_step: Executes a training step.
        validation_step: Executes a validation step.
        training_epoch_end: Computes training loss at the end of an epoch.
        validation_epoch_end: Computes validation loss at the end of an epoch.
        epoch_end: Prints epoch-level information.
        evaluate: Evaluates the model on the validation set. 
                  It computes the validation loss and returns it.
        do_training: Performs training on the training set. 
                     It executes the training step for each batch and computes the training loss.
        fit: Orchestrates the training process. 
             It iterates through epochs, performs training and validation, and saves the model if specified.
        predict: Generates predictions on a given data loader.
    """
    def __init__(self, config, device=None, *args, **kwargs):

        super().__init__()
        self.n_inputs = kwargs.get('n_inputs', 1282)
        self.n_outputs = kwargs.get('n_outputs', 11459)
        self.device = device if device else torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.num_hidden_layers = config.num_hidden_layers
        self.max_node = config.max_node
        self.same_num_nodes = config.same_num_nodes
        self.node_divition_factor = config.node_divition_factor
        self.activation_layer = config.activation_layer
        
        if self.num_hidden_layers > 0:
            print(f'You are using a model with {self.num_hidden_layers} hidden layers')
            # Intermediate layers
            for i in range(self.num_hidden_layers):
                if i==0:
                    layer = nn.Linear(self.n_inputs, self.max_node)
                else:
                    if self.same_num_nodes and ((i+1) < self.num_hidden_layers):                                                    
                        layer = nn.Linear(layer.out_features, self.max_node)
                    else:
                        layer = nn.Linear(layer.out_features, round(layer.out_features/self.node_divition_factor))    
                self.add_module(f'linear{i}', layer)
                bn_layer = nn.BatchNorm1d(layer.out_features)
                self.add_module(f'bn{i}', bn_layer)
                activation_name = self.get_activation_module(self.activation_layer)
                self.add_module(f'a_function{i}', activation_name)
            final_hidden_layer_input_size = layer.out_features
        else:
            # When num_hidden_layers is 0, final hidden layer input size is the input size
            print('You are using a model with zero hidden layers')
            final_hidden_layer_input_size = self.n_inputs
        # Final Layer
        self.final_layer = nn.Linear(final_hidden_layer_input_size, self.n_outputs) 
        self.add_module('final_layer', self.final_layer)

    def get_activation_module(self, activation_name):
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

    def training_step(self, batch, optimizer):
        inputs, targets, gen_expr = batch
        inputs, targets, gen_expr = inputs.to(self.device), targets.to(self.device), gen_expr.to(self.device)
        out = self(inputs, gen_expr)   # Generate predictions
        loss = F.mse_loss(out, targets)    # Calculate loss
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        return {'loss': loss.detach()}

    def validation_step(self, batch):
        inputs, targets, gb = batch
        inputs, targets, gb = inputs.to(self.device), targets.to(self.device), gb.to(self.device)
        with torch.no_grad():
            out = self(inputs, gb)             # Generate predictions
            loss_val = F.mse_loss(out, targets)    # Calculate loss
        return {'val_loss': loss_val.detach()}

    def training_epoch_end(self, outputs):
        batch_losses = [x['loss'] for x in outputs]
        epoch_loss = torch.stack(batch_losses).mean()   # Combine losses
        return {'loss': epoch_loss.item()}

    def validation_epoch_end(self, outputs):
        batch_losses = [x['val_loss'] for x in outputs]
        epoch_loss = torch.stack(batch_losses).mean()   # Combine losses
        return {'val_loss': epoch_loss.item()}

    #def epoch_end(self, epoch, loss, result):
    #    print("Epoch [{}], loss: {:.4f}, val_loss: {:.4f}".format(epoch, loss['loss'], result['val_loss']))  

    def evaluate(self, val_loader):
        outputs = [self.validation_step(batch) for batch in val_loader]
        return self.validation_epoch_end(outputs)

    def do_training(self, train_loader, optimizer):
        outputs = [self.training_step(batch, optimizer) for batch in train_loader]
        return self.training_epoch_end(outputs)

    def fit(self, epochs, train_loader, val_loader, optimizer, save_model=True, path=None, model_name=None):
        self.train() # set model to train mode
        history_train = []
        history_val = []

        for epoch in tqdm(range(epochs), desc="Training", unit="epoch"):
            # Training Phase
            train_epoch_end = self.do_training(train_loader, optimizer)
            history_train.append(train_epoch_end)
            # Validation phase
            val_epoch_end = self.evaluate(val_loader)
            history_val.append(val_epoch_end)
            #self.epoch_end(epoch, train_epoch_end, val_epoch_end)
            
            if epoch % 10 == 0:  # Show tqdm progress every 10 epochs
                tqdm.write(f'Epoch {epoch}/{epochs} - Training Loss: {train_epoch_end}, Validation Loss: {val_epoch_end}')

        if save_model:
            torch.save(self.state_dict(), os.path.join(path, model_name)) 
            print(f'Model saved succesfully in {path}')   
        return history_train, history_val

    def predict(self, data_loader):
        y_true = []
        y_pred = []
        self.eval()  # Set the model to evaluation mode
        with torch.no_grad():
            #self = self.to(self.device)
            for batch in data_loader:
                inputs, targets, gb = batch
                inputs, targets, gb = inputs.to(self.device), targets.to(self.device), gb.to(self.device)
                out = self(inputs, gb)  # Generate predictions
                y_true.append(targets.cpu().numpy())
                y_pred.append(out.detach().cpu().numpy())
        y_true = np.concatenate(y_true).flatten()
        y_pred = np.concatenate(y_pred).flatten()
        return y_pred, y_true