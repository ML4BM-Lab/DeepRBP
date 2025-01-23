# src/deeprbp/training_module/train_model.py

import torch
import numpy as np
from tqdm import tqdm

class TrainPredictor:
    """
    A class to train and evaluate a machine learning model for predicting transcript abundance.

    This class encapsulates the training and validation processes for the `PredictorModel`, which
    is designed to predict transcript levels based on RNA-binding protein (RBP) and gene expression data.
    It handles the entire training loop, including epoch management, loss calculation, and model evaluation.

    Args:
        model (PredictorModel): An instance of the PredictorModel to be trained and evaluated.
        config (dict): Configuration dictionary containing settings for training, such as device (CPU/GPU) and logging options.
        input_features (tuple): Tuple of input feature names used for model training (default: ('scaled_rbp_expr_log2p_tpm', 'gn_expr_each_iso_tpm')).
        output_features (tuple): Tuple of output feature names that the model will predict (default: ('trans_expr_log2p_tpm_df',)).
    """
    def __init__(self, model, config, input_features=('scaled_rbp_expr_log2p_tpm', 'gn_expr_each_iso_tpm'), output_features=('trans_expr_log2p_tpm',)):
        """
        Initializes the TrainPredictor class with a model, configuration, and feature specifications.
        """
        self.model = model
        self.config = config
        self.device = torch.device('cuda' if config['cuda'] and torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
        self.input_features = input_features  # Features to be used as inputs
        self.output_features = output_features  # Features to be predicted

    def prepare_batch(self, batch):
        """
        Extracts inputs and targets from the given batch.

        This auxiliary method converts input and target features from the batch into tensors,
        moving them to the specified device (CPU or GPU) and ensuring they are of type float.

        Parameters:
            batch (dict): A dictionary containing the batch of data.

        Returns:
            tuple: A tuple containing:
                - rbp_expr (Tensor): Input tensor for RNA-binding protein expression.
                - targets (Tensor): Target tensor for gene expression predictions.
                - gen_expr (Tensor): Input tensor for gene expression data.
        """
        inputs = [batch[feature].to(self.device).float() for feature in self.input_features]
        targets = [batch[feature].to(self.device).float() for feature in self.output_features]
        rbp_expr, gen_expr = inputs
        targets = torch.stack(targets).squeeze(0)
        return rbp_expr, targets, gen_expr
    
    def train_one_epoch(self, train_loader):
        """
        Trains the model for one epoch using the provided DataLoader.

        This method iterates through the training data, performing a training step for each batch
        and collecting the losses.

        Parameters:
            train_loader (DataLoader): DataLoader containing the training dataset.

        Returns:
            float: The average loss for the training epoch.
        """
        self.model.train() # Set the model to training mode
        losses = []
        for batch in train_loader:
            rbp_expr, targets, gen_expr = self.prepare_batch(batch)
            loss = self.model.train_step(rbp_expr, targets, gen_expr)  # Training step
            losses.append(loss)
        return torch.stack(losses).mean().item()  # Average loss over the epoch
    
    def validate_one_epoch(self, val_loader):
        """
        Validates the model for one epoch using the provided DataLoader.

        This method evaluates the model on the validation dataset, collecting validation losses
        without updating model parameters.

        Parameters:
            val_loader (DataLoader): DataLoader containing the validation dataset.

        Returns:
            float: The average validation loss for the epoch.
        """
        self.model.eval() # Set the model to evaluation mode
        val_losses = []
        with torch.no_grad():
            for batch in val_loader:
                rbp_expr, targets, gen_expr = self.prepare_batch(batch)
                val_loss = self.model.validate_step(rbp_expr, targets, gen_expr)  # Validation step
                val_losses.append(val_loss)
        return torch.stack(val_losses).mean().item() # Average validation loss
    
    def fit(self, train_loader, val_loader, epochs):
        """
        Trains the model for a specified number of epochs.

        This method orchestrates the training and validation processes, logging progress and returning
        the training and validation loss history.

        Parameters:
            train_loader (DataLoader): DataLoader for the training dataset.
            val_loader (DataLoader): DataLoader for the validation dataset.
            epochs (int): The number of epochs to train the model.

        Returns:
            tuple: A tuple containing:
                - list: History of training losses for each epoch.
                - list: History of validation losses for each epoch.
        """
        train_history = []
        val_history = []

        for epoch in tqdm(range(epochs), desc="Training", unit="epoch"):
            # Training Phase
            train_loss = self.train_one_epoch(train_loader)
            train_history.append(train_loss)
            # Validation Phase
            val_loss = self.validate_one_epoch(val_loader)
            val_history.append(val_loss)
            # Display progress
            if epoch % self.config['print_every'] == 0:
                tqdm.write(f'Epoch {epoch}/{epochs} - Training Loss: {train_loss:.4f} 📉, 'f'Validation Loss: {val_loss:.4f} 📉')
        return train_history, val_history
    
    def generate_predictions(self, data_loader):
        """
        Generates predictions from the model using the provided DataLoader.

        This method evaluates the model on the provided dataset, collecting both the true values
        and the model predictions.

        Parameters:
            data_loader (DataLoader): DataLoader for the dataset on which to generate predictions.
        """
        true_values = []
        predictions = []
        self.model.eval()  # Set the model to evaluation mode

        with torch.no_grad():
            for batch in data_loader:
                rbp_expr, targets, gen_expr = self.prepare_batch(batch)
                out = self.model(rbp_expr, gen_expr) # Generate predictions
                true_values.append(targets.cpu().numpy())
                predictions.append(out.detach().cpu().numpy())

        true_values = np.concatenate(true_values).flatten()
        concatenated_predictions = np.concatenate(predictions)
        flattened_predictions = concatenated_predictions.flatten()
        return flattened_predictions, true_values, concatenated_predictions
