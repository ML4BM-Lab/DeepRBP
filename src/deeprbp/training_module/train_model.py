# src/deeprbp/training_module/train_model.py

import os
import torch
import numpy as np
from tqdm import tqdm
import optuna

from .model import PredictorModel

class TrainPredictor:
    """
    A class to train and evaluate a machine learning model for predicting transcript abundance in log2(TPM+1)

    This class encapsulates the training and validation processes for the `PredictorModel`, which
    is designed to predict transcript levels based on RNA-binding protein (RBP) and gene expression data.
    It handles the entire training loop, including epoch management, loss calculation, and model evaluation.

    Args:
        model (PredictorModel): An instance of the PredictorModel to be trained and evaluated.
        config (dict): Configuration object containing settings for training, such as device (CPU/GPU) and logging options.
        input_features (tuple, optional): Tuple of input feature names used for model training. These features represent
                                           the RNA-binding protein expression and gene expression data extracted from
                                           the DataLoader batches (default: ('scaled_rbp_expr_log2p_tpm', 
                                           'gn_expr_each_iso_tpm')).
        output_features (tuple, optional): Tuple of output feature names that the model will predict. These features
                                            represent the target transcript abundance values that the model aims to
                                            predict, derived from the DataLoader batches (default: ('trans_expr_log2p_tpm',)).
    """
    def __init__(self, model, config, input_features=None, output_features=None):
        """
        Initializes the TrainPredictor class with a model, configuration, and feature specifications.
        """
        self.config = config
        self.device = torch.device('cuda' if self.config.get('cuda') and torch.cuda.is_available() else 'cpu')
        print(f"🖥️  Using device: {self.device} for training.")
        self.model = model
        self.is_trained = False 
        # Set default input and output features if not provided
        self.input_features = input_features if input_features is not None else ('scaled_rbp_expr_log2p_tpm', 'gn_expr_each_iso_tpm')  # Feature keys to be used as inputs
        self.output_features = output_features if output_features is not None else ('trans_expr_log2p_tpm',) # Features keys to be predicted
        # Sent the model to the device
        self.model.to(self.device)
        print(f"✅ Model has been moved to: {next(self.model.parameters()).device}")
        # Set random seed for reproducibility
        seed = self.config.get('seed', 42)
        torch.manual_seed(seed)
        np.random.seed(seed)
        if self.device.type == 'cuda':
            torch.cuda.manual_seed(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False  
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
        #print(f"📦 Batch data moved to device: {inputs[0].device}, targets: {targets[0].device}")
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
            assert rbp_expr.device == self.device, "rbp_expr is not on the correct device."
            assert targets.device == self.device, "targets is not on the correct device."
            assert gen_expr.device == self.device, "gen_expr is not on the correct device."
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
    def fit(self, train_loader, val_loader, epochs, path_save_results=None, optuna_trial=None):
        """
        Trains the model for a specified number of epochs.

        This method orchestrates the training and validation processes, logging progress and returning
        the training and validation loss history. It also manages the saving and loading of the best model
        based on validation loss if configured to do so.

        Parameters:
            train_loader (DataLoader): DataLoader for the training dataset.
            val_loader (DataLoader): DataLoader for the validation dataset.
            epochs (int): The number of epochs to train the model.
            path_save_results (str, optional): Directory path where the best model will be saved. 
                                                If None, the best model will not be saved.
            optuna_trial (optuna.Trial, optional): The Optuna trial object for reporting metrics and pruning.

        Returns:
            tuple: A tuple containing:
                - list: History of training losses for each epoch.
                - list: History of validation losses for each epoch.
                - PredictorModel: The trained model instance.

        Raises:
            ValueError: If the path to save the model is invalid or if the training process encounters issues.

        Notes:
            - The best model is saved only if the configuration parameter `save_best_model` is set to True.
            - If `save_best_model` is True and `path_save_results` is provided, the best model weights are saved whenever the validation loss improves.
            - After training, if `save_best_model` is True and a valid path is provided, the best model is loaded back into the instance for further use.
            - Ensure that the model is trained before calling `generate_predictions`, as predictions made without training may not be reliable.
        """
        # Check if saving the best model is configured
        save_best_model = self.config.get('save_best_model', False)
        # Raise an error if attempting to save the best model without a specified path
        if save_best_model and path_save_results is None:
            raise ValueError("When 'save_best_model' is True, 'path_save_results' must be specified to save the model.")
        train_history = []
        val_history = []
        print_every = self.config.get('print_every', 1)
        best_val_mse = float('inf')  # Initialize with infinity
        print(f"\nStarting training for {epochs} epochs... 🚀")
        for epoch in tqdm(range(epochs), desc="Training", unit="epoch"):
            # Training Phase
            train_loss = self.train_one_epoch(train_loader)
            train_history.append(train_loss)
            # Validation Phase
            val_loss = self.validate_one_epoch(val_loader)
            val_history.append(val_loss)
            # Report the validation loss to Optuna if optuna_trial is provided
            if optuna_trial is not None:
                optuna_trial.report(val_loss, epoch)  # Report the current validation loss
                if optuna_trial.should_prune():  # Check if the trial should be pruned
                     print(f"🚫 Pruning the current trial due to poor performance at epoch {epoch + 1}.")   
                     raise optuna.TrialPruned()
            # Display progress
            if epoch % print_every == 0:
                tqdm.write(f'Epoch {epoch}/{epochs} - Training Loss: {train_loss:.4f}, '
                           f'Validation Loss: {val_loss:.4f}')
            # Save the best model if configured to do so
            if save_best_model and val_loss < best_val_mse:
                best_val_mse = val_loss
                if path_save_results:
                    self.model.save_model(path_save_results, 'best_model.pt')
                    print(f'💾 Best model saved at epoch {epoch + 1} with validation loss: {best_val_mse:.4f}')  # Print when saving the model
        # Load the best model after training if configured to do so
        if save_best_model and path_save_results:
            print("🔄 Loading the best model...")
            self.model = PredictorModel.load_model(os.path.join(path_save_results, 'best_model.pt'), self.config, 
                                                   input_size=self.model.input_size, output_size=self.model.output_size)
        self.is_trained = True
        print("\n🏁 Training completed!")
        return train_history, val_history # not commited yet:, self.model
    def generate_predictions(self, data_loader):
        """
        Generates log2(tpm+1) predictions from the model using the provided DataLoader.

        This method evaluates the model on the provided dataset, collecting both the true values
        and the model predictions.

        Parameters:
            data_loader (DataLoader): DataLoader for the dataset on which to generate predictions.

        Returns:
            tuple: A tuple containing:
                - np.ndarray: batched predictions
                - np.ndarray: batched labels

        Raises:
            Warning: If the model has not been trained yet, indicating that predictions may not be reliable.
        """
        if not self.is_trained:
            print("⚠️ Warning: The model has not been trained yet! Predictions may not be reliable.")
        true_values = []
        predictions = []
        self.model.eval()  # Set the model to evaluation mode
        with torch.no_grad():
            for batch in data_loader:
                rbp_expr, targets, gen_expr = self.prepare_batch(batch)
                assert rbp_expr.device == self.device, "rbp_expr is not on the correct device."
                assert gen_expr.device == self.device, "gen_expr is not on the correct device."
                out = self.model(rbp_expr, gen_expr) # Generate predictions
                true_values.append(targets.cpu().numpy())
                predictions.append(out.detach().cpu().numpy())
        true_values = np.concatenate(true_values)
        predictions = np.concatenate(predictions)
        return predictions, true_values
    
        

