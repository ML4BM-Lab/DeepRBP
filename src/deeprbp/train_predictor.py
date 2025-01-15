
import torch
import numpy as np
from tqdm import tqdm

class TrainPredictor:
    def __init__(self, model, config, input_features=('rbp_expr', 'gene_expr'), output_features=('trans_expr',)):
        self.model = model
        self.config = config
        self.device = torch.device('cuda' if config['cuda'] and torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
        self.input_features = input_features  # Features que serán usadas como entrada
        self.output_features = output_features  # Features que serán usadas como salida

    def prepare_batch(self, batch):
        """
        Auxiliary method to extract the inputs and target from the batch.
        """
        inputs = [batch[feature].to(self.device).float() for feature in self.input_features]
        targets = [batch[feature].to(self.device).float() for feature in self.output_features]
        rbp_expr, gen_expr = inputs
        targets = torch.stack(targets).squeeze(0)
        return rbp_expr, targets, gen_expr
    
    def train_one_epoch(self, train_loader):
        """
        Train the model for one epoch using the provided DataLoader.
        """
        self.model.train()  
        losses = []
        for batch in train_loader:
            rbp_expr, targets, gen_expr = self.prepare_batch(batch)
            loss = self.model.train_step(rbp_expr, targets, gen_expr)  # Usamos las variables definidas directamente
            losses.append(loss)
        return torch.stack(losses).mean().item()  # Promedio de pérdidas
    
    def validate_one_epoch(self, val_loader):
        self.model.eval()  
        val_losses = []
        with torch.no_grad():
            for batch in val_loader:
                rbp_expr, targets, gen_expr = self.prepare_batch(batch)
                val_loss = self.model.validate_step(rbp_expr, targets, gen_expr)  # Validation step
                val_losses.append(val_loss)
        return torch.stack(val_losses).mean().item()  # Promedio de pérdidas
    
    def fit(self, train_loader, val_loader, epochs):
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
                tqdm.write(f'Epoch {epoch}/{epochs} - Training Loss: {train_loss:.4f}, 'f'Validation Loss: {val_loss:.4f}')
        return train_history, val_history
    
    def generate_predictions(self, data_loader):
        """
        Genera las predicciones del modelo sobre un DataLoader de datos.
        """
        true_values = []
        predictions = []
        self.model.eval()  # Ponemos el modelo en modo evaluación
        with torch.no_grad():
            for batch in data_loader:
                rbp_expr, targets, gen_expr = self.prepare_batch(batch)
                out = self.model(rbp_expr, gen_expr) # Generar las predicciones
                true_values.append(targets.cpu().numpy())
                predictions.append(out.detach().cpu().numpy())
        true_values = np.concatenate(true_values).flatten()
        concatenated_predictions = np.concatenate(predictions)
        flattened_predictions = concatenated_predictions.flatten()
        return flattened_predictions, true_values, concatenated_predictions

    # def save_model(self, output_dir, model_name):
    #     """Saves the model to the specified directory."""
    #     # Save the model
    #     os.makedirs(output_dir, exist_ok=True)
    #     torch.save(self.model.state_dict(), os.path.join(output_dir, model_name))
    #     print(f'Model saved successfully in {output_dir}')
        