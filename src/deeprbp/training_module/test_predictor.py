# src/deeprbp/training_module/test_predictor.py

import argparse
from .pipeline import DeepRBPredictorPipeline   

def parse_args():   
    parser = argparse.ArgumentParser(description='Evaluate the DeepRBP predictor using a test data.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file of the test data.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results.')
    parser.add_argument('--trained_files_dir', type=str, required=True, help='Directory to the trained model and scaler.')
    return parser.parse_args()
    
def main():
    args = parse_args()
    # Create an instance of the DeepRBPredictorPipeline with the provided configuration and output directory
    pipeline = DeepRBPredictorPipeline(args.config_path, args.output_dir)
    set_name='test'
    # Import the test data using the pipeline
    test_data = pipeline.import_data()
    # Load the scaler used during training from the specified directory
    pipeline.load_scaler(args.trained_files_dir)
    # Scale the test data using the loaded scaler
    scaled_test_data = pipeline.scale_data(test_data)
    # Get the DataLoader for the scaled test data
    test_loader = pipeline.get_loaders(scaled_test_data)[0]
    # Set up the model trainer with the test DataLoader and the path to the model weights
    pipeline.setup_model_trainer(test_loader, path_to_weights=f'{args.trained_files_dir}/best_model.pt')
    # Evaluate the model on the test dataset
    pipeline.evaluate_model(test_loader, set_name)
    # Evaluate the model's performance per category in the test dataset
    pipeline.eval_model_per_category(scaled_test_data, set_name)  

if __name__ == "__main__":
    main()

# CHANGE THIS WHOLE CODE
# mira q ahora va a ser asi: 
 # Si se necesita escalar los datos, realiza la operación aquí
        # scaled_test_data = self.scale_data(test_data)  # Método que deberías implementar para escalar

        # # Crear el TensorDataset para los datos de test
        # test_dataset = self.create_tensor_dataset(scaled_test_data)

        # # Crear el DataLoader para los datos de test
        # test_loader = self.create_data_loader(test_dataset, batch_size=self.val_batch_size, shuffle=False, drop_last=False)
        
       # Test on final datasets

# Evaluate general
# (inserta codigo aquí bro)
#trainer.test(model, sparseGO_data)