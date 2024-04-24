## Model Directory (`./model`)

The `./model` directory contains the necessary code and utilities to create your own model using processed data.

### Directory Structure

- **DLModelClass**
  - Contains the `DeepRBP` class, the predictive model with specialized Deep Learning functions.
  
- **utils**
  - A utility folder containing essential Python files:
    - `Config.py`: Manages and defines the configuration for the deep learning model.
    - `Generate_Datasets.py`: Processes raw data from datasets (TCGA, GTeX, and real kds) to create inputs for the model.
    - `Utils.py`: Provides helper functions and utilities for common or complex tasks across the project.
    - `Plots.py`: Contains functions for data visualization, creating charts to interpret and present results.

### Model Output

The trained model, along with its configuration and the scaler used, will be saved in `./model/output`.

### Creating the Model

You can create the model by following the provided Jupyter notebook tutorial `.ipynb` or using the shell script `run_DeepRBP_predictor.sh`. This shell script invokes `main_predictor.py` to train the DeepRBP model with the specified configurations.

### Example Model

We've included an example trained model in the `./model/output/example` directory. This model was trained with samples from TCGA of all tumor types, utilizing two hidden layers with 1024 and 128 nodes, a batch size of 128, and 1,000 epochs.





