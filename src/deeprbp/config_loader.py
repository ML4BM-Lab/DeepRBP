# /scratch/jsanchoz/DeepRBP/src/deeprbp/config_loader.py

import yaml
from datetime import datetime
import random
import string
from dataclasses import dataclass, field
from typing import Dict, List, Any

@dataclass
class Config:
    source_name: str
    data_paths: Dict[str, str]
    model: Dict[str, any]
    training: Dict[str, any]
    sample_category: str
    select_samples: List[str]
    output_dir: str = field(default_factory=str)
    seed: int = 0
    plot_results: bool = False

    def __getitem__(self, key: str) -> Any:
        #First, we try to get the value of the class attributes.
        if hasattr(self, key):
            return getattr(self, key)
        #Si no está en los atributos, intentamos acceder a los valores dentro de data_paths o cualquier otro diccionario.
        if key in self.data_paths:
            return self.data_paths[key]
        if key in self.model:
            return self.model[key]
        if key in self.training:
            return self.training[key]
        raise KeyError(f"[Config] Key '{key}' not found in the configuration.")

def generate_unique_id(config: Config, suffix_length: int = 3) -> str:
    """Generates a unique identifier for the configuration based on the source name and selected tumor types."""
    timestamp = datetime.now().strftime("%Y-%m-%d")
    source_train = config.source_name
    tumor_types = '-'.join([t.split('_')[0] for t in config.select_samples])
    # Generate a random alphanumeric suffix to ensure uniqueness
    suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=suffix_length))
    return f"{source_train}_{tumor_types}_{timestamp}_{suffix}"

def load_config(yaml_path: str) -> Config:
    """Loads the configuration from a YAML file and returns a Config object."""
    print('[load_config] Loading configuration... Let\'s make some magic! 💫')

    try:
        with open(yaml_path, 'r') as file:
            config_dict = yaml.safe_load(file)
        
        # Ensure all required configuration keys are present
        required_keys = ['source_name', 'data_paths', 'model', 'training', 'sample_category', 'select_samples', 'seed', 'plot_results']
        for key in required_keys:
            if key not in config_dict:
                raise KeyError(f"[load_config] Missing required key: {key} in the configuration file.")

        # Create the Config object
        config = Config(**config_dict)

        # If output_dir is not specified, generate a unique path
        if not config.output_dir.strip():
            print('[load_config] output_dir is empty, generating a unique path... 🎯')
            unique_id = generate_unique_id(config)
            config.output_dir = f'/scratch/jsanchoz/DeepRBP/output/results/analysis/{unique_id}/train_prediction_model'
        return config

    except FileNotFoundError:
        print(f"[load_config] Error: The YAML file '{yaml_path}' was not found.")
        raise
    except yaml.YAMLError as e:
        print(f"[load_config] Error parsing the YAML file: {e}")
        raise
    except Exception as e:
        print(f"[load_config] An unexpected error occurred: {e}")
        raise

