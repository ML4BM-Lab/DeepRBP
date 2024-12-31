# /scratch/jsanchoz/DeepRBP/src/deeprbp/config_loader.py
import yaml
from datetime import datetime
import random
import string
from dataclasses import dataclass, field
from typing import Dict, List, Any, Tuple

@dataclass
class Config:
    """
    Represents the configuration loaded from a YAML file.
    """
    source_name: str
    data_paths: Dict[str, str]
    sample_category: str
    select_samples: List[str]
    model: Dict[str, Any]
    training: Dict[str, Any] = field(default_factory=dict)
    explainability: Dict[str, Any] = field(default_factory=dict)
    output_dir: str = field(default_factory=str)
    seed: int = 0
    plot_results: bool = False

    def __getitem__(self, key: str) -> Any:
        """ Access configuration values using dictionary-style indexing. """
        return self.get(key)
    
    def get(self, key: str, default=None) -> Any:
        """ Retrieve a value for a given key, or raise an exception if key does not exist. """
        value = self[key]
        if value is None and default is None:
            raise KeyError(f"Key '{key}' not found in configuration.")
        return value

class ConfigParser:
    def __init__(self, config_path: str):
        self.config_path = config_path
        self.config = self._load_config()

    def _load_config(self) -> Config:
        """ Loads the YAML configuration file into a Config object. """
        try:
            with open(self.config_path, "r") as file:
                config_dict = yaml.safe_load(file)

            # Default values for missing keys
            required_defaults = {
                'source_name': 'default_source',
                'data_paths': {},
                'model': {
                    'input_size': 1348,
                    'output_size': 11459,
                    'num_hidden_layers': 2,
                    'max_node': 1024,
                    'uniform_nodes': True,
                    'node_shrink_factor': 2,
                    'activation_func': "relu",
                    'learning_rate': 0.001,
                    'optimizer_name': "adamW",
                    'cuda': True
                },

                'training': {
                    'epochs': 100,
                    'batch_size': 128,
                    'print_every': 10,
                    'train_test_split': True,
                    'train_val_split': True,
                    'test_fraction': 0.2,
                    'val_fraction': 0.15
                },

                'explainability': {  # Nueva sección para la explicación
                    'trained_model_path': '',
                    'model_file': '',
                    'scaler_path': '',
                    'explanation_method': 'DeepLIFT',
                    'reference_data': 'knockdown_reference',
                    'batch_reduction_method': 't-statistic',
                    'gene_collapse_method': 'max_absolute_value',
                    'postar_matrix_path': '',
                    'postar_file': ''
                },

                'sample_category': 'detailed_category',
                'select_samples': ['all'],  # Valor por defecto si no se especifica
                'output_dir': '',
                'seed': 0,
                'plot_results': False
            }
            
            # Fill missing keys with defaults
            for key, default_value in required_defaults.items():
                config_dict.setdefault(key, default_value)

            if not config_dict['output_dir']:
                unique_id = self._generate_unique_id(config_dict)
                config_dict['output_dir'] = f'/scratch/jsanchoz/DeepRBP/output/results/analysis/{unique_id}/train_prediction_model'
            
            # Merge any existing model and training configurations with the defaults
            config_dict['model'] = {**required_defaults['model'], **config_dict.get('model', {})}
            config_dict['training'] = {**required_defaults['training'], **config_dict.get('training', {})}
            # Create the Config object
            return Config(**config_dict)
        
        except FileNotFoundError:
            raise FileNotFoundError(f"[ConfigParser] The YAML file '{self.config_path}' was not found.")
        except yaml.YAMLError as e:
            raise ValueError(f"[ConfigParser] Error parsing the YAML file: {e}")
        except Exception as e:
            raise RuntimeError(f"[ConfigParser] An unexpected error occurred: {e}")
        
    def _generate_unique_id(self, config_dict: Dict) -> str:
        """Generates a unique identifier for the configuration based on the source name and selected tumor types."""
        timestamp = datetime.now().strftime("%Y-%m-%d")
        source_train = config_dict.get("source_name", "default_source")
        tumor_types = '-'.join([t.split('_')[0] for t in config_dict.get('select_samples', ['no_samples'])])
        suffix = ''.join(random.choices(string.digits, k=3))
        return f"{source_train}_{tumor_types}_{timestamp}_{suffix}"
    
    def get_base_config(self) -> Dict[str, Any]:
        """
        Retrieves the base configuration (data and sample selection).
        
        Returns:
        - dict: The full configuration (including data_paths, sample_category, select_samples, etc.)
        """
        return {
            "source_name": self.config.source_name,
            "data_paths": self.config.data_paths,
            "sample_category": self.config.sample_category,
            "select_samples": self.config.select_samples,
            "output_dir": self.config.output_dir,
            "seed": self.config.seed,
            "plot_results": self.config.plot_results,
        }
    
    def get_model_training_config(self) -> Dict[str, Any]:
        return {
            **self.config.model,  # Merge model configuration
            **self.config.training  # Merge training configuration
        }
    
    def get_explainability_config(self) -> Dict[str, Any]:
        return self.config.explainability 
