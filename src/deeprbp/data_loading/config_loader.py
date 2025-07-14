# src/deeprbp/data_loading/config_loader.py

import os
import yaml

class ConfigParser:
    def __init__(self, config_path=None, **kwargs):
        """
        Initializes the ConfigParser with the path to the YAML configuration file or keyword arguments.

        Parameters:
        - config_path (str): Optional path to the YAML configuration file.
        - **kwargs: Key-value pairs to use as configuration data.
        """
        self.config_path = config_path
        if config_path:
            self.config_data = self._load_config()
        else:
            self.config_data = kwargs  # use provided kwargs directly 
    def _load_config(self):
        """
        Loads the configuration data from the YAML file.

        Returns:
        - dict: The configuration data as a dictionary.
        """
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(f"The configuration file {self.config_path} does not exist.")
        with open(self.config_path, 'r') as file:
            config_data = yaml.safe_load(file)
        return config_data
    def get(self, key, default=None):
        """
        Retrieves the value for a given key from the configuration data.

        Parameters:
        - key (str): The key to retrieve from the configuration data.
        - default: The default value to return if the key does not exist.

        Returns:
        - The value associated with the key, or the default value if not found.
        """
        return self.config_data.get(key, default)
    def update(self, key, value):
        """
        Updates a specific key in the configuration data.

        Parameters:
        - key (str): The key to update.
        - value: The new value to set for the key.
        """
        self.config_data[key] = value
    def __str__(self):
        """
        Returns a string representation of the configuration data.

        Returns:
        - str: A formatted string of the configuration data.
        """
        return yaml.dump(self.config_data, default_flow_style=False)