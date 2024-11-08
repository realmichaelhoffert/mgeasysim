import os
import yaml
from datetime import datetime

class Config:
    def __init__(self):
        # Path to the configuration file

        config_path = os.path.join(os.path.dirname(__file__), "config.yaml")
        self.config_path = config_path
        self._config_data = self._load_config(config_path)

    def _clear(self):
        self._config_data = {'info':{'software':'mgeasysim'}}
        self._save_config(self.config_path)
    
    def _load_config(self, config_path):
        """Load configuration file
        """
        if not os.path.exists(config_path):
            with open(config_path, 'w') as handle:
                yaml.dump({'info':{'software':'mgeasysim'}}, handle)

        with open(config_path, "r") as file:
                return yaml.safe_load(file)
        
    def _save_config(self, config_path):
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, "w") as file:
            yaml.dump(self._config_data, file)
        
    def get(self, section, option, default=None):
        return self._config_data.get(section, {}).get(option, default)
    
    def set(self, section, option, value):
        # Update in-memory configuration
        if section not in self._config_data:
            self._config_data[section] = {}
        self._config_data[section][option] = value
        
        # Write to user configuration file
        self._save_config(self.config_path)
    
# Singleton pattern to ensure one instance of Config across the app
config = Config()