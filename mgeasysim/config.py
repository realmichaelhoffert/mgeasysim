import os
import yaml

class Config:
    def __init__(self):
        # Path to the configuration file
        config_path = os.path.join(os.path.dirname(__file__), "config.yaml")
        self._config_data = self._load_config(config_path)

    def _load_config(self, config_path):
        with open(config_path, "r") as file:
            return yaml.safe_load(file)

    def get(self, section, option, default=None):
        return self._config_data.get(section, {}).get(option, default)


# Singleton pattern to ensure one instance of Config across the app
config = Config()