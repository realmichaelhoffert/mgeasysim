import os
import yaml
import glob
from pprint import pprint, pformat
import warnings

def load_configuration(path : str, write_new = False):

    """Load configuration file
    returns Config: empty if file doesn't exist and write_new = True
    error if write_new = False
    otherwise tries to load
    """
    cleaned_path = os.path.abspath(path)
    
    if os.access(cleaned_path, os.W_OK):
        with open(path, "r") as file:
            data = yaml.safe_load(file)

        if data is None or data == {}:
            warnings.warn(f"Loaded empty config at {cleaned_path}")
            data = {}
                    
        c = Config(path)
        c._set_data(data)
        
    else:
        if write_new:
            print(f'writing new config to {path}')
            c = Config(path)
            c._set_data({})
            c._save_config(c.config_path)
            
        else:
            raise FileNotFoundError(f'File not found or not writable: {path}')

    return c

def preview_configuration(path : str):
    with open(path, "r") as file:
        data = yaml.safe_load(file)
    print(f'Peek at config {path}:')
    pprint(data)

class Config:
    def __init__(self, fp):
        # Path to the configuration file
        
        # if os.path.exists(os.path.abspath(fp)):
            # config_path = os.path.join(os.path.abspath(fp)), "config.yaml")
        self.config_path = fp
            # self._config_data = self._load_config(config_path)

    # def _clear(self):
    #     self._config_data = {'info':{'software':'mgeasysim'}}
    #     self._save_config(self.config_path)

    def __str__(self):
        if hasattr(self, '_config_data'):
            pathstr = f'mgeasysim.config.Config at {self.config_path}\n'
            return  pathstr + pformat(self._config_data)
        else:
            return pathstr + '{}'
        
    def _set_data(self, data_load):
        self._config_data = data_load
        
    def _save_config(self, config_path):
        # os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, "w") as file:
            yaml.dump(self._config_data, file)
        
    def get(self, section, option):
        if section not in self._config_data:
            raise KeyError(f"Section '{section}' not found in configuration.")
        
        if option not in self._config_data[section]:
            raise KeyError(f"Option '{option}' not found in section '{section}'.")
    
        return self._config_data[section][option]
    
    def set(self, section, option, value):
        # Update in-memory configuration
        if section not in self._config_data:
            self._config_data[section] = {}
        self._config_data[section][option] = value
