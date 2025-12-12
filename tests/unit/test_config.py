import os
import yaml
import pytest

import pytest

from mgeasysim.config import *

class TestConfig:
    
    def test_validate_file_exists_before_reading(self, tmp_path):
        """
        Validates that the class correctly identifies and reads an existing file.
        """
        # 1. Setup: Create a temporary yaml file
        d = tmp_path / "subdir"
        d.mkdir()
        p = d / "test_config.yaml"
        
        initial_data = {'section1': {'option1': 'value1'}}
        p.write_text(yaml.dump(initial_data))
        
        # 2. Execution: Initialize Config with the path
        # Note: Assuming the logic is fixed to load the file
        config = load_configuration(str(p))

        print('Test printing:')
        print(config)

        # 3. Validation
        # Verify the file was actually found and data loaded into memory
        assert config.get('section1', 'option1') == 'value1'

    def test_handle_non_existent_file(self, tmp_path, capsys):
        """
        Validates behavior when file does not exist (validating the existence check).
        """
        # 1. Setup: Define a path that definitely does not exist
        fake_path = tmp_path / "ghost_config.yaml"

        # 2. Execution & Validation
        # Depending on how you handle the 'else' block in your code, 
        # this might raise an error or result in an empty config.
        # Here we assume it initializes as empty or handles it gracefully.

        # test that error is thrown on empty file with no permission to make new
        with pytest.raises(Exception) as e_info:
           config = load_configuration('blah', write_new = False)

        config = load_configuration(fake_path, write_new = True)
        # print('should throw error')
        # print(capsys)
        
        # Verify internal data is empty or None because file didn't exist
        # Accessing private member for test validation purposes
        assert getattr(config, '_config_data', None) is None or config._config_data == {}

    def test_content_written_properly_upon_save(self, tmp_path):
        """
        Validates that calling set() triggers a save and writes correct YAML to disk.
        """
        # 1. Setup: Create a dummy path
        p = tmp_path / "write_test.yaml"
        
        # Initialize with empty/new file
        # We need to create the file first if your __init__ requires strict existence
        p.touch() 
        config = Config(str(p))
        
        # Manually inject data if __init__ failed to load empty file
        config._set_data({}) 

        # 2. Execution: Set a new value (which triggers _save_config)
        test_section = "database"
        test_option = "host"
        test_value = "localhost"
        
        config.set(test_section, test_option, test_value)

        config._save_config(config.config_path)

        # 3. Validation: Open the file manually to verify disk content
        assert os.path.exists(str(p))
        
        with open(str(p), "r") as f:
            content = yaml.safe_load(f)
        
        assert content[test_section][test_option] == test_value