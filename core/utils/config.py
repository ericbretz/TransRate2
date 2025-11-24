import os
import yaml
from pathlib             import Path
from core.utils.printout import PrintOut

class ConfigManager:    
    def __init__(self, highlight_color, background_color, version='0.1.0'):
        self.highlight_color  = highlight_color
        self.background_color = background_color
        self.version          = version
        self.config_dir       = Path.home() / '.transrate2'
        self.printClass       = PrintOut('', self.highlight_color, self.background_color)
        self.printout         = self.printClass.printout
        
    def set_quiet(self, quiet: bool):
        self.printClass.set_quiet(quiet)
    
    def set_nocolor(self, nocolor: bool):
        self.printClass.set_nocolor(nocolor)
        
    def get_defaults_dict(self):
        return {
            # Basic
            'output_dir'              : os.getcwd(),
            'threads'                 : 4,
            'clutter'                 : False,
            
            # Assembly
            'assembly'                : None,
            'left'                    : None,
            'right'                   : None,
            'reference'               : None,
            'bam'                     : None,
            
            # Aligner
            'bowtie2'                 : True,
            'hisat2'                  : False,
            
            # Options
            'quiet'                   : False,
            'debug'                   : False,
            'nocolor'                 : False,
        }
    
    def get_parameter_types(self):
        return {
            # Basic
            'output_dir'              : 'string', 
            'threads'                 : 'integer',
            'clutter'                 : 'boolean',
            
            # Assembly
            'assembly'                : 'string',
            'left'                    : 'string',
            'right'                   : 'string',
            'reference'               : 'string',
            'bam'                     : 'string',
            
            # Aligner
            'bowtie2'                 : 'boolean',
            'hisat2'                  : 'boolean',
            
            # Options
            'quiet'                   : 'boolean',
            'debug'                   : 'boolean',
            'nocolor'                 : 'boolean',
        }
    
    def load_config(self, config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _extract_values(self, config):
        values = {}
        
        for section_name, section in config.items():
            if section_name == 'transrate2_config':
                continue
                
            if isinstance(section, dict):
                if 'description' in section:
                    for param_name, param_data in section.items():
                        if param_name != 'description' and isinstance(param_data, dict):
                            if 'value' in param_data:
                                values[param_name] = param_data['value']
                else:
                    for key, value in section.items():
                        if key != 'transrate2_config':
                            values[key] = value
        return values
    
    def save_config(self, config, output_path):
        self.config_dir.mkdir(exist_ok=True)
        
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False, indent=2)
    
    def create_config(self, output_path):
        config = {
            'transrate2_config': {
                'version': self.version,
                'description': 'TransRate2 Configuration File'
            },
            **self.get_defaults_dict()
        }
        self.save_config(config, output_path)
        out_name = str(output_path)
        out_name = out_name if len(out_name) < 32 else '...' + out_name[-29:]
        self.printout('info', f"Config created at: {out_name}")
    
    def validate_config(self, config):
        parameter_types = self.get_parameter_types()
        defaults        = self.get_defaults_dict()
        
        for key, value in config.items():
            if key == 'transrate2_config':
                continue
                
            if key not in defaults:
                self.printout('error', f"Unknown config parameter: {key}")
                continue
            
            expected_type = parameter_types.get(key, 'any')
            if expected_type == 'integer' and not isinstance(value, int):
                self.printout('error', f"{key} must be an integer, got {type(value).__name__}")
                return False
            elif expected_type == 'float' and not isinstance(value, (int, float)):
                self.printout('error', f"{key} must be a number, got {type(value).__name__}")
                return False
            elif expected_type == 'boolean' and not isinstance(value, bool):
                self.printout('error', f"{key} must be a boolean, got {type(value).__name__}")
                return False
            elif expected_type == 'string' and value is not None and not isinstance(value, str):
                self.printout('error', f"{key} must be a string, got {type(value).__name__}")
                return False
        
        return True

if __name__ == "__main__":
    config_manager = ConfigManager('\033[94m', '\033[44m', '0.1.0')
    config_manager.create_config(Path('./test_config.yaml'))
