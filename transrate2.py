'''
██████┐██████┐  █████┐ ███┐  ██┐ ██████┐██████┐  █████┐ ██████┐███████┐██████┐ 
└─██┌─┘██┌──██┐██┌──██┐████┐ ██│██┌────┘██┌──██┐██┌──██┐└─██┌─┘██┌────┘└────██┐
  ▓▓│  ▓▓▓▓▓▓┌┘▓▓▓▓▓▓▓│▓▓┌▓▓┐▓▓│└▓▓▓▓▓┐ ▓▓▓▓▓▓┌┘▓▓▓▓▓▓▓│  ▓▓│  ▓▓▓▓▓┐    ▓▓▓┌─┘
  ▒▒│  ▒▒┌──▒▒┐▒▒┌──▒▒│▒▒│└▒▒▒▒│ └───▒▒┐▒▒┌──▒▒┐▒▒┌──▒▒│  ▒▒│  ▒▒┌──┘  ▒▒┌──┘  
  ░░│  ░░│  ░░│░░│  ░░│░░│ └░░░│░░░░░░┌┘░░│  ░░│░░│  ░░│  ░░│  ░░░░░░░┐░░░░░░░┐
  └─┘  └─┘  └─┘└─┘  └─┘└─┘  └──┘└─────┘ └─┘  └─┘└─┘  └─┘  └─┘  └──────┘└──────┘
'''

import os
import sys
import argparse
from glob                import glob
from pathlib             import Path
from core.utils.logo     import print_logo, print_help, print_args
from core.datahub        import DataHub
from core.utils.deps     import Deps
from core.utils.printout import PrintOut
from core.utils.config   import ConfigManager

'''
Skip flags:
m: Skip Aligner
a: Skip Assembly
c: Skip Contig
r: Skip Reference
s: Skip Summary
'''

MAJOR   = 2                                                            
MINOR   = 9                                                            
PATCH   = 3

VERSION = f"{MAJOR}.{MINOR}.{PATCH}"                                   

class TransRate2:
    def __init__(self):
        from core.utils.logo import get_transrate_colors
        yellow, green, red = get_transrate_colors()
        self.highlight_color, self.background_color = yellow[1], yellow[0]
        
        self.printClass        = PrintOut(1, self.highlight_color, self.background_color)                                                                      
        self.printout          = self.printClass.printout                                                         
        self.config_manager    = ConfigManager(self.highlight_color, self.background_color, VERSION)                                                        
        self.is_gz             = False
        self.quiet             = False
        self.nocolor           = False                                                                           

        self.original_defaults                      = self.config_manager.get_defaults_dict().copy()                                   
        self.original_defaults['config']            = None
        self.original_defaults['config_create']     = False
        self.original_defaults['log']               = 3
        self.original_defaults['help']              = False
        self.original_defaults['version']           = False

    def parser(self, config_values=None):
        parser = argparse.ArgumentParser(description="TransRate2", add_help=False)

        self.defaults_dict                    = self.config_manager.get_defaults_dict()
        self.defaults_dict['config']          = None
        self.defaults_dict['config_create']   = False
        self.defaults_dict['log']             = 3
        self.defaults_dict['help']            = False
        self.defaults_dict['version']         = False
        self.defaults_dict['SAVE']            = False
        
        if config_values:
            for key, value in config_values.items():
                if key in self.defaults_dict:
                    self.defaults_dict[key] = value

        def error(message):
            quiet_mode = '--quiet' in sys.argv or '-q' in sys.argv
            if not quiet_mode:
                print_help(self.highlight_color, self.defaults_dict)
                message       = message.replace('unknown option:', '')                                   
                message       = message.strip()                                                          
                chunks        = []                                                                       
                current_chunk = ""                                                                       
                words         = message.split()                                                          
                for word in words:
                    if len(current_chunk) + len(word) + 1 <= 58:
                        if current_chunk:
                            current_chunk += " " + word
                        else:
                            current_chunk = word
                    else:
                        if current_chunk:
                            chunks.append(current_chunk)
                        current_chunk = word
                if current_chunk:
                    chunks.append(current_chunk)
                for chunk in chunks:
                    self.printout('error', chunk)
                self.printout('error', 'Try --help for more information.')
            sys.exit(2)
        parser.error = error

        # Basic
        parser.add_argument("--output-dir", "-o",   type=str,    help="Output directory",                           default=self.defaults_dict['output_dir'])
        parser.add_argument("--threads",    "-t",   type=int,    help="Number of threads for parallel processing",  default=self.defaults_dict['threads'])
        parser.add_argument("--clutter",    "-c",                help="Remove intermediate files after completion", default=self.defaults_dict['clutter'],          action="store_true")

        # Assembly
        parser.add_argument("--assembly",   "-a",   type=str,    help="Assembly file(s), comma separated",          default=self.defaults_dict['assembly'])
        parser.add_argument("--left",       "-l",   type=str,    help="Left reads file",                            default=self.defaults_dict['left'])
        parser.add_argument("--right",      "-r",   type=str,    help="Right reads file",                           default=self.defaults_dict['right'])
        parser.add_argument("--reference",  "-f",   type=str,    help="Reference file",                             default=self.defaults_dict['reference'])
        parser.add_argument("--bam",        "-x",   type=str,    help="BAM file",                                   default=self.defaults_dict['bam'])

        # Aligner
        parser.add_argument("--bowtie2",    "-b",                help="Use Bowtie2 aligner", action="store_true")
        parser.add_argument("--hisat2",     "-s",                help="Use Hisat2 aligner",                        action="store_true")
        
        # Config
        parser.add_argument("--config",             type=str,    help="Path to configuration file",                 default=None,                                   nargs="?", const=True)
        parser.add_argument("--config-create",      type=str,    help="Create a configuration template",            default=False,                                  nargs="?", const="config.yaml")
        parser.add_argument("--config-save",        type=str,    help="Save current arguments to config file",      default=False,                                  nargs="?", const="config.yaml")

        # Standard
        parser.add_argument("--version",    "-v",                help="Print version",                                                                              action="store_true")
        parser.add_argument("--log",                type=int,    help=argparse.SUPPRESS,                            default=self.defaults_dict['log'],              choices=[0, 1, 2, 3, 4])
        parser.add_argument("--help",       "-h",                help=argparse.SUPPRESS,                            default=self.defaults_dict['help'],             action="store_true")
        parser.add_argument("--quiet",      "-q",                help="Suppress terminal output",                   default=self.defaults_dict['quiet'],            action="store_true")
        parser.add_argument("--debug",                           help="Enable debug mode",                          default=self.defaults_dict['debug'],            action="store_true")
        parser.add_argument("--nocolor",                         help="Disable colored terminal output",            default=self.defaults_dict['nocolor'],          action="store_true")

        # DEV
        parser.add_argument("--mode_multi",                      help=argparse.SUPPRESS,                            default=False,                                  action="store_true")
        
        return parser.parse_args()
    
    
    def _validate_files(self, args):
        if not any([args.assembly, args.left, args.right, args.reference, args.bam]):
            self.printout('error', 'At least one input file must be provided')
            sys.exit(1)
        
        if args.bam:
            if not (args.left or args.right):
                self.printout('error', 'BAM file requires at least one reads file (--left or --right)')
                sys.exit(1)

        if args.assembly:
            assembly_files = [f.strip() for f in args.assembly.split(',')]
            for assembly_file in assembly_files:
                if not os.path.exists(assembly_file):
                    self.printout('error', f'Assembly file does not exist: {assembly_file}')
                    sys.exit(1)
                
                file_extension = os.path.splitext(assembly_file)[1]
                if file_extension == '.gz':
                    base_name      = os.path.splitext(assembly_file)[0]                                                
                    file_extension = os.path.splitext(base_name)[1] + file_extension                                   
                
                if file_extension not in ['.fasta', '.fa', '.fna', '.faa', '.fa.gz', '.fna.gz', '.faa.gz']:
                    self.printout('error', f'Invalid assembly file extension: {file_extension}')
                    sys.exit(1)

        files_to_check = [
            (args.left,     ['.fastq', '.fq', '.fastq.gz', '.fq.gz']),
            (args.right,    ['.fastq', '.fq', '.fastq.gz', '.fq.gz']),
            (args.reference,['.fasta', '.fa', '.fna', '.fa.gz', '.fna.gz', '.fa.gz']),
            (args.bam,      ['.bam'])
        ]
        
        for file_path, expected_extensions in files_to_check:
            if file_path and not os.path.exists(file_path):
                self.printout('error', f'File does not exist: {file_path}')
                sys.exit(1)
            
            if file_path:
                file_extension = os.path.splitext(file_path)[1]
                
                if file_extension == '.gz':
                    base_name      = os.path.splitext(file_path)[0]                                                    
                    file_extension = os.path.splitext(base_name)[1] + file_extension                                   
                    self.is_gz     = True                                                                              
                
                if file_extension not in expected_extensions:
                    self.printout('error', f'Invalid file extension: {file_extension}')
                    sys.exit(1)

    def _mode_check(self, args):
        if args.bam and not args.assembly:
            # BAM-only mode: BAM file without assembly, limited analysis (no assembly stats)
            args.mode = 4
        elif args.bam and args.assembly and not (args.left or args.right):
            # BAM + Assembly but no reads
            args.mode = 3
        elif args.assembly and not (args.left or args.right):
            # Assembly only mode
            args.mode = 0
        elif (args.assembly or args.bam) and (args.left or args.right) and (bool(args.left) != bool(args.right)):
            # Assembly/BAM + single-end reads mode
            args.mode = 1
        elif (args.assembly or args.bam) and (args.left or args.right) and (bool(args.left) == bool(args.right)):
            # Assembly/BAM + paired-end reads mode
            args.mode = 2

        if args.bowtie2:
            args.aligner = 'bowtie2'
            args.hisat2  = False
        elif args.hisat2: 
            args.aligner = 'hisat2'
            args.bowtie2 = False
        else: 
            args.aligner = 'bowtie2'
            args.bowtie2 = True
            args.hisat2  = False

    def _run_multi_mode(self, args, assembly_files, passed_args):
        args.mode_multi   = True                                            
        original_assembly = args.assembly                                   

        all_results       = []                                              
        
        for i, assembly_file in enumerate(assembly_files, 1):
            try:
                self.printout('metric', f'Processing assembly {i}/{len(assembly_files)}: {os.path.basename(assembly_file)}')
                
                args.assembly = assembly_file
                
                dataHub = DataHub(args)
                dataHub.run()
                
                results = dataHub.get_assembly_results()
                all_results.append(results)
                
            except Exception as e:
                self.printout('error', f'Failed to process assembly {assembly_file}: {str(e)}')
                continue                                                                   

    def run(self):
        if len(sys.argv) == 1:
            if not self.quiet:
                print_logo(VERSION, nocolor=self.nocolor)
            print_help(self.highlight_color, self.original_defaults, nocolor=self.nocolor)
            sys.exit()
            
        args = self.parser()
        
        self.quiet = getattr(args, 'quiet', False)
        self.printClass.set_quiet(self.quiet)
        self.config_manager.set_quiet(self.quiet)
        
        self.nocolor = getattr(args, 'nocolor', False)
        if self.nocolor:
            self.highlight_color  = ''
            self.background_color = ''
            self.printClass.set_nocolor(True)
            self.config_manager.set_nocolor(True)
        
        display_only_flags = ['nocolor', 'debug']
        has_display_flags  = any(getattr(args, flag, False) for flag in display_only_flags if hasattr(args, flag))
        has_input_files    = any([args.assembly, args.left, args.right, args.reference, args.bam])
        has_other_commands = any([
            args.help, args.version, getattr(args, 'config_create', False), 
            getattr(args, 'config_save', False), getattr(args, 'config', None) is not None
        ])
        
        has_only_display_flags = has_display_flags and not has_input_files and not has_other_commands
        
        if has_only_display_flags:
            if not self.quiet:
                print_logo(VERSION, nocolor=self.nocolor)
                print_help(self.highlight_color, self.original_defaults, nocolor=self.nocolor)
            sys.exit()
        
        if not self.quiet:
            print_logo(VERSION, nocolor=self.nocolor)
        
        args.mode_multi = False                                   
        args.mode       = 0                                       

        if args.help:
            print_help(self.highlight_color, self.original_defaults, self.quiet, self.nocolor)
            sys.exit()
        if args.version:
            if not self.quiet:
                print(VERSION)
            sys.exit()
        
        if getattr(args, 'config_create', False):
            config_filename = args.config_create if args.config_create != "config.yaml" else "config.yaml"
            if not config_filename.endswith('.yaml'):
                config_filename += '.yaml'
            config_path = Path(args.output_dir) / config_filename if args.output_dir else Path(config_filename)
            self.config_manager.create_config(config_path)
            sys.exit(0)

        if args.output_dir:
            Path(args.output_dir).mkdir(parents=True, exist_ok=True)

        config_values = None
        if args.config is not None:
            if args.config is True:
                self.printout('error', 'No configuration file path provided with --config')
                self.printout('error', 'Usage: --config /path/to/config.yaml')
                self.printout('error', 'Or use --config-create to create a new config file')
                sys.exit(1)
            else:
                config_path = Path(args.config)
                if config_path.exists():
                    config_values = self.config_manager.load_config(config_path)
                else:
                    self.printout('error', f'Config file not found: {config_path}')
                    sys.exit(1)
        
        if config_values:
            args = self.parser(config_values)

        config_save_used = getattr(args, 'config_save', False)
        if config_save_used:
            config_dict = {}
            for key, value in args.__dict__.items():
                if key in self.config_manager.get_defaults_dict():
                    config_dict[key] = value

            config_filename = args.config_save if args.config_save != "config.yaml" else "config.yaml"
            if not config_filename.endswith('.yaml'):
                config_filename += '.yaml'
            config_path = Path(args.output_dir) / config_filename if args.output_dir else Path(config_filename)
            self.config_manager.save_config(config_dict, config_path)
            out_name = str(config_path.absolute()) if len(str(config_path.absolute())) < 32 else '...' + str(config_path.absolute())[-29:]
            self.printout('info', f"Current arguments saved to {out_name}")
        
        args.highlight_color  = self.highlight_color
        args.background_color = self.background_color
        
        self._validate_files(args)
        self._mode_check(args)
        
        passed_args = {k: v for k, v in args.__dict__.items()
                       if k in self.original_defaults and v != self.original_defaults[k] and v != '-'}
        
        if hasattr(args, 'aligner'):
            passed_args['aligner'] = args.aligner
        
        print_args(args, passed_args)
        deps = Deps(args.log, self.highlight_color, self.background_color, nocolor=self.nocolor)
        deps.check_deps(quiet=self.quiet)
        
        if args.assembly:
            assembly_files = [f.strip() for f in args.assembly.split(',')]
            if len(assembly_files) > 1:
                self._run_multi_mode(args, assembly_files, passed_args)
            else:
                dataHub = DataHub(args)
                dataHub.run()
        else:
            dataHub = DataHub(args)
            dataHub.run() 

if __name__ == "__main__":
    main = TransRate2()
    main.run()