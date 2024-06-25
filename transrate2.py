import __init__
import sys
import os
import argparse
import psutil
from modules.run import MAIN
from modules.terminal.reqs import req_check
from modules.terminal.logo import Logo
from modules.terminal.helpoptions import Help
from modules.terminal.path import PathCheck


if __name__ == "__main__":
    class Parser:
        def __init__(self):
            #### Info ####
            self.version     = __init__.__version__
            self.author      = __init__.__author__
            self.description = __init__.__description__

            #### Settings ####
            self.quiet       = False
            self.skip        = False
            self.plot        = False
            self.clutter     = False
            self.threads     = 1

            #### Args ####
            self.args = None

            #### Prints ####
            self.logo        = Logo(self.version).logo
            self.helpblocks  = Help()

            #### Terminal Colors ####
            self.colors = {
                'red'   : '\033[0;31m',
                'green' : '\033[0;32m',
                'yellow': '\033[0;33m',
                'blue'  : '\033[0;34m',
                'purple': '\033[0;35m',
                'cyan'  : '\033[0;36m',
                'white' : '\033[0;37m',
                'reset' : '\033[0m'
            }
            self.term_colors = [self.colors['yellow'], self.colors['red'], self.colors['green']]
            self.color       = self.colors['yellow']
            self.reset  = self.colors['reset']

            #### Assembly Info ####
            self.assembly_count = 0
            self.assembly_run   = 1
            self.assembly_name  = ''
            self.assembly_base  = ''
            self.assembly_multi = False
            self.assembly_list  = []

            #### Read Info ####
            self.singleend   = False
            self.single      = None
            self.pairedend   = False
            self.solo        = False

            #### Aligner Info ####
            self.aligner     = None

            #### Output Info ####
            self.output      = None

            #### Initiate ####
            # Parse arguments
            self.parse()
            # Check requirements
            self.check_requirements()
            # Check for multiple assemblies
            self.assembly_counter()
            # Set color
            self.color_iter()
            # Check for multiple aligners
            self.aligners()
            # Check for valid paths
            self.path_check()
            # Check for valid threads
            self.thread_check()
            # Run
            self.run()

        def parse(self):
            ############ Arguments ##################################
            # --assembly, -a: Assembly file                         #
            # --left, -l: Left reads file                           #
            # --right, -r: Right reads file                         #
            # --reference, -f: Reference file                       #
            # --output, -o: Output directory                        #
            # --threads, -t: Number of threads                      #
            # --star, -s: Use STAR aligner (default)                #
            # --snap, -p: Use Snap aligner                          #
            # --bowtie2, -b: Use Bowtie2 aligner                    #
            # --clutter, -c: Remove clutter from output directory   #
            # --help, -h: Display this help message                 #
            # --quiet, -q: Supress terminal output                  #
            # --plot, -P: Create .png plots of results              #
            # --skip, -k: Skip to transrate                         #
            # --version, -v: Display version                        #
            #########################################################

            parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter, description='Transrate2')
            parser.add_argument('--assembly', '-a', type=str , help='Assembly file')
            parser.add_argument('--left', '-l', type=str , help='Left reads file')
            parser.add_argument('--right', '-r', type=str , help='Right reads file')
            parser.add_argument('--reference', '-f', type=str , help='Reference file')
            parser.add_argument('--output', '-o', type=str , help='Output directory')
            parser.add_argument('--threads', '-t', type=int , help='Number of threads')
            parser.add_argument('--star', '-s', action='store_true', default=False, help='Use STAR aligner (default)')
            parser.add_argument('--snap', '-p', action='store_true', help='Use Snap aligner')
            parser.add_argument('--bowtie2', '-b', action='store_true', help='Use Bowtie2 aligner')
            parser.add_argument('--clutter', '-c', action='store_true', help='Remove clutter from output directory')
            parser.add_argument('--help', '-h', action='store_true', help='Display this help message')
            parser.add_argument('--quiet', '-q', action='store_true', help='Supress terminal output')
            parser.add_argument('--plot', '-P', action='store_true', help='Create .png plots of results')
            parser.add_argument('--skip', '-k', action='store_true', help='Skip to transrate')
            parser.add_argument('--version', '-v', action='store_true', help='Display version')

            self.args = parser.parse_args()

            # Print help message
            if self.args.help or len(sys.argv) < 2:
                print(f'{self.logo}')
                print(f'{self.helpblocks}')
                sys.exit()

            # Print version
            if self.args.version:
                print(f'{__init__.__version__}')
                sys.exit()

            # Supress output
            if not self.args.quiet:
                print(f'{self.logo}')

            # Settings
            self.threads = int(self.args.threads) if self.args.threads else self.threads
            self.quiet   = self.args.quiet
            self.plot    = self.args.plot
            self.skip    = self.args.skip
            self.clutter = self.args.clutter

        # Check for required programs
        def check_requirements(self):
            star = True if (self.args.star + self.args.bowtie2 + self.args.snap) == 0 else True if self.args.star else False
            req_check(star, self.args.bowtie2, self.args.snap)

        # Count assemblies in run
        def assembly_counter(self):
            if not self.args.assembly:
                return
            assembly = self.args.assembly.split(',')
            self.assembly_count = len(assembly)
            if self.assembly_count > 1:
                self.assembly_multi = True
            self.assembly_list  = assembly

        # Name first assembly run
        def assembly_namer(self, assembly):
            self.assembly_name = os.path.basename(assembly)
            self.assembly_base = os.path.splitext(self.assembly_name)[0]
            return self.assembly_base

        def color_iter(self):
            self.color = self.term_colors[(self.assembly_run % 3) - 1]

        # Check for multiple aligners
        def aligners(self):
            if self.args.star + self.args.bowtie2 + self.args.snap > 1:
                print(f'{self.colors["red"]}  ┌{"─" * 27}{self.reset}    Aligner Error   {self.colors["red"]}{"─" * 27}┐{self.reset}')
                print(f'{self.colors["red"]}  │ {self.reset} Please choose only one aligner to use.{self.colors["red"]}{" " * 34}│{self.reset}')
                print(f'{self.colors["red"]}  │  •{self.reset} STAR    (-s) is the default aligner.{self.colors["red"]}{" " * 34}│{self.reset}')
                print(f'{self.colors["red"]}  │  •{self.reset} Snap    (-p){" " * 58}{self.colors["red"]}│{self.reset}')
                print(f'{self.colors["red"]}  │  •{self.reset} Bowtie2 (-b){" " * 58}{self.colors["red"]}│{self.reset}')
                print(f'{self.colors["red"]}  │ {self.reset} Or use (--help) for more options.{self.colors["red"]}{" " * 39}│{self.reset}')
                print(f'{self.colors["red"]}  └{"─" * 74}┘{self.reset}')
                sys.exit(1)
            # Set aligner
            if self.args.bowtie2:
                self.aligner = 'bowtie2'
            elif self.args.snap:
                self.aligner = 'snap'
            else:
                self.aligner = 'star'

        # Check for valid paths
        def path_check(self):
            path  = PathCheck(self.assembly_list, self.args.left, self.args.right, self.args.output, self.args.reference, self.args)
            error = path.check()
            if error:
                error = f'...{error[-61:]}' if len(error) > 64 else error
                print(f'{self.colors["red"]}  ┌{"─" * 29}{self.reset}{"Path Error":^16}{self.colors["red"]}{"─" * 29}┐{self.reset}')
                print(f'{self.colors["red"]}  │ {self.reset} The following path does not exist:{self.colors["red"]}{" " * 38}│{self.reset}')
                print(f'{self.colors["red"]}  │ {self.reset} {error}{" " * (72 - len(error))}{self.colors["red"]}│{self.reset}')
                print(f'{self.colors["red"]}  └{"─" * 74}┘{self.reset}')
                sys.exit(1)

            # Check if single, paired, or solo
            if bool(self.args.left) != bool(self.args.right):
                self.singleend = True
                self.single = self.args.left if self.args.left else self.args.right
            elif bool(self.args.left) + bool(self.args.right) != 0:
                self.pairedend = True
            else:
                self.solo = True
            # Set output
            if not self.args.output:
                self.args.output = os.getcwd()
                self.output = os.getcwd()
            else:
                self.output = self.args.output

        # Check for valid threads
        def thread_check(self):
            # Limit threads to available cores
            self.threads = max(min(self.threads, psutil.cpu_count(logical=True)), 4)    

        def parameter_block(self):
            assembly = f'...{"".join(self.assembly_name[-61:])}' if len(self.assembly_name) > 64 else self.assembly_name
            left = ''
            right = ''
            reads = ''
            reference = ''
            output = ''

            if self.singleend:
                reads = f'...{"".join(os.path.basename(self.single)[-61:])}' if len(self.single) > 64 else os.path.basename(self.single)
            elif not self.solo:
                left = f'...{"".join(os.path.basename(self.args.left)[-61:])}' if len(self.args.left) > 64 else os.path.basename(self.args.left)
                right = f'...{"".join(self.args.right[-61:])}' if len(self.args.right) > 64 else os.path.basename(self.args.right)
            if self.args.reference:
                reference = f'...{"".join(os.path.basename(self.args.reference)[-61:])}' if len(self.args.reference) > 64 else os.path.basename(self.args.reference)
            output  = f'...{"".join(self.args.output[-61:])}' if len(self.output) > 64 else self.output
            aligner = self.aligner.capitalize()
            threads = self.threads
            plots   = "True" if self.args.plot else None
            run = f'    ({self.assembly_run}/{self.assembly_count})' if self.assembly_count > 1 else ''
            clutter = "True" if self.clutter else None
            param_dict = {
                'Assembly': assembly + run,
                'Left': left,
                'Right': right,
                'Reads': reads,
                'Reference': reference,
                'Output': output,
                'Aligner': aligner,
                'Threads': threads,
                'Plots': plots,
                'Clutter': clutter
            }
            print(f'{self.color}  ┌{"─" * 74}┐{self.reset}')
            for k,v in param_dict.items():
                if v:
                    print(f'{self.color}  │ {self.reset} {k:<10} {v:<61}{self.color}│{self.reset}')
            print(f'{self.color}  └{"─" * 74}┘{self.reset}')

        def run(self):
            t2 = MAIN()
            for i,a in enumerate(self.assembly_list):
                if i == 0:
                    self.blank_dct = t2.__dict__.copy()
                
                #### Settings ####
                t2.__dict__       = self.blank_dct.copy()
                t2.TR_DCT         = {}
                t2.QUIET          = self.quiet
                t2.SKIP           = self.skip
                t2.CLUTTER        = self.clutter
                t2.THREADS        = self.threads
                t2.COLOR          = self.color

                #### Analyses ####
                # Assembly
                t2.ASSEMBLY       = a
                t2.ASSEMBLY_COUNT = self.assembly_count
                t2.ASSEMBLY_MULTI = self.assembly_multi
                t2.ASSEMBLY_RUN   = self.assembly_run
                t2.ASSEMBLY_NAME  = self.assembly_namer(a)

                # Reads
                t2.LEFT           = self.args.left
                t2.RIGHT          = self.args.right
                t2.READMODE       = 1 if self.singleend else 2 if self.pairedend else 0
                t2.READS          = True if self.singleend or self.pairedend else False
                
                # Plotting
                t2.PLOT           = self.args.plot

                # Reference
                t2.REFERENCE      = True if self.args.reference else False
                t2.REFERENCE_FILE = self.args.reference

                # Output
                t2.OUTPUT         = self.args.output

                # Aligner
                t2.ALIGNER        = self.aligner

                #### Run ####
                self.parameter_block() if not self.quiet else None
                t2.run()

                #### Iterate ####
                self.assembly_run += 1
                self.color_iter()
    term = Parser()
