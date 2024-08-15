import psutil

class Help:
    def __init__(self):
        self.assembly_block: dict = {}
        self.aligner_block : dict = {}
        self.other_block   : dict = {}
        self.mode_block    : dict = {}
        self.logical       : str  = str(psutil.cpu_count(logical=True))
        self.physical      : str  = str(psutil.cpu_count(logical=False))
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
        self.color = self.colors['yellow']
        self.reset = self.colors['reset']
        
        self.option_blocks()


    def __str__(self):
        block = self.block_build()
        return block

    def option_blocks(self):
        self.assembly_block = {'Assembly': {
                'Assembly'        : ['--assembly', '-a', 'Path to assembly file (FASTA)'],
                'Left Reads'      : ['--left', '-l', 'Path to left reads file (FASTQ)'],
                'Right Reads'     : ['--right', '-r', 'Path to right reads file (FASTQ)'],
                'Reference'       : ['--reference', '-f', 'Path to reference file (FASTA)'],
                'Output Directory': ['--outdir', '-o', 'Path to output directory'],
                'Threads'         : ['--threads', '-t', 'Number of threads to use'],
            }}
        
        self.aligner_block = {'Aligners': {
                'STAR'            : ['--star', '-s', 'Use STAR aligner (default)'],
                'Snap'            : ['--snap', '-p', 'Use Snap aligner'],
                'Bowtie2'         : ['--bowtie2', '-b', 'Use Bowtie2 aligner'],
            }}
        
        self.other_block = {'Others': {
                'Plot'            : ['--plot', '-P', 'Create .png plots of results'],
                'Clutter'         : ['--clutter', '-c', 'Remove intermediate files'],
                'Quiet'           : ['--quiet', '-q', 'Supress terminal output'],
                'Fix'             : ['--fix', '-x', 'Fixes rare bug where Transrate2 hangs'],
                'Help'            : ['--help', '-h', 'Display this help message'],
                'Version'         : ['--version', '-v', 'Display version'],
            }}
        
        self.mode_block = {
                'Assembly' : ['-a', 'Run assembly analysis only.'],
                'Reads'    : ['-a -l -r', 'Run assembly with paired-end reads analysis'],
                'All'      : ['-a -l -r -f', 'Run assembly with paired-end reads and reference'],
                'Single'   : ['-a -l', 'Run assembly with single-end reads analysis'],
                'Reference': ['-a -f', 'Run assembly with reference analysis'],
            }
        return
    
    def top_bars(self, title: str, color: str = 'yellow'):
        bar_left  = f'{self.color}  ┌{"─" * 29}{self.reset}'
        bar_right = f'{self.color}{"─" * 29}┐{self.reset}'
        bar_title = f'{self.reset}{title:^16}{self.reset}'
        return f'{bar_left}{bar_title}{bar_right}'
    
    def bottom_bars(self, color: str = 'yellow'):
        bar_bottom= f'{self.color}  └{"─" * 74}┘{self.reset}'
        return bar_bottom
        
    def block_build(self):
        titles         = ['Help Options', 'Mode Types', '#Threads']
        top_bar        = [self.top_bars(t) for t in titles]
        bottom_bar     = self.bottom_bars() + '\n'
        block_assembly = self.mid_block(self.assembly_block['Assembly'])
        block_aligner  = self.mid_block(self.aligner_block['Aligners'])
        block_other    = self.mid_block(self.other_block['Others'])
        block_mode     = self.mid_block(self.mode_block)
        block_threads  = self.thread_block()
        mid_bar        = f'{self.color}  │ {" " * 73}│{self.reset}'
        title_assembly = self.mid_title('Assembly')
        title_aligner  = self.mid_title('Aligners')
        title_other    = self.mid_title('Others')

        help = f'{top_bar[0]}\n{mid_bar}\n'
        assembly = f'{title_assembly}\n{mid_bar}\n{block_assembly}\n{mid_bar}\n'
        aligner = f'{title_aligner}\n{mid_bar}\n{block_aligner}\n{mid_bar}\n'
        other = f'{title_other}\n{mid_bar}\n{block_other}\n{mid_bar}\n'

        mode = f'{top_bar[1]}\n{mid_bar}\n{block_mode}\n{mid_bar}\n'

        threads = f'{top_bar[2]}\n{mid_bar}\n{block_threads}\n{mid_bar}\n'

        return f'{help}{assembly}{aligner}{other}{bottom_bar}{mode}{bottom_bar}{threads}{bottom_bar}'
        
    def mid_title(self, title: str):
        left  = f'{self.color}  ├{"┄"*29}{self.reset}'
        right = f'{self.color}{"┄"*29}┤{self.reset}'
        mid   = f'{title:^16}'
        return f'{left}{mid}{right}'
    
    def mid_block(self, block: dict):
        out       : list = []
        out_string: str  = ''
        for v in block.values():
            if len(v) == 3:
                v0_len = len(v[0])
                v2_len = len(v[2])
                out.append(f'{self.color}  │{self.reset} {v[0]}{" " * (19-v0_len)}{v[1]}{" " * 8}{v[2]}{" " * (44-v2_len)}{self.color}│{self.reset}')
            elif len(v) == 2:
                v0_len = len(v[0])
                v1_len = len(v[1])
                out.append(f'{self.color}  │{self.reset} {v[0]}{" " * (20-v0_len)}{v[1]}{" " * (53-v1_len)}{self.color}│{self.reset}')
        for o in out:
            out_string += f'{o}\n'
        return out_string.rstrip('\n')
    
    def thread_block(self):
        title = f'{self.color}  │{self.reset} Threads:{" " * 65}{self.color}│{self.reset}'
        log_length = len(self.logical)
        phy_length = len(self.physical)
        phy_str = f'{" " * (4 - phy_length)}{self.physical}'
        log_str = f'{" " * (4 - log_length)}{self.logical}'
        logical  = f'{self.color}  │{self.reset}   Logical: {log_str}{" " * 58}{self.color}│{self.reset}'
        physical = f'{self.color}  │{self.reset}   Physical:{phy_str}{" " * 58}{self.color}│{self.reset}'
        return f'{title}\n{logical}\n{physical}'
    
    def color_set(self, color: str):
        self.color = self.colors[color]
        return

