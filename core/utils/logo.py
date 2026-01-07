import os

COLORS = {
    'red'    : {'bg': '\033[41m', 'fg': '\033[31m'},
    'green'  : {'bg': '\033[42m', 'fg': '\033[32m'},
    'yellow' : {'bg': '\033[43m', 'fg': '\033[33m'},
    'blue'   : {'bg': '\033[44m', 'fg': '\033[34m'},
    'magenta': {'bg': '\033[45m', 'fg': '\033[35m'},
    'cyan'   : {'bg': '\033[46m', 'fg': '\033[36m'},
    'white'  : {'bg': '\033[47m', 'fg': '\033[37m'},
    'reset'  : '\033[0m'
}

def get_transrate_colors():
    return [
        (COLORS['yellow']['bg'], COLORS['yellow']['fg']),
        (COLORS['green']['bg'] , COLORS['green']['fg']),
        (COLORS['red']['bg']   , COLORS['red']['fg'])
    ]

def draw_box(lines, hcolor, nocolor=False):
    width  = 78
    reset  = '' if nocolor else COLORS['reset']
    hcolor = '' if nocolor else hcolor
    top    = f'{hcolor}╭{"─" * width}╮{reset}'
    bottom = f'{hcolor}╰{"─" * width}╯{reset}'
    print(top)
    for line in lines:
        print(f'{hcolor}│ {reset}{line:<{width-1}}{hcolor}│{reset}')
    print(bottom)

class Logo:
    def __init__(self, version='0.1.0', nocolor=False):
        self.version = version
        self.nocolor = nocolor
        if nocolor:
            self.colors = {
                'red'    : '',
                'green'  : '',
                'yellow' : '',
                'blue'   : '',
                'magenta': '',
                'cyan'   : '',
                'white'  : '',
                'end'    : ''
            }
        else:
            self.colors = {
                'red'    : COLORS['red']['fg'],
                'green'  : COLORS['green']['fg'],
                'yellow' : COLORS['yellow']['fg'],
                'blue'   : COLORS['blue']['fg'],
                'magenta': COLORS['magenta']['fg'],
                'cyan'   : COLORS['cyan']['fg'],
                'white'  : COLORS['white']['fg'],
                'end'    : COLORS['reset']
            }
        self.trans  = [
                    '██████┐██████┐  █████┐ ███┐  ██┐ ██████┐',
                    '└─██┌─┘██┌──██┐██┌──██┐████┐ ██│██┌────┘',
                    '  ▓▓│  ▓▓▓▓▓▓┌┘▓▓▓▓▓▓▓│▓▓┌▓▓┐▓▓│└▓▓▓▓▓┐ ',
                    '  ▒▒│  ▒▒┌──▒▒┐▒▒┌──▒▒│▒▒│└▒▒▒▒│ └───▒▒┐',
                    '  ░░│  ░░│  ░░│░░│  ░░│░░│ └░░░│░░░░░░┌┘',
                    '  └─┘  └─┘  └─┘└─┘  └─┘└─┘  └──┘└─────┘ '
                    ]
        self.rate  =  [
                    '██████┐  █████┐ ██████┐███████┐',
                    '██┌──██┐██┌──██┐└─██┌─┘██┌────┘',
                    '▓▓▓▓▓▓┌┘▓▓▓▓▓▓▓│  ▓▓│  ▓▓▓▓▓┐  ',
                    '▒▒┌──▒▒┐▒▒┌──▒▒│  ▒▒│  ▒▒┌──┘  ',
                    '░░│  ░░│░░│  ░░│  ░░│  ░░░░░░░┐',
                    '└─┘  └─┘└─┘  └─┘  └─┘  └──────┘'
                    ]
        self.two   =  [
                    '██████┐ ',
                    '└────██┐',
                    '  ▓▓▓┌─┘',
                    '▒▒┌──┘  ',
                    '░░░░░░░┐',
                    '└──────┘'
                    ]
        self.border      = '└┘┌┐─│'
        self.shading     = '░▓◠'
        self.info        = ['             Quality analysis for de-novo transcriptome assemblies             ',
                            '             ░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ ░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ ░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░      EC Bretz',
                           f'{self.version:>80}']

        self.top_border  = '┌' + '─' * 78 + '┐'
        self.side_border = '│'
        self.bot_border  = '└' + '─' * 78 + '┘'

    def print_bottom(self):
        shading_chars = self.shading
        white         = self.colors['white']
        green         = self.colors['green']
        yellow        = self.colors['yellow']
        red           = self.colors['red']
        end           = self.colors['end']
        
        for line_idx, line in enumerate(self.info):
            if line_idx == 1:
                sections = []
                start    = 0
                while start < len(line):
                    start_pos = line.find('░', start)
                    if start_pos == -1:
                        break
                    end_pos = line.find('░', start_pos + 1)
                    if end_pos == -1:
                        break
                    sections.append((start_pos, end_pos + 1))
                    start = end_pos + 1
                
                section_colors  = [green, yellow, red]
                line_parts      = []
                current_color   = None
                current_segment = []
                
                for i, char in enumerate(line):
                    char_color = white 
                    
                    if char in shading_chars:
                        for section_idx, (section_start, section_end) in enumerate(sections):
                            if section_start <= i < section_end:
                                char_color = section_colors[section_idx]
                                break
                    
                    if char_color != current_color:
                        if current_segment:
                            line_parts.append(f"{current_color}{''.join(current_segment)}{end}")
                            current_segment = []
                        current_color = char_color
                    
                    current_segment.append(char)
                
                if current_segment:
                    line_parts.append(f"{current_color}{''.join(current_segment)}{end}")
                
                print(''.join(line_parts))
            else: 
                print(f"{white}{line}{end}")

    def print_top(self):
        border_chars = self.border
        white        = self.colors['white']
        green        = self.colors['green']
        yellow       = self.colors['yellow']
        red          = self.colors['red']
        end          = self.colors['end']
        
        for i in range(6):
            line_parts = []
            
            for section, color in [(self.trans[i], green), (self.rate[i], yellow), (self.two[i], red)]:
                current_color = None
                current_segment = []
                
                for char in section:
                    char_color = white if char in border_chars else color
                    
                    if char_color != current_color:
                        if current_segment:
                            line_parts.append(f"{current_color}{''.join(current_segment)}{end}")
                            current_segment = []
                        current_color = char_color
                    
                    current_segment.append(char)
                
                if current_segment:
                    line_parts.append(f"{current_color}{''.join(current_segment)}{end}")
            
            print(''.join(line_parts))

    def print_logo(self):
        self.print_top()
        self.print_bottom()
        return

def print_logo(version, quiet=False, nocolor=False):
    if not quiet:
        version_str = f'{version}'
        
        logo = Logo(version_str, nocolor=nocolor)
        logo.print_logo()
        
        yellow, green, red = get_transrate_colors()
        
        reset = '' if nocolor else COLORS['reset']
        if not nocolor:
            title_line = f'{yellow[0]}   {reset}{green[0]}  {reset}{red[0]}{version_str.center(70)}{reset}{green[0]}  {reset}{yellow[0]}   {reset}'
            # print(title_line)
    
    if nocolor:
        return '', ''
    else:
        yellow, green, red = get_transrate_colors()
        return yellow[1], yellow[0]

def print_help(hcolor, defaults, quiet=False, nocolor=False):
    if quiet:
        return
    help_lines = [
        'BASIC:',
        f'--output-dir            -o    DIR     Output directory',
        f'--threads               -t    INT     Number of threads',
        f'--clutter               -c    BOOL    Remove intermediate files',
        '',
        'ASSEMBLY:',
        f'--assembly              -a    FILE    Assembly file(s), comma separated',
        f'--left                  -l    FILE    Left reads file',
        f'--right                 -r    FILE    Right reads file',
        f'--reference             -f    FILE    Reference file',
        f'--bam                   -x    FILE    BAM file',
        f'',
        f'ALIGNER:',
        f'--bowtie2               -b    BOOL    Use Bowtie2 aligner',
        f'--hisat2                -s    BOOL    Use Hisat2 aligner',
        f'',
        f'OPTIONS:',
        f'--quiet                 -q    BOOL    Suppress terminal output',
        f'--debug                       BOOL    Enable debug mode',
        f'--nocolor                     BOOL    Disable color in terminal output',
        f'',
        f'CONFIG:',
        f'--config                      PATH    Path to configuration file',
        f'--config-create [NAME]        BOOL    Create configuration template',
        f'--config-save [NAME]          BOOL    Save current args to config file',
        f'',
        f'STANDARD:',
        f'--version               -v    BOOL    Print version',
        f'--help                  -h    BOOL    Show this help message',
    ]
    draw_box(help_lines, hcolor, nocolor=nocolor)

def print_args(args, passed_args):
    if getattr(args, 'quiet', False):
        return
    nocolor = getattr(args, 'nocolor', False)
    hcolor  = '' if nocolor else getattr(args, 'highlight_color', COLORS['blue']['fg'])
    
    assembly_display = 'None'
    if args.assembly:
        assembly_files = [f.strip() for f in args.assembly.split(',')]
        if len(assembly_files) > 1:
            assembly_display = assembly_files
        else:
            assembly_display = args.assembly
    
    is_single_end = (bool(args.left) != bool(args.right))
    
    assembly_section = {
        'assembly'                  : [f'Assembly file(s):', assembly_display],
        'reference'                 : [f'Reference:', args.reference if args.reference else 'None'],
        'bam'                       : [f'BAM file:', args.bam if args.bam else 'None'],
    }
    
    if is_single_end:
        reads_file                  = args.left if args.left else args.right
        reads_key                   = 'left' if args.left else 'right'
        assembly_section[reads_key] = [f'Reads:', reads_file if reads_file else 'None']
    else:
        assembly_section['left']    = [f'Left reads:', args.left if args.left else 'None']
        assembly_section['right']   = [f'Right reads:', args.right if args.right else 'None']
    
    arg_mappings = {
        'BASIC:': {
            'output_dir'                : [f'Output directory:', args.output_dir if args.output_dir else 'Same as input'],
            'threads'                   : [f'Threads:', args.threads],
            'clutter'                   : [f'Remove intermediate files:', args.clutter],
        },
        'ASSEMBLY:':                       assembly_section,
        'ALIGNER:': {
            'aligner'                   : [f'Aligner:', getattr(args, 'aligner', 'bowtie2').title()],
        },
        'OPTIONS:': {
            'quiet'                     : [f'Quiet mode:', args.quiet],
            'debug'                     : [f'Debug mode:', args.debug],
            'log'                       : [f'Log level:', args.log],
        }
    }
    
    def truncate_value(value, max_len=41):
        value_str = str(value)
        if len(value_str) > max_len:
            return '...' + value_str[-(max_len-3):]
        return value_str
    
    args_list = []
    for stage, arg_dict in arg_mappings.items():
        stage_args = []
        for k, v in arg_dict.items():
            if k in passed_args:
                if isinstance(v[1], list):
                    truncated_first = truncate_value(v[1][0])
                    stage_args.append(f'  {v[0]:<33} {truncated_first}')
                    for item in v[1][1:]:
                        truncated_item = truncate_value(item)
                        stage_args.append(f'  {"":>33} {truncated_item}')
                else:
                    truncated_val = truncate_value(v[1])
                    stage_args.append(f'  {v[0]:<33} {truncated_val}')
        
        if stage_args:
            args_list.append(stage)
            args_list.extend(stage_args)
    
    if args_list:
        draw_box(args_list, hcolor, nocolor=nocolor)
    else:
        single_line = ['Using all default arguments']
        draw_box(single_line, hcolor, nocolor=nocolor)

if __name__ == "__main__":
    logo = Logo(version='0.1.0')
    logo.print_logo()
    logo.print_args({'test': 'value', 'another_arg': 123, 'flag': True})
