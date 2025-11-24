import os
import sys
import time
from typing import Any, Dict, List, Union

class PrintOut:
    def __init__(self, level, hc=None, bc=None):
        self.level            = level
        self.hc               = hc or '\033[34m'
        self.bc               = bc or '\033[44m'
        self.quiet            = False
        self.nocolor          = False
        self.current_stage    = None
        self.stage_start_time = None
        self.stage_metrics    = {}

        self._init_styles()

    def _init_styles(self):
        if self.nocolor:
            self.styles = {
                'title': {
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
                'subtitle': {
                    'color'   : '\033[7m',
                    'width'   : 80,
                    'char'    : ' '
                },
                'metric': {
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
                'info': {
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
                'progress': {
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
                'error': {
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
                'debug': {
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
                'success': {    
                    'color'   : '',
                    'width'   : 80,
                    'char'    : ' '
                },
            }
        else:
            self.styles = {
            'title': {
                'color'   : self.bc,
                'width'   : 80,
                'char'    : ' '
            },
            'subtitle': {
                'color'   : self.bc,
                'width'   : 80,
                'char'    : ' '
            },
            'metric': {
                'color'   : self.bc,
                'width'   : 80,
                'char'    : ' '
            },
            'info': {
                'color'   : '\033[43m',
                'width'   : 80,
                'char'    : ' '
            },
            'progress': {
                'color'   : '\033[47m',
                'width'   : 80,
                'char'    : ' '
            },
            'error': {
                'color'   : '\033[41m',
                'width'   : 80,
                'char'    : ' '
            },
            'debug': {
                'color'   : '\033[44m',
                'width'   : 80,
                'char'    : ' '
            },
            'success': {    
                'color'   : '\033[42m',
                'width'   : 80,
                'char'    : ' '
            },
        }
 
        self.key_translate = {
        }
        self.map_assembly = {
            'assembly'          : ('Assembly', str),
            'nSeqs'             : ('# Seqs', int),
            'bases'             : ('Bases', int),
            'smallest'          : ('Smallest', int),
            'largest'           : ('Largest', int),
            'meanLength'        : ('Mean Length', int),
            'medianLength'      : ('Median Length', int),
            'stdLength'         : ('Std Length', int),
            'nUnder200'         : ('# Under 200', int),
            'nOver1k'           : ('# Over 1k', int),
            'nOver10k'          : ('# Over 10k', int),
            'nWithOrf'          : ('# With Orf', int),
            'meanOrfPercent'    : ('Mean Orf Percent', float),
            'n90'               : ('n90', int),
            'n70'               : ('n70', int),
            'n50'               : ('n50', int),
            'n30'               : ('n30', int),
            'n10'               : ('n10', int),
            'gcCount'           : ('GC Count', int),
            'pGC'               : ('p GC', float),
            'basesN'            : ('# Bases N', int),
            'pN'                : ('p N', float),
            'fragments'         : ('# Fragments', int),
            'fragmentsMapped'   : ('# Fragments Mapped', int),
            'bothMapped'        : ('# Both Mapped', int),
            'pFragmentsMapped'  : ('p Fragments Mapped', float),
            'softclipped'       : ('# Softclipped', int),
            'pSoftclipped'      : ('p Softclipped', float),
            'goodMappings'      : ('# Good Mappings', int),
            'pGoodMappings'     : ('p Good Mappings', float),
            'badMappings'       : ('# Bad Mappings', int),
            'potentialBridges'  : ('# Potential Bridges', int),
            'basesUncovered'    : ('# Bases Uncovered', int),
            'pBasesUncovered'   : ('p Bases Uncovered', float),
            'contigsUncovBase'  : ('# Contigs Bases Uncovered', int),
            'pContigsUncovbase' : ('p Contigs Bases Uncovered', float),
            'contigsUncovered'  : ('# Contigs Uncovered', int),
            'pContigsUncovered' : ('p Contigs Uncovered', float),
            'contigsLowcovered' : ('# Contigs Lowcovered', int),
            'pContigsLowcovered': ('p Contigs Lowcovered', float),
            'contigsSegmented'  : ('# Contigs Segmented', int),
            'pContigsSegmented' : ('p Contigs Segmented', float),
            'goodContigs'       : ('# Good Contigs', int),
            'badContigs'        : ('# Bad Contigs', int),
            'cutoff'            : ('Cutoff', float),
            'weighted'          : ('Weighted', float),
            'optimalScore'      : ('Optimal Score', float),
            'score'             : ('Score', float),

            # Reference statistics
            'CRBBhits'          : ('CRBB Hits', int),
            'nContigsWithCRBB'  : ('# Contigs w/ CRBB', int),
            'pContigsWithCRBB'  : ('p Contigs w/ CRBB', float),
            'nRefsWithCRBB'     : ('# Refs w/ CRBB', int),
            'pRefsWithCRBB'     : ('% Refs w/ CRBB', float),
            'rbhPerReference'   : ('RBH per Reference', float),
            'cov25'             : ('Coverage ≥25%', int),
            'pCov25'            : ('p Coverage ≥25%', float),
            'cov50'             : ('Coverage ≥50%', int),
            'pCov50'            : ('p Coverage ≥50%', float),
            'cov75'             : ('Coverage ≥75%', int),
            'pCov75'            : ('p Coverage ≥75%', float),
            'cov85'             : ('Coverage ≥85%', int),
            'pCov85'            : ('p Coverage ≥85%', float),
            'cov95'             : ('Coverage ≥95%', int),
            'pCov95'            : ('p Coverage ≥95%', float),
            'referenceCoverage' : ('Reference Coverage', float),
            }

        self.map_contig = {
        #     'name'              : 'Name',
        #     'length'            : 'Length',
        #     'fragments'         : '# Fragments',
        #     'gcCount'           : 'GC Count',
        #     'pGC'               : 'p GC',
        #     'basesUncovered'    : '# Bases Uncovered',
        #     'pBasesCovered'     : 'p Bases Covered',
        #     'bridges'           : '# Bridges',
        #     'bothMapped'        : '# Both Mapped',
        #     'properPair'        : '# Proper Pair',
        #     'good'              : '# Good',
        #     'pGood'             : 'p Good',
        #     'orfLength'         : 'Orf Length',
        #     'pNotSegmented'     : 'p Not Segmented',
        #     'pSeqTrue'          : 'p Seq True',
        #     'softclipped'       : '# Softclipped',
        #     'pSoftclipped'      : 'p Softclipped',
        #     'effLength'         : 'Eff Length',
        #     'effCount'          : '# Eff Count',
        #     'tpm'               : 'TPM',
        #     'coverage'          : 'Coverage',
        #     'sCnuc'             : 'sCnuc',
        #     'sCcov'             : 'sCcov',
        #     'sCord'             : 'sCord',
        #     'sCseg'             : 'sCseg',
        #     'score'             : 'Score',
        }

    def fmt_key(self, key: str) -> str: 
        if key in self.key_translate:
            return self.key_translate[key]
        else:
            return key
    
    def fmt_lst_cnt(self, lst: List[Any]) -> str:
        return f"{len(lst)}"

    def fmt_dict(self, dct: Dict[Any, Any]) -> List[tuple]:
        dct_list = []
        for k, v in dct.items():
            if isinstance(v, list):
                dct_list.append((k, self.fmt_lst_cnt(v)))
            elif isinstance(v, dict):
                nested_items = self.fmt_dict(v)
                for nk, nv in nested_items:
                    dct_list.append((nk, nv))
            else:
                dct_list.append((k, str(v)))
        return dct_list
    
    def fmt_str(self, string: str, value: bool = False) -> str:
        if value is None:
            return string
        elif len(string) >= 20 and not value:
            return '...' + string[-17:]
        elif len(string) >= 56 and value:
            return '...' + string[-53:]
        else:
            return string

    def check_type(self, content: Any) -> None:
        if isinstance(content, dict):
            content = self.fmt_dict(content)
        elif isinstance(content, list):
            if content and isinstance(content[0], tuple):
                content = content
            else:
                content = self.fmt_lst_cnt(content)
        elif isinstance(content, str):
            content = content
        elif isinstance(content, int):
            content = str(content)
        else:
            raise ValueError(f"Invalid content type: {type(content)}")
        return content
    
    def check_style(self, style: str, content: Any) -> None:
        if style == 'title':
            self.p_title(content)
        elif style == 'subtitle':
            self.p_subtitle(content)
        elif style == 'metric':
            self.p_metric(content)
        elif style == 'info':
            self.p_info(content)
        elif style == 'error':
            self.p_error(content)
        elif style == 'debug':
            self.p_debug(content)
        elif style == 'success':
            self.p_success(content)
        elif style == 'progress':
            self.p_progress(content)
        elif style == 'final':
            self.p_final(content)
        else:
            raise ValueError(f"Invalid style: {style}")
        return content

    def p_title(self, title: str,) -> None:
        color  = self.styles['title']['color']
        width  = self.styles['title']['width']
        char   = self.styles['title']['char']
        string = f'{title:^{width}}'
        reset  = '\033[0m'
        print(f"{color}{string}{reset}")

    def p_subtitle(self, subtitle: str) -> None:
        color  = self.styles['subtitle']['color']
        width  = self.styles['subtitle']['width']
        char   = self.styles['subtitle']['char']
        string = f'{subtitle:^{width}}'
        reset  = '\033[0m'
        print(f"{color}{string}{reset}")

    def p_metric(self, metric: str) -> None:
        color    = self.styles['metric']['color']
        width    = self.styles['metric']['width']
        char     = self.styles['metric']['char']
        if isinstance(metric, list):
            for line in metric:
                if isinstance(line, tuple) and len(line) == 2:
                    key, value     = line
                    metric_str     = self.fmt_str(f'    {self.fmt_key(key)}', value=False)
                    metric_str     = f'{metric_str:<20}│'
                    value_str      = self.fmt_str(str(value), value=True)
                    value_str      = f'{value_str:<{width - 21}}'
                    reset          = '\033[0m'
                    print(f"{metric_str}{reset} {value_str}")
                else:
                    reset =  '\033[0m'
                    print(f"{color}{str(line):^{width}}{reset}")
        elif isinstance(metric, dict):
            for key, value in metric.items():
                if key not in self.key_translate:
                    continue
                metric_str = self.fmt_str(f'    {self.fmt_key(key)}', value=False)
                metric_str = f'{metric_str:<20}│'
                value_str  = self.fmt_str(str(value), value=True)
                value_str  = f'{value_str:<{width - 21}}'
                reset      = '\033[0m'
                print(f"{metric_str}{reset} {value_str}")
        else:
            metric_str = self.fmt_str(self.fmt_key(str(metric)), value=None)
            white_bg   = '' if self.nocolor else '\033[47m'
            reset      = '\033[0m'
            print(f"{white_bg}{metric_str:^{width}}{reset}")

    def p_info(self, info: str) -> None:
        color  = self.styles['info']['color']
        width  = self.styles['info']['width']
        char   = self.styles['info']['char']
        title  = f'{"INFO":<20}'
        string = f'{info:<{width - 20}}'
        reset  = '\033[0m'
        print(f"{color}{title}{reset} {string}")

    def p_error(self, error: str) -> None:
        color  = self.styles['error']['color']
        width  = self.styles['error']['width']
        char   = self.styles['error']['char']
        title  = f'{"ERROR":<20}'
        string = f'{error:<{width - 20}}'
        reset  = '\033[0m'
        print(f"{color}{title}{reset} {string}")

    def p_debug(self, debug: str) -> None:
        color  = self.styles['debug']['color']
        width  = self.styles['debug']['width']
        char   = self.styles['debug']['char']
        title  = f'{"DEBUG":<20}'
        string = f'{debug:<{width - 20}}'
        reset  = '\033[0m'
        print(f"{color}{title}{reset} {string}")

    def p_success(self, success: str) -> None:
        color  = self.styles['success']['color']
        width  = self.styles['success']['width']
        char   = self.styles['success']['char']
        title  = f'{"SUCCESS":<20}'
        string = f'{success:<{width - 20}}'
        reset  = '\033[0m'
        print(f"{color}{title}{reset} {string}")

    def p_progress(self, progress: str) -> None:
        color    = self.styles['progress']['color']
        width    = self.styles['progress']['width']
        white_bg = '' if self.nocolor else '\033[47m'
        reset    = '\033[0m'
        print(f"{white_bg}{progress:^{width}}{reset}", end='\r', flush=True)
    
    def start_progress_stage(self, stage_name: str) -> None:
        self.current_stage    = stage_name
        self.stage_start_time = time.time()
        self.stage_metrics    = {}
    
    def start_section(self, section_name: str) -> None:
        self.current_section = section_name
        self.section_metrics = {}
    
    def add_stage_to_section(self, stage_metrics: Dict[str, Any]) -> None:
        if hasattr(self, 'section_metrics'):
            self.section_metrics.update(stage_metrics)
    
    def complete_section(self) -> None:
        if hasattr(self, 'section_metrics') and self.section_metrics:
            # print()
            for key, value in self.section_metrics.items():
                if key in self.map_assembly:
                    key_format = self.map_assembly[key][0]
                    key_type   = self.map_assembly[key][1]
                    try:
                        if key_type == int and isinstance(value, (str, float)):
                            converted_value = int(float(value))
                        elif key_type == float and isinstance(value, (str, int)):
                            converted_value = float(value)
                        else:
                            converted_value = key_type(value)
                        formatted_value = self._format_metric_value(converted_value)
                    except (ValueError, TypeError):
                        formatted_value = self._format_metric_value(value)
                elif key in self.map_contig:
                    key_format      = self.map_contig[key]
                    formatted_value = self._format_metric_value(value)
                else:
                    key_format      = key
                    formatted_value = str(value)
                metric_line = f"  {key_format:<25} {formatted_value}"
                print(metric_line[:80])
            print()
    
    def update_progress_bar(self, current: int, total: int, step_name: str = "", metrics: Dict[str, Any] = None) -> None:
        if self.quiet:
            return
            
        if metrics:
            self.stage_metrics.update(metrics)
        
        progress     = (current / total) * 100 if total > 0 else 0
        elapsed      = time.time() - self.stage_start_time if self.stage_start_time else 0
        
        bar_width    = 32
        filled_width = int(bar_width * progress / 100)
        bar          = '█' * filled_width + '░' * (bar_width - filled_width)
        
        progress_str = f"{progress:5.1f}%"
        step_display = step_name[:20] if step_name else "Processing"
        time_str     = f"{elapsed:4.1f}s"
        
        final_line   = f"\r[{bar}] {progress_str} │ {step_display:^20} │ {time_str:>13}"
        
        print(final_line, end='', flush=True)
    
    def complete_progress_stage(self, final_metrics: Dict[str, Any] = None) -> None:
        if self.quiet:
            return
            
        print()
        
        if final_metrics:
            self.stage_metrics.update(final_metrics)
            
        if hasattr(self, 'section_metrics') and self.stage_metrics:
            self.section_metrics.update(self.stage_metrics)
        
        self.current_stage    = None
        self.stage_start_time = None
        self.stage_metrics    = {}
    
    def _format_metric_value(self, value: Any) -> str:
        if isinstance(value, float):
            # if value > 1000000:
            #     return f"{value/1000000:.2f}M"
            # elif value > 1000:
            #     return f"{value/1000:.2f}k"
            # else:
            #     return f"{value:.3f}"
            return f"{value:.2f}"
        elif isinstance(value, int):
            # if value > 1000000:
                # return f"{value//1000000}M"
            # elif value > 1000:
                # return f"{value//1000}k"
            # else:
                # return str(value)
            return str(value)
        else:
            return str(value)[:15]

    def p_final(self, files: dict) -> None:
        def draw_box(lines, hcolor):
            width  = 78
            reset  = '' if self.nocolor else '\033[0m'
            top    = f'{hcolor}╭{"─" * width}╮{reset}'
            bottom = f'{hcolor}╰{"─" * width}╯{reset}'
            print(top)
            for line in lines:
                print(f'{hcolor}│ {reset}{line:<{width-1}}{hcolor}│{reset}')
            print(bottom)
        out_list = []
        for dir, file in files:
            s = f'{dir:<18}│ {self.fmt_str(str(file), value=True)}'
            out_list.append(s)
        draw_box(out_list, self.hc)

    def set_quiet(self, quiet: bool) -> None:
        self.quiet = quiet
    
    def set_nocolor(self, nocolor: bool) -> None:
        self.nocolor = nocolor
        if self.nocolor:
            self.hc = ''
            self.bc = ''
        else:
            self.hc = self.hc or '\033[34m'
            self.bc = self.bc or '\033[44m'
        self._init_styles()
    
    def printout(self, style: str, content: Any) -> None:
        if self.quiet:
            return
        content = self.check_type(content)
        self.check_style(style, content)

if __name__ == "__main__":
    import time
    test_dict = {
        'assembly'    : '/path/to/assembly.fasta',
        'threads'     : 8,
        'contig_count': 15243,
        'n50'         : 1524,
        'total_length': 45234523,
        'gc_content'  : 45.2,
        'mapping_rate': 89.5,
        'coverage'    : 12.4,
    }
    
    printout = PrintOut(level='info', hc='\033[94m', bc='\033[44m')
    
    printout.printout('title',      'TransRate2 Analysis')
    printout.printout('subtitle',   'Assembly Quality Assessment')
    printout.printout('metric',     'Assembly Statistics')
    printout.printout('metric',      test_dict)
    printout.printout('info',       'Analysis completed successfully')
    printout.printout('success',    'All quality checks passed')
    
    for i in range(10):
        printout.printout('progress', f"Processing contigs: {i*10}/100%")
        time.sleep(0.1)
    print()
