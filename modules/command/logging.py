import os
import time

class Logging:
    def __init__(self, main):
        #### Logging ####
        self.LOG_PATH       = None  # Path to log files
        self.LOGDCT         = {}    # Dictionary of process start and end times
        self.STDOUT         = None  # Path to stdout log file
        self.STDERR         = None  # Path to stderr log file
        self.RUNTIME        = time.strftime('%Y%m%d%H%M%S') # Time of run start
        
        #### Log files ####
        self.LOGFILES     = {
            'snap_paired'   : [f'{self.RUNTIME}_snap_stdout.log',             f'{self.RUNTIME}_snap_stderr.log'],
            'snap_index'    : [f'{self.RUNTIME}_snap_index_stdout.log',       f'{self.RUNTIME}_snap_index_stderr.log'],
            'bowtie2_index' : [f'{self.RUNTIME}_bowtie2_index_stdout.log',    f'{self.RUNTIME}_bowtie2_index_stderr.log'],
            'bowtie2'       : [f'{self.RUNTIME}_bowtie2_stdout.log',          f'{self.RUNTIME}_bowtie2_stderr.log'],
            'star_index'    : [f'{self.RUNTIME}_star_index_stdout.log',       f'{self.RUNTIME}_star_index_stderr.log'],
            'star'          : [f'{self.RUNTIME}_star_stdout.log',             f'{self.RUNTIME}_star_stderr.log'],
            'salmon'        : [f'{self.RUNTIME}_salmon_stdout.log',           f'{self.RUNTIME}_salmon_stderr.log'],
            'samtools_sort' : [f'{self.RUNTIME}_samtools_stdout.log',         f'{self.RUNTIME}_samtools_stderr.log'],
            'samtools_index': [f'{self.RUNTIME}_samtools_index_stdout.log',   f'{self.RUNTIME}_samtools_index_stderr.log'],
            'transrate'     : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'assembly'      : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'assembly_solo' : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'reference'     : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'file'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'frag'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'seq'           : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'good'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'base'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'sgmt'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'plots'         : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
        }

    def log_time(self, main, process, type):
        # Create dictionary of process start and end times
        if type == 'start':
            self.LOGDCT[process] = {}
            self.LOGDCT[process]['start'] = time.perf_counter()
        elif type == 'end':
            self.LOGDCT[process]['end'] = time.perf_counter()
        else:
            pass

    def log_set(self, main, process):
        # Set log file paths
        self.STDOUT = os.path.join(main.LOG_PATH, self.LOGFILES[process][0])
        self.STDERR = os.path.join(main.LOG_PATH, self.LOGFILES[process][1])

    def log_write(self, main, stdout, stderr):
        # Write stdout and stderr to log files
        with open(self.STDOUT, 'a') as stdout_f:
            stdout_f.write(stdout.decode('utf-8'))
        with open(self.STDERR, 'a') as stderr_f:
            stderr_f.write(stderr.decode('utf-8'))