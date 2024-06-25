import time
import os
import subprocess
import pysam
import sys
import threading
import pandas as pd
from modules.analyses import file,frag,seq,base,sgmt,good,assembly,reference,plot
from modules.command.logging import Logging
from modules.command.fileset import FileSet
from modules.command.output import Output
from modules.command.tr2 import Tr2
from modules.command.bt2 import BT2
from modules.command.star import STAR
from modules.command.snap import SNAP
from modules.command.salmon import Salmon
from modules.command.samtools import Samtools
from modules.terminal.stats import Stats
import warnings
import shutil

warnings.filterwarnings('ignore')

class MAIN:
    def __init__(self):
        #### Settings ####
        self.CLUTTER        = False # Remove intermediate files
        self.COLOR          = None  # Color for terminal output
        self.FILESET        = FileSet() # Set up file paths and names
        self.QUIET          = False # Suppress all messages
        self.SKIP           = False # Skip to TransRate
        self.THREADS        = 1     # Number of threads to use

        #### Analyses ####
        self.PLOT           = False # Run Plot
        self.READS          = False # Run Reads if not solo
        self.REFERENCE      = False # Run Reference

        #### Input ####
        # Assembly
        self.ASSEMBLY       = None  # Path to assembly .fasta file
        self.ASSEMBLY_COUNT = 0     # Total number of assemblies
        self.ASSEMBLY_MULTI = False # Multiple assemblies bool
        self.ASSEMBLY_NAME  = ''    # Assembly name
        self.ASSEMBLY_ROW   = 0     # Number of rows in assembly .fasta file
        self.ASSEMBLY_RUN   = 0     # Current assembly run
        
        # Reads
        self.LEFT           = None  # Path to left reads .fastq file
        self.READMODE       = None  # 0: Assembly Only, 1: Single-end 2: Paired-end
        self.RIGHT          = None  # Path to right reads .fastq file
        self.SINGLE         = None
        
        # Reference
        self.REFERENCE_FILE = None  # Path to reference .fasta file

        # Output
        self.OUTPUT         = None  # Path to output directory
        self.TRANSRATE_PATH = None  # Path to Transrate output directory

        # Aligner
        self.ALIGNER        = None  # Aligner to use (bowtie2, star, snap)
        self.READ_COUNT     = None  # Read Count
        self.ALIGNER_INDEX  = None  # Path to aligner index
        self.ALIGNER_PATH   = None  # Path to aligner
        
        # Salmon
        self.SALMON_PATH    = None  # Path to salmon
        self.SALMON_QUANT   = None  # Path to salmon quant.sf file
        
        # TransRate
        self.FUNCTIONS      = [file.file,frag.frag,seq.seq,good.good,base.base,sgmt.sgmt,'csv'] # TransRate functions
        self.TR_DCT         = {}    # TransRate dictionary
        self.PLOT           = False # Plot results
        # Logging
        self.LOG            = Logging(self) # Logging class
        self.LOG_PATH       = None  # Path to log files
        self.ERROR          = False # Error flag
        self.STOPPED        = False # Stop flag

        # File Info
        self.ASSEMBLY_FILE  = ''    # Path to assembly .csv file
        self.BAM_ALIGNER    = None  # Path to bam file
        self.BAM_SALMON     = None  # Path to salmon file
        self.BAM_SORTED     = None  # Path to sorted bam file
        self.CONTIG_FILE    = ''    # Path to contig .csv file

        # Terminal Output
        self.TERMINALCLASS  = None  # Terminal class
        self.OUT_THREAD     = None  # Terminal thread
        self.PRINTED        = True
        self.STAGE          = None  # Current stage
        self.TEND           = None  # TransRate2 Final Print
        self.TIMES          = {}    # Time dictionary
        self.TPRINTED       = True  # Finished printing TransRate2 stage
        self.TSTAGE         = None  # Transrate stage
        self.TSTART         = None  # TransRate2 started
        self.TSTOP          = None  # TransRate2 stopped

    def run(self):
        self.FILESET.run(self)      # Set file names based on assembly count
        
        if self.READMODE == 1:
            self.SINGLE = self.LEFT if self.LEFT else self.RIGHT
        
        if not self.QUIET:
            self.OUT_THREAD = threading.Thread(target=self.terminal)
            self.OUT_THREAD.start()

        self.LOG_PATH = os.path.join(self.OUTPUT, 'logs')

        if not os.path.exists(self.LOG_PATH):
            os.makedirs(self.LOG_PATH)

        if not self.SKIP:
            if self.READMODE != 0:
                # Run aligner (STAR, SNAP, Bowtie2)
                self.aligner_run()
                # Run salmon
                self.salmon_run()
                # Run Samtools
                self.samtools_run()
        # Run TransRate2
        self.tr2_run()
        # Run Assembly Stats
        self.assembly_run()
        # Run Plot
        if self.PLOT:
            self.plots_run()

        if self.REFERENCE:
            self.reference_run()

        if self.CLUTTER:
            self.clutter_run()

        self.STAGE = 'Finished'

        while not self.PRINTED:
            pass

        self.OUT_THREAD.join()
        
        self.stats_run()

        return
    
    def aligner_run(self):
        if self.ALIGNER == 'star':
            star = STAR(self)

            self.STAGE   = 'STAR Index'
            self.STARTED = True
            self.TIMES[self.STAGE] = time.perf_counter()

            star.star_index(self)

            self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
            self.PRINTED = False

            self.STAGE   = 'STAR'
            self.STARTED = True
            self.TIMES[self.STAGE] = time.perf_counter()

            star.star(self)

            self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
            self.PRINTED = False

            return
        
        elif self.ALIGNER == 'snap':
            snap = SNAP(self)

            self.STAGE   = 'Snap Index'
            self.STARTED = True
            self.TIMES[self.STAGE] = time.perf_counter()

            snap.snap_index(self)

            self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
            self.PRINTED = False

            self.STAGE   = 'Snap'
            self.STARTED = True
            self.TIMES[self.STAGE] = time.perf_counter()

            snap.snap(self)

            self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
            self.PRINTED = False
            return
        
        elif self.ALIGNER == 'bowtie2':
            bt2 = BT2(self)

            self.STAGE   = 'Bowtie2 Index'
            self.STARTED = True
            self.TIMES[self.STAGE] = time.perf_counter()

            bt2.bowtie2_index(self)

            self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
            self.PRINTED = False

            self.STAGE   = 'Bowtie2'
            self.STARTED = True
            self.TIMES[self.STAGE] = time.perf_counter()

            bt2.bowtie2(self)

            self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
            self.PRINTED = False
            return
        
        else:
            self.OUT_THREAD.join()
            # print('\n')
            # print('Aligner not found')
            # sys.exit(3)
            self.LOG.error_out(self, 'Aligner', 'Aligner not found')

    def salmon_run(self):
        salmon = Salmon(self)

        self.STAGE   = 'Salmon'
        self.TIMES[self.STAGE] = time.perf_counter()

        salmon.salmon(self)

        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

    def samtools_run(self):
        samtools = Samtools(self)

        self.STAGE   = 'Samtools Sort'
        self.TIMES[self.STAGE] = time.perf_counter()

        samtools.samtools_sort(self)

        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

        self.STAGE   = 'Samtools Index'
        self.TIMES[self.STAGE] = time.perf_counter()

        samtools.samtools_index(self)

        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

    def tr2_run(self):
        tr2 = Tr2(self)

        self.STAGE   = 'TransRate2'
        self.TIMES[self.STAGE] = time.perf_counter()

        tr2.run(self)

        self.PRINTED = False

    def assembly_run(self):
        self.STAGE   = 'Assembly Stats'
        self.TIMES[self.STAGE] = time.perf_counter()

        assembly_class = assembly.Assembly(self)
        assembly_class.run()

        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

    def plots_run(self):
        self.STAGE   = 'Plotting Results'
        self.TIMES[self.STAGE] = time.perf_counter()

        plots_class = plot.plots(self)
        plots_class.run(self)

        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

    def reference_run(self):
        self.STAGE   = 'Reference'
        self.TIMES[self.STAGE] = time.perf_counter()

        reference_class = reference.REFERENCE(self)
        reference_class.run()

        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

    def clutter_run(self):
        self.STAGE   = 'Clutter Removal'
        self.TIMES[self.STAGE] = time.perf_counter()

        delfolders = ['star', 'bowtie2', 'snap', 'salmon']
        for root, dirs, files in os.walk(self.OUTPUT):
            for file in files:
                if not file.endswith('.csv') and 'logs' not in root and not file.endswith('.fa') and not file.endswith('.png'):
                    os.remove(os.path.join(root, file))
            for dir in dirs:
                if any(folder in dir for folder in delfolders) or os.path.getsize(os.path.join(root, dir)) == 64:
                    shutil.rmtree(os.path.join(root, dir))


        self.TIMES[self.STAGE] = time.perf_counter() - self.TIMES[self.STAGE]
        self.PRINTED = False

    def stats_run(self):

        stats_class = Stats(self)
        stats_class.stats(self)

    def terminal(self):
        self.TERMINALCLASS  = Output(self)
        self.STARTED = True
        while True:
            self.TERMINALCLASS.run(self)
            time.sleep(0.1) # How often should the terminal update



