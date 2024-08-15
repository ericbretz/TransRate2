import os
import subprocess
import sys
import time
import numpy as np
import pandas as pd

class Tr2:
    def __init__(self, main):
        self.TSTART = False
    def run(self, main):
        functions = main.FUNCTIONS

        for func in functions: # Run TransRate2 functions
            if func == 'csv':
                func = type('', (), {})()
                func.__name__ = 'csv'

            if main.READMODE == 1: # if Single-end, do not run good or segmentation analyses
                if func.__name__ in ['good', 'sgmt']:
                    continue

            main.TSTAGE = func.__name__ # Set the current TransRate2 stage

            while not main.TPRINTED and not main.QUIET: # Wait for terminal to print previous stage
                pass

            main.TIMES[func.__name__] = time.perf_counter() # Start timer
            
            # print(f'\n\n{func.__name__}\n\n')
            if func.__name__ == 'csv':
                self.csv_run(main)
            else:
                threads = main.THREADS if not main.FIX else 1
                main.TR_DCT = func(main.BAM_SORTED, main.TR_DCT, threads, main.READMODE)

            main.TIMES[func.__name__] = time.perf_counter() - main.TIMES[func.__name__] # Stop timer
            main.TPRINTED = False # Wait for terminal to print current stage

        main.TIMES[main.STAGE] = time.perf_counter() - main.TIMES[main.STAGE]
        main.TSTOP = True # Tell output that TransRate has finished
        
        while not main.TEND:
            pass

        return
    
    def csv_run(self, main):
        if main.READMODE != 0:
            columns = ['name', 'p_seqtrue', 'bridges', 'length', 'fragments',
                        'both_mapped', 'properpair', 'good', 'basesuncovered', 'p_notsegmented'] if main.READMODE == 2 else \
                        ['name', 'p_seqtrue', 'length', 'fragments', 'basesuncovered']
            csvout = pd.DataFrame(columns=columns)

            for key, value in main.TR_DCT.items():
                csvout.loc[len(csvout)] = value['stats']

            csvout.to_csv(main.TR_CSV, index=False)


