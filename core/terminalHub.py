import pandas as pd
import sys
import os
import time

class TerminalMain:
    def __init__(self, debug, mode, aligner):
        self.debug       = debug
        self.mode        = mode
        self.tools       = [aligner, 'salmon', 'samtools'] if self.mode > 0 else ['assembly', 'csv']
        self.toolTimes   = []
        self.toolLine    = 0
        self.toolPrinted = False
        self.toolTime    = 0

    def main(self):
        pass

    def runInfo(self, mode, assembly, reference, left, right, bam, aligner, output, threads, clutter, debug):
        self.args = {
            'mode'      : 'Paired' if mode == 2 else 'Single' if mode == 1 else 'Assembly',
            'assembly'  : assembly,
            'reference' : reference,
            'left'      : left,
            'right'     : right,
            'bam'       : bam,
            'aligner'   : aligner,
            'output'    : output,
            'threads'   : threads,
            'clutter'   : clutter,
            'debug'     : debug
        }

        infoPaths = ['assembly', 'left', 'right', 'bam', 'output']
        runBanner = f'\033[7m{" " * 38}Info{" " * 38}\033[0m'
        print(runBanner)

        for path in infoPaths:
            if self.args.get(path):
                print(f"    {path.capitalize():<10} {self.args[path][-60:]}")
        remaining_keys = [
            key for key in self.args
            if key not in infoPaths and self.args[key] not in [False, None, ''] and not pd.isna(self.args[key])
        ]
        num_full_rows = len(remaining_keys) // 4
        remainder = len(remaining_keys) % 4
        underline = lambda text: f"\033[4m{text:<19}\033[0m"

        for i in range(num_full_rows):
            keys_with_values = [
                key for key in remaining_keys[i*4:(i+1)*4]
                if self.args[key] not in [False, None, ''] and not pd.isna(self.args[key])
            ]
            key_row = " ".join(f"{underline(key.capitalize()):<20}" for key in keys_with_values)
            print(key_row)

            value_row = " ".join(f"{str(self.args[key]).capitalize()[:19]:<19}" for key in keys_with_values)
            print(value_row)

        if remainder:
            keys_with_values = [
                key for key in remaining_keys[-remainder:]
                if self.args[key] not in [False, None, ''] and not pd.isna(self.args[key])
            ]
            key_row = " ".join(f"{underline(key.capitalize()):<20}" for key in keys_with_values)
            print(key_row)

            value_row = " ".join(f"{str(self.args[key]).capitalize()[:19]:<19}" for key in keys_with_values)
            print(value_row)
                  
    def processBanner(self):
        print(f'\033[7m{" " * 33}Preprocessing {" " * 33}\033[0m')
        
    def assemblyBanner(self):
        print(f'\033[7m{" " * 36}Assembly{" " * 36}\033[0m')

    def contigBanner(self):
        self.toolPrint('samtools')
        print(f'\n\033[7m{" " * 37}Contig{" " * 37}\033[0m')

    def toolPrint(self, tool):
        max_width = 80
        num_tools = len(self.tools)
        if not self.toolPrinted:
            if num_tools > 0:
                column_width = max_width  // num_tools
                tool_row     = ""

                for t in self.tools:
                    t = t.capitalize() if t != 'csv' else 'CSV'
                    tool_row += f"{t:^{column_width}}"
                print(tool_row)
                self.toolPrinted = True
                self.toolTime    = time.perf_counter()
        else:
            self.toolTimes.append(time.perf_counter() - self.toolTime)
            self.toolTime = time.perf_counter()
            column_width  = max_width // num_tools
            time_row = ""
            for t in self.toolTimes:
                time_str = f"{t:.2f}s"
                time_row += f"{time_str:^{column_width}}"
            print(time_row, end='\r')


    def assemblyPrint(self, assemblyDct):
        self.assemblyDct = pd.read_csv(assemblyDct)
        assemblyStats = f'\033[7m{" " * 33}Assembly Stats{" " * 33}\033[0m'
        print(assemblyStats)
        keys = self.assemblyDct.columns
        values = self.assemblyDct.iloc[0]

        num_full_rows = len(keys) // 4
        remainder     = len(keys) % 4
        underline     = lambda text: f"\033[4m{text:<19}\033[0m"

        for i in range(num_full_rows):
            key_row = " ".join(f"{underline(key):<20}" for key in keys[i*4:(i+1)*4])
            print(key_row)

            value_row = " ".join(f"{str(value)[-19:]:<19}" for value in values[i*4:(i+1)*4])
            print(value_row)
            print()

        if remainder:
            key_row = " ".join(f"{underline(key):<20}" for key in keys[-remainder:])
            print(key_row)

            value_row = " ".join(f"{str(value)[-19:]:<19}" for value in values[-remainder:])
            print(value_row)

class QuietMain:
    def __init__(self, debug, assembly):
        pass

    def runInfo(self, mode, assembly, reference, left, right, bam, aligner, output, threads, clutter, debug):
        pass

    def toolBanner(self):
        pass

    def toolPrint(self, tool):
        pass

    def assemblyPrint(self):
        pass

