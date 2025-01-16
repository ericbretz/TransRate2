import os
import subprocess
import sys

class Salmon:
    def __init__(self, main):
        return
    
    def quant(self, main):
        salmon_cmd = ['salmon', 'quant', '--libType', 'A', '--alignments', main.alignerBam, '--targets', main.assembly, '--noErrorModel', f'--threads={main.threads}', '--sampleOut', '--sampleUnaligned', '--output', main.salmonDir]
        salmon_run = subprocess.Popen(salmon_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = salmon_run.communicate()
        returncode = salmon_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        if os.path.exists(os.path.join(main.salmonDir, 'postSample.bam')):
            os.rename(os.path.join(main.salmonDir, 'postSample.bam'), main.salmonBam)
        return returncode