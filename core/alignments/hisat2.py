import os
import subprocess
import sys
class Hisat2:
    def __init__(self, main):
        return
    
    def index(self, main):
        hisat2_index_cmd = ['hisat2-build', main.assembly, os.path.join(main.alignerIndex, main.assemblyName)]
        hisat2_index_run = subprocess.Popen(hisat2_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = hisat2_index_run.communicate()
        returncode = hisat2_index_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        return returncode
    
    def align(self, main):
        prefix = os.path.join(main.alignerDir , main.assemblyName + '_')
        index = os.path.join(main.alignerIndex, main.assemblyName)
        hisat2_cmd = ['hisat2', '-p', str(main.threads), '-q', '-x', index]
        hisat2_cmd.extend(['-U', main.single, '-S', prefix + 'aligned.sam'] if main.mode == 1 else ['-1', main.left, '-2', main.right, '-S', prefix + 'aligned.sam'])
        hisat2_run = subprocess.Popen(hisat2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = hisat2_run.communicate()
        returncode = hisat2_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        return returncode