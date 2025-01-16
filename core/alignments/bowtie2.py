import os
import subprocess
import sys

class Bowtie2:
    def __init__(self, main):
        return
    
    def index(self, main):
        bt2_index_cmd = ['bowtie2-build', '--threads', str(main.threads), main.assembly, main.alignerIndex]
        bt2_index_run = subprocess.Popen(bt2_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = bt2_index_run.communicate()
        returncode = bt2_index_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        return returncode
    
    def align(self, main):
        bt2_cmd = ['bowtie2', '--threads', f'{main.threads}', '--very-sensitive', '--phred33', '--no-unal', '--no-mixed', '--no-discordant', '--rdg', '1000,1000', '--rfg', '1000,1000', '-x', main.alignerIndex]
        bt2_cmd.extend(['-U', main.single, '-S', main.alignerBam] if main.mode == 1 else ['-1', main.left, '-2', main.right, '-S', main.alignerBam])
        bt2_run = subprocess.Popen(bt2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = bt2_run.communicate()
        returncode = bt2_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        return returncode