import os
import subprocess
import sys

class Samtools:
    def __init__(self, main):
        return
    
    def index(self, main):
        samtools_index_cmd = ['samtools', 'index', '-b', f'-@{main.threads}', main.sortedBam]
        samtools_index_run = subprocess.Popen(samtools_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = samtools_index_run.communicate()
        returncode = samtools_index_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        return returncode
    
    def sort(self, main):
        samtools_sort_cmd = ['samtools', 'sort', f'-@{main.threads}', main.salmonBam, '-o', main.sortedBam]
        samtools_sort_run = subprocess.Popen(samtools_sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = samtools_sort_run.communicate()
        returncode = samtools_sort_run.returncode
        if returncode != 0:
            print(f"\033[91m\nError: {stderr.decode('utf-8')}\033[0m")
            sys.exit(1)
        return returncode