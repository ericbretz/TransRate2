import os
import sys
import subprocess

class Dep:
    def __init__(self):
        self.packages = ['pandas', 'pysam', 'numpy', 'argparse', 'psutil']

    def depCheck(self):
        def tool_installed(tool):
            check = subprocess.Popen(['which', tool], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout,stderr = check.communicate()
            return len(stdout.decode('utf-8')) > 0
        
        for package in self.packages:
            try:
                __import__(package)
            except ImportError:
                print("Checking Dependencies...", end="\r")
                if os.system('pip install ' + package + ' > /dev/null 2>&1') != 0:
                    sys.exit("Error: " + package + " not installed. Exiting.")

        missing = []
        tools          = {
            'salmon'      : 'https://github.com/COMBINE-lab/salmon',
            'samtools'    : 'https://github.com/samtools/samtools',
            'diamond'     : 'https://github.com/bbuchfink/diamond',
            'hisat2'      : 'https://github.com/DaehwanKimLab/hisat2',
            'bowtie2'     : 'https://github.com/BenLangmead/bowtie2'
            }

        for t,l in tools.items():
            if not tool_installed(t):
                missing.append((t,l))

        return missing