import os
import sys
import shutil
from core.utils.printout import PrintOut

class Deps:
    def __init__(self, log_level, hcolor=None, bcolor=None, nocolor=False):
        self.hcolor     = '' if nocolor else (hcolor or '\033[34m')
        self.bcolor     = '' if nocolor else (bcolor or '\033[44m')
        self.printClass = PrintOut(log_level, self.hcolor, self.bcolor)
        if nocolor:
            self.printClass.set_nocolor(True)
        self.printout   = self.printClass.printout

        self.deps = {
            'bowtie2'   : 'bowtie2',
            'hisat2'    : 'hisat2',
            'samtools'  : 'samtools',
            'salmon'    : 'salmon',
        }

        self.optional_deps = {
            'blast'     : 'blastn',
            'diamond'   : 'diamond',
        }

    def check_deps(self, check_optional=False, quiet=False):
        if quiet:
            self.printClass.set_quiet(True)
        self.printout('subtitle', 'Checking dependencies')
        missing_deps = []
        found_deps   = {}
        
        for dep, cmd in self.deps.items():
            if not shutil.which(cmd):
                missing_deps.append(dep)
            else:
                found_deps[dep] = 'found'
        
        if check_optional:
            for dep, cmd in self.optional_deps.items():
                if not shutil.which(cmd):
                    found_deps[dep] = 'missing (optional)'
                else:
                    found_deps[dep] = 'found'
        
        if found_deps:
            self.printout('metric', found_deps)
            print()
        
        if missing_deps:
            self.printout('error', f'Missing required dependencies: {", ".join(missing_deps)}')
            self.printout('error', 'Please install missing dependencies before running TransRate2')
            sys.exit(1)
        # else:
        #     self.printout('success', 'All required dependencies found')

    def check_python_packages(self):
        required_packages = [
            'numpy',
            'pandas', 
            'pysam',
            'biopython',
            'pyyaml'
        ]
        
        missing_packages = []
        found_packages = {}
        
        for package in required_packages:
            try:
                __import__(package)
                found_packages[package] = 'found'
            except ImportError:
                missing_packages.append(package)
                found_packages[package] = 'missing'
        
        if found_packages:
            self.printout('subtitle', 'Python Package Dependencies')
            self.printout('metric', found_packages)
        
        if missing_packages:
            self.printout('error', f'Missing Python packages: {", ".join(missing_packages)}')
            self.printout('error', 'Install with: pip install ' + ' '.join(missing_packages))
            return False
        return True

if __name__ == '__main__':
    deps = Deps(3, '\033[94m', '\033[44m')
    deps.check_deps(check_optional=True)
    deps.check_python_packages()
