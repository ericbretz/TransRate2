from util.deps import Dep
dep = Dep()
dep.depCheck()
import os
import sys
import argparse
import psutil
from core.coreHub import TransRate
from util.help import HelpText

class Transrate2:
    def __init__(self):
        self.version     = '2.6.5'
        self.author      = 'EC. Bretz'
        self.description = 'Quality analysis for de-novo transcriptome assemblies'
        self.t2          = TransRate()

    def parseArgs(self):
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter, description='Transrate2')
        parser.add_argument('--assembly',  '-a', type=str,            help='Assembly file(s), comma seperated')
        parser.add_argument('--left',      '-l', type=str,            help='Left reads file')
        parser.add_argument('--right',     '-r', type=str,            help='Right reads file')
        parser.add_argument('--reference', '-f', type=str,            help=argparse.SUPPRESS)
        parser.add_argument('--output',    '-o', type=str,            help='Output directory')
        parser.add_argument('--threads',   '-t', type=int,            help='Number of threads', default=1)
        parser.add_argument('--bam',       '-x', type=str,            help=argparse.SUPPRESS)
        parser.add_argument('--hisat2',    '-s', action='store_true', help='Use Hisat2 aligner')
        parser.add_argument('--bowtie2',   '-b', action='store_true', help='Use Bowtie2 aligner')
        parser.add_argument('--clutter',   '-c', action='store_true', help='Remove clutter from output directory')
        parser.add_argument('--help',      '-h', action='store_true', help='Display this help message')
        parser.add_argument('--quiet',     '-q', action='store_true', help='Supress terminal output')
        parser.add_argument('--version',   '-v', action='store_true', help=argparse.SUPPRESS)
        parser.add_argument('--debug',     '-d', action='store_true', help=argparse.SUPPRESS)
        self.args = parser.parse_args()

        if self.args.help or len(sys.argv) < 2:
            helpText = HelpText(self.version)
            helpText.run()
            sys.exit()
        else:
            HelpText(self.version).printLogo() if not self.args.quiet else None

        if self.args.version:
            print(f'{self.version}')
            sys.exit()

    def configuration(self):
        def validFile(file, extensions):
            return any(file.endswith(ext) for ext in extensions) and os.path.exists(file)

        file_extensions = {
            0: ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'],
            1: ['.fq', '.fastq', '.fq.gz', '.fastq.gz'],
            2: ['.bam', '.sam']
        }

        files_to_check = [
            (0, self.args.assembly, "assembly"),
            (0, self.args.reference, "reference"),
            (1, self.args.left, "read"),
            (1, self.args.right, "read"),
            (2, self.args.bam, "bam")
        ]

        for file_type, file, file_desc in files_to_check:
            if file and any(not validFile(f, file_extensions[file_type]) for f in file.split(',')):
                sys.exit(f"Error: {file} contains an invalid {file_desc} file")

        aligner_map = {
            True: 'hisat2',
            False: 'bowtie2'
        }
        self.t2.aligner = aligner_map[self.args.hisat2]

        mode_map = {
            self.args.bam: ('BAM', self.args.bam),
            self.args.left and self.args.right: (2, (self.args.left, self.args.right)),
            self.args.left or self.args.right: (1, self.args.left or self.args.right),
            True: (0, None)
        }
        self.t2.mode, mode_value = next((mode, value) for condition, (mode, value) in mode_map.items() if condition)
        if self.t2.mode == 'BAM':
            self.t2.bam = mode_value
        elif self.t2.mode == 2:
            self.t2.left, self.t2.right = mode_value
        elif self.t2.mode == 1:
            self.t2.single = mode_value

        self.assemblies       = self.args.assembly.split(',')
        self.t2.assemblyTotal = len(self.assemblies)
        self.t2.reference     = self.args.reference or ''
        self.t2.output        = self.args.output or 'transrate2'
        self.t2.threads       = max(min(self.args.threads, psutil.cpu_count(logical=True)), 4)
        self.t2.clutter       = self.args.clutter
        self.t2.quiet         = self.args.quiet
        self.t2.debug         = self.args.debug
        self.t2.mode          = bool(self.args.left) + bool(self.args.right)


    def run(self):
        try:
            for assembly in self.assemblies:
                self.t2.assemblyCount += 1
                self.t2.assembly = assembly
                self.t2.run()
        except KeyboardInterrupt:
            print("")
            

main = Transrate2()
main.parseArgs()
main.configuration()
main.run()
