from core.fileSystem import FileSetup
from core.alignments.hisat2 import Hisat2
from core.alignments.bowtie2 import Bowtie2
from core.alignments.salmon import Salmon
from core.alignments.samtools import Samtools
from core.contigHub import ContigHub
from core.assemblyHub import AssemblyHub
from core.csvProcess import CSV
from core.terminalHub import TerminalMain, QuietMain
import pandas as pd
import shutil
import glob
import os
import sys

class TransRate:
    def __init__(self):
        # Aligner
        self.aligner       :str  = ''
        # Bam Files
        self.alignerBam    :str  = ''
        self.salmonBam     :str  = ''
        self.sortedBam     :str  = ''
        self.bam           :str  = ''
        # Index
        self.alignerIndex  :str  = ''
        self.salmonIndex   :str  = ''
        # Assembly
        self.assembly      :str  = ''
        self.assemblyName  :str  = ''
        self.assemblyTotal :int  = 0
        self.assemblyCount :int  = 0
        # Reads
        self.left          :str  = ''
        self.right         :str  = ''
        self.single        :str  = ''
        self.readCount     :int  = 0
        self.refList       :list = []
        # Reference
        self.reference     :str  = ''
        # Output
        self.output        :str  = ''
        # Threads
        self.threads       :int  = 1
        # Mode
        self.mode          :int  = 0
        # Bool Flags
        self.clutter       :bool = False
        self.quiet         :bool = False
        self.debug         :bool = False
        # Directories
        self.transrateDir  :str  = ''
        self.alignerDir    :str  = ''
        self.logDir        :str  = ''
        self.salmonDir     :str  = ''
        # CSVs
        self.aHeaders      :list = []
        self.cHeaders      :list = []
        self.trCsv         :str  = ''
        self.contigCSV     :str  = ''
        self.assemblyCSV   :str  = ''
        self.scoreOptCSV   :str  = ''
        self.goodContig    :str  = ''
        self.badContig     :str  = ''
        # Paths
        self.salmonQuant   :str  = ''
        # Dictionaries
        self.trDct         :dict = {}
        self.basesDct      :dict = {}
        self.assemblyDct   :dict = {}
        # DataFrames
        self.contigDF      :pd.DataFrame = pd.DataFrame()
        self.assemblyDF    :pd.DataFrame = pd.DataFrame()

    def fileSetup(self):
        fs = FileSetup(self)
        fs.run(self)

    def run(self):
        try:
            self.fileSetup()
            self.main()
        except KeyboardInterrupt:
            print("")
        finally:
            pass

    def main(self):
        term = TerminalMain(self.debug, self.mode, self.aligner) if not self.quiet else QuietMain(self.debug, self.mode)
        term.runInfo(self.mode, self.assembly, self.reference, self.left, self.right, self.bam, self.aligner, self.output, self.threads, self.clutter, self.debug)
        term.processBanner()
        csv = CSV(self)
        modeFuncs = {
            'aligner' : self.runAligner,
            'salmon'  : self.runSalmon,
            'samtools': self.runSamtools,
            'contigs' : self.runContigs,
        }

        if self.mode:
            for task in modeFuncs:
                if not self.bam:
                    if task == 'aligner':
                        term.toolPrint('hisat2' if self.aligner == 'hisat2' else 'bowtie2')
                    else:
                        term.toolPrint(task) if task != 'contigs' else term.contigBanner()
                    modeFuncs[task]()
                else:
                    term.toolPrint(task)
                    self.alignerBam = self.bam
                    modeFuncs[task]()

        term.assemblyBanner()
        # term.toolPrint('assembly')
        self.runAssembly()
        # term.toolPrint('csv')
        self.runCSV(csv)
        term.assemblyPrint(self.assemblyCSV)
        if self.clutter:
            self.runClutter()

    def runAligner(self):
        aligner = Hisat2(self) if self.aligner == 'hisat2' else Bowtie2(self)
        aligner.index(self)
        aligner.align(self)

    def runSalmon(self):
        Salmon(self).quant(self)

    def runSamtools(self):
        samtools = Samtools(self)
        samtools.sort(self)
        samtools.index(self)

    def runContigs(self):
        ContigHub(self.__dict__).run(self.__dict__)

    def runAssembly(self):
        AssemblyHub(self.__dict__).run(self.__dict__)

    def runCSV(self, csv):
        if self.mode > 0:
            csv.contigCSV(self)
        csv.assemblyCSV(self)

    def runClutter(self):
        for file in glob.glob(os.path.join(self.transrateDir, '*.csv')) + glob.glob(os.path.join(self.transrateDir, '*.fa')):
            shutil.move(file, os.path.join(self.output, os.path.basename(file)))

        for directory in [self.transrateDir, self.alignerDir, self.salmonDir]:
            if os.path.exists(directory):
                shutil.rmtree(directory)



