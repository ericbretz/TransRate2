import numpy as np
from pysam import FastaFile

class NStats:
    
    def mainRun(self, args):
        assembly_data                       = args
        bases                               = assembly_data['assemblyDF']['bases'].iloc[0]
        fasta                               = FastaFile(assembly_data['assembly'])
        lengths                             = fasta.lengths
        baseDct                             = {'assembly': {'n90': 0, 
                                                            'n70': 0, 
                                                            'n50': 0, 
                                                            'n30': 0, 
                                                            'n10': 0}}

        self.b90                            = [int(bases * 0.9), False]
        self.b70                            = [int(bases * 0.7), False]
        self.b50                            = [int(bases * 0.5), False]
        self.b30                            = [int(bases * 0.3), False]
        self.b10                            = [int(bases * 0.1), False]

        baseDct['assembly']['smallest']     = min(lengths)
        baseDct['assembly']['largest']      = max(lengths)
        baseDct['assembly']['medianLength'] = int(np.median(lengths))
        baseDct['assembly']['meanLength']   = int(np.mean(lengths))
        baseDct['assembly']['stdLength']    = int(np.std(lengths))
        baseDct['assembly']['nSeqs']        = len(lengths)

        self.sortedLengths = sorted(lengths, reverse=True)

        def nCalc():
            total = 0
            for length in self.sortedLengths:
                total += length
                if total >= self.b90[0] and not self.b90[1]:
                    baseDct['assembly']['n90'] = length
                    self.b90[1] = True
                if total >= self.b70[0] and not self.b70[1]:
                    baseDct['assembly']['n70'] = length
                    self.b70[1] = True
                if total >= self.b50[0] and not self.b50[1]:
                    baseDct['assembly']['n50'] = length
                    self.b50[1] = True
                if total >= self.b30[0] and not self.b30[1]:
                    baseDct['assembly']['n30'] = length
                    self.b30[1] = True
                if total >= self.b10[0] and not self.b10[1]:
                    baseDct['assembly']['n10'] = length
                    self.b10[1] = True
            return
        
        nCalc()
        return baseDct