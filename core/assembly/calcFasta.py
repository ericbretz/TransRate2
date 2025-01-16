from pysam import FastaFile

class CalcFasta:
    def __init__(self):
        return
    
    def mainRun(self, args):
        main         = args
        fasta        = FastaFile(main['assembly'])
        fastaLengths = fasta.lengths
        count        = len(fastaLengths)
        meanLen      = sum(fastaLengths) / count
        nBases       = main['assemblyDF']['basesN'].iloc[0]
        bases        = main['assemblyDF']['bases'].iloc[0]
        gcCount      = main['assemblyDF']['gcCount'].iloc[0]
        orfLenSum    = main['assemblyDF']['orfLenSum'].iloc[0]
        iterDct      = {'assembly': {'pGC'            : 0, 
                                     'meanOrfPercent' : 0,
                                     'fragments'      : 0,
                                     'pN'             : 0
                                   }}
        
        iterDct['assembly']['pN']             = nBases / bases
        iterDct['assembly']['meanOrfPercent'] = (300 * orfLenSum) / (count * meanLen)
        iterDct['assembly']['pGC']            = gcCount / bases
        iterDct['assembly']['fragments']      = main['readCount']
        return iterDct