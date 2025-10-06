from pysam import FastaFile

class CalcFasta:
    def __init__(self):
        return
    
    def mainRun(self, args):
        assembly_data = args
        fasta         = FastaFile(assembly_data['assembly'])
        fastaLengths  = fasta.lengths
        count         = len(fastaLengths)
        meanLen       = sum(fastaLengths) / count
        nBases        = assembly_data['assemblyDF']['basesN'].iloc[0]
        bases         = assembly_data['assemblyDF']['bases'].iloc[0]
        gcCount       = assembly_data['assemblyDF']['gcCount'].iloc[0]
        orfLenSum     = assembly_data['assemblyDF']['orfLenSum'].iloc[0]
        iterDct       = {'assembly': {'pGC'            : 0, 
                                      'meanOrfPercent' : 0,
                                      'fragments'      : 0,
                                      'pN'             : 0
                                    }}
        
        iterDct['assembly']['pN']             = nBases / bases
        iterDct['assembly']['meanOrfPercent'] = (300 * orfLenSum) / (count * meanLen)
        iterDct['assembly']['pGC']            = gcCount / bases
        iterDct['assembly']['fragments']      = assembly_data['readCount']
        return iterDct