from pysam import FastaFile

class IterFasta:
    
    def __init__(self):
        return
    
    def mainRun(self, args):
        assembly_data = args[0]
        i             = args[1]
        fasta         = FastaFile(assembly_data['assembly'])
        refs          = fasta.references[i::assembly_data['threads']]
        iterDct = {'assembly':  {'nUnder200': 0,
                                 'nOver1k'  : 0,
                                 'nOver10k' : 0,
                                 'nWithOrf' : 0,
                                 'orfLenSum': 0,
                                 'bases'    : 0,
                                 'basesN'   : 0,
                                 'gcCount'  : 0,
                                 'refCount' : len(refs)
                                },
                    'contigs': {
                        ref:    {'orfLength': 0, 
                                'pGC'       : 0,
                                'gcCount'   : 0,
                                'basesN'    : 0} 
                                for ref in refs}
                                }
        for ref in refs:
            seq                                  = fasta.fetch(ref)
            orfLength, orfLong                   = self.orfLongest(seq)
            iterDct['contigs'][ref]['orfLength'] = orfLength
            iterDct['assembly']['nWithOrf']     += orfLong
            iterDct['assembly']['orfLenSum']    += orfLength
            iterDct['assembly']['nUnder200']    += 1 if len(seq) < 200 else 0
            iterDct['assembly']['nOver1k']      += 1 if len(seq) > 1000 else 0
            iterDct['assembly']['nOver10k']     += 1 if len(seq) > 10000 else 0

            gc, n                                = self.baseCounts(seq)
            iterDct['contigs'][ref]['pGC']       = gc / len(seq)
            iterDct['assembly']['bases']        += len(seq)
            iterDct['assembly']['basesN']       += n
            iterDct['assembly']['gcCount']      += gc
            iterDct['contigs'][ref]['gcCount']   = gc
            iterDct['contigs'][ref]['basesN']    = n

        return iterDct
        
    def orfLongest(self, seq):
        longest  = 0
        len_list = [0, 0, 0]
        sl       = len(seq)
        for i in range(sl - 2):
            if seq[i:i+3] == 'atg':
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 1
            elif seq[i:i+3] in ['tag', 'taa', 'tga']:
                longest = max(longest, len_list[i % 3])
                len_list[i % 3] = -1
            else:
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 0
        longest  = max(longest, max(len_list))
        len_list = [0, 0, 0]
        for i in range(sl - 1, 1, -1):
            if seq[i-2:i+1] == 'cat':
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 1
            elif seq[i-2:i+1] in ['cta', 'tta', 'tct']:
                longest = max(longest, len_list[i % 3])
                len_list[i % 3] = -1
            else:
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 0

        orfLength = max(longest, max(len_list))
        orfLong   = 1 if orfLength > 149 else 0
        return orfLength, orfLong
    
    def baseCounts(self, seq):
        gcCount = seq.upper().count('G') + seq.upper().count('C')
        nBases  = seq.upper().count('N')
        return gcCount, nBases