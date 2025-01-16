import pandas as pd

class SalmonStats:
    def mainRun(self, args):
        main        = args
        length_dict = dict(zip(main['contigDF']['name'], main['contigDF']['length']))
        salmonDct   = {'assembly': 
                            {'contigsLowcovered' : 0, 
                            'contigsUncovered'   : 0,
                            },
                      'contigs': {
                            ref: {'tpm'          : 0, 
                                  'effLength'    : 0, 
                                  'effCount'     : 0} 
                            for ref in main['refList']}
                      }
        salmonFile = pd.read_csv(main['salmonQuant'], sep='\t', header=0).set_index('Name').fillna(0)
        
        effLength = salmonFile['EffectiveLength'].to_numpy()
        effCount  = salmonFile['NumReads'].to_numpy()
        tpm       = salmonFile['TPM'].to_numpy()
        
        valid_indices = [i for i, ref in enumerate(salmonFile.index) if ref in length_dict]
        
        for i in valid_indices:
            ref = salmonFile.index[i]
            coverage = salmonDct['contigs'][ref]['coverage'] = int(effCount[i] * length_dict[ref] / effLength[i])
            salmonDct['assembly']['contigsLowcovered']      += int(0 < coverage < 10)
            salmonDct['assembly']['contigsUncovered']       += int(coverage == 0)
            salmonDct['contigs'][ref]['tpm']                 = tpm[i]
            salmonDct['contigs'][ref]['effLength']           = effLength[i]
            salmonDct['contigs'][ref]['effCount']            = effCount[i]
        
        return salmonDct