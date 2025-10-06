import warnings
warnings.filterwarnings('ignore', category=UserWarning)
import pandas as pd

class SalmonStats:
    def mainRun(self, args):
        assembly_data = args
        length_dict   = dict(zip(assembly_data['contigDF']['name'], assembly_data['contigDF']['length']))
        salmonDct     = {'assembly': 
                              {'contigsLowcovered' : 0, 
                              'contigsUncovered'   : 0,
                              },
                        'contigs': {
                              ref: {'tpm'          : 0, 
                                    'effLength'    : 0, 
                                    'effCount'     : 0} 
                              for ref in assembly_data['refList']}
                        }
        salmonFile = pd.read_csv(assembly_data['salmonQuant'], sep='\t', header=0).set_index('Name').fillna(0)
        
        effLength = salmonFile['EffectiveLength'].to_numpy()
        effCount  = salmonFile['NumReads'].to_numpy()
        tpm       = salmonFile['TPM'].to_numpy()
        
        for i, salmon_ref in enumerate(salmonFile.index):
            matching_contig = None
            for contig_name in length_dict.keys():
                if salmon_ref.startswith(contig_name) or contig_name in salmon_ref:
                    matching_contig = contig_name
                    break
            
            if matching_contig:
                coverage                                            = int(effCount[i] * length_dict[matching_contig] / effLength[i])
                salmonDct['assembly']['contigsLowcovered']         += int(0 < coverage < 10)
                salmonDct['assembly']['contigsUncovered']          += int(coverage == 0)
                salmonDct['contigs'][matching_contig]['tpm']        = tpm[i]
                salmonDct['contigs'][matching_contig]['effLength']  = effLength[i]
                salmonDct['contigs'][matching_contig]['effCount']   = effCount[i]
                salmonDct['contigs'][matching_contig]['coverage']   = coverage
        
        return salmonDct