import os
import pandas as pd
import json
from pysam import FastaFile

class FileSetup:
    
    def __init__(self, main):
        self.mode         = main.mode
        self.aligner      = main.aligner
        self.output       = main.output


        main.assemblyName = os.path.basename(main.assembly).split('.')[0]
        main.alignerDir   = os.path.join(main.output, f'{main.aligner}_{main.assemblyName}')
        main.transrateDir = os.path.join(main.output, f'transrate2_{main.assemblyName}')
        main.logDir       = os.path.join(main.output, f'logs_{main.assemblyName}')
        main.salmonDir    = os.path.join(main.output, f'salmon_{main.assemblyName}')
        main.alignerIndex = os.path.join(main.alignerDir, 'index')
        main.salmonIndex  = os.path.join(main.salmonDir, 'index')
    
    def createDirs(self, main):
        for d in [main.transrateDir, main.logDir]:
            os.makedirs(d, exist_ok=True)

        if self.mode > 0:
            for d in [main.alignerDir, main.alignerIndex, main.salmonDir, main.salmonIndex]:
                os.makedirs(d, exist_ok=True)
        
    def nameSetup(self, main):
        main.alignerBam   = os.path.join(main.alignerDir, main.assemblyName + '_aligned.sam' if main.aligner == 'hisat2' else main.assemblyName + '.bam')
        main.salmonBam    = os.path.join(main.salmonDir, main.assemblyName + '.postSample' + '.bam')
        main.sortedBam    = os.path.join(main.salmonDir, main.assemblyName + '.postSample.sorted.bam')

        main.salmonQuant  = os.path.join(main.salmonDir, 'quant.sf')
        main.trCsv        = os.path.join(main.transrateDir, main.assemblyName + '.transrate2.csv')
        main.contigCSV    = os.path.join(main.transrateDir, main.assemblyName + '.contigs.csv')
        main.assemblyCSV  = os.path.join(main.transrateDir, 'assembly.csv')
        main.scoreOptCSV  = os.path.join(main.transrateDir, 'assembly_score_optimisation.csv')
        main.goodContig   = os.path.join(main.transrateDir, 'good.' + main.assemblyName + '.fa')
        main.badContig    = os.path.join(main.transrateDir, 'bad.' + main.assemblyName + '.fa')
        main.assemblyTmp  = self.assemblyTmp(main)
        main.toolTmp      = self.toolTmp(main)
        
    def assemblyTmp(self, main):
        path = os.path.join(main.transrateDir, main.assemblyName + '.assembly.json')
        if os.path.exists(path):
            os.remove(path)
        
        # Create a dictionary with headers as keys and empty values
        empty_data = {header: None for header in self.aHeaders(main)}
        
        # Write the empty data to a JSON file
        with open(path, 'w') as f:
            json.dump(empty_data, f, indent=4)
        
        return path
    
    def toolTmp(self, main):
        path = os.path.join(main.transrateDir, main.assemblyName + '.tool.tmp')
        if os.path.exists(path):
            os.remove(path)
        return path

    def cHeaders(self, main):
        base = [
            'name', 'length', 'gcCount', 'pGC', 'orfLength', 'fragments', 'softclipped', 
            'pSoftclipped', 'basesUncovered', 'pBasesCovered', 'pSeqTrue', 
            'effLength', 'effCount', 'tpm', 'coverage', 'sCnuc', 'sCcov'
        ]
        
        paired = [
            'bridges', 'properPair', 'good', 'pGood', 'score', 'pNotSegmented', 
            'sCord', 'sCseg', 'bothMapped'
        ]
        
        combined_headers = base + paired if main.mode == 2 else base
        
        final_order = [
            'name', 'length', 'fragments', 'gcCount', 'pGC', 'basesUncovered','pBasesCovered',
            'bridges', 'bothMapped', 'properPair', 'good', 'pGood', 'orfLength', 'pNotSegmented',
            'pSeqTrue', 'softclipped', 'pSoftclipped',
            'effLength', 'effCount', 'tpm', 'coverage', 'sCnuc',
            'sCcov', 'sCord', 'sCseg', 'score'
        ]
        
        cHeaders = [header for header in final_order if header in combined_headers]
        return cHeaders
    
    def aHeaders(self, main):
        base = [
            'assembly', 'nSeqs', 'bases', 'smallest', 'largest', 'basesN', 'meanLength',
            'medianLength', 'stdLength', 'nUnder200', 'nOver1k', 'nOver10k', 
            'nWithOrf', 'meanOrfPercent', 'n90', 'n70', 'n50', 'n30', 'n10',
            'gcCount', 'pGC', 'basesN', 'pN'
        ]

        single = [
            'fragments', 'fragmentsMapped', 'pFragmentsMapped', 'softclipped', 
            'pSoftclipped', 'basesUncovered', 'pBasesUncovered', 
            'contigsUncovBase', 'pContigsUncovbase', 'contigsUncovered', 
            'pContigsUncovered', 'contigsLowcovered', 'pContigsLowcovered',
        ]

        paired = [
            'goodMappings', 'pGoodMappings', 'badMappings', 'potentialBridges', 
            'contigsSegmented', 'pContigsSegmented', 'score', 
            'optimalScore', 'cutoff', 'weighted', 'goodContigs', 'badContigs',
        ]

        final_order = [
            'assembly', 'nSeqs', 'bases', 'smallest', 'largest', 'meanLength',
            'medianLength', 'stdLength', 'nUnder200', 'nOver1k', 'nOver10k', 
            'nWithOrf', 'meanOrfPercent', 'n90', 'n70', 'n50', 'n30', 'n10',
            'gcCount', 'pGC', 'basesN', 'pN', 'fragments', 'fragmentsMapped', 
            'pFragmentsMapped', 'softclipped', 'pSoftclipped', 
            'goodMappings', 'pGoodMappings', 'badMappings', 'potentialBridges', 
            'basesUncovered', 'pBasesUncovered', 'contigsUncovBase', 'pContigsUncovbase', 
            'contigsUncovered', 'pContigsUncovered', 'contigsLowcovered', 'pContigsLowcovered', 
            'contigsSegmented', 'pContigsSegmented','goodContigs', 'badContigs',
            'cutoff', 'weighted', 'score', 'optimalScore', 
        ]

        if main.mode == 2:
            combined_headers = base + single + paired
        elif main.mode == 1:
            combined_headers = base + single
        else:
            combined_headers = base

        aHeaders = [header for header in final_order if header in combined_headers]
        return aHeaders

    def rowSetup(self, main):
        fasta = FastaFile(main.assembly)
        refs  = fasta.references
        if main.mode > 0:
            main.cHeaders = self.cHeaders(main)
            contig_data = [{'name': ref, 'length': fasta.get_reference_length(ref)} for ref in refs]
            main.contigDF = pd.DataFrame(contig_data, columns=main.cHeaders)
        else:
            pd.DataFrame()

    def argsSetup(self, main):
        main.args = {'args': {}, 
                     'tools': ['salmon', 'samtools'],
                     'assemblyTmp': main.assemblyTmp
                     }

        if main.mode > 0 and main.aligner:
            main.args['tools'].append(main.aligner)

        def truncate(value):
            return '...' + str(value)[-60:] if len(str(value)) > 63 else str(value)

        main.args['args'] = {
            'Assembly': truncate(main.assembly),
            'Output': truncate(main.output),
            'Mode': {2: 'Paired', 1: 'Single'}.get(main.mode, 'Assembly Only'),
            'Left': truncate(main.left) if main.left else None,
            'Right': truncate(main.right) if main.right else None,
            'Aligner': str(main.aligner) if main.mode > 0 else None,
            'Clutter': str(main.clutter) if main.clutter else None,
            'Threads': str(main.threads) if main.threads else None,
            'Reference': truncate(main.reference) if main.reference else None,
            'Single': str(main.single) if main.single else None,
        }
        main.args['args'] = {k: v for k, v in main.args['args'].items() if v is not None}

        

    def run(self, main):
        self.createDirs(main)
        self.nameSetup(main)
        self.rowSetup(main) if main.mode > 0 else None
        self.argsSetup(main)
        main.aHeaders = self.aHeaders(main)
        assemblyVal = {'assembly': main.assemblyName}
        main.assemblyDF = pd.DataFrame(columns=main.aHeaders)
        main.assemblyDF.loc[0] = assemblyVal
        return 