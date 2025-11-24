import os
import numpy as np
import pandas as pd

class CSV:
    def __init__(self, printout_class):
        self.printout_class = printout_class
        
        self.intVars = [
            'nSeqs'            , 'bases'           , 'smallest'         , 'largest'         , 'meanLength'      ,
            'medianLength'     , 'stdLength'       , 'nUnder200'        , 'nOver1k'         , 'nOver10k'        ,
            'nWithOrf'         , 'n90'             , 'n70'              , 'n50'             , 'n30'             ,
            'n10'              , 'gcCount'         , 'basesN'           , 'fragments'       , 'fragmentsMapped' ,
            'softclipped'      , 'goodMappings'    , 'badMappings'      , 'potentialBridges', 'basesUncovered'  ,
            'contigsUncovBase' , 'contigsUncovered', 'contigsLowcovered', 'contigsSegmented', 'goodContigs'     ,
            'badContigs'       , 'length'          , 'bothMapped'       , 'properPair'      , 'good'            ,
            'orfLength'        , 'effCount'        , 'coverage'         , 'CRBBhits'        , 'nContigsWithCRBB',
            'nRefsWithCRBB'    , 'cov25'           , 'cov50'            , 'cov75'           , 'cov85'           ,
            'cov95'            
        ]
        
        self.floatVars = [
            'meanOrfPercent'   , 'pGC'             , 'proportionN'      , 'pFragmentsMapped' , 'pSoftclipped'      ,
            'pGoodMappings'    , 'pBasesUncovered' , 'pContigsUncovbase', 'pContigsUncovered', 'pContigsLowcovered',
            'pContigsSegmented', 'cutoff'          , 'weighted'         , 'score'            , 'optimalScore'      ,
            'pBasesCovered'    , 'pGood'           , 'pNotSegmented'    , 'pSeqTrue'         , 'effLength'         ,
            'tpm'              , 'sCnuc'           , 'sCcov'            , 'sCord'            , 'sCseg'             ,
            'pN'               , 'pContigsWithCRBB', 'rbhPerReference'  , 'pRefsWithCRBB'    , 'pCov25'            ,
            'pCov50'           , 'pCov75'          , 'pCov85'           , 'pCov95'           , 'referenceCoverage'
        ]

    def contigCSV(self, assembly_data):
        contigDF  = assembly_data['contigDF'][assembly_data['cHeaders']].copy()
        
        for col in self.intVars:
            if col in contigDF.columns:
                contigDF[col] = contigDF[col].fillna(0)
                contigDF[col] = contigDF[col].replace([np.inf, -np.inf], 0)
        
        dtype_map = {col: 'int' for col in self.intVars if col in contigDF.columns}
        dtype_map.update({col: 'float' for col in self.floatVars if col in contigDF.columns})
        
        contigDF  = contigDF.astype(dtype_map)
        
        float_columns_to_round           = [col for col in self.floatVars if col in contigDF.columns]
        contigDF[float_columns_to_round] = contigDF[float_columns_to_round].round(4)
        
        contigDF.to_csv(assembly_data['contigCSV'], index=False)
        contig_str = f'...{assembly_data["contigCSV"][-57:]}' if len(assembly_data['contigCSV']) > 60 else assembly_data['contigCSV']
        print(f"{'Contig CSV:':<19} {contig_str}")
        return
    
    def assemblyCSV(self, assembly_data):
        assemblyDF             = assembly_data['assemblyDF'].copy()
        assemblyDF['assembly'] = os.path.basename(assembly_data['assembly'])
        assemblyDF             = assemblyDF[assembly_data['aHeaders']]
        
        for col in self.intVars:
            if col in assemblyDF.columns:
                assemblyDF[col] = assemblyDF[col].fillna(0)
                assemblyDF[col] = assemblyDF[col].replace([np.inf, -np.inf], 0)

        dtype_map = {col: 'int' for col in self.intVars if col in assemblyDF.columns}
        dtype_map.update({col: 'float' for col in self.floatVars if col in assemblyDF.columns})
        
        assemblyDF = assemblyDF.astype(dtype_map)
        
        float_columns_to_round             = [col for col in self.floatVars if col in assemblyDF.columns]
        assemblyDF[float_columns_to_round] = assemblyDF[float_columns_to_round].round(4)
        
        assemblyCSV_path = os.path.join(assembly_data['dict_dir']['results'], 'assembly.csv')
        assemblyDF.to_csv(assemblyCSV_path, index=False)

        assembly_str = f'...{assemblyCSV_path[-57:]}' if len(assemblyCSV_path) > 60 else assemblyCSV_path
        print(f"{'Assembly CSV:':<19} {assembly_str}")
        return
