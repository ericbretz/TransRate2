import os

class CSV:
    def __init__(self, main):
        self.intVars = [
            'nSeqs'            , 'bases'           , 'smallest'         , 'largest'         , 'meanLength'      ,
            'medianLength'     , 'stdLength'       , 'nUnder200'        , 'nOver1k'         , 'nOver10k'        ,
            'nWithOrf'         , 'n90'             , 'n70'              , 'n50'             , 'n30'             ,
            'n10'              , 'gcCount'         , 'basesN'           , 'fragments'       , 'fragmentsMapped' ,
            'softclipped'      , 'goodMappings'    , 'badMappings'      , 'potentialBridges', 'basesUncovered'  ,
            'contigsUncovBase' , 'contigsUncovered', 'contigsLowcovered', 'contigsSegmented', 'goodContigs'     ,
            'badContigs'       , 'length'          , 'bothMapped'       , 'properPair'      , 'good'            ,
            'orfLength'        , 'effCount'        , 'coverage'         ,
        ]
        
        self.floatVars = [
            'meanOrfPercent'   , 'pGC'            , 'proportionN'      , 'pFragmentsMapped' , 'pSoftclipped'      ,
            'pGoodMappings'    , 'pBasesUncovered', 'pContigsUncovbase', 'pContigsUncovered', 'pContigsLowcovered',
            'pContigsSegmented', 'cutoff'         , 'weighted'         , 'score'            , 'optimalScore'      ,
            'pBasesCovered'    , 'pGood'          , 'pNotSegmented'    , 'pSeqTrue'         , 'effLength'         ,
            'tpm'              , 'sCnuc'          , 'sCcov'            , 'sCord'            , 'sCseg'             ,
            'pN'
        ]

        return
    
    def contigCSV(self, main):
        main.contigDF = main.contigDF[main.cHeaders]
        
        dtype_map = {col: 'int' for col in self.intVars if col in main.contigDF.columns}
        dtype_map.update({col: 'float' for col in self.floatVars if col in main.contigDF.columns})
        
        main.contigDF = main.contigDF.astype(dtype_map)
        
        float_columns_to_round = [col for col in self.floatVars if col in main.contigDF.columns]
        main.contigDF[float_columns_to_round] = main.contigDF[float_columns_to_round].round(4)
        
        main.contigDF.to_csv(main.contigCSV, index=False)
        return
    
    def assemblyCSV(self, main):
        main.assemblyDF['assembly'] = os.path.basename(main.assembly)
        main.assemblyDF = main.assemblyDF[main.aHeaders]

        dtype_map = {col: 'int' for col in self.intVars if col in main.assemblyDF.columns}
        dtype_map.update({col: 'float' for col in self.floatVars if col in main.assemblyDF.columns})
        
        main.assemblyDF = main.assemblyDF.astype(dtype_map)
        
        float_columns_to_round = [col for col in self.floatVars if col in main.assemblyDF.columns]
        main.assemblyDF[float_columns_to_round] = main.assemblyDF[float_columns_to_round].round(4)
        
        main.assemblyDF.to_csv(main.assemblyCSV, index=False)
        return