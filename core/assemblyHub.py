from multiprocessing import Pool
from core.assembly import iterFasta
from core.assembly import calcFasta
from core.assembly import nStats
from core.assembly import salmonStats
from core.assembly import iterContig
from core.assembly import calcContig
from core.assembly import scoreCalc
from core.assembly import goodContig
import time
import pandas as pd

class AssemblyHub:
    def __init__(self, mainDct):
        pass
    
    def run(self, mainDct):

        metric_classes = {
            "Salmon Stats": (salmonStats.SalmonStats(), mainDct['mode'] > 0, False, 'Salmon Stats'),
            "Fasta Iter"  : (iterFasta.IterFasta()    , True               , True,  ''),
            "Fasta Calc"  : (calcFasta.CalcFasta()    , True               , False, 'Assembly Quality'),
            "N Stats"     : (nStats.NStats()          , True               , False, ''),
            "Contig Iter" : (iterContig.IterContig()  , mainDct['mode'] > 0, True,  ''),
            "Contig Calc" : (calcContig.CalcContig()  , mainDct['mode'] > 0, False, 'Contig Stats'),
            "Score Calc"  : (scoreCalc.ScoreCalc()    , mainDct['mode'] > 0, False, 'Assembly Score'),
            "Good Contig" : (goodContig.goodContig()  , mainDct['mode'] > 0, False, 'Good/Bad Contigs')
        }
        timer = time.perf_counter()
        for name, (metric_class, condition, threaded, banner) in metric_classes.items():
            if condition:
                timer = time.perf_counter()
                timeStarted = time.strftime("%H:%M:%S", time.gmtime())
                if banner:
                    print(f'{timeStarted:<30}{banner:<30}', end='\r')
                self.metricProcess(name, mainDct, metric_class, threaded=threaded)
                if banner:
                    print(f'{timeStarted:<30}{banner:<30}{time.perf_counter() - timer:>19.2f}s')
            if mainDct['debug']:
                print(f'{name} total took {time.perf_counter() - timer:.2f} seconds')

    def metricProcess(self, name, mainDct, function, threaded=False):
        timer = time.perf_counter()
        if threaded:
            with Pool(mainDct['threads']) as p:
                results = p.map(function.mainRun, [[mainDct, i] for i in range(mainDct['threads'])])
        else:
            results = [function.mainRun(mainDct)]
        if mainDct['debug']:
            print(f'{name} processtook {time.perf_counter() - timer:.2f} seconds')

        self.resultsProcess(mainDct, results)


    def resultsProcess(self, mainDct, results):
        timer = time.perf_counter()
        for result in results:
            for dfType, dfDct in result.items():
                if dfType == 'assembly':
                    self.assemblyProcess(mainDct, dfDct)
                elif dfType == 'contigs':
                    self.contigProcess(mainDct, dfDct) 
        self.save_assemblyDF_to_json(mainDct['assemblyTmp'], mainDct['assemblyDF'], mainDct)
        if mainDct['debug']:
            print(f'    resultsProcess took {time.perf_counter() - timer:.2f} seconds')

    def assemblyProcess(self, mainDct, dfDct):
        timer = time.perf_counter()
        if mainDct['assemblyDF'].empty:
            mainDct['assemblyDF'].loc[0] = None
        
        idx = mainDct['assemblyDF'].index[0]
        new_cols = set(dfDct.keys()) - set(mainDct['assemblyDF'].columns)
        if new_cols:
            for col in new_cols:
                mainDct['assemblyDF'][col] = None
        for key, value in dfDct.items():
            if pd.isna(mainDct['assemblyDF'].loc[idx, key]):
                mainDct['assemblyDF'].loc[idx, key] = 0
            mainDct['assemblyDF'].loc[idx, key] += value
        if mainDct['debug']:
            print(f'    assemblyProcess took {time.perf_counter() - timer:.2f} seconds')

    def save_assemblyDF_to_json(self, path, df, mainDct):
        timer = time.perf_counter()
        filtered_df = df.loc[:, df.columns.isin(mainDct['aHeaders'])]
        
        with open(path, 'w') as f:
            filtered_df.to_json(f, orient='records', lines=True)
        if mainDct['debug']:
            print(f'    save_assemblyDF_to_json took {time.perf_counter() - timer:.2f} seconds')
        
    def contigProcess(self, mainDct, dfDct):
        timer = time.perf_counter()
        if mainDct['mode'] == 0:
            return
        if not dfDct:
            return
        
        innerDict = next(iter(dfDct.values()), {})
        newCols = set(innerDict.keys()) - set(mainDct['contigDF'].columns)
        
        for col in newCols:
            mainDct['contigDF'][col] = None
        
        dfUpdate = pd.DataFrame.from_dict(dfDct, orient='index')
        dfUpdate.index.name = 'name'
        
        mainDct['contigDF'] = mainDct['contigDF'].merge(dfUpdate, on='name', how='left', suffixes=('', '_new'))
        
        for col in dfUpdate.columns:
            mainDct['contigDF'][col] = mainDct['contigDF'][col].combine_first(mainDct['contigDF'][f'{col}_new'])
            mainDct['contigDF'].drop(columns=[f'{col}_new'], inplace=True)

        if mainDct['debug']:
            print(f'    contigProcess took {time.perf_counter() - timer:.2f} seconds')