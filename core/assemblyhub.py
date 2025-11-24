import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
from multiprocessing import Pool
from core.assembly import iterFasta
from core.assembly import calcFasta
from core.assembly import nStats
from core.assembly import salmonStats
from core.assembly import iterContig
from core.assembly import calcContig
from core.assembly import scoreCalc
from core.assembly import goodContig
# from core.assembly import reference
    
import pandas as pd

class AssemblyHub:
    def __init__(self, assembly_data, printout_class):
        self.printout_class = printout_class
        self.printout       = printout_class.printout
        pass
    
    def run(self, assembly_data):
        metric_classes = {
            "Analysing FASTA"      : (iterFasta.IterFasta()    , True                , True),
            "FASTA Statistics"     : (calcFasta.CalcFasta()    , True                , False),
            "N-Statistics"         : (nStats.NStats()          , True                , False),
            "Salmon Quantification": (salmonStats.SalmonStats(), assembly_data['mode'] > 0 , False),
            "Analysing Contigs"    : (iterContig.IterContig()  , assembly_data['mode'] > 0 , True),
            "Contig Calculations"  : (calcContig.CalcContig()  , assembly_data['mode'] > 0 , False),
            "Assembly Scoring"     : (scoreCalc.ScoreCalc()    , assembly_data['mode'] > 0 , False),
            "Quality Assessment"   : (goodContig.goodContig()  , assembly_data['mode'] > 0 , False),
            # "Reference Analysis"   : (reference.Reference()    , assembly_data['reference'], False)
        }
        
        self.printout_class.start_section("Assembly Stats")
        self.printout_class.start_progress_stage("Assembly Analysis")
        
        active_metrics = [(name, metric_class, threaded) for name, (metric_class, condition, threaded) in metric_classes.items() if condition]
        total_tasks    = 0

        for name, metric_class, threaded in active_metrics:
            if threaded:
                total_tasks += assembly_data['threads']
            else:
                total_tasks += 1
        
        completed_tasks = 0
        
        for i, (name, metric_class, threaded) in enumerate(active_metrics):
            progress = (completed_tasks / total_tasks) * 100
            self.printout_class.update_progress_bar(int(progress), 100, name)
            
            if threaded:
                stage_metrics    = self.metricProcessWithProgress(name, assembly_data, metric_class, completed_tasks, total_tasks)
                completed_tasks += assembly_data['threads']
            else:
                stage_metrics    = self.metricProcess(name, assembly_data, metric_class, threaded=False)
                completed_tasks += 1
                progress         = (completed_tasks / total_tasks) * 100
                self.printout_class.update_progress_bar(int(progress), 100, name)
            
            self.printout_class.add_stage_to_section(stage_metrics)
        
        self.printout_class.complete_progress_stage()
        self.printout_class.complete_section()
        
        return

    def metricProcess(self, name, assembly_data, function, threaded=False):
        if threaded:
            with Pool(assembly_data['threads']) as p:
                results = p.map(function.mainRun, [[assembly_data, i] for i in range(assembly_data['threads'])])
        else:
            results = [function.mainRun(assembly_data)]

        stage_metrics = self._get_stage_metrics(name, results, assembly_data)
        
        self.resultsProcess(assembly_data, results)
        
        return stage_metrics
    
    def metricProcessWithProgress(self, name, assembly_data, function, completed_tasks, total_tasks):
        import time
        
        num_threads = assembly_data['threads']
        
        with Pool(num_threads) as pool:
            async_results = []
            for i in range(num_threads):
                async_result = pool.apply_async(function.mainRun, ([assembly_data, i],))
                async_results.append(async_result)
            
            completed_in_stage = 0
            while completed_in_stage < num_threads:
                time.sleep(0.05)
                
                new_completed = sum(1 for result in async_results if result.ready())
                if new_completed > completed_in_stage:
                    completed_in_stage      = new_completed
                    current_total_completed = completed_tasks + completed_in_stage
                    progress                = (current_total_completed / total_tasks) * 100
                    self.printout_class.update_progress_bar(int(progress), 100, name)
            
            results = [result.get() for result in async_results]
        
        stage_metrics = self._get_stage_metrics(name, results, assembly_data)
        self.resultsProcess(assembly_data, results)
        
        return stage_metrics

    def resultsProcess(self, assembly_data, results):
        for result in results:
            for dfType, dfDct in result.items():
                if dfType == 'assembly':
                    self.assemblyProcess(assembly_data, dfDct)
                elif dfType == 'contigs':
                    self.contigProcess(assembly_data, dfDct) 
        self.save_assemblyDF_to_json(assembly_data['assemblyTmp'], assembly_data['assemblyDF'], assembly_data)

    def assemblyProcess(self, assembly_data, dfDct):
        if assembly_data['assemblyDF'].empty:
            assembly_data['assemblyDF'].loc[0] = None
        
        idx      = assembly_data['assemblyDF'].index[0]
        new_cols = set(dfDct.keys()) - set(assembly_data['assemblyDF'].columns)
        if new_cols:
            for col in new_cols:
                assembly_data['assemblyDF'][col] = None
        for key, value in dfDct.items():
            if pd.isna(assembly_data['assemblyDF'].loc[idx, key]):
                assembly_data['assemblyDF'].loc[idx, key] = 0
            assembly_data['assemblyDF'].loc[idx, key] += value

    def save_assemblyDF_to_json(self, path, df, assembly_data):
        filtered_df = df.loc[:, df.columns.isin(assembly_data['aHeaders'])]
        
        with open(path, 'w') as f:
            filtered_df.to_json(f, orient='records', lines=True)
        
    def contigProcess(self, assembly_data, dfDct):
        if assembly_data['mode'] == 0:
            return
        if not dfDct:
            return
        
        innerDict = next(iter(dfDct.values()), {})
        newCols   = set(innerDict.keys()) - set(assembly_data['contigDF'].columns)
        
        for col in newCols:
            assembly_data['contigDF'][col] = None
        
        dfUpdate            = pd.DataFrame.from_dict(dfDct, orient='index')
        dfUpdate.index.name = 'name'
        
        assembly_data['contigDF'] = assembly_data['contigDF'].merge(dfUpdate, on='name', how='left', suffixes=('', '_new'))
        
        for col in dfUpdate.columns:
            assembly_data['contigDF'][col] = assembly_data['contigDF'][col].combine_first(assembly_data['contigDF'][f'{col}_new'])
            assembly_data['contigDF'].drop(columns=[f'{col}_new'], inplace=True)
    
    def _get_stage_metrics(self, stage_name, results, assembly_data):
        available_metrics = {}
        
        for key in self.printout_class.map_assembly.keys():
            if key in assembly_data['assemblyDF'].columns:
                value = assembly_data['assemblyDF'][key].iloc[0]
                if pd.isna(value):
                    continue
                    
                if key in [
                    'nSeqs'           , 'bases'            , 'smallest'        , 'largest'         , 'meanLength'       , 'medianLength'    , 
                    'stdLength'       , 'nUnder200'        , 'nOver1k'         , 'nOver10k'        , 'nWithOrf'         , 'n90'             , 
                    'n70'             , 'n50'              , 'n30'             , 'n10'             , 'gcCount'          , 'basesN'          ,
                    'fragments'       , 'fragmentsMapped'  , 'bothMapped'      , 'softclipped'     , 'goodMappings'     , 'contigsUncovered',
                    'badMappings'     , 'potentialBridges' , 'basesUncovered'  , 'contigsUncovBase', 'contigsLowcovered', 'contigsSegmented', 
                    'goodContigs'     , 'badContigs'
                          ]:
                    available_metrics[key] = self._format_large_number(value)

                elif key in [
                    'pGC'              , 'meanOrfPercent'    , 'pN'               , 'pFragmentsMapped', 
                    'pSoftclipped'     , 'pGoodMappings'     , 'pBasesUncovered'  , 'pContigsUncovbase', 
                    'pContigsUncovered', 'pContigsLowcovered', 'pContigsSegmented'
                    ]:
                    # available_metrics[key] = f"{value * 100:.2f}%"
                    available_metrics[key] = f"{value:.2f}"
                elif key in ['score', 'optimalScore', 'cutoff', 'weighted']:
                    available_metrics[key] = f"{value:.2f}"
                else:
                    available_metrics[key] = str(value)
        
        return available_metrics
    
    def _format_large_number(self, num):
        # if num >= 1000000:
        #     return f"{num/1000000:.1f}M"
        # elif num >= 1000:
        #     return f"{num/1000:.1f}k"
        # else:
        #     return str(num)
        return str(num)