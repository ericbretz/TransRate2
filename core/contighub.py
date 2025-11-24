import warnings
warnings.filterwarnings('ignore')
from core.contig import base, frag, good, seqs, sgmt, file
from multiprocessing import Pool
import pandas as pd
    

class ContigHub:
    def __init__(self, contig_data, printout_class):
        self.printout_class = printout_class
        self.printout       = printout_class.printout
        self.metrics        = [base, frag, good, seqs, sgmt] if contig_data['mode'] == 2 else [base, frag, seqs]
        self.mask           = contig_data['contigDF']['name'].values
        self.bannerDct      = {
            'base' : 'Base Stats',
            'frag' : 'Fragmentation',
            'good' : 'Good/Bad Reads',
            'seqs' : 'Read Coverage',
            'sgmt' : 'Contig Univariance',
        }
        return

    def run(self, contig_data):
        file.main(contig_data)
        
        self.printout_class.start_section("Contig Stats")
        self.printout_class.start_progress_stage("Contig Analysis")
        
        total_metrics   = len(self.metrics)
        total_threads   = contig_data['threads']
        total_tasks     = total_metrics * total_threads
        completed_tasks = 0
        
        for i, metric in enumerate(self.metrics):
            metric_name = self.bannerDct[metric.__name__.split('.')[-1]]
            progress    = (completed_tasks / total_tasks) * 100
            self.printout_class.update_progress_bar(int(progress), 100, metric_name)
            
            completed_tasks += self.metricPoolWithProgress(contig_data, metric, completed_tasks, total_tasks, metric_name)
            stage_metrics    = self._get_stage_metrics(metric, contig_data)
            self.printout_class.add_stage_to_section(stage_metrics)
        
        self.printout_class.complete_progress_stage()
        self.printout_class.complete_section()
        
        return
                
    def metricPool(self, contig_data, metric):
        with Pool(contig_data['threads']) as p:
            results = p.map(metric.mainRun, [[contig_data, i] for i in range(contig_data['threads'])])

        updates = []
        for result in results:
            for ref, categories in result.items():
                if 'other' in categories:
                    contig_data['basesDct'][ref]['bases'] = categories['other']['bases']
                updates.append((ref, categories))
        
        if updates:
            dfUpdate = pd.DataFrame.from_dict(dict(updates), orient='index')
            dfUpdate.index.name = 'name'
            
            contig_data['contigDF'] = contig_data['contigDF'].merge(dfUpdate, on='name', how='left', suffixes=('', '_new'))
            
            for col in dfUpdate.columns:
                new_col_name = f'{col}_new'
                if new_col_name in contig_data['contigDF'].columns:
                    contig_data['contigDF'][col] = contig_data['contigDF'][col].combine_first(contig_data['contigDF'][new_col_name])
                    contig_data['contigDF'].drop(columns=[new_col_name], inplace=True)
        return
    
    def metricPoolWithProgress(self, contig_data, metric, completed_tasks, total_tasks, metric_name):
        import multiprocessing as mp
        import time
        
        num_threads = contig_data['threads']
        
        with Pool(num_threads) as pool:
            async_results = []
            for i in range(num_threads):
                async_result = pool.apply_async(metric.mainRun, ([contig_data, i],))
                async_results.append(async_result)
            
            completed_in_stage = 0
            while completed_in_stage < num_threads:
                time.sleep(0.05)
                
                new_completed = sum(1 for result in async_results if result.ready())
                if new_completed > completed_in_stage:
                    completed_in_stage      = new_completed
                    current_total_completed = completed_tasks + completed_in_stage
                    progress                = (current_total_completed / total_tasks) * 100
                    self.printout_class.update_progress_bar(int(progress), 100, metric_name)
            
            results = [result.get() for result in async_results]
        
        updates = []
        for result in results:
            for ref, categories in result.items():
                if 'other' in categories:
                    contig_data['basesDct'][ref]['bases'] = categories['other']['bases']
                updates.append((ref, categories))
        
        if updates:
            dfUpdate                = pd.DataFrame.from_dict(dict(updates), orient='index')
            dfUpdate.index.name     = 'name'
            contig_data['contigDF'] = contig_data['contigDF'].merge(dfUpdate, on='name', how='left', suffixes=('', '_new'))
            
            for col in dfUpdate.columns:
                new_col_name = f'{col}_new'
                if new_col_name in contig_data['contigDF'].columns:
                    contig_data['contigDF'][col] = contig_data['contigDF'][col].combine_first(contig_data['contigDF'][new_col_name])
                    contig_data['contigDF'].drop(columns=[new_col_name], inplace=True)
        
        return num_threads
    
    def _get_stage_metrics(self, metric, contig_data):
        metric_type = metric.__name__.split('.')[-1]
        df = contig_data['contigDF']
        
        available_metrics = {}
        
        for col_key in self.printout_class.map_contig.keys():
            if col_key in df.columns:
                if col_key == 'name':
                    continue
                elif col_key in [
                            'length'    , 'fragments' , 'gcCount', 'basesUncovered', 'bridges'    , 
                            'bothMapped', 'properPair', 'good'   , 'orfLength'     , 'softclipped', 
                            'effLength' , 'effCount'
                               ]:
                    available_metrics[col_key] = self._format_large_number(df[col_key].sum())
                elif col_key in ['pGC', 'pBasesCovered', 'pGood', 'pNotSegmented', 'pSeqTrue', 'pSoftclipped']:
                    available_metrics[col_key] = f"{df[col_key].mean() * 100:.2f}%"
                elif col_key in ['tpm', 'coverage', 'sCnuc', 'sCcov', 'sCord', 'sCseg', 'score']:
                    available_metrics[col_key] = f"{df[col_key].mean():.3f}"
        
        return available_metrics
    
    def _format_large_number(self, num):
        if pd.isna(num):
            return "N/A"
        if num >= 1000000:
            return f"{num/1000000:.1f}M"
        elif num >= 1000:
            return f"{num/1000:.1f}k"
        else:
            return str(int(num))

                    