from core.contigs import base, file, frag, good, seqs, sgmt
from multiprocessing import Pool
import time
import pandas as pd
    

class ContigHub:
    def __init__(self, mainDct):
        self.metrics = [base, frag, good, seqs, sgmt] if mainDct['mode'] == 2 else [base, frag, seqs]
        self.mask    = mainDct['contigDF']['name'].values
        self.bannerDct = {
            'base' : 'Base Stats',
            'frag' : 'Fragmentation',
            'good' : 'Good/Bad Reads',
            'seqs' : 'Read Coverage',
            'sgmt' : 'Contig Univariance',
        }
        return
    
    def run(self, mainDct):
        file.main(mainDct)
        for metric in self.metrics:
            timer = time.perf_counter()
            timeStarted = time.strftime("%H:%M:%S", time.gmtime())
            banner = f'{self.bannerDct[metric.__name__.split(".")[-1]]:<30}'
            print(f'{timeStarted:<30}{banner}', end='\r')
            self.metricPool(mainDct, metric)
            print(f'{timeStarted:<30}{banner}{time.perf_counter() - timer:>19.2f}s')
            if mainDct['debug']:
                print(f'{metric.__name__} took {time.perf_counter() - timer:.2f} seconds')
        return
                
    def metricPool(self, mainDct, metric):
        timer = time.perf_counter()
        with Pool(mainDct['threads']) as p:
            results = p.map(metric.mainRun, [[mainDct, i] for i in range(mainDct['threads'])])
        if mainDct['debug']:
            print(f'metricPool took {time.perf_counter() - timer:.2f} seconds')
        
        timer = time.perf_counter()
        updates = []
        for result in results:
            for ref, categories in result.items():
                if 'other' in categories:
                    mainDct['basesDct'][ref]['bases'] = categories['other']['bases']
                updates.append((ref, categories))
        
        if updates:
            dfUpdate = pd.DataFrame.from_dict(dict(updates), orient='index')
            dfUpdate.index.name = 'name'
            
            mainDct['contigDF'] = mainDct['contigDF'].merge(dfUpdate, on='name', how='left', suffixes=('', '_new'))
            
            for col in dfUpdate.columns:
                new_col_name = f'{col}_new'
                if new_col_name in mainDct['contigDF'].columns:
                    mainDct['contigDF'][col] = mainDct['contigDF'][col].combine_first(mainDct['contigDF'][new_col_name])
                    mainDct['contigDF'].drop(columns=[new_col_name], inplace=True)
                else:
                    if mainDct['debug']:
                        print(f"Warning: Column '{new_col_name}' not found in contigDF after merge.")
        
        if mainDct['debug']:
            print(f'results took {time.perf_counter() - timer:.2f} seconds')
        return

                    