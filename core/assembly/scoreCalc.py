import warnings
warnings.filterwarnings('ignore', category=UserWarning)
import numpy as np
import pandas as pd
import math

class ScoreCalc:
    def __init__(self):
        return
    
    def mainRun(self, args):
        assembly_data = args
        self.scoreDct = {'assembly': {
                                    'score'       : 0,
                                    'optimalScore': 0,
                                    'cutoff'      : 0,
                                    'weighted'    : 0
                                    },
                         'contigs': {ref: {
                                    'score'       : 0,
                                    'sCnuc'       : 0,
                                    'sCcov'       : 0,
                                    'sCord'       : 0,
                                    'sCseg'       : 0
                                    } 
                                    for ref in assembly_data['refList']}}
        contigDF = assembly_data['contigDF']
        
        if assembly_data['mode'] == 2:
            self.good      = contigDF.set_index('name')['good'].to_dict()
            self.goodTotal = sum(self.good.values())
            self.scores(assembly_data, contigDF)
            self.optimal_score(assembly_data)
            self.weighted_score(assembly_data, contigDF)
        else:
            contigDF.apply(lambda row: self.update_contig_scores(row), axis=1)
        
        return self.scoreDct
    
    def update_contig_scores(self, row):
        ref_name = row['name']
        self.scoreDct['contigs'][ref_name]['sCnuc'] = self.sCnuc(row)
        self.scoreDct['contigs'][ref_name]['sCcov'] = self.sCcov(row)
    
    def scores(self, assembly_data, contigDF):
        scoreList = contigDF.apply(lambda row: self.calculate_scores(row), axis=1).tolist()
        self.rawScore(scoreList, assembly_data)
        return scoreList

    def calculate_scores(self, row):
        sCnuc = self.sCnuc(row)
        sCcov = self.sCcov(row)
        sCord = self.sCord(row)
        sCseg = self.sCseg(row)
        prod  = sCnuc * sCcov * sCord * sCseg
        s     = max(prod, 1e-2)
        ref   = row['name']
        self.scoreDct['contigs'][ref].update({'sCnuc': sCnuc, 'sCcov': sCcov, 'sCord': sCord, 'sCseg': sCseg, 'score': s})
        return s

    def rawScore(self, scoreList, assembly_data):
        harm_mean = self.harmMean(scoreList)
        if assembly_data['readCount'] > 0:
            read_support = self.goodTotal / assembly_data['readCount']
            self.scoreDct['assembly']['score'] = harm_mean * read_support
        else:
            self.scoreDct['assembly']['score'] = harm_mean
    
    def harmMean(self, score):
        f = np.array([s for s in score if s > 0])
        return (len(f) / np.sum(1.0/f)) if f.size > 0 else 0

    def sCnuc(self, row):
        return max(row['pSeqTrue'], 1e-2)
    
    def sCcov(self, row):
        return max(row['pBasesCovered'], 1e-2)
    
    def sCord(self, row):
        return max(row['pGood'], 1e-2)
    
    def sCseg(self, row):
        return max(row['pNotSegmented'], 1e-2)
    
    def optimal_score(self, assembly_data):
        opt_scores   = [self.scoreDct['contigs'][ref]['score'] for ref in self.scoreDct['contigs']]
        opt_good     = self.goodTotal
        cutoffscore  = {}
        score_sorted = sorted(((self.scoreDct['contigs'][ref]['score'], self.good[ref]) for ref in self.scoreDct['contigs']), key=lambda x: x[0])
        
        for score, good in score_sorted:
            opt_scores.remove(score)
            opt_good -= good
            opt_count = len(opt_scores)
            
            if opt_count > 0 and assembly_data['readCount'] > 0:
                harm_mean    = self.harmMean(opt_scores)
                read_support = opt_good / assembly_data['readCount']
                new_score    = harm_mean * read_support
            else:
                new_score = 0
            
            cutoffscore[score] = new_score
        
        optimal, cutoff = max(cutoffscore.items(), key=lambda x: x[1])
        self.scoreDct['assembly'].update({'optimalScore': optimal, 'cutoff': cutoff})
        pd.DataFrame(list(cutoffscore.items()), columns=['cutoff', 'score']).round(5).to_csv(assembly_data['scoreOptCSV'], index=False)
    
    def weighted_score(self, assembly_data, contigDF):
        wscore = contigDF.apply(lambda row: self.scoreDct['contigs'][row['name']]['score'] * row['tpm'], axis=1)
        pWeighted = wscore.mean() if not wscore.empty else 0
        self.scoreDct['assembly']['weighted'] = pWeighted * (self.goodTotal / assembly_data['readCount'])