import pandas as pd
import os
import math
import numpy as np
import sys
import time

class Assembly:
    def __init__(self, main):
        #### Setup ####
        self.CSV      = main.TR_CSV
        self.CONTIG_F = main.CONTIG_FILE
        self.FA       = main.ASSEMBLY
        self.SNAP     = main.READ_COUNT
        self.OUTDIR   = main.OUTPUT
        self.BASE     = main.ASSEMBLY_NAME
        self.DCT      = main.TR_DCT
        self.HEADERS  = []
        self.CHEADERS = []
        self.MULTI    = main.ASSEMBLY_MULTI
        self.SINGLE   = True if main.READMODE == 1 else False
        self.SOLO     = True if not main.READMODE else False
        self.STAR     = True if main.ALIGNER == 'star' else False
        self.BT2      = True if main.ALIGNER == 'bt2' else False
        self.SQ       = main.SALMON_QUANT
        # self.SALMONF  = pd.read_csv(os.path.join(main.SALMON_PATH,'quant.sf'), sep='\t')


        #### Dataframes ####
        self.OUTPUT_DF = None
        self.CONTIG_DF = None

        #### Metrics ####
        self.fragments   = main.READ_COUNT

        #### Iter / ORF ####
        self.frame       = ''
        self.bases       = 0
        self.gc          = 0
        self.read_length = 0
        self.lengths     = []
        self.first       = True
        self.seqs        = 0
        self.n_count     = 0
        self.n_prop      = 0
        self.n_under_200 = 0
        self.n_over_1k   = 0
        self.n_over_10k  = 0
        self.n_with_orf  = 0
        self.orf_len_sum = 0
        self.lengths     = []

        #### Assembly Stats ####
        self.n_prop         = 0
        self.gc             = 0
        self.smallest       = 0
        self.largest        = 0
        self.mean           = 0
        self.median         = 0
        self.std            = 0
        self.lengths_sorted = []
        self.base_90        = 0
        self.base_70        = 0
        self.base_50        = 0
        self.base_30        = 0
        self.base_10        = 0
        self.n90            = 0
        self.n70            = 0
        self.n50            = 0
        self.n30            = 0
        self.n10            = 0
        self.b90            = False
        self.b70            = False
        self.b50            = False
        self.b30            = False
        self.b10            = False
        self.base_total     = 0

        #### Contig Stats ####
        self.prop_gc_contig      = 0
        self.orf_length_contig   = 0
        self.p_good_contig       = 0
        self.p_bases_covered_contig = 0
        self.good_dct_contigs = {}
        self.total_length = 0
        self.total_softclipped = 0
        self.total_softclipped_length = 0
        self.p_softclipped = 0.0
        self.p_softclipped_length = 0.0

        #### Salmon Stats ####
        self.SALMONF            = ''
        self.contigs_uncovered  = 0
        self.contigs_lowcovered = 0

        #### Assembly Stats ####
        self.assembly             = ''
        self.fragments_mapped     = 0
        self.p_fragments_mapped   = 0
        self.good_mappings        = 0
        self.p_good_mappings      = 0
        self.bad_mappings         = 0
        self.potential_bridges    = 0
        self.bases_uncovered      = 0
        self.p_bases_uncovered    = 0
        self.p_contigs_uncovered  = 0
        self.p_contigs_lowcovered = 0
        self.contigs_segmented    = 0
        self.p_contigs_segmented  = 0
        self.contig_uncovbase     = 0
        self.p_contig_uncovbase   = 0
        self.score               = []
        self.goodcontig           = 0
        self.goodcontiglist       = []

    def run(self):
        self.BASE = ''.join(os.path.basename(self.FA).split('.')[:-1])
        self.headers()
        self.starsnap() if not self.SOLO else None
        self.file_iter()
        self.basestats()
        self.nstats()
        self.salmonstats() if not self.SOLO else None
        self.assemblystats() if not self.SOLO else None
        self.score_run() if not self.SOLO and not self.SINGLE else None
        self.sc_metrics() if not self.SOLO else None
        self.softclip() if not self.SOLO else None
        self.assembly_output()
        self.contig_output() if not self.SOLO else None

        while not os.path.exists(self.outfile):
            pass
        while not os.path.exists(self.CONTIG_F) and not self.SOLO:
            pass

    def headers(self):
        if self.SINGLE:
            self.HEADERS = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                            'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                            'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10', 
                            'gc', 'bases_n', 'proportion_n', 'fragments', 'fragments_mapped', 
                            'p_fragments_mapped', 'softclipped', 'p_softclipped', 'p_softclipped_length', 'bases_uncovered', 'p_bases_uncovered', 
                            'contigs_uncovbase', 'p_contigs_uncovbase', 'contigs_uncovered', 
                            'p_contigs_uncovered', 'contigs_lowcovered', 'p_contigs_lowcovered', 
                            'sCnuc_Harmonic', 'sCcov_Harmonic', 'sCnuc_Geometric', 'sCcov_Geometric']
        elif self.SOLO:
            self.HEADERS = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                            'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                            'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10',
                            'gc', 'bases_n', 'proportion_n']
        else:
            self.HEADERS = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                            'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                            'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10', 
                            'gc', 'bases_n', 'proportion_n', 'fragments', 'fragments_mapped', 
                            'p_fragments_mapped', 'softclipped', 'p_softclipped',
                            'p_softclipped_length', 'good_mappings', 'p_good_mapping', 'bad_mappings', 
                            'potential_bridges', 'bases_uncovered', 'p_bases_uncovered', 
                            'contigs_uncovbase', 'p_contigs_uncovbase', 'contigs_uncovered', 
                            'p_contigs_uncovered', 'contigs_lowcovered', 'p_contigs_lowcovered', 
                            'contigs_segmented', 'p_contigs_segmented', 'geometric_score', 'harmonic_score', 'optimal_score', 
                            'cutoff', 'weighted', 'goodcontig', 'sCnuc_Harmonic', 'sCcov_Harmonic',
                            'sCord_Harmonic', 'sCseg_Harmonic', 'sCnuc_Geometric', 'sCcov_Geometric',
                            'sCord_Geometric', 'sCseg_Geometric']
            
        if not self.SOLO:
            self.CHEADERS = [
            'name', 'length', 'prop_gc', 'orf_length', 'bridges', 'properpair',
            'fragments', 'both_mapped', 'good', 'p_good', 'basesuncovered', 
            'p_bases_covered', 'p_seqtrue', 'softclipped', 'p_softclipped', 'p_softclipped_length', 'score', 'p_notsegmented', 'eff_length', 
            'eff_count', 'tpm', 'coverage', 'sCnuc', 'sCcov', 'sCord', 'sCseg'
            ]
            self.CONTIG_DF = pd.read_csv(self.CSV, sep=',')
            missing_columns = list(set(self.CHEADERS) - set(self.CONTIG_DF.columns))
            for column in missing_columns:
                self.CONTIG_DF[column] = 0
        elif not self.SINGLE:
            self.CHEADERS = [
            'name', 'length', 'prop_gc', 'orf_length',
            'fragments', 'softclipped', 'p_softclipped', 'p_softclipped_length', 
            'both_mapped', 'basesuncovered', 'p_bases_covered', 'p_seqtrue', 
            'eff_length', 'eff_count', 'tpm', 'coverage', 'sCnuc', 'sCcov',
            ]

        self.OUTPUT_DF = pd.DataFrame(columns=self.HEADERS) 

    def starsnap(self):
        return
        # if self.STAR:
        #     for line in open(self.SNAP, 'r'):
        #         if 'Number of input reads' in line:
        #             self.fragments = int(line.split('	')[1])
        # elif self.BT2:
        #     with open(self.SNAP, 'r') as file:
        #         lines = file.readlines()
        #         self.fragments = int((lines[0].split(' ')[0]))
        # else:
        #     with open(self.SNAP, 'r') as file:
        #         lines = file.readlines()
        #         if self.SINGLE:
        #             self.fragments = int(lines[3].strip().split(' ')[0].replace(',', ''))
        #         else:
        #             self.fragments = int(int(lines[3].strip().split(' ')[0].replace(',', '')) / 2)

    def longest_orf(self):
        longest  = 0
        len_list = [0, 0, 0]
        sl       = len(self.frame)
        str      = self.frame
        for i in range(sl - 2):
            if str[i:i+3] == 'atg':
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 1
            elif str[i:i+3] in ['tag', 'taa', 'tga']:
                longest = max(longest, len_list[i % 3])
                len_list[i % 3] = -1
            else:
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 0
        longest  = max(longest, max(len_list))
        len_list = [0, 0, 0]
        for i in range(sl - 1, 1, -1):
            if str[i-2:i+1] == 'cat':
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 1
            elif str[i-2:i+1] in ['cta', 'tta', 'tct']:
                longest = max(longest, len_list[i % 3])
                len_list[i % 3] = -1
            else:
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 0
        return max(longest, max(len_list))
    
    def orf(self):
        self.lengths.append(self.read_length)
        self.n_under_200       += self.read_length < 200
        self.n_over_1k         += self.read_length > 1000
        self.n_over_10k        += self.read_length > 10000
        orf                     = self.longest_orf()
        self.orf_length_contig  = orf
        self.n_with_orf        += orf > 149
        self.orf_len_sum       += orf
        self.read_length        = 0

    def file_iter(self):
        with open(self.FA, 'r') as file:
            name = ''
            for line in file:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    self.seqs += 1
                    if not self.first:
                        self.orf()
                        self.good_dct_contigs[name] = {'orf_len': 0, 'prop_gc': 0}
                        self.prop_gc_contig = round((self.frame.lower().count('g') + self.frame.lower().count('c')) / len(self.frame), 4)
                        self.good_dct_contigs[name]['prop_gc'] = self.prop_gc_contig
                        self.good_dct_contigs[name]['orf_len'] = self.orf_length_contig
                        self.frame        = ''
                    name = line.split(' ')[0].strip('>')
                    self.first = False
                    continue
                self.read_length += len(line)
                self.bases       += len(line)
                self.gc          += line.lower().count('g') + line.lower().count('c')
                self.n_count     += line.lower().count('n')
                self.frame       += line

            self.orf()
            self.good_dct_contigs[name] = {'orf_len': 0, 'prop_gc': 0}
            self.good_dct_contigs[name]['orf_len'] = self.orf_length_contig
            self.prop_gc_contig = round((self.frame.lower().count('g') + self.frame.lower().count('c')) / len(self.frame), 4)
            self.good_dct_contigs[name]['prop_gc'] = self.prop_gc_contig

    def basestats(self):
        self.n_prop           = round(self.n_count/self.bases, 4)
        self.gc               = round(self.gc/self.bases, 4)
        self.smallest         = min(self.lengths)
        self.largest          = max(self.lengths)
        self.mean             = round(sum(self.lengths)/len(self.lengths), 4)
        self.median           = self.lengths[int(len(self.lengths)/4)]
        self.std              = round((sum([(x - self.mean)**2 for x in self.lengths]) / len(self.lengths))**0.5, 4)
        self.lengths_sorted   = sorted(self.lengths)
        self.base_90          = int(self.bases * 0.9)
        self.base_70          = int(self.bases * 0.7)
        self.base_50          = int(self.bases * 0.5)
        self.base_30          = int(self.bases * 0.3)
        self.base_10          = int(self.bases * 0.1)
        self.mean_orf_percent = round((300 * self.orf_len_sum) / (self.seqs * self.mean), 4)

    def nstats(self):
        for length in self.lengths_sorted[::-1]:
            self.base_total += length
            if self.base_total > self.base_90 and self.b90 == False:
                self.n90 = length
                self.b90 = True
            if self.base_total > self.base_70 and self.b70 == False:
                self.n70 = length
                self.b70 = True
            if self.base_total > self.base_50 and self.b50 == False:
                self.n50 = length
                self.b50 = True
            if self.base_total > self.base_30 and self.b30 == False:
                self.n30 = length
                self.b30 = True
            if self.base_total > self.base_10 and self.b10 == False:
                self.n10 = length
                self.b10 = True

    def salmonstats(self):
        self.SALMONF                          = pd.read_csv(self.SQ, sep='\t')
        self.CONTIG_DF['eff_length']          = self.SALMONF['EffectiveLength'].fillna(0)
        self.CONTIG_DF['eff_count']           = self.SALMONF['NumReads'].fillna(0)
        self.CONTIG_DF['tpm']                 = self.SALMONF['TPM'].fillna(0).astype(float)
        self.CONTIG_DF['eff_length']          = self.CONTIG_DF['eff_length'].where(self.CONTIG_DF['eff_length'] != 0, 0)
        self.CONTIG_DF['eff_count']           = self.CONTIG_DF['eff_count'].where(self.CONTIG_DF['eff_count'] != 0, 0)
        self.CONTIG_DF['tpm']                 = self.CONTIG_DF['tpm'].where(self.CONTIG_DF['tpm'] != 0, 0)
        self.CONTIG_DF['eff_length']          = self.CONTIG_DF['eff_length'].where(self.SALMONF['Length'] != 0, 0)
        coverage                              = self.CONTIG_DF['eff_count'] * self.CONTIG_DF['length'] / self.CONTIG_DF['eff_length']
        self.contigs_uncovered               += (coverage < 1).sum()
        self.contigs_lowcovered              += (coverage < 10).sum()

    def assemblystats(self):
        self.assembly                         = self.BASE
        self.fragments_mapped                 = self.CONTIG_DF['fragments'].sum()
        self.p_fragments_mapped               = round(self.fragments_mapped / self.fragments, 4)
        self.bases_uncovered                  = self.CONTIG_DF['basesuncovered'].sum()
        self.p_bases_uncovered                = round(self.bases_uncovered / self.bases, 4)
        self.p_contigs_uncovered              = round(self.contigs_uncovered / self.CONTIG_DF.shape[0], 4)
        self.p_contigs_lowcovered             = round(self.contigs_lowcovered / self.CONTIG_DF.shape[0], 4)
        self.contig_uncovbase                 = round((self.CONTIG_DF['basesuncovered'] > 0).sum(), 4)
        self.p_contig_uncovbase               = round(self.contig_uncovbase / self.CONTIG_DF.shape[0], 4)
        if not self.SINGLE:
            self.good_mappings                = self.CONTIG_DF['good'].sum()
            self.p_good_mappings              = round(self.good_mappings / self.fragments, 4)
            self.bad_mappings                 = round(self.fragments_mapped - self.good_mappings, 4)
            self.potential_bridges            = (self.CONTIG_DF['bridges'] > 0).sum()
            self.contigs_segmented            = round((self.CONTIG_DF['p_notsegmented'] < 0.5).sum(), 4)
            self.p_contigs_segmented          = round(self.contigs_segmented / self.CONTIG_DF.shape[0], 4)
    
    def softclip(self):
        for index,row in self.CONTIG_DF.iterrows():
            self.CONTIG_DF.loc[index, 'softclipped']          = self.DCT[row['name']]['seq']['softclipped']
            self.CONTIG_DF.loc[index, 'p_softclipped']        = self.DCT[row['name']]['seq']['p_softclipped']
            self.CONTIG_DF.loc[index, 'p_softclipped_length'] = self.DCT[row['name']]['seq']['p_softclipped_length']
            self.total_softclipped                           += 1 if self.DCT[row['name']]['seq']['softclipped'] > 0 else 0
            self.total_softclipped_length                    += self.DCT[row['name']]['seq']['softclipped_length']
        try:
            self.p_softclipped                                = self.total_softclipped / self.fragments
            self.p_softclipped_length                         = self.total_softclipped_length / self.bases
        except:
            pass
        return

    def scores(self):
        score = []
        for index, row in self.CONTIG_DF.iterrows():
            prod = np.max([(row['length'] - row['basesuncovered'])/row['length'], 0.01]) * \
                    np.max([row['p_notsegmented'], 0.01]) * \
                    np.max([row['good'] / row['fragments'] if row['fragments'] > 1 else 0.01, 0.01]) * \
                    np.max([row['p_seqtrue'], 0.01])
            s = np.max([prod, 0.01])
            score.append(s)
        return score
    
    def raw_score(self, harmonic=False):
        contig_score = self.geomean(self.score) if not harmonic else self.harmmean(self.score)
        return contig_score * (self.good_mappings / self.fragments)
    
    def geomean(self, x):
        return math.exp(sum([math.log(max(v, 0.01)) for v in x]) / len(x))
    
    def harmmean(self, x):
        return len(x) / sum([1.0 / v if v != 0 else 1.0 / 0.01 for v in x])

    def optimal_score(self):
        opt_product   = 0
        opt_good      = 0
        opt_product   = self.CONTIG_DF['score'].apply(math.log).sum()
        opt_good      = self.CONTIG_DF['good'].sum()
        opt_count     = len(self.CONTIG_DF)
        cutoffscore   = {}
        score_sorted  = self.CONTIG_DF[['score', 'good']].sort_values(by='score')
        for index, row in score_sorted.iterrows():
            opt_product -= math.log(row['score'])
            opt_good    -= row['good']
            opt_count   -= 1 if opt_count > 1 else 0
            new_score    = math.exp(opt_product / opt_count) * (opt_good / self.fragments)
            cutoffscore[row['score']] = new_score
        optimal = 0
        cutoff  = 0
        for c, score in cutoffscore.items():
            if score > optimal:
                optimal = score
                cutoff = c
        cutoffscore_df = pd.DataFrame(list(cutoffscore.items()), columns=['cutoff', 'score'])
        cutoffscore_df = cutoffscore_df.round(5)
        scorefile      = 'assembly_score_optimisation.csv' if not self.MULTI else f'{self.BASE}_assembly_score_optimisation.csv'
        cuttofffile    = os.path.join(self.OUTDIR, 'transrate', scorefile)
        cutoffscore_df.to_csv(cuttofffile, index=False)
        return [round(optimal,6), round(cutoff,6)]
    
    def weighted_score(self):
        wscore = []
        for index, row in self.CONTIG_DF.iterrows():
            wscore.append(row['score'] * row['tpm'])
        contig_weighted = sum(wscore) / len(wscore)
        return round(contig_weighted * (self.good_mappings / self.fragments), 4)
    
    def score_run(self):
        self.score = self.scores()
        self.CONTIG_DF['score'] = [round(score, 4) for score in self.score]
        self.geo_raw_score_val  = self.raw_score(False)
        self.har_raw_score_val  = self.raw_score(True)
        self.optimal_score_val  = self.optimal_score() 

    def sc_metrics(self):
        for index, row in self.CONTIG_DF.iterrows():
            try:    
                self.CONTIG_DF.loc[index, 'p_good']          = round(row['good'] / row['fragments'], 4)
            except:
                self.CONTIG_DF.loc[index, 'p_good']          = 0.0
            try:
                self.CONTIG_DF.loc[index, 'p_bases_covered'] = round((row['length'] - row['basesuncovered']) / row['length'], 4)
            except:
                self.CONTIG_DF.loc[index, 'p_bases_covered'] = 0.0
            try:
                self.CONTIG_DF.loc[index, 'sCnuc']           = row['p_seqtrue']
            except:
                self.CONTIG_DF.loc[index, 'sCnuc']           = 0.0
            try:
                self.CONTIG_DF.loc[index, 'sCcov']           = round((row['length'] - row['basesuncovered']) / row['length'], 4)
            except:
                self.CONTIG_DF.loc[index, 'sCcov']           = 0.0
            try:
                self.CONTIG_DF.loc[index, 'sCseg']           = round(row['p_notsegmented'], 4)
            except:
                self.CONTIG_DF.loc[index, 'sCseg']           = 0.0
            try:
                self.CONTIG_DF.loc[index, 'coverage']        = round(row['eff_count'] * self.DCT[row['name']]['seq']['max_read_length'] / row['eff_length'], 4)
            except:
                self.CONTIG_DF.loc[index, 'coverage']        = 0.0
            try:
                self.CONTIG_DF.loc[index, 'orf_length']      = self.good_dct_contigs[row['name']]['orf_len']
            except:
                self.CONTIG_DF.loc[index, 'orf_length']      = 0
            try:
                self.CONTIG_DF.loc[index, 'prop_gc']         = self.good_dct_contigs[row['name']]['prop_gc']
            except:
                self.CONTIG_DF.loc[index, 'prop_gc']         = 0
            if not self.SINGLE:
                try:
                    self.CONTIG_DF.loc[index, 'sCord']       = round(row['good'] / row['fragments'], 4)
                except:
                    self.CONTIG_DF.loc[index, 'sCord']       = 0.0

    def good_contigs(self):
        for index, row in self.CONTIG_DF.iterrows():
            if self.CONTIG_DF.loc[index, 'score'] > self.optimal_score_val[1]:
                self.goodcontig += 1
                self.goodcontiglist.append(self.CONTIG_DF.loc[index, 'name'])
        try:
            with open(os.path.join(self.OUTDIR, 'transrate', f'good.{self.assembly}.fa'), 'w') as good:
                good.write('')
            with open(os.path.join(self.OUTDIR, 'transrate', f'bad.{self.assembly}.fa'), 'w') as bad:
                bad.write('')
        except:
            pass
        with open(self.FA, 'r') as file:
            fa = file.readlines()
            lineindex = 0
            for i, line in enumerate(fa):
                frame = []
                if line.startswith('>'):
                    if line.strip('\n').strip('>') in self.goodcontiglist:
                        frame.append(line.strip('\n'))
                        while True:
                            lineindex += 1
                            if lineindex == len(fa) or fa[lineindex].startswith('>'):
                                break
                            frame.append(fa[lineindex].upper().strip('\n'))
                        with open(os.path.join(self.OUTDIR, 'transrate', f'good.{self.assembly}.fa'), 'a') as good:
                            good.write('\n'.join(frame))
                            good.write('\n')
                    else:
                        frame.append(line.strip('\n'))
                        while True:
                            lineindex += 1
                            if lineindex == len(fa) or fa[lineindex].startswith('>'):
                                break
                            frame.append(fa[lineindex].upper().strip('\n'))
                        with open(os.path.join(self.OUTDIR, 'transrate', f'bad.{self.assembly}.fa'), 'a') as bad:
                            bad.write('\n'.join(frame))
                            bad.write('\n')


        if 'goodcontig' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'goodcontig'] = self.goodcontig

    def assembly_output(self):
        if 'assembly' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'assembly']             = self.BASE
        if 'n_seqs' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_seqs']               = self.seqs
        if 'smallest' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'smallest']             = self.smallest
        if 'largest' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'largest']              = self.largest
        if 'n_bases' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_bases']              = self.bases
        if 'mean_len' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'mean_len']             = round(self.mean, 4)
        if 'median_len' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'median_len']           = round(self.median, 4)
        if 'std_len' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'std_len']              = round(self.std, 4)
        if 'n_under_200' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_under_200']          = self.n_under_200
        if 'n_over_1k' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_over_1k']            = self.n_over_1k
        if 'n_over_10k' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_over_10k']           = self.n_over_10k
        if 'n_with_orf' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_with_orf']           = self.n_with_orf
        if 'mean_orf_percent' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'mean_orf_percent']     = round(self.mean_orf_percent, 4)
        if 'n90' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n90']                  = self.n90
        if 'n70' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n70']                  = self.n70
        if 'n50' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n50']                  = self.n50
        if 'n30' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n30']                  = self.n30
        if 'n10' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n10']                  = self.n10
        if 'gc' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'gc']                   = self.gc
        if 'bases_n' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'bases_n']              = self.n_count
        if 'proportion_n' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'proportion_n']         = round(self.n_prop, 4)
        if 'fragments' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'fragments']            = self.fragments
        if 'fragments_mapped' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'fragments_mapped']     = self.fragments_mapped
        if 'softclipped' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'softclipped']          = self.total_softclipped
        if 'p_softclipped' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_softclipped']        = round(self.p_softclipped, 4)
        if 'p_softclipped_length' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_softclipped_length'] = round(self.p_softclipped_length, 4)
        if 'p_fragments_mapped' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_fragments_mapped']   = round(self.p_fragments_mapped, 4)
        if 'good_mappings' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'good_mappings']        = self.good_mappings
        if 'p_good_mapping' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_good_mapping']       = round(self.p_good_mappings, 4)
        if 'bad_mappings' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'bad_mappings']         = self.bad_mappings
        if 'potential_bridges' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'potential_bridges']    = self.potential_bridges
        if 'bases_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'bases_uncovered']      = self.bases_uncovered
        if 'p_bases_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_bases_uncovered']    = round(self.p_bases_uncovered, 4)
        if 'contigs_uncovbase' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_uncovbase']    = self.contig_uncovbase
        if 'p_contigs_uncovbase' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_uncovbase']  = round(self.p_contig_uncovbase, 4)
        if 'contigs_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_uncovered']    = self.contigs_uncovered
        if 'p_contigs_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_uncovered']  = round(self.p_contigs_uncovered, 4)
        if 'contigs_lowcovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_lowcovered']   = self.contigs_lowcovered
        if 'p_contigs_lowcovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_lowcovered'] = round(self.p_contigs_lowcovered, 4)
        if 'contigs_segmented' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_segmented']    = self.contigs_segmented
        if 'p_contigs_segmented' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_segmented']  = round(self.p_contigs_segmented, 4)
        if 'geometric_score' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'geometric_score']      = round(self.geo_raw_score_val, 4)
        if 'harmonic_score' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'harmonic_score']       = round(self.har_raw_score_val, 4)
        if 'optimal_score' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'optimal_score']        = round(self.optimal_score_val[0], 4)
        if 'cutoff' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'cutoff']               = round(self.optimal_score_val[1], 4)
        if 'weighted' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'weighted']             = round(self.weighted_score(), 4)
        if 'sCnuc_Harmonic' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCnuc_Harmonic']       = round(self.harmmean(self.CONTIG_DF['sCnuc']), 4)
        if 'sCcov_Harmonic' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCcov_Harmonic']       = round(self.harmmean(self.CONTIG_DF['sCcov']), 4)
        if 'sCord_Harmonic' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCord_Harmonic']       = round(self.harmmean(self.CONTIG_DF['sCord']), 4)
        if 'sCseg_Harmonic' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCseg_Harmonic']       = round(self.harmmean(self.CONTIG_DF['sCseg']), 4)
        if 'sCnuc_Geometric' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCnuc_Geometric']      = round(self.geomean(self.CONTIG_DF['sCnuc']), 4)
        if 'sCcov_Geometric' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCcov_Geometric']      = round(self.geomean(self.CONTIG_DF['sCcov']), 4)
        if 'sCord_Geometric' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCord_Geometric']      = round(self.geomean(self.CONTIG_DF['sCord']), 4)
        if 'sCseg_Geometric' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCseg_Geometric']      = round(self.geomean(self.CONTIG_DF['sCseg']), 4)



        self.good_contigs() if not self.SOLO and not self.SINGLE else None

        self.outfile    = os.path.join(self.OUTDIR, 'transrate', 'assembly.csv')
        if os.path.exists(self.outfile) and self.MULTI:
            self.OUTPUT_DF.to_csv(self.outfile, index=False, mode='a', header=False)
        else:
            self.OUTPUT_DF.to_csv(self.outfile, index=False)

    def contig_output(self):
        # cfile = 'contigs.csv' if not self.MULTI else f'{self.assembly}.contigs.csv'
        self.CONTIG_DF = self.CONTIG_DF[self.CHEADERS]
        self.CONTIG_DF.to_csv(self.CONTIG_F, index=False)
