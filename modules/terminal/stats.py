import os
import pandas as pd
import sys

class Stats:
    def __init__(self, main):
        self.assemblies = pd.read_csv(main.ASSEMBLY_FILE).iloc[[-1,-1]].reset_index(drop=True)
        self.contigs = pd.read_csv(main.CONTIG_FILE) if main.READMODE != 0 else None

    def stats(self, main):
        clear = ' ' * 66
        print(clear)
        
        def contig_stats(self, main):
            contiglabel = f'{main.COLOR}  ┌{"─" * 28}\033[m  Contig Metrics  {main.COLOR}{"─" * 28}┐\033[m'
            contig = {
                    'n_seqs'              : '# Seqs',
                    'smallest'            : 'Smallest',
                    'largest'             : 'Largest',
                    'n_bases'             : '# Bases',
                    'mean_len'            : 'Mean Len',
                    'median_len'          : 'Median Len',
                    'std_len'             : 'Std Len',
                    'n_under_200'         : '# < 200',
                    'n_over_1k'           : '# > 1k',
                    'n_over_10k'          : '# > 10k',
                    'n_with_orf'          : '# ORF',
                    'mean_orf_percent'    : 'Mean ORF %',
                    'n90'                 : 'N90',
                    'n70'                 : 'N70',
                    'n50'                 : 'N50',
                    'n30'                 : 'N30',
                    'n10'                 : 'N10',
                    'gc'                  : 'GC',
                    'bases_n'             : '# N',
                    'proportion_n'        : 'p̂ N',
                    }
            print(contiglabel) # Print top label bar
            for k,v in self.assemblies.items(): # Print each contig metric
                try:
                    if 'p̂' in k or 'p̂' in contig[k]:
                        print(f' {main.COLOR}  •\033[m {contig[k]:<23}: {v[0]}')
                    else:
                        print(f' {main.COLOR}  •\033[m {contig[k]:<22}: {v[0]}')
                except:
                    pass

            print(f'{main.COLOR}  └{"─" * 74}┘\033[m') # Print bottom bar

        def readmap_stats(self, main):
            readmaplabel = f'{main.COLOR}  ┌{"─" * 28}\033[m   Read Mapping   {main.COLOR}{"─" * 28}┐\033[m'
            readmap = {
                        'fragments'           : '# Fragments',
                        'fragments_mapped'    : '# Fragments Mapped',
                        'p_fragments_mapped'  : 'p̂ Fragments Mapped',
                        'good_mappings'       : '# Good Mappings',
                        'p_good_mapping'      : 'p̂ Good Mappings',
                        'bad_mappings'        : '# Bad Mappings',
                        'potential_bridges'   : '# Potential Bridges',
                        'bases_uncovered'     : '# Bases Uncovered',
                        'p_bases_uncovered'   : 'p̂ Bases Uncovered',
                        'contigs_uncovbase'   : '# Contigs Uncovbase',
                        'p_contigs_uncovbase' : 'p̂ Contigs Uncovbase',
                        'contigs_uncovered'   : '# Contigs Uncovered',
                        'p_contigs_uncovered' : 'p̂ Contigs Uncovered',
                        'contigs_lowcovered'  : '# Contigs Lowcovered',
                        'p_contigs_lowcovered': 'p̂ Contigs Lowcovered',
                        'contigs_segmented'   : '# Contigs Segmented',
                        'p_contigs_segmented' : 'p̂ Contigs Segmented',
                        'softclipped'         : '# Softclipped',
                        'p_softclipped'       : 'p̂ Softclipped',
                        'p_softclipped_length': 'p̂ Softclipped Length',
                        }
            print(readmaplabel) # Print top label bar
            for k,v in self.assemblies.items(): # Print each read mapping metric
                try:
                    if 'p̂' in k or 'p̂' in readmap[k]:
                        print(f' {main.COLOR}  •\033[m {readmap[k]:<23}: {v[0]}')
                    else:
                        print(f' {main.COLOR}  •\033[m {readmap[k]:<22}: {v[0]}')
                except:
                    pass
            print(f'{main.COLOR}  └{"─" * 74}┘\033[m') # Print bottom bar
            
        def score_stats(self, main):
            scorelabel = f'{main.COLOR}  ┌{"─" * 28}\033[m  Quality Scores  {main.COLOR}{"─" * 28}┐\033[m'
            score = {
                        'optimal_score'       : 'Optimal Score',
                        'cutoff'              : 'Cutoff',
                        'weighted'            : 'Weighted',
                        'goodcontig'          : 'Good Contig',
                        'badcontig'           : 'Bad Contig',
            }
            if any(key in self.assemblies.keys() for key in score.keys()):
                print(scorelabel) # Print top label bar
                for k,v in self.assemblies.items(): # Print each score metric
                    try:
                        if type(v[0]) == float:
                            v[0] = f'{v[0]:.4f}'
                        if 'p̂' in k or 'p̂' in score[k]:
                            print(f' {main.COLOR}  •\033[m {score[k]:<23}: {v[0]}')
                        else:
                            print(f' {main.COLOR}  •\033[m {score[k]:<22}: {v[0]}')
                    except:
                        pass
                print(f'{main.COLOR}  └{"─" * 74}┘\033[m') # Print bottom bar

            scdict = {
                'score_h': ['Harmonic Score',   self.assemblies['harmonic_score'].values[0]],
                'score_g': ['Geometric Score',  self.assemblies['geometric_score'].values[0]],
                'scnuc_h': ['Avg. sCnuc',       self.assemblies['sCnuc_Harmonic'].values[0]],
                'sccov_h': ['Avg. sCcov',       self.assemblies['sCcov_Harmonic'].values[0]],
                'scord_h': ['Avg. sCord',       self.assemblies['sCord_Harmonic'].values[0]],
                'scseg_h': ['Avg. sCseg',       self.assemblies['sCseg_Harmonic'].values[0]],
                'scnuc_g': ['Avg. sCnuc',       self.assemblies['sCnuc_Geometric'].values[0]],
                'sccov_g': ['Avg. sCcov',       self.assemblies['sCcov_Geometric'].values[0]],
                'scord_g': ['Avg. sCord',       self.assemblies['sCord_Geometric'].values[0]],
                'scseg_g': ['Avg. sCseg',       self.assemblies['sCseg_Geometric'].values[0]],
            } if main.READMODE == 2 else {
                # 'score_h' : ['Harmonic Score',  self.assemblies['harmonic_score'].values[0]],
                # 'score_g' : ['Geometric Score', self.assemblies['geometric_score'].values[0]],
                'scnuc_h' : ['Avg. sCnuc',      self.assemblies['sCnuc_Harmonic'].values[0]],
                'sccov_h' : ['Avg. sCcov',      self.assemblies['sCcov_Harmonic'].values[0]],
                'scnuc_g' : ['Avg. sCnuc',      self.assemblies['sCnuc_Geometric'].values[0]],
                'sccov_g' : ['Avg. sCcov',      self.assemblies['sCcov_Geometric'].values[0]],
            } if main.READMODE == 1 else None
            hlabel = f'{main.COLOR}  ┌{"─" * 28}\033[m Harmonic  Scores {main.COLOR}{"─" * 28}┐\033[m'
            print(hlabel) # Print top label bar
            for k,v in scdict.items(): # Print each score metric
                if not '_h' in k:
                    continue
                print(f' {main.COLOR}  •\033[m {v[0]:<22}: {round(v[1],4):.4f}')
            print(f'{main.COLOR}  └{"─" * 74}┘\033[m') # Print bottom bar
            
            glabel = f'{main.COLOR}  ┌{"─" * 28}\033[m Geometric Scores {main.COLOR}{"─" * 28}┐\033[m'
            print(glabel) # Print top label bar
            for k,v in scdict.items(): # Print each score metric
                if not '_g' in k:
                    continue
                print(f' {main.COLOR}  •\033[m {v[0]:<22}: {round(v[1],4):.4f}')
            print(f'{main.COLOR}  └{"─" * 74}┘\033[m') # Print bottom bar


        def ref_stats(self, main):
            reflabel = f'{main.COLOR}  ┌{"─" * 25}\033[m  Reference Statistics  {main.COLOR}{"─" * 25}┐\033[m'
            clear = ' ' * 80
            assemblies = pd.read_csv(main.ASSEMBLY_FILE).to_dict()

            refs = {
                'CRBB_hits' : 'CRBB Hits',
                'n_contigs_with_CRBB' : '# Contigs with CRBB',
                'p_contigs_with_CRBB' : 'p̂ Contigs with CRBB',
                'rbh_per_reference' : 'RBH per Reference',
                'n_refs_with_CRBB' : '# Refs with CRBB',
                'p_refs_with_CRBB' : 'p̂ Refs with CRBB',
                'cov25' : 'Coverage 25',
                'p_cov25' : 'p̂ Coverage 25',
                'cov50' : 'Coverage 50',
                'p_cov50' : 'p̂ Coverage 50',
                'cov75' : 'Coverage 75',
                'p_cov75' : 'p̂ Coverage 75',
                'cov85' : 'Coverage 85',
                'p_cov85' : 'p̂ Coverage 85',
                'cov95' : 'Coverage 95',
                'p_cov95' : 'p̂ Coverage 95',
                'reference_coverage' : 'Reference Coverage',
            }
            print(reflabel) # Print top label bar
            for k,v in self.assemblies.items(): # Print each reference metric
                try:
                    if type(v[0]) == float:
                        v[0] = f'{v[0]:.4f}'
                    if 'p̂' in k or 'p̂' in refs[k]:
                        print(f' {main.COLOR}  •\033[m {refs[k]:<23}: {v[0]}')
                    else:
                        print(f' {main.COLOR}  •\033[m {refs[k]:<22}: {v[0]}')
                except:
                    pass
            print(f'{main.COLOR}  └{"─" * 74}┘\033[m') # Print bottom bar

        contig_stats(self, main) 
        readmap_stats(self, main) if main.READMODE != 0 else None
        score_stats(self, main) if main.READMODE != 0 else None
        ref_stats(self, main) if main.REFERENCE else None
