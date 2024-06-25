import os
import pandas as pd
import sys
import plotnine as p9
import warnings
import numpy as np

class plots:
    def __init__(self, main):
        self.assembly      = pd.read_csv(main.ASSEMBLY_FILE, sep=',')
        self.assembly      = self.assembly.tail(1)
        self.contigs       = pd.read_csv(main.CONTIG_FILE, sep=',') if main.READMODE != 0 else None
        self.assembly_name = main.ASSEMBLY_NAME
        self.output        = main.TRANSRATE_PATH + '/plots'
        self.colors        = ['#2b9eb3', '#f8333c', '#fcab10', '#44af69', '#6a4c94']

    def run(self, main):
        if main.READMODE == 0:
            return
        self.plot_assembly_stats(main)
        if main.READMODE == 2:
            self.plot_score()
            self.plot_score_vs_scnuc()
            self.plot_score_vs_sccov()
            self.plot_combined_violins(main)
            self.plot_score_vs_scstats(main)
            self.plot_score_vs_scord()
            self.plot_score_vs_scseg()

    def plot_score(self):
        p = (p9.ggplot(self.contigs, p9.aes(x='score')) +
            p9.geom_histogram(bins=100, fill=self.colors[4], boundary=0) +
            p9.scale_x_continuous(breaks=[i/10 for i in range(11)]) +
            p9.labs(title=self.assembly_name + ' - Score Distribution') +
            p9.xlab('Contig Score') +
            p9.ylab('Count') +
            p9.scale_x_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.scale_y_continuous(limits=[0, 4500], breaks=[i * 500 for i in range(10)]) +
            p9.geom_text(p9.aes(label='stat(count)'), stat='bin', bins=100, size=8, va='bottom', format_string=' {:.0f}', angle=90, boundary=0) +
            p9.theme_light() +
            p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))
        self.print_plot(p, 'score.png')

    def plot_score_vs_scnuc(self):
        p = (p9.ggplot(self.contigs, p9.aes(x='score', y='sCnuc')) +
            p9.geom_point(size=0.5, color=self.colors[0]) +
            p9.labs(title=self.assembly_name + ' - sCnuc', subtitle='The proportion of nucleotides in the mapped reads that are the same as those in the assembled contig.') +
            p9.xlab('Contig Score') +
            p9.theme_light() +
            p9.scale_x_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.scale_y_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))
        self.print_plot(p, 'score_vs_scnuc.png')

    def plot_score_vs_sccov(self):
        p = (p9.ggplot(self.contigs, p9.aes(x='score', y='sCcov')) +
            p9.geom_point(size=0.5, color=self.colors[1]) +
            p9.labs(title=self.assembly_name + ' - sCcov', subtitle='The proportion of nucleotides in the contig that have supporting read data.') +
            p9.xlab('Contig Score') +
            p9.theme_light() +
            p9.scale_x_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.scale_y_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))
        self.print_plot(p, 'score_vs_sccov.png')

    def plot_score_vs_scord(self):
        p = (p9.ggplot(self.contigs, p9.aes(x='score', y='sCord')) +
            p9.geom_point(size=0.5, color=self.colors[2]) +
            p9.labs(title=self.assembly_name + ' - sCord', subtitle='The extent to which the order of the bases in the contig are correct by analyzing the pairing information in the mapped reads.') +
            p9.xlab('Contig Score') +
            p9.theme_light() +
            p9.scale_x_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.scale_y_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))
        self.print_plot(p, 'score_vs_scord.png')

    def plot_score_vs_scseg(self):
        p = (p9.ggplot(self.contigs, p9.aes(x='score', y='sCseg')) +
            p9.geom_point(size=0.5, color=self.colors[3]) +
            p9.labs(title=self.assembly_name + ' - sCseg', subtitle='The probability that the coverage depth of the transcript is univariate.') +
            p9.xlab('Contig Score') +
            p9.theme_light() +
            p9.scale_x_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.scale_y_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
            p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))
        self.print_plot(p, 'score_vs_scseg.png')
    
    def plot_score_vs_scstats(self, main):
        scores = pd.DataFrame({'score': self.contigs['score']})

        sCnuc = pd.DataFrame({'sCnuc': self.contigs['sCnuc'].clip(lower=0.01) / 
                                (self.contigs['sCnuc'].clip(lower=0.01) + 
                                self.contigs['sCcov'].clip(lower=0.01) + 
                                self.contigs['sCseg'].clip(lower=0.01) if main.READMODE == 2 else 0 + 
                                self.contigs['sCord'].clip(lower=0.01) if main.READMODE == 2 else 0)})
        sCcov = pd.DataFrame({'sCcov': self.contigs['sCcov'].clip(lower=0.01) /
                                (self.contigs['sCnuc'].clip(lower=0.01) + 
                                self.contigs['sCcov'].clip(lower=0.01) + 
                                self.contigs['sCseg'].clip(lower=0.01) if main.READMODE == 2 else 0 + 
                                self.contigs['sCord'].clip(lower=0.01) if main.READMODE == 2 else 0)})
        sCord = pd.DataFrame({'sCord': self.contigs['sCord'].clip(lower=0.01) /
                                (self.contigs['sCnuc'].clip(lower=0.01) + 
                                self.contigs['sCcov'].clip(lower=0.01)  + 
                                self.contigs['sCseg'].clip(lower=0.01)  + 
                                self.contigs['sCord'].clip(lower=0.01))})
        sCseg = pd.DataFrame({'sCseg': self.contigs['sCseg'].clip(lower=0.01) /
                                (self.contigs['sCnuc'].clip(lower=0.01) + 
                                self.contigs['sCcov'].clip(lower=0.01)  + 
                                self.contigs['sCseg'].clip(lower=0.01)  + 
                                self.contigs['sCord'].clip(lower=0.01))})
        data = [scores, sCnuc, sCcov, sCord, sCseg] if main.READMODE == 2 else [scores, sCnuc, sCcov]
        stats = pd.concat(data, axis=1)
        stats = pd.melt(stats, id_vars='score', var_name='stat', value_name='value')
        p = (p9.ggplot(stats, p9.aes(x='score', y='value')) +
                p9.geom_col(p9.aes(fill='stat'), position='fill', width=0.01) +
                p9.labs(title=self.assembly_name + ' - Scores by Proportion of sC Stats') +
                p9.xlab('Contig Score') +
                p9.ylab('Proportion of Score') +
                p9.scale_x_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
                p9.scale_y_continuous(breaks=[i/20 for i in range(21)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
                p9.scale_fill_manual(values=self.colors, guide=p9.guide_legend(title='Stat')) +
                p9.theme_light() +
                p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))

        self.print_plot(p, 'score_vs_scstats.png')

    def plot_assembly_stats(self, main):
        stats = ['p_softclipped', 'p_fragments_mapped', 
                 'p_good_mapping', 'p_bases_uncovered', 
                 'p_contigs_uncovbase', 'p_contigs_uncovered', 
                 'p_contigs_lowcovered', 'p_contigs_segmented',
                 'geometric_score', 'harmonic_score', 'optimal_score',
                 'sCnuc_Harmonic', 'sCcov_Harmonic', 'sCord_Harmonic', 'sCseg_Harmonic',
                 'sCnuc_Geometric', 'sCcov_Geometric', 'sCord_Geometric', 'sCseg_Geometric'] if main.READMODE == 2 else \
                ['p_softclipped', 'p_fragments_mapped',
                 'p_bases_uncovered', 'p_contigs_uncovbase', 
                 'p_contigs_uncovered', 'p_contigs_lowcovered',
                 'sCnuc_Harmonic', 'sCcov_Harmonic',
                 'sCnuc_Geometric', 'sCcov_Geometric']

        assembly_stats = self.assembly[stats]
        assembly_stats = pd.melt(assembly_stats, var_name='stat', value_name='value')

        p = (p9.ggplot(assembly_stats, p9.aes(x='stat', y='value')) +
                p9.geom_col(fill=self.colors[4], size=0.5) +
                p9.scale_y_continuous(breaks=[i/10 for i in range(11)], limits=[0, 1], minor_breaks=[i/100 for i in range(101)]) +
                p9.geom_text(p9.aes(label='value', y='value'), size=8, va='bottom') +
                p9.labs(title=self.assembly_name + ' - Assembly Stats') +
                p9.xlab('Assembly Statistic') +
                p9.ylab('Value') +
                p9.theme_light() +
                p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)))
        
        self.print_plot(p, 'assembly_stats.png')
    
    def plot_combined_violins(self, main):
        if main.READMODE == 2:
            data = self.contigs[['sCnuc', 'sCcov', 'sCord', 'sCseg']]
        elif main.READMODE == 1:
            data = self.contigs[['sCnuc', 'sCcov']]
        p = (p9.ggplot(pd.melt(data, var_name='stat', value_name='value'), p9.aes(x='stat', y='value')) +
            p9.geom_violin(p9.aes(fill='stat')) +
            p9.labs(title=self.assembly_name + ' - sCstats Violin Plots') +
            p9.xlab('sCstat') +
            p9.ylab('Density of Values') +
            p9.theme_light() +
            p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1)) +
            p9.scale_fill_manual(values=self.colors, guide=p9.guide_legend(title='Stat')))

        self.print_plot(p, 'scstats_violins.png')

    def print_plot(self, p, name='plot.png'):
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        p.save(os.path.join(self.output, f'{self.assembly_name}_{name}'), dpi=300, height=10, width=20)
