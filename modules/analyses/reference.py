from modules.analyses.hit import Hit
import sys
import subprocess
import os
import math
import csv
import pandas as pd

class REFERENCE:

    def __init__(self, main):        
        # file paths
        self.q2p_blast      = ''
        self.p2q_blast      = ''
        self.fasta_f        = main.ASSEMBLY
        self.prot_f         = main.REFERENCE_FILE
        self.reference      = {}
        # seq_count
        self.seqs           = 0
        self.query_dict     = {}
        self.p_seqs         = 0
        # load_outputs
        self.target_dict    = {}
        self.target         = ''
        self.query          = ''
        self.query_results  = {}
        self.target_results = {}
        self.q_count        = 0
        self.t_count        = 0
        # find_reciprocals
        self.reciprocals    = {}
        self.missed         = {}
        self.evalues        = []
        self.longest        = 0
        # find_secondaries
        self.length_hash    = {}
        self.fitting        = {}
        # comparitive metrics
        self.comp_stats     = {}
        self.assembly       = {}
        self.contig_hits    = 0
        self.has_crb        = {}

        self.assemblycsv    = main.ASSEMBLY_FILE
        self.output         = main.TRANSRATE_PATH
        self.fa_translated  = ''
        self.threads        = main.THREADS
        self.multi          = main.ASSEMBLY_MULTI
        self.rh_name        = ''.join(os.path.basename(self.fasta_f).split('.')[:-1])

    #### RUN ####
    def run(self):
        self.diamond()
        self.seq_count()
        self.load_outputs()
        self.find_reciprocals()
        self.find_secondaries()
        self.reference_hits()
        self.rbh()
        self.ref_coverage()
        self.write_output()


    #### CRBBLAST ####

    def seq_count(self):
        with open(self.fasta_f) as f:
            seq = ''
            for line in f:
                if line.startswith('>'):
                    try:
                        self.query_results[name]['seq'] = seq
                        self.query_results[name]['prot'] = self.prot_guess(seq)
                    except:
                        pass
                    seq = ''
                    line = line.strip('>').strip('\n').split(' ')[0]
                    name = line
                    if line not in self.query_results:
                        self.query_results[line] = {'prot': False, 'hits': [], 'seq': ''}
                    self.seqs += 1
                else:
                    seq += line.strip('\n')
                try:
                    self.query_results[name]['seq'] = seq
                    self.query_results[name]['prot'] = self.prot_guess(seq)
                except:
                    pass
        with open(self.prot_f) as f:
            seq = ''
            for line in f:
                if line.startswith('>'):
                    try:
                        self.target_results[name]['seq'] = seq
                        self.target_results[name]['prot'] = self.prot_guess(seq)
                    except Exception as e:
                        pass
                    seq = ''
                    line = line.strip('>').strip('\n').split(' ')[0]
                    name = line
                    if line not in self.target_results:
                        self.target_results[line] = {'prot': False, 'hits': [], 'seq': ''}
                    self.p_seqs += 1
                else:
                    seq += line.strip('\n')
                try:
                    self.target_results[name]['seq'] = seq
                    self.target_results[name]['prot'] = self.prot_guess(seq)
                except Exception as e:
                    pass
        return self.seqs

    def prot_guess(self, seq):
        seq = seq.upper()
        length  = len(seq)
        a_count = seq.count('A')
        t_count = seq.count('T')
        g_count = seq.count('G')
        c_count = seq.count('C')
        total_count = a_count + t_count + g_count + c_count
        if total_count / length >= 0.9:
            return False
        else:
            return True

    def load_outputs(self):
        with open(self.q2p_blast, encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                h = Hit(line, self.query, self.target)
                hit = h.lst()
                if line[0] not in self.query_dict:
                    self.query_dict[line[0]] = []
                self.query_dict[line[0]].append(h.lst())
                if hit['query'] not in self.query_results:
                    self.query_results[hit['query']] = []
                self.query_results[hit['query']]['hits'].append(hit)
                self.q_count += 1

        with open(self.p2q_blast, encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                h = Hit(line, self.target, self.query)
                hit = h.lst()
                if line[0] not in self.target_dict:
                    self.target_dict[line[0]] = []
                self.target_dict[line[0]].append(h.lst())
                if hit['query'] not in self.target_results:
                    self.target_results[hit['query']]['hits'] = []
                if hit not in self.target_results[hit['query']]['hits']:
                    self.target_results[hit['query']]['hits'].append(hit)
                self.t_count += 1
        
        return [self.q_count, self.t_count]

    def find_reciprocals(self):
        self.reciprocals = {}
        self.missed = {}
        self.evalues = []
        self.longest = 0
        hits = 0
        for query_id, list_of_hits in self.query_results.items():
            for query_index, target_hit in enumerate(list_of_hits['hits']):
                if target_hit['target'] in self.target_results:
                    list_of_hits_2 = self.target_results[target_hit['target']]['hits']
                    for target_index, query_hit2 in enumerate(list_of_hits_2):
                        if query_index == 0 and target_index == 0 and query_id == query_hit2['target']:
                            e = float(target_hit['evalue'])
                            if e == 0:
                                e = 1e-200
                            e = -math.log10(e)
                            if query_id not in self.reciprocals:
                                self.reciprocals[query_id] = []
                            self.reciprocals[query_id].append(target_hit)
                            hits += 1
                            if target_hit['alnlen'] > self.longest:
                                self.longest = target_hit['alnlen']
                            self.evalues.append({'e': e, 'length': target_hit['alnlen']})
                        elif query_id == query_hit2['target']:
                            if query_id not in self.missed:
                                self.missed[query_id] = []
                            self.missed[query_id].append(target_hit)
        return hits

    def find_secondaries(self):
        with open(f'{self.output}/evalues_data', 'w') as io:
            for h in self.evalues:
                self.length_hash.setdefault(h['length'], []).append(h)
                io.write(f"{h['length']}\t{h['e']}\n")

            for center in range(10, self.longest + 1):
                e = 0
                count = 0
                s = int(center * 0.1)
                s = s if s >= 5 else 5
                for side in range(-s, s + 1):
                    if center + side in self.length_hash:
                        for point in self.length_hash[center + side]:
                            e += point['e']
                            count += 1
                if count > 0:
                    mean = e / count
                    if center - 1 in self.fitting:
                        if self.fitting[center - 1] > mean:
                            self.fitting[center] = self.fitting[center - 1]
                        else:
                            self.fitting[center] = mean
                    else:
                        self.fitting[center] = mean

            with open(f'{self.output}/fitting_data', "w") as io:
                for center, mean in self.fitting.items():
                    io.write(f"{center}\t{mean}\n")

            hits = 0
            for id, lst in self.missed.items():
                for hit in lst:
                    l = int(hit['alnlen'])
                    e = hit['evalue']
                    e = -math.log10(1e-200) if e == 0 else -math.log10(e)
                    if l in self.fitting:
                        if e >= self.fitting[l]:
                            if id not in self.reciprocals:
                                self.reciprocals[id] = [hit]
                                hits += 1
            return hits
        
    def size(self):
        hits = 0
        for k,v in self.reciprocals.items():
            for hit in v:
                hits += 1
        return hits
    
    def write_output(self):
        for k,v in self.comp_stats.items():
            if type(v) == float:
                self.comp_stats[k] = round(v, 4)
            else:
                self.comp_stats[k] = int(v)
        s = ""
        rhfile = f'{self.output}/reciprocal_hits.csv' if not self.multi else f'{self.output}/{self.rh_name}_reciprocal_hits.csv'
        with open(rhfile, "w") as f:
            if self.reciprocals is not None:
                s += "query\ttarget\tid\talnlen\tevalue\tbitscore\t"
                s += "qstart..qend\ttstart..tend\tqlen\ttlen\n"
                for id, hits in self.reciprocals.items():
                    for hit in hits:
                        s += f"{hit['query']}\t{hit['target']}\t{hit['id']}\t{hit['alnlen']}\t{hit['evalue']}\t{hit['bitscore']}\t{hit['qstart']}..{hit['qend']}\t{hit['tstart']}..{hit['tend']}\t{hit['qlen']}\t{hit['tlen']}\n"
                f.write(s)
        # if self.assemblycsv:
        df = pd.DataFrame([self.comp_stats])
        existing_df = pd.read_csv(self.assemblycsv)
        existing_bottom = existing_df.iloc[-1]
        if existing_df.shape[0] > 1:
            existing_bottom[df.columns] = df.values[0]
            existing_df.iloc[-1] = existing_bottom
        else:
            existing_df = pd.concat([existing_df, df], axis=1)
        existing_df.to_csv(self.assemblycsv, index=False)
        os.remove(self.fa_translated)
        
    #### COMPARITIVE METRICS ####
    
    def reference_hits(self):
        for name, contig in self.target_results.items():
            self.reference.setdefault(name, {'hits': [], 'crbb': False, 'prot': self.target_results[name]['prot'], 'seq': self.target_results[name]['seq']})
        for query_id, lst in self.reciprocals.items():
            for hit in lst:
                if hit['target'] not in self.reference:
                    pass
                self.reference[hit['target']]['hits'].append(hit)
        self.comp_stats['CRBB_hits']           = self.seqs
        self.comp_stats['n_contigs_with_CRBB'] = len(self.reciprocals)
        self.comp_stats['p_contigs_with_CRBB'] = round(len(self.reciprocals) / self.seqs, 6)


    def rbh(self):
        n_refs_with_recip = 0
        total_crbb_hits = 0
        for x in self.reference:
            try:
                if len(self.reference[x]['hits']) > 0:
                    n_refs_with_recip += 1
                    self.reference[x]['crbb'] = True
            except:
                pass
            try:
                total_crbb_hits += len(self.reference[x]['hits'])
            except:
                pass
        self.comp_stats['rbh_per_reference'] = round(total_crbb_hits / self.p_seqs, 6)
        self.comp_stats['n_refs_with_CRBB']  = n_refs_with_recip
        self.comp_stats['p_refs_with_CRBB']  = round(n_refs_with_recip / self.p_seqs, 6)

    def ref_coverage(self):
        coverage_thresholds = [0.25, 0.5, 0.75, 0.85, 0.95]
        coverage_totals     = [0, 0, 0, 0, 0]
        total_coverage      = 0
        total_length        = 0

        for k,v in self.reference.items():
            if v['prot'] and v['crbb']:
                for hit in v['hits']:
                    coverage  = (3 * hit['alnlen']) - (3 * hit['mismatches']) - (3 * hit['gaps'])
                    coverage /= hit['tlen']
                    hit['coverage'] = coverage
            elif not v['prot'] and v['crbb']:
                for hit in v['hits']:
                    coverage  = hit['alnlen'] - hit['mismatches'] - hit['gaps']
                    coverage /= hit['tlen']
                    hit['coverage'] = coverage

        for k,v in self.reference.items():
            if v['prot'] and v['crbb']:
                covered = [False] * (len(v['seq']) * 3)
            else:
                covered = [False] * len(v['seq'])
            for i, hit in enumerate(v['hits']):
                if v['prot']:
                    if hit['qstart'] % 3 == 0:
                        tstart = 3 * hit['tstart'] - 4
                        tend   = 3 * hit['tend']
                    elif hit['qstart'] % 3 == 1:
                        tstart = 3 * hit['tstart'] - 2
                        tend   = 3 * hit['tend']
                    elif hit['qstart'] % 3 == 2:
                        tstart = 3 * hit['tstart'] - 3
                        tend   = 3 * hit['tend'] - 1
                    if hit['qlen'] % 3 == 1:
                        tend  += 1
                    elif hit['qlen'] % 3 == 2:
                        tend  += 2
                else:
                    tstart = hit['tstart']
                    tend   = hit['tend']
                for b in range(tstart, tend + 1):
                    covered[b - 1] = True
            coverage = sum(covered)
            ref_p = coverage / len(covered)
            v['reference_coverage'] = ref_p
            for index, n in enumerate(coverage_thresholds):
                if ref_p >= n:
                    coverage_totals[index] += 1

            total_coverage += coverage
            total_length   += len(v['seq']) * 3 if v['prot'] else len(v['seq'])

        # calculate proportion of ref sequences with coverage over thresholds
        for i, p in enumerate(coverage_thresholds):
            self.comp_stats[f"cov{(100 * p):.0f}"]   = coverage_totals[i]
            self.comp_stats[f"p_cov{(100 * p):.0f}"] = round(coverage_totals[i] / len(self.reference), 6)

        self.comp_stats['reference_coverage'] = round(total_coverage / total_length, 6)

    def translator(self):
        self.fa_translated = self.output + '/' + os.path.basename(self.fasta_f).split('.')[0] + '_translated.fa'
        if os.path.exists(self.fa_translated):
            os.remove(self.fa_translated)
        else:
            os.mknod(self.fa_translated)

        codons = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }

        with open(self.fasta_f) as f:
            seq = []
            prot = ''
            final = ''
            for line in f:
                if line.startswith('>'):
                    seq = []
                    prot = ''
                    final = ''
                    line = line.strip('\n').split(' ')[0]
                    name = line + '\n'
                    with open(self.fa_translated, 'a') as f_out:
                        f_out.write(name)
                    line = next(f)
                    while not line.startswith('>'):
                        seq.append(line.strip('\n'))
                        try:
                            line = next(f)
                        except:
                            break
                    for x in seq:
                        for i in range(0, len(x), 3):
                            codon = x[i:i+3]
                            prot += codons[codon]
                        final += prot
                        prot = ''
                    with open(self.fa_translated, 'a') as f_out:
                        for i in range(0, len(final), 80):
                            f_out.write(final[i:i+80] + '\n')

    def diamond(self):

        self.p2q_db = self.output + '/' + os.path.basename(self.prot_f).split('.')[0] + "_" + os.path.basename(self.fasta_f).split('.')[0] + '_db.dmnd'
        self.q2p_db = self.output + '/' + os.path.basename(self.fasta_f).split('.')[0] + "_" + os.path.basename(self.prot_f).split('.')[0] + '_db.dmnd'

        self.p2q_blast = self.p2q_db.strip('.dmnd') + '.blast'
        self.q2p_blast = self.q2p_db.strip('.dmnd') + '.blast'

        while self.translator():
            pass

        prot_db  = ['diamond', 'makedb', '--threads', str(self.threads), '--in', self.prot_f, '--db', self.p2q_db]
        query_db = ['diamond', 'makedb', '--threads', str(self.threads), '--in', self.fa_translated, '--db', self.q2p_db]
        blastx   = ['diamond', 'blastx', '--threads', str(self.threads), '-b1', '--db', self.p2q_db, '--out', self.q2p_blast, '--query', self.fasta_f, '--evalue', '0.00001', '--max-target-seqs', '50', '--masking', '0', '--ultra-sensitive', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
        blastp   = ['diamond', 'blastp', '--threads', str(self.threads), '-b1', '--db', self.q2p_db, '--out', self.p2q_blast, '--query', self.prot_f, '--evalue', '0.00001', '--max-target-seqs', '50', '--masking', '0', '--ultra-sensitive', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

        with open(f'{self.output}/make_prot_db.log', 'w') as f:
            make_prot_db = subprocess.Popen(prot_db, stdout=f, stderr=f, preexec_fn=os.setsid)
            make_prot_db_pid = make_prot_db.pid
            make_prot_db.communicate()

        with open(f'{self.output}/make_query_db.log', 'w') as f:
            make_query_db = subprocess.Popen(query_db, stdout=f, stderr=f, preexec_fn=os.setsid)
            make_query_db_pid = make_query_db.pid
            make_query_db.communicate()

        with open(f'{self.output}/blastx.log', 'w') as f:
            blastx = subprocess.Popen(blastx, stdout=f, stderr=f, preexec_fn=os.setsid)
            blastx_pid = blastx.pid
            blastx.communicate()

        with open(f'{self.output}/blastp.log', 'w') as f:
            blastp = subprocess.Popen(blastp, stdout=f, stderr=f, preexec_fn=os.setsid)
            blastp_pid = blastp.pid
            blastp.communicate()







