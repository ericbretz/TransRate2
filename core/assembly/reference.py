import warnings
warnings.filterwarnings('ignore', category=UserWarning)
import os
import subprocess
import sys
import math
import pandas as pd
import json

class Reference:

    def __init__(self):
        return
    
    def mainRun(self, args):
        assembly_data         = args
        multi        :bool    = assembly_data.get('multi', False)
        rh_name      :str     = assembly_data['assembly_base_name']
        
        if not rh_name or rh_name.strip() == '':
            rh_name = os.path.splitext(os.path.basename(assembly_data['assembly']))[0]
        
        rh_name = "".join(c for c in rh_name if c.isalnum() or c in "._-")
        
        assemblycsv         :str = assembly_data.get('assemblycsv', os.path.join(assembly_data['dict_dir']['results'], 'assembly.csv'))
        assemblyPath        :str = assembly_data['assembly']
        referencePath       :str = assembly_data['reference']
        referenceOutput     :str = assembly_data['dict_dir']['temp_reference']
        threads             :int = assembly_data['threads']
        translatedPath      :str = os.path.join(referenceOutput, os.path.basename(assemblyPath).split('.')[0] + '_translated.fa')
        
        if not self.prot_guess(referencePath):
            translatedReferencePath = os.path.join(referenceOutput, os.path.basename(referencePath).split('.')[0] + '_translated.fa')
            self.translator(referencePath, translatedReferencePath)
            referencePath = translatedReferencePath

        assemblyTranslated  :str = self.translator(assemblyPath, translatedPath)
        p2q_name            :str = os.path.basename(referencePath).split('.')[0] + "_" + os.path.basename(assemblyPath).split('.')[0]
        q2p_name            :str = os.path.basename(assemblyPath).split('.')[0] + "_" + os.path.basename(referencePath).split('.')[0]
        p2q_db_path         :str = os.path.join(referenceOutput, p2q_name + '_db.dmnd')
        q2p_db_path         :str = os.path.join(referenceOutput, q2p_name + '_db.dmnd')
        p2q_db_blast        :str = os.path.join(referenceOutput, p2q_name + '_db.blast')
        q2p_db_blast        :str = os.path.join(referenceOutput, q2p_name + '_db.blast')

        self.diamond(p2q_db_path, q2p_db_path, p2q_db_blast, q2p_db_blast, assemblyPath, translatedPath, referencePath, threads, referenceOutput)
        query_results, target_results, n_seqs , p_seqs           = self.seq_count(assemblyPath, referencePath)
        query_dict   , target_dict   , q_count, t_count          = self.diamond_results(q2p_db_blast, p2q_db_blast, query_results, target_results)
        reciprocals  , missed        , evalues, longest, hits    = self.find_reciprocal_hits(query_results, target_results)
        reciprocals  , hits          , missed , evalues, longest = self.find_secondaries(evalues, referenceOutput, longest, missed, reciprocals)
        reference    , comp_stats                                = self.reference_hits(reciprocals, target_results, n_seqs)
        reference    , comp_stats_updated                        = self.rbh(reference, p_seqs)
        comp_stats.update(comp_stats_updated)
        reference    , comp_stats_updated                        = self.ref_coverage(reference)
        comp_stats.update(comp_stats_updated)
        self.write_output(comp_stats, assemblycsv, referenceOutput, multi, rh_name, reciprocals)
        comp_stats = {'assembly': comp_stats}
        return comp_stats

    def diamond(self, p2q_db_path, q2p_db_path, p2q_db_blast, q2p_db_blast, assemblyPath, translatedPath, referencePath, threads, referenceOutput):

        prot_db  = ['diamond', 'makedb', '--threads', str(threads), '--in', referencePath, '--db', p2q_db_path]
        query_db = ['diamond', 'makedb', '--threads', str(threads), '--in', translatedPath, '--db', q2p_db_path]
        blastx   = ['diamond', 'blastx', '--threads', str(threads), '-b1', '--db', p2q_db_path, '--out', q2p_db_blast, '--query', assemblyPath, '--evalue', '0.00001', '--max-target-seqs', '50', '--masking', '0', '--ultra-sensitive', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
        blastp   = ['diamond', 'blastp', '--threads', str(threads), '-b1', '--db', q2p_db_path, '--out', p2q_db_blast, '--query', referencePath, '--evalue', '0.00001', '--max-target-seqs', '50', '--masking', '0', '--ultra-sensitive', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

        with open(f'{referenceOutput}/make_prot_db.log', 'w') as f:
            make_prot_db     = subprocess.Popen(prot_db, stdout=f, stderr=f, preexec_fn=os.setsid)
            make_prot_db_pid = make_prot_db.pid
            make_prot_db.communicate()

        with open(f'{referenceOutput}/make_query_db.log', 'w') as f:
            make_query_db     = subprocess.Popen(query_db, stdout=f, stderr=f, preexec_fn=os.setsid)
            make_query_db_pid = make_query_db.pid
            make_query_db.communicate()
         

        with open(f'{referenceOutput}/blastx.log', 'w') as f:
            blastx     = subprocess.Popen(blastx, stdout=f, stderr=f, preexec_fn=os.setsid)
            blastx_pid = blastx.pid
            blastx.communicate()

        with open(f'{referenceOutput}/blastp.log', 'w') as f:
            blastp     = subprocess.Popen(blastp, stdout=f, stderr=f, preexec_fn=os.setsid)
            blastp_pid = blastp.pid
            blastp.communicate()

    def translator(self, assemblyPath, translatedPath):
        if os.path.exists(translatedPath):
            os.remove(translatedPath)
        else:
            os.mknod(translatedPath)

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

        def translate_sequence(sequence):
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3]
                yield codons.get(codon, 'X')

        with open(translatedPath, 'w') as f_out:
            with open(assemblyPath) as f:
                lines = ''
                for line in f:
                    if line.startswith('>'):
                        if lines:
                            f_out.write(''.join(translate_sequence(lines)) + '\n')
                        f_out.write(line)
                        lines = ''
                    else:
                        lines += line.strip()

                if lines:
                    f_out.write(''.join(translate_sequence(lines)) + '\n')

    def seq_count(self, assemblyPath, referencePath):
        query_results  = {}
        target_results = {}
        n_seqs         = [0]
        p_seqs         = [0]
        for type, file, seqs in [(query_results, assemblyPath, n_seqs), (target_results, referencePath, p_seqs)]:
            with open(file) as f:
                seq   = ''
                lines = f.readlines()
                for line in lines:
                    if line.startswith('>'):
                        seqs[0] += 1
                        if seq:
                            type[name]['seq']        = seq
                            type[name]['prot_guess'] = self.prot_guess(seq)
                        seq  = ''
                        name = line[1:].split()[0]
                        type[name] = {}
                    else:
                        seq += line.strip('\n')

                if seq:
                    type[name]['seq']        = seq
                    type[name]['prot_guess'] = self.prot_guess(seq)
        return query_results, target_results, n_seqs[0], p_seqs[0]

    def prot_guess(self, seq):
        if os.path.isfile(seq):
            with open(seq) as f:
                lines = f.readlines()
                seq = ''
                for line in lines:
                    if not line.startswith('>'):
                        seq += line.strip()
                    elif seq:
                        break
        seq = seq.upper()
        total_count = sum(seq.count(nuc) for nuc in 'ATGC')
        return total_count / len(seq) < 0.9
    
    def diamond_results(self, q2p_db_blast, p2q_db_blast, query_results, target_results):
        q_count = 0
        t_count = 0

        for db, dict, count in [(q2p_db_blast, query_results, q_count), (p2q_db_blast, target_results, t_count)]:
            with open(db, encoding='utf-8') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('#'):
                        continue
                    else:
                        line = line.strip().split('\t')
                        hit = self.hit(line)
                        if hit['query'] not in dict:
                            dict[hit['query']] = {'hits': []}
                        elif 'hits' not in dict[hit['query']]:
                            dict[hit['query']]['hits'] = []
                        dict[hit['query']]['hits'].append(hit)
                        count += 1
        return query_results, target_results, q_count, t_count

    def hit(self, line):
        return {
            'query'     : line[0].split()[0],
            'target'    : line[1].split()[0],
            'id'        : line[2],
            'alnlen'    : int(line[3]),
            'mismatches': int(line[4]),
            'gaps'      : int(line[5]),
            'qstart'    : int(line[6]),
            'qend'      : int(line[7]),
            'tstart'    : int(line[8]),
            'tend'      : int(line[9]),
            'evalue'    : float(line[10]),
            'bitscore'  : float(line[11]),
            'qlen'      : float(line[12]),
            'tlen'      : float(line[13])
        }

    def find_reciprocal_hits(self, query_results, target_results):
        reciprocals = {}
        missed      = {}
        evalues     = []
        longest     = 0
        hits        = 0
        for query_id, query_hits in query_results.items():
            if 'hits' not in query_hits:
                continue
            for query_index, target_hit in enumerate(query_hits['hits']):
                target_id = target_hit['target']
                if target_id in target_results and 'hits' in target_results[target_id]:
                    target_hits = target_results[target_id]['hits']
                    for target_index, query_hit2 in enumerate(target_hits):
                        if query_index == 0 and target_index == 0 and query_id == query_hit2['target']:
                            e = float(target_hit['evalue'])
                            if e == 0:
                                e = 1e-200
                            e = -math.log10(e)
                            if query_id not in reciprocals:
                                reciprocals[query_id] = []
                            reciprocals[query_id].append(target_hit)
                            hits += 1
                            if target_hit['alnlen'] > longest:
                                longest = target_hit['alnlen']
                            evalues.append({'e': e, 'length': target_hit['alnlen']})
                        elif query_id == query_hit2['target']:
                            if query_id not in missed:
                                missed[query_id] = []
                            missed[query_id].append(target_hit)
        return reciprocals, missed, evalues, longest, hits

    def find_secondaries(self, evalues, referenceOutput, longest, missed, reciprocals):
        length_dict = {}
        fitting     = {}
        with open(f'{referenceOutput}/evalues_data', 'w') as f:
            for h in evalues:
                length_dict.setdefault(h['length'], []).append(h)
                f.write(f'{h["length"]}\t{h["e"]}\n')

        for center in range(10, longest + 1):
            e     = 0
            count = 0
            s     = int(center*0.1)
            s     = s if s >= 5 else 5
            for side in range(-s, s +1):
                if center + side in length_dict:
                    for point in length_dict[center + side]:
                        e     += point['e']
                        count += 1
            if count:
                mean = e / count
                if center - 1 in fitting:
                    if fitting[center - 1] > mean:
                        fitting[center] = fitting[center - 1]
                    else:
                        fitting[center] = mean
                else:
                    fitting[center] = mean
        
        with open(f'{referenceOutput}/fitting_data', 'w') as f:
            for center, e in sorted(fitting.items()):
                f.write(f'{center}\t{e}\n')
        
        hits = 0
        for id, lst in missed.items():
            for hit in lst:
                l = int(hit['alnlen'])
                e = hit['evalue']
                e = -math.log10(1e-200) if e == 0 else -math.log10(e)
                if l in fitting:
                    if e >= fitting[l]:
                        if id not in reciprocals:
                            reciprocals[id] = []
                        reciprocals[id].append(hit)
                        hits += 1
        return reciprocals, hits, missed, evalues, longest
    
    def reference_hits(self, reciprocals, target_results, n_seqs):
        reference  = {}
        comp_stats = {}
        for name, contig in target_results.items():
            reference.setdefault(name, {'hits': [], 'crbb': False, 'prot': contig['prot_guess'], 'seq': contig['seq']})
        for query_id, lst in reciprocals.items():
            for hit in lst:
                if hit['query'] in reference:
                    reference[hit['query']]['hits'].append(hit)
                    comp_stats['CRBBhits'] = n_seqs
                    comp_stats['nContigsWithCRBB'] = len(reciprocals)
                    comp_stats['pContigsWithCRBB'] = comp_stats['nContigsWithCRBB'] / len(reference) if len(reference) > 0 else 0
        return reference, comp_stats

    def rbh(self, reference, p_seqs):
        n_refs_with_recip = 0
        total_crbb_hits   = 0
        comp_stats        = {}
        for x in reference:
            if 'hits' in reference[x] and len(reference[x]['hits']) > 0:
                n_refs_with_recip += 1
                reference[x]['crbb'] = True
            if 'hits' in reference[x]:
                total_crbb_hits += len(reference[x]['hits'])
        comp_stats['rbhPerReference'] = round(total_crbb_hits / p_seqs, 6) if p_seqs > 0 else 0
        comp_stats['nRefsWithCRBB']   = n_refs_with_recip
        comp_stats['pRefsWithCRBB']   = round(n_refs_with_recip / p_seqs, 6) if p_seqs > 0 else 0
        return reference, comp_stats
    
    def ref_coverage(self, reference):
        coverage_thresholds = [0.25, 0.5, 0.75, 0.85, 0.95]
        coverage_totals     = [0, 0, 0, 0, 0]
        total_coverage      = 0
        total_length        = 0
        comp_stats          = {}

        for k,v in reference.items():
            if v['prot'] and v['crbb']:
                for hit in v['hits']:
                    coverage  = (3 * hit['alnlen']) - (3 * hit['mismatches']) - (3 * hit['gaps'])
                    coverage /= hit['tlen'] if hit['tlen'] > 0 else 1
                    hit['coverage'] = coverage
            elif not v['prot'] and v['crbb']:
                for hit in v['hits']:
                    coverage  = hit['alnlen'] - hit['mismatches'] - hit['gaps']
                    coverage /= hit['tlen'] if hit['tlen'] > 0 else 1
                    hit['coverage'] = coverage

        for k,v in reference.items():
            if v['crbb']:
                if v['prot']:
                    covered = [False] * ((len(v['seq']) - 1) * 3)
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
                        if 0 <= b - 1 < len(covered):
                            covered[b - 1] = True
                coverage = sum(covered)
                ref_p    = coverage / len(covered) if len(covered) > 0 else 0
                v['reference_coverage'] = ref_p
                for index, n in enumerate(coverage_thresholds):
                    if ref_p >= n:
                        coverage_totals[index] += 1

                total_coverage += coverage
                total_length   += len(v['seq']) * 3 if v['prot'] else len(v['seq'])

        for i, p in enumerate(coverage_thresholds):
            comp_stats[f"cov{(100 * p):.0f}"]  = coverage_totals[i]
            comp_stats[f"pCov{(100 * p):.0f}"] = round(coverage_totals[i] / len(reference), 6)

        comp_stats['referenceCoverage'] = round(total_coverage / total_length, 6) if total_length > 0 else 0
        return reference, comp_stats
    
    def write_output(self, comp_stats, assemblycsv, output, multi, rh_name, reciprocals):
        for k,v in comp_stats.items():
            if type(v) == float:
                comp_stats[k] = round(v, 4)
            else:
                comp_stats[k] = int(v)
        
        s = ""
        if multi:
            safe_rh_name = rh_name if rh_name else "assembly"
            rhfile       = f'{output}/{safe_rh_name}_reciprocal_hits.csv'
        else:
            rhfile       = f'{output}/reciprocal_hits.csv'
            
        os.makedirs(output, exist_ok=True)
            
        with open(rhfile, "w") as f:
            if reciprocals is not None and len(reciprocals) > 0:
                s += "query\ttarget\tid\talnlen\tevalue\tbitscore\t"
                s += "qstart..qend\ttstart..tend\tqlen\ttlen\n"
                for id, hits in reciprocals.items():
                    for hit in hits:
                        s += f"{hit['query']}\t{hit['target']}\t{hit['id']}\t{hit['alnlen']}\t{hit['evalue']}\t{hit['bitscore']}\t{hit['qstart']}..{hit['qend']}\t{hit['tstart']}..{hit['tend']}\t{hit['qlen']}\t{hit['tlen']}\n"
            f.write(s)
        
        if assemblycsv and assemblycsv.strip():
            csv_dir = os.path.dirname(assemblycsv)
            if csv_dir:
                os.makedirs(csv_dir, exist_ok=True)
            
            try:
                df = pd.DataFrame([comp_stats])
                
                if os.path.exists(assemblycsv):
                    existing_df = pd.read_csv(assemblycsv)
                    
                    if len(existing_df) > 0:
                        existing_bottom = existing_df.iloc[-1].copy()
                        
                        for col in df.columns:
                            if col not in existing_df.columns:
                                existing_df[col] = None
                        
                        for col in df.columns:
                            try:
                                existing_bottom[col] = df.iloc[0][col]
                            except (ValueError, TypeError) as ve:
                                existing_bottom[col] = str(df.iloc[0][col])
                        
                        existing_df.iloc[-1] = existing_bottom
                    else:
                        existing_df = df
                        
                    existing_df.to_csv(assemblycsv, index=False)
                else:
                    df.to_csv(assemblycsv, index=False)
                    
            except FileNotFoundError as fnf:
                print(f"Warning: Could not access CSV path {assemblycsv}: {fnf}")
            except PermissionError as pe:
                print(f"Warning: Permission denied accessing {assemblycsv}: {pe}")
            except Exception as e:
                print(f"Warning: Could not update assembly CSV {assemblycsv}: {e}")
                pass