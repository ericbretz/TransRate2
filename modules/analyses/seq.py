import pysam
import sys
from multiprocessing import Pool, set_start_method

def process_task(args):
    split = args[0]
    bamfile = args[1]
    threads = args[2]
    bam = pysam.AlignmentFile(bamfile, 'rb')
    refs = bam.references[split::threads]
    tmp_dct = {ref: {'seq': {
                    'true'                : 0,
                    'count'               : 0,
                    'softclipped'         : 0,
                    'softclipped_length'  : 0,
                    'total_length'        : 0,
                    'total_reads'         : 0,
                    'p_softclipped'       : 0.0,
                    'p_softclipped_length': 0.0
                }, 'stats': {
                    'p_seqtrue': 0}} for ref in refs}

    cigdict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}

    for ref in refs:
        total_softclipped = 0
        total_reads = 0
        total_length = 0
        total_softclipped_length = 0
        softclipped = 0
        max_read_length = 0
        count = 0
        for fetch in bam.fetch(reference=ref):
            count += 1
            if count > 5000:
              break
            max_read_length = max(max_read_length, fetch.query_length)
            try:
                total_reads += 1
                name_r  = fetch.reference_name
                len_r   = fetch.reference_length
                nm      = fetch.get_tag('NM', with_value_type = True)
                scale   = (len_r - 35) / len_r
                tmp_dct[name_r]['seq']['count']      += 1
                tmp_dct[name_r]['seq']['true']       += (((len_r - nm[0]) / len_r) - scale) * (1 / (1 - scale))
                seqtrue = round(tmp_dct[name_r]['seq']['true'] / tmp_dct[name_r]['seq']['count'], 4)
                tmp_dct[name_r]['stats']['p_seqtrue'] = seqtrue
            except:
                pass
            try:
                cigar         = fetch.cigartuples
                f_start       = fetch.query_alignment_start
                f_end         = fetch.query_alignment_end
                f_length      = fetch.query_length
                total_length += f_length

                if cigar:
                    for cig in cigar:
                        letter = cigdict[cig[0]]
                        if letter == 'D':
                            continue
                        elif letter == 'S':
                            softclipped = 1
                            total_softclipped += 1
                            total_softclipped_length += cig[1]
            except:
                pass

        try:
            tmp_dct[ref]['seq']['softclipped']        = total_softclipped
            tmp_dct[ref]['seq']['softclipped_length'] = total_softclipped_length
            tmp_dct[ref]['seq']['total_length']       = total_length
            tmp_dct[ref]['seq']['total_reads']        = total_reads
            tmp_dct[ref]['seq']['max_read_length']    = max_read_length
            if total_reads:
                tmp_dct[ref]['seq']['p_softclipped'] = float(round(total_softclipped / total_reads, 4))
            if total_length:
                tmp_dct[ref]['seq']['p_softclipped_length'] = float(round(total_softclipped_length / total_length, 4))
            
        except:
            continue

    return tmp_dct

def seq(bamfile, bamdct, threads, single):

    try:
        set_start_method('spawn')
    except:
        pass
    with Pool(threads) as p:
        results = p.map(process_task, [[i, bamfile, threads] for i in range(threads)])

    for result in results:
        for key in result:
            if key in bamdct:
                for subkey in result[key]:
                    if isinstance(bamdct[key][subkey], dict):
                        bamdct[key][subkey].update(result[key][subkey])
                    else:
                        bamdct[key][subkey] += result[key][subkey]
            else:
                bamdct[key] = result[key]

    for k in bamdct:
        if bamdct[k]['stats']['p_seqtrue'] == 0:
            bamdct[k]['stats']['p_seqtrue'] = 1
    return bamdct