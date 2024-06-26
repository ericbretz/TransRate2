import pysam
from multiprocessing import Pool, set_start_method
import sys

def process_task(args):
    split = args[0]
    bamfile = args[1]
    threads = args[2]
    single = args[3]
    bam = pysam.AlignmentFile(bamfile, 'rb')
    refs = bam.references[split::threads]
    tmp_dct = {ref: {'stats': {'fragments': 0, 'properpair': 0, 'both_mapped': 0, 'bridges': 0}} for ref in refs}

    for ref in refs:
        for fetch in bam.fetch(reference=ref):
            try:
                name_r = fetch.reference_name
                if single == 2:
                    if (fetch.is_read1 and fetch.is_mapped) or (fetch.is_read2 and not fetch.mate_is_mapped):
                        tmp_dct[name_r]['stats']['fragments'] += 1
                    if fetch.is_read1 and fetch.is_paired:
                        if fetch.is_proper_pair:
                            tmp_dct[name_r]['stats']['properpair'] += 1
                    if fetch.is_read1 and fetch.is_mapped and fetch.mate_is_mapped:
                        tmp_dct[name_r]['stats']['both_mapped'] += 1
                    if fetch.reference_id != fetch.next_reference_id and fetch.is_read1:
                        tmp_dct[name_r]['stats']['bridges'] += 1
                else:
                    if fetch.is_mapped:
                        tmp_dct[name_r]['stats']['fragments'] += 1
            except:
                pass

    return tmp_dct

def frag(bamfile, bamdct, threads, single):
    try:
        set_start_method('spawn')
    except:
        pass
    with Pool(threads) as p:
        results = p.map(process_task, [[i, bamfile, threads, single] for i in range(threads)])

    p.terminate()

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
    return bamdct