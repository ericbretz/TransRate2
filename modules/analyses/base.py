import numpy as np
import pysam
from multiprocessing import Pool, set_start_method

def process_task(args):
    split = args[0]
    bamfile = args[1]
    threads = args[2]
    bamdct = args[3]
    bam = pysam.AlignmentFile(bamfile, 'rb',)
    refs = bam.references[split::threads]
    tmp_dct = {ref: {'stats': {'basesuncovered': 0}, 'coverage': {'bases': []}} for ref in refs}
    for ref in refs:
        length = bamdct[ref]['stats']['length']
        cov = bam.count_coverage(ref, stop=length, quality_threshold=0)
        cov_sum = np.sum(cov, 0)
        uncover = length - np.count_nonzero(cov_sum)
        tmp_dct[ref]['coverage']['bases'] = cov_sum
        tmp_dct[ref]['stats']['basesuncovered'] = uncover

    return tmp_dct

def base(bamfile, bamdct, threads, single):
    try:
        set_start_method('spawn')
    except:
        pass
    with Pool(threads) as p:
        results = p.map(process_task, [[i, bamfile, threads, bamdct] for i in range(threads)])

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
