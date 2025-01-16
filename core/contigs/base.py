import numpy as np
from pysam import AlignmentFile

def mainRun(args):
    main    = args[0]
    i       = args[1]
    bam     = AlignmentFile(main['sortedBam'], 'rb')  
    refs    = bam.references[i::main['threads']]
    baseDct = {ref: {'basesUncovered': 0,
                     'other': {'bases': []}} for ref in refs}

    for ref in refs:
        length                                  = bam.get_reference_length(ref)
        cov                                     = bam.count_coverage(ref, stop=length, quality_threshold=0)
        cov_sum                                 = np.sum(cov, 0)
        baseDct[ref]['basesUncovered']          = length - np.count_nonzero(cov_sum)
        baseDct[ref]['other']['bases']          = cov_sum.tolist()
    return baseDct  
