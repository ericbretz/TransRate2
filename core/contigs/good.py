import numpy as np
from pysam import AlignmentFile

def realDist(ref, bam):
    mateDct = {}
    count = 0
    m = 0
    s = 0

    for fetch in bam.fetch(reference=ref):
        query_name = fetch.query_name
        if query_name not in mateDct:
            mateDct[query_name] = {}
        mateDct[query_name]['f' if fetch.is_forward else 'r'] = (fetch.query_alignment_start, fetch.query_length)

    mateDct = {k: v for k, v in mateDct.items() if 'f' in v and 'r' in v}

    for alignments in mateDct.values():
        pos1, len1 = alignments['f']
        pos2, len2 = alignments['r']

        if pos1 >= 0 and pos2 >= 0:
            fragment = abs(pos1 - pos2) + (len1 if pos1 > pos2 else len2)
            count += 1
            delta = fragment - m
            m += delta / count
            s += delta * (fragment - m)

    if count:
        mean = m
        s = np.sqrt(s / count)
        rD = int((3 * s) + mean)
        return rD, mateDct
    else:
        return 0, mateDct

def mainRun(args):
    main = args[0]
    i = args[1]
    bam = AlignmentFile(main['sortedBam'], 'rb')
    refs = bam.references[i::main['threads']]
    goodDct = {ref: {'good': 0, 'pGood': 0} for ref in refs}

    for ref in refs:
        realisticDistance, mateDct = realDist(ref, bam)
        if not realisticDistance:
            continue

        for alignments in mateDct.values():
            pos1, len1 = alignments['f']
            pos2, len2 = alignments['r']
            if abs(pos1 - pos2) + (len1 if pos1 > pos2 else len2) <= realisticDistance:
                goodDct[ref]['good'] += 1

        fragments = main['contigDF'][main['contigDF']['name'] == ref]['fragments'].iloc[0]
        goodDct[ref]['pGood'] = goodDct[ref]['good'] / fragments if fragments else 0

    return goodDct