import numpy as np
from pysam import AlignmentFile

def calculate_fragment_distance(pos1, pos2, end1, end2):
    leftmost_start = min(pos1, pos2)
    rightmost_end = max(end1, end2)
    
    return rightmost_end - leftmost_start

def collect_fragment_distances(ref, bam):
    mateDct   = {}
    distances = []

    for fetch in bam.fetch(reference=ref):
        if fetch.is_unmapped or fetch.reference_start is None or fetch.reference_end is None:
            continue
        
        if not fetch.mate_is_unmapped and fetch.next_reference_id != fetch.reference_id:
            continue
            
        query_name = fetch.query_name
        if query_name not in mateDct:
            mateDct[query_name] = {}

        mateDct[query_name]['f' if fetch.is_forward else 'r'] = (
            fetch.reference_start, 
            fetch.reference_end
        )

    mateDct = {k: v for k, v in mateDct.items() if 'f' in v and 'r' in v}

    for query_name, alignments in mateDct.items():
        pos1, end1 = alignments['f']
        pos2, end2 = alignments['r']

        if pos1 >= 0 and pos2 >= 0 and end1 > pos1 and end2 > pos2:
            fragment = calculate_fragment_distance(pos1, pos2, end1, end2)
            if fragment > 0:
                distances.append(fragment)

    return distances, mateDct

def calculate_global_threshold(bam, refs):
    all_distances = []
    
    for ref in refs:
        distances, _ = collect_fragment_distances(ref, bam)
        all_distances.extend(distances)
    
    if len(all_distances) == 0:
        return 0
    
    all_distances = np.array(all_distances)
    median        = np.median(all_distances)
    mad           = np.median(np.abs(all_distances - median))
    
    threshold = int(median + 2 * mad)
    
    return threshold

def mainRun(args):
    contig_data = args[0]
    i           = args[1]
    bam         = AlignmentFile(contig_data['dict_file']['samtools_bam'], 'rb')
    refs        = bam.references[i::contig_data['threads']]
    goodDct     = {ref: {'good': 0, 'pGood': 0} for ref in refs}

    global_threshold = calculate_global_threshold(bam, refs)
    
    if not global_threshold:
        return goodDct

    for ref in refs:
        distances, mateDct = collect_fragment_distances(ref, bam)
        
        if len(distances) == 0:
            continue

        if contig_data['mode'] == 2:
            for fetch in bam.fetch(reference=ref):
                if fetch.is_read1 and fetch.is_mapped and fetch.mate_is_mapped:
                    query_name = fetch.query_name
                    if (query_name in mateDct and 
                        'f' in mateDct[query_name] and 
                        'r' in mateDct[query_name]):
                        
                        pos1, end1 = mateDct[query_name]['f']
                        pos2, end2 = mateDct[query_name]['r']
                        
                        fragment_distance = calculate_fragment_distance(pos1, pos2, end1, end2)
                        
                        if fragment_distance <= global_threshold:
                            goodDct[ref]['good'] += 2
        else:
            for query_name, alignments in mateDct.items():
                pos1, end1 = alignments['f']
                pos2, end2 = alignments['r']
                
                fragment_distance = calculate_fragment_distance(pos1, pos2, end1, end2)
                
                if fragment_distance <= global_threshold:
                    goodDct[ref]['good'] += 1

        fragments             = contig_data['contigDF'][contig_data['contigDF']['name'] == ref]['fragments'].iloc[0]
        goodDct[ref]['pGood'] = goodDct[ref]['good'] / fragments if fragments else 0

    return goodDct
