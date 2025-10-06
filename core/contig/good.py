import numpy as np
from pysam import AlignmentFile

def is_splice_read(cigar):
    if cigar is None:
        return False
    
    for op, length in cigar:
        if op == 3:
            return True
    return False

def get_intron_length(cigar):
    if cigar is None:
        return 0
    
    intron_length = 0
    for op, length in cigar:
        if op == 3:
            intron_length += length
    return intron_length

def calculate_appropriate_distance(pos1, len1, pos2, len2, is_splice_pair, intron_len1, intron_len2):
    genomic_distance = abs(pos1 - pos2) + (len1 if pos1 > pos2 else len2)
    
    if is_splice_pair and (intron_len1 > 0 or intron_len2 > 0):
        transcript_distance = genomic_distance - intron_len1 - intron_len2
        min_distance        = len1 + len2
        if transcript_distance < min_distance:
            transcript_distance = min_distance
            
        return transcript_distance
    else:
        return genomic_distance

def collect_fragment_distances(ref, bam):
    """Collect fragment distances and mate pair information for a contig"""
    mateDct     = {}
    splice_info = {}
    distances   = []

    for fetch in bam.fetch(reference=ref):
        query_name = fetch.query_name
        if query_name not in mateDct:
            mateDct[query_name]     = {}
            splice_info[query_name] = {}
        
        mateDct[query_name]['f' if fetch.is_forward else 'r'] = (
            fetch.query_alignment_start, 
            fetch.query_alignment_length
        )
        
        splice_info[query_name]['f' if fetch.is_forward else 'r'] = {
            'is_splice' : is_splice_read(fetch.cigar),
            'intron_len': get_intron_length(fetch.cigar)
        }

    mateDct = {k: v for k, v in mateDct.items() if 'f' in v and 'r' in v}

    for query_name, alignments in mateDct.items():
        if query_name not in splice_info:
            continue
            
        pos1, len1 = alignments['f']
        pos2, len2 = alignments['r']

        if pos1 >= 0 and pos2 >= 0:
            f_splice = splice_info[query_name].get('f', {'is_splice': False, 'intron_len': 0})
            r_splice = splice_info[query_name].get('r', {'is_splice': False, 'intron_len': 0})
            
            is_splice_pair = f_splice['is_splice'] or r_splice['is_splice']
            
            fragment = calculate_appropriate_distance(
                pos1, len1, pos2, len2, 
                is_splice_pair,
                f_splice['intron_len'], 
                r_splice['intron_len']
            )
            
            distances.append(fragment)

    mateDct['_splice_info'] = splice_info
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

        splice_info = mateDct.get('_splice_info', {})

        if contig_data['mode'] == 2:
            for fetch in bam.fetch(reference=ref):
                if fetch.is_read1 and fetch.is_mapped and fetch.mate_is_mapped:
                    query_name = fetch.query_name
                    if (query_name in mateDct and 
                        'f' in mateDct[query_name] and 
                        'r' in mateDct[query_name]):
                        
                        pos1, len1 = mateDct[query_name]['f']
                        pos2, len2 = mateDct[query_name]['r']
                        
                        if query_name in splice_info:
                            f_splice = splice_info[query_name].get('f', {'is_splice': False, 'intron_len': 0})
                            r_splice = splice_info[query_name].get('r', {'is_splice': False, 'intron_len': 0})
                            
                            is_splice_pair = f_splice['is_splice'] or r_splice['is_splice']
                            
                            fragment_distance = calculate_appropriate_distance(
                                pos1, len1, pos2, len2,
                                is_splice_pair,
                                f_splice['intron_len'],
                                r_splice['intron_len']
                            )
                        else:
                            fragment_distance = abs(pos1 - pos2) + (len1 if pos1 > pos2 else len2)
                        
                        if fragment_distance <= global_threshold:
                            goodDct[ref]['good'] += 2
        else:
            for query_name, alignments in mateDct.items():
                if query_name == '_splice_info':
                    continue
                    
                pos1, len1 = alignments['f']
                pos2, len2 = alignments['r']
                
                if query_name in splice_info:
                    f_splice = splice_info[query_name].get('f', {'is_splice': False, 'intron_len': 0})
                    r_splice = splice_info[query_name].get('r', {'is_splice': False, 'intron_len': 0})
                    
                    is_splice_pair = f_splice['is_splice'] or r_splice['is_splice']
                    
                    fragment_distance = calculate_appropriate_distance(
                        pos1, len1, pos2, len2,
                        is_splice_pair,
                        f_splice['intron_len'],
                        r_splice['intron_len']
                    )
                else:
                    fragment_distance = abs(pos1 - pos2) + (len1 if pos1 > pos2 else len2)
                
                if fragment_distance <= global_threshold:
                    goodDct[ref]['good'] += 1

        fragments             = contig_data['contigDF'][contig_data['contigDF']['name'] == ref]['fragments'].iloc[0]
        goodDct[ref]['pGood'] = goodDct[ref]['good'] / fragments if fragments else 0

    return goodDct
