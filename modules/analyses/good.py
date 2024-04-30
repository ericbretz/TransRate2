import math
import pysam
from multiprocessing import Pool


def process_task(args):
    split = args[0]
    bamfile = args[1]
    threads = args[2]
    realistic_distance = args[3]
    bam = pysam.AlignmentFile(bamfile, 'rb')
    refs = bam.references[split::threads]
    tmp_dct = {ref: {'stats': {'good': 0}, 'good' : {'goodlst' : []}} for ref in refs}

    for ref in refs:
        for fetch in bam.fetch(reference=ref):
            if not fetch.is_mapped:
                continue

            name_r = fetch.reference_name

            if not (fetch.is_read1 and fetch.mate_is_mapped):
                continue

            if name_r != fetch.next_reference_name:
                continue

            ldist = max(fetch.reference_start-fetch.next_reference_start,
                        fetch.next_reference_start-fetch.reference_start)
                        
            if ldist > realistic_distance:
                continue

            is_reversed = fetch.is_reverse
            is_mate_reversed = fetch.mate_is_reverse

            if not is_reversed and is_mate_reversed:

                if fetch.reference_start < fetch.next_reference_start:
                    tmp_dct[name_r]['stats']['good'] += 1
                    if fetch.reference_name not in tmp_dct[name_r]['good']['goodlst']:
                        tmp_dct[name_r]['good']['goodlst'].append(fetch.reference_name)

            elif is_reversed and not is_mate_reversed:

                if fetch.next_reference_start < fetch.reference_start:
                    tmp_dct[name_r]['stats']['good'] += 1
                    if fetch.reference_name not in tmp_dct[name_r]['good']['goodlst']:
                        tmp_dct[name_r]['good']['goodlst'].append(fetch.reference_name)
                    
    return tmp_dct

def good(bamfile, bamdct, threads, single):
    count    =  0
    name     =  ''
    prev     =  ''
    len1     = -1
    len2     = -1
    pos1     = -1
    pos2     = -1
    mean     =  0
    s        =  0
    m        = -1
    realistic_distance = 0

    with pysam.AlignmentFile(bamfile, 'rb') as bam:

        for fetch in bam.fetch():
            if count >= 10000:
                break
            if name != "":
                prev = name
                pos2 = pos1
                len2 = len1

            if not fetch.is_secondary and not fetch.is_supplementary:
                name = fetch.query_name
                pos1 = fetch.reference_start
                len1 = fetch.query_length

                if prev == name:
                    if pos1 >= 0 and pos2 >= 0:
                        is_reversed = fetch.is_reverse
                        is_mate_reversed = fetch.mate_is_reverse

                        if not is_reversed and is_mate_reversed:
                            if pos1 > fetch.next_reference_start:
                                continue
                        elif is_reversed and not is_mate_reversed:
                            if fetch.next_reference_start > pos1:
                                continue

                        if pos1 > pos2:
                            fragment = pos1 - pos2 + len1
                        else:
                            fragment = pos2 - pos1 + len2
                        if count > 0:
                            mn = m + (fragment - m) / count
                            s = s + ((fragment - m) * (fragment - mn))
                            m = mn
                        else:
                            m = fragment
                        mean += fragment
                        count += 1

        mean = mean / count
        s = math.sqrt(s / (count - 1))
        realistic_distance = int((3 * s) + mean)

        with Pool(threads) as p:
            results = p.map(process_task, [[i, bamfile, threads, realistic_distance] for i in range(threads)])

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