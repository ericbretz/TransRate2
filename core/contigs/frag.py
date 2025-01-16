from pysam import AlignmentFile

def mainRun(args):
    main = args[0]
    i = args[1]
    bam = AlignmentFile(main['sortedBam'], 'rb')
    refs = bam.references[i::main['threads']]
    fragDct = {ref: {'fragments': 0, 
                     'properPair': 0,
                     'bothMapped': 0,
                     'bridges': 0} for ref in refs}

    for ref in refs:
        for fetch in bam.fetch(reference=ref):
            try:
                name = fetch.reference_name
                frag_info = fragDct[name]
                is_read1 = fetch.is_read1
                is_mapped = fetch.is_mapped
                mate_is_mapped = fetch.mate_is_mapped
                is_proper_pair = fetch.is_proper_pair
                reference_id = fetch.reference_id
                next_reference_id = fetch.next_reference_id

                if main['mode'] == 2:
                    if (is_read1 and is_mapped) or (fetch.is_read2 and not mate_is_mapped):
                        frag_info['fragments'] += 1
                    if is_read1 and fetch.is_paired and is_proper_pair:
                        frag_info['properPair'] += 1
                    if is_read1 and is_mapped and mate_is_mapped:
                        frag_info['bothMapped'] += 1
                    if reference_id != next_reference_id and is_read1:
                        frag_info['bridges'] += 1
                else:
                    if is_mapped:
                        frag_info['fragments'] += 1
            except KeyError:
                pass
    return fragDct