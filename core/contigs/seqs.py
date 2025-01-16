from pysam import AlignmentFile

def mainRun(args):
    main = args[0]
    i    = args[1]
    bam  = AlignmentFile(main['sortedBam'], 'rb')
    refs = bam.references[i::main['threads']]
    seqDct = {ref: {'true'               : 0,
                    'seqCount'           : 0,
                    'softclipped'        : 0,
                    'softclippedLength'  : 0,
                    'softTotalLength'    : 0,
                    'totalReads'         : 0,
                    'pSoftclipped'       : 0.0,
                    'pSeqTrue'           : 0} for ref in refs}

    cigdict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}

    for ref in refs:
        softclipped              = 0
        total_reads              = 0
        total_length             = 0
        softclippedLength        = 0
        max_read_length          = 0
        for fetch in bam.fetch(reference=ref):
            max_read_length = max(max_read_length, fetch.query_length)
            try:
                # Could not figure out the point of scale from TR1
                total_reads                          += 1
                name_r                                = fetch.reference_name
                len_r                                 = fetch.reference_length
                nm                                    = fetch.get_tag('NM', with_value_type = True)
                seqDct[name_r]['seqCount']           += 1
                # scale                                 = (len_r - 35) / len_r
                # seqDct[name_r]['true']               += (((len_r - nm[0]) / len_r) - scale) * (1 / (1 - scale))
                seqDct[name_r]['true']               += (len_r - nm[0]) / len_r
                seqDct[name_r]['pSeqTrue']            = seqDct[name_r]['true'] / seqDct[name_r]['seqCount']
                
            except:
                pass
            try:
                cigar         = fetch.cigartuples
                f_length      = fetch.query_length
                total_length += f_length

                if cigar:
                    for cig in cigar:
                        letter = cigdict[cig[0]]
                        if letter == 'D':
                            continue
                        elif letter == 'S':
                            softclipped += 1
                            softclippedLength += cig[1]
            except:
                pass

        try:
            seqDct[ref]['softclipped']       = softclipped
            seqDct[ref]['softclippedLength'] = softclippedLength
            seqDct[ref]['softTotalLength']   = total_length
            seqDct[ref]['totalReads']        = total_reads
            seqDct[ref]['max_read_length']   = max_read_length
            if total_reads:
                seqDct[ref]['pSoftclipped']  = float(softclipped / total_reads)
            
        except:
            continue

    return seqDct