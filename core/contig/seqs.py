from pysam import AlignmentFile

def calculate_cigar_accuracy(cigar):
    if not cigar:
        return 0.0
    
    aligned_bases   = 0
    explicit_errors = 0
    
    for op, length in cigar:
        if op == 0:
            aligned_bases   += length
        elif op == 1:
            explicit_errors += length
        elif op == 2: 
            explicit_errors += length
        elif op == 8:
            explicit_errors += length
            aligned_bases   += length
        elif op == 7:
            aligned_bases   += length
    
    if aligned_bases == 0:
        return 0.0
    accuracy = max(0.0, (aligned_bases - explicit_errors) / aligned_bases)
    
    return min(1.0, accuracy)

def mainRun(args):
    contig_data = args[0]
    i           = args[1]
    bam         = AlignmentFile(contig_data['dict_file']['samtools_bam'], 'rb')
    refs        = bam.references[i::contig_data['threads']]
    
    seqDct = {ref: {'true'             : 0,
                    'seqCount'         : 0,
                    'softclipped'      : 0,
                    'softclippedLength': 0,
                    'softTotalLength'  : 0,
                    'totalReads'       : 0,
                    'pSoftclipped'     : 0.0,
                    'pSeqTrue'         : 0} for ref in refs}

    cigdict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}

    for ref in refs:
        softclipped       = 0
        total_reads       = 0
        total_length      = 0
        softclippedLength = 0
        max_read_length   = 0
        
        for fetch in bam.fetch(reference=ref):
            max_read_length = max(max_read_length, fetch.query_length)
            
            try:
                total_reads    += 1
                name_r          = fetch.reference_name
                cigar_accuracy  = calculate_cigar_accuracy(fetch.cigar)
                
                if cigar_accuracy > 0:
                    seqDct[name_r]['seqCount'] += 1
                    seqDct[name_r]['true']     += cigar_accuracy
                    seqDct[name_r]['pSeqTrue']  = seqDct[name_r]['true'] / seqDct[name_r]['seqCount']
                
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
                            softclipped       += 1
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
