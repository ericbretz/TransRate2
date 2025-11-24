from pysam import AlignmentFile

def mainRun(args):
    contig_data = args[0]
    i           = args[1]
    bam         = AlignmentFile(contig_data['dict_file']['samtools_bam'], 'rb')
    refs        = bam.references[i::contig_data['threads']]
    baseDct     = {ref: {'basesUncovered': 0, 'other': {'bases': []}} for ref in refs}

    for ref in refs:
        length   = bam.get_reference_length(ref)
        coverage = [0] * length
        
        for read in bam.fetch(reference=ref):
            if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                continue

            start = read.reference_start
            end   = read.reference_end
            
            for pos in range(start, min(end, length)):
                coverage[pos] += 1
        
        baseDct[ref]['basesUncovered'] = coverage.count(0)
        baseDct[ref]['other']['bases'] = coverage
        
    return baseDct
