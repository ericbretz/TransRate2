from pysam import AlignmentFile

def main(contig_data):
    bam_file = contig_data['dict_file']['samtools_bam']
    
    with AlignmentFile(bam_file, 'rb') as bamfile:
        refList = bamfile.references
        contig_data['refCount']  = len(refList)
        contig_data['refList']   = refList
        contig_data['readCount'] = bamfile.count()
        for ref in refList:
            contig_data['basesDct'][ref] = {
                'bases': [],
            }
    return