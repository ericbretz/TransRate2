from pysam import AlignmentFile

def main(mainDct):
    with AlignmentFile(mainDct['sortedBam'], 'rb') as bamfile:
        refList = bamfile.references
        mainDct['refCount']  = len(refList)
        mainDct['refList']   = refList
        mainDct['readCount'] = bamfile.count()
        for ref in refList:
            mainDct['basesDct'][ref] = {
                'bases': [],
            }
    return