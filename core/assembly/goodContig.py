from pysam import FastaFile

class goodContig:
    def __init__(self):
        self.goodContigs     = 0
        self.badContigs      = 0
        self.goodContiglist  = []
        self.goodDct         = {'assembly': {'goodContig': 0, 'badContig': 0}}
        return
    
    def mainRun(self, args):
        assembly_data = args
        self.goodContig(assembly_data)
        return self.goodDct
    
    def goodContig(self, assembly_data):
        cutoff = assembly_data['assemblyDF']['cutoff'].iloc[0]
        fasta  = FastaFile(assembly_data['assembly'])

        self.goodContiglist = assembly_data['contigDF'][assembly_data['contigDF']['score'] > cutoff]['name'].values.tolist()
        self.goodContigs    = len(self.goodContiglist)
        self.badContigs     = len(assembly_data['contigDF']) - self.goodContigs
        
        def write_fasta(file_path, contigs):
            with open(file_path, 'w') as file:
                for contig in contigs:
                    file.write(f'>{contig}\n{fasta.fetch(contig)}\n')
            with open(file_path, 'r+') as file:
                lines = file.readlines()
                if lines and lines[-1] == '\n':
                    file.seek(0)
                    file.writelines(lines[:-1])
                    file.truncate()

        write_fasta(assembly_data['goodContig'], self.goodContiglist)
        bad_contigs = set(fasta.references) - set(self.goodContiglist)
        write_fasta(assembly_data['badContig'], bad_contigs)
        
        self.goodDct['assembly']['goodContigs'] = self.goodContigs
        self.goodDct['assembly']['badContigs']  = self.badContigs
        return