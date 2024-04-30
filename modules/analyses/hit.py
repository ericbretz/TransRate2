class Hit:
    def __init__(self, lsti, query, reference):
        # if len(lsti) < 14:
        #     raise RuntimeError("unexpected number of columns")
        self.lsti       = lsti
        self.query      = lsti[0].split()[0]
        self.target     = lsti[1].split()[0]
        self.id         = lsti[2]
        self.alnlen     = int(lsti[3])
        self.mismatches = int(lsti[4])
        self.gaps       = int(lsti[5])
        self.qstart     = int(lsti[6])
        self.qend       = int(lsti[7])
        self.tstart     = int(lsti[8])
        self.tend       = int(lsti[9])
        self.evalue     = float(lsti[10])
        self.bitscore   = float(lsti[11])
        self.qlen       = float(lsti[12])
        self.tlen       = float(lsti[13])

    def lst(self):
        s = {'query': self.query, 'target': self.target, 'id': self.id, 'alnlen': self.alnlen, 'mismatches': self.mismatches,
                        'gaps': self.gaps, 'qstart': self.qstart, 'qend': self.qend, 'tstart': self.tstart, 'tend': self.tend,
                        'evalue': self.evalue, 'bitscore': self.bitscore, 'qlen': self.qlen, 'tlen': self.tlen}
        return s
