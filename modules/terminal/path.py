import os
class PathCheck:
    def __init__(self, assembly=None, left=None, right=None, output=None, reference=None, bam=None):
        self.assembly  = assembly
        self.left      = left
        self.right     = right
        self.reads     = [left, right] if left and right else [left] if left else [right] if right else None
        self.reference = reference
        self.output    = output
        self.bam       = bam

    def check(self):
        if self.assembly:
            for a in self.assembly:
                if not os.path.exists(a):
                    return a

        if self.reads:
            for r in self.reads:
                if not os.path.exists(r):
                    return r
                
        if self.reference:
            if not os.path.exists(self.reference):
                return self.reference
            
        if self.output:
            if not os.path.exists(self.output):
                os.makedirs(self.output)
        return None