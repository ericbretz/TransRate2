import math
import sys
import pysam
from multiprocessing import Pool
import sys

def process_task(args):
    split = args[0]
    bamfile = args[1]
    threads = args[2]
    bamdct = args[3]
    bam = pysam.AlignmentFile(bamfile, 'rb')
    refs = bam.references[split::threads]
    tmp_dct = {ref: {'stats': {'p_notsegmented': 0}} for ref in refs}

    for ref in refs:
        uncovered = bamdct[ref]['coverage']['bases']
        pnotseg = setPNotSegmented(len(uncovered), uncovered, 0.9)
        tmp_dct[ref]['stats']['p_notsegmented'] = round(pnotseg, 4)

    return tmp_dct

def sgmt(bamfile, bamdct, threads, single):
    with Pool(threads) as p:
        results = p.map(process_task, [[i, bamfile, threads, bamdct] for i in range(threads)])

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

def setPNotSegmented(ref_length, coverage, nullprior):
    states = [0]*30
    bin_width = (ref_length // 30) + 1  # ceil
    counter = 0
    pos = 0
    total = 0
    for i in range(ref_length):
        total += coverage[i]
        counter += 1
        if counter == bin_width or i == ref_length - 1:
            total = max(1, total)
            states[pos] = min(24, max(0, int(math.log2(total / counter))))
            pos += 1
            counter = 0
            total = 0
    segmenter = Segmenter(states, nullprior)
    return segmenter.prob_k_given_R(0)  # prob of 0 change points (1 segment)

class Segmenter:
    def __init__(self, seq, nullprior=0.7):
        self.seq = seq
        self.nullprior = nullprior
        self.states_ = [0]*24
        self.pRk_ = [-1.0, -1.0]
        self.maxk = 1
        self.load_states()

    def load_states(self):
        for i in range(len(self.seq)):
            if self.seq[i] < 24:
                self.states_[self.seq[i]] += 1
            else:
                self.states_[23] += 1

    def marginal_likelihood_R(self):
        sum_ = 0.0
        for k in range(2):
            sum_ += self.prior_k(k) * self.prob_R_given_k(k)
        return sum_

    def prior_k(self, k):
        if k == 0:
            return self.nullprior
        else:
            return (1.0 - self.nullprior)

    def prob_R_given_k(self, k):
        result = self.pRk_[k]
        if result == -1.0:
            if k == 0:
                result = self.prob_R_given_zero_k()
            elif k == 1:
                result = self.prob_R_given_unit_k()
            self.pRk_[k] = result
        return result

    def prob_R_given_zero_k(self):
        if self.pRk_[0] != -1.0:
            return self.pRk_[0]
        result = self.prob_R_given_k_rhs(self.states_, len(self.seq))
        self.pRk_[0] = result
        return result

    def prob_R_given_unit_k(self):
        if self.pRk_[1] != -1.0:
            return self.pRk_[1]
        lstates = self.states_.copy()
        total = len(self.seq)
        pmat = [[0.0]*total for _ in range(total)]
        for i in range(total):
            segstates = lstates.copy()
            for j in range(total-1, i-1, -1):
                segresult = self.prob_R_given_k_rhs(segstates, j+1-i)
                pmat[i][j] = segresult
                segstates[self.seq[j]] -= 1
            lstates[self.seq[i]] -= 1
        pA = max((1.0 / (total - 1)), sys.float_info.min)
        result = 0.0
        for i in range(total - 1):
            left = pmat[0][i]
            right = pmat[i+1][total-1]
            product = left * right * pA
            result += product
        self.pRk_[1] = result
        return result

    def prob_R_given_k_rhs(self, states, length):
        nstates = len(states)
        l = math.gamma(nstates)
        upper = 1.0
        for p in states:
            upper *= math.gamma(p + 1)
        lower = math.gamma(length + nstates)
        return l * (upper / lower)

    def prob_k_given_R(self, k):
        pRk = self.prob_R_given_k(k)
        pk = self.prior_k(k)
        mlR = self.marginal_likelihood_R()
        return pRk * pk / mlR
