#!/usr/bin/env python

"""
Add quant to bed

BED: 
chr2L   62010   62016   chr2L:62010-62016       255     +       copia   5143    1.0000  28      full

quant_txt:
Geneid  Chr     Start   End     Strand  Length  results/align/minimap2//ONT_sample-1.bam
chr2L:62010-62016       chr2L   62011   62016   +       6       133
"""


import os
import sys


def load_fc(x):
    """
    load featureCounts txt
    """
    d = {}
    with open(x) as r:
        for line in r:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            d[p[0]] = p[-1] 
    return d


def add_quant(bed, quant):
    """
    Add quant to bed file
    """
    q = load_fc(quant)
    with open(bed) as r:
        for line in r:
            p = line.strip().split('\t')
            c = q.get(p[3], 0) # count
            try:
                pct = '{:.4f}'.format(int(p[9])/float(c))
            except:
                pct = '1000' # 
            p.insert(10, str(c)) # quant
            p.insert(11, pct) # pct
            print('\t'.join(p))


def main():
    if len(sys.argv) < 3:
        print('quant_bed.py {bed} {quant} > out.bed')
        sys.exit(1)
        
    add_quant(sys.argv[1], sys.argv[2])
    

if __name__ == '__main__':
    main()