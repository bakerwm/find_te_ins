#!/usr/bin/env python

"""
convert bed8 to fasta
bed8 is the output of extract_ins.py

bed8: bed6 + {length} + {seq}
"""


import os
import sys

def bed2fa(x, out=None):
    """
    bed8: (last column store sequence)
    bed: column-4 => name
    bed: last column => seq
    """
    if isinstance(x, str):
        if out is None:
            out = x.replace('.bed', '.fa')
        # run
        with open(x) as r, open(out, 'wt') as w:
            for line in r:
                p = line.strip().split('\t')
                w.write('>{}\n{}\n'.format(p[3], p[-1]))
                
def main():
    if len(sys.argv) < 2:
        print('Usage: bed2fa.py in.bed [out.fa]')
        sys.exit(1)
    in_bed = sys.argv[1]
    out_fa = sys.argv[2] if len(sys.argv) > 2 else None
    bed2fa(in_bed, out_fa)
    

if __name__ == '__main__':
    main()
    
# EOF