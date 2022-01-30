#!/usr/bin/env python

"""
Extract Insertions from Alignment file (BAM)
"""


import os
import sys
import argparse
import pysam
from hiseq.utils.utils import update_obj, log


class Cigar(object):
    """
    x : int, 0-8
    
    cigar table:
    
    BAM  Op  query  ref
    0    M   yes    yes
    1    I   yes    no
    2    D   no     yes
    3    N   no     yes
    4    S   yes    no 
    5    H   no     no 
    6    P   no     no 
    7    =   yes    yes
    8    X   yes    yes
    
    # convert CIGAR BAM to Op
    """
    def __init__(self, x):
        self.x = self.to_str(x)
        self.on_ref = self.x in 'MDN=X'
        self.on_query = self.x in 'MIS=X'
        
        
    def to_str(self, x):
        """
        MIDNSHP=X 
        012345678
        """
        out = None
        if isinstance(x, str):
            if x in 'MIDNSHP=X':
                out = x
        elif isinstance(x, int):
            if x in range(9):
                out = 'MIDNSHP=X'[x]
        return out


class GetIns(object):
    """
    Fetch insertions from alignment, pysam.AlignedSegment
    """
    def __init__(self, x, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.x = x
        # default arguments
        args = {
            'min_size': 300,
            'max_size': 14000
        }
        self = update_obj(self, args, force=False)
        
        
    def get_ins(self):
        """
        a is CIGAR (tuple format) 
        """
        c_ref = self.x.reference_start
        c_query = 0
        out = [] # output
        for a,b in self.x.cigar:
            if a == 1 and b > self.min_size and b < self.max_size:
                # insertion: on query, not on ref
                # export TE insertion
                ## coordinates on ref
                strand = '-' if read.is_reverse else '+'
                bed = [
                    read.reference_name, 
                    c_ref, 
                    c_ref + 1, 
                    '{}:{}'.format(read.reference_name, c_ref),
                    255, 
                    strand
                ]
                ## sequences on query
                q_s = c_query
                q_e = q_s + b
                if strand == '+':
                    ins = read.query_sequence[q_s:q_e]
                else:
                    s = read.query_sequence[::-1]
                    ins = s[q_s:q_e]
                    ins = s[::-1]
                ## output
                bed = bed + [b, ins] 
                bed = list(map(str, bed))
                out.append('\t'.join(bed))
                ## update 
                c_query += b
            else:
                cigar = Cigar(a)
                if cigar.on_ref:
                    c_ref += b
                if cigar.on_query:
                    c_query += b
        return out
    

## arguments
def get_args():
    """ 
    Parsing arguments
    """
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ python cnr_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -g dm6',
        '2. for specific index',
        '$ python cnr_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te',
    ])  
    parser = argparse.ArgumentParser(
        prog='extract_ins.py',
        description='extract insertions from Alignment file: BAM',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('bam', required=True,
        help='BAM alignment file')
    parser.add_argument('-m', '--min-size', dest='min_size', type=int,
        default=300, help='minimal size of insertion, default [300]')
    parser.add_argument('-M', '--max-size', dest='max_size', type=int,
        default=14000, help='maximal size of insertion, default [14000]')
    return parser


def main():
    args = vars(get_args().parse_args())
    # run
    sam = pysam.AlignmentFile(args['bam'])
    for read in sam:
        ins = GetIns(read).get_ins()
        if len(ins) > 0:
            print('\n'.join(ins))


if __name__ == '__main__':
    main()