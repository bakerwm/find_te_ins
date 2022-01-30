#!/usr/bin/env python

"""
Extract Insertions from Alignment file (BAM)
"""

import os
import sys
import pysam 
import argparse


def load_fa_idx(x):
    """
    fa length
    column-1: id
    column-2: length 
    """
    d = {}
    with open(x) as r:
        for line in r:
            p = line.strip().split('\t')
            d[p[0]] = p[1]
    return d


###### extract TE annotation 
def anno_te(x, ref_idx, cutoff=0.6):
    """
    Annotation of TEs
    read map to TE consensus
    
    group   covered by length
    full    >cutoff (60%?)
    p5      <cutoff (left)
    p3      <cutoff (right)
    """
    # load fa idx
    ref = load_fa_idx(ref_idx)
    # load alignment
    sam = pysam.AlignmentFile(x)
    out = []
    read = next(sam)
    for read in sam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        # strand
        strand = '-' if read.is_reverse else '+'
        # query: alignment
        q_name = read.query_name
        q_len = read.query_length
        q_align_len = read.query_alignment_length
        q_pct = '{:.4f}'.format(int(q_align_len)/int(q_len))
        q_align_pos = '{}\t{}'.format(read.query_alignment_start, read.query_alignment_end)
        # reference: alignment
        ref_name = read.reference_name
        ref_len = ref.get(ref_name, 1)
        ref_align_len = read.reference_length
        ref_pct = '{:.4f}'.format(int(ref_align_len)/int(ref_len))
        ref_align_pos = '{}\t{}'.format(read.reference_start, read.reference_end)
        # output format
        # query_name: read_id:start:end,chr:start:end
        # expand 
        r1, r2 = q_name.split(',')
        p1 = r1.split(':') + r2.split(':')
        p2 = [
            read.query_alignment_start, 
            read.query_alignment_end,
            q_pct,
            ref_name,
            read.reference_start, 
            read.reference_end,
            ref_len,
            ref_pct,
            strand,
        ]
        p = p1 + p2 # combine
        # classification
        # full,p5,p3
        if float(ref_pct) > cutoff:
            tag = 'full'
        else:
            # tag = 'part'
            tag = 'p5' if int(read.reference_end) < int(ref_len)/2 else 'p3'
        p.append(tag)
        # to str
        p = list(map(str, p))
        out.append('\t'.join(p))
    return out


def get_args():
    """ 
    Parsing arguments
    """
    example = '\n'.join([
        'Examples:',
        '$ python cnr_r1.py bam > out_ins.bed'
    ])  
    parser = argparse.ArgumentParser(
        prog='extract_ins.py',
        description='extract insertions from Alignment file: BAM',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('bam', help='BAM alignment file')
    parser.add_argument('-x', '--ref-idx', dest='ref_idx', required=True,
        help='Index file of the ref.fa (Transposon)')
    parser.add_argument('-c', '--cutoff', dest='cutoff', type=float,
        default=0.6, help='cutoff to define the Transposon insertion, 0-1, default: [0.6]')
    return parser


def main():
    args = get_args().parse_args()
    # run
    out = anno_te(args.bam, args.ref_idx, args.cutoff)
    print('\n'.join(out))
            
            
if __name__ == '__main__':
    main()