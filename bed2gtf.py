#!/usr/bin/env python3

import os, sys

if len(sys.argv) < 2:
	sys.exit('bed2gtf.py <in.bed>')

bed = sys.argv[1]

n = 0
with open(bed) as fi:
    for line in fi:
        tab = line.strip().split('\t')
        if len(tab) < 3:
            continue
        n += 1 #
        # bed3 or bed6
        if len(tab) > 3:
            name = tab[3]
            strand = tab[5]
        else:
            name = 'g{:06d}'.format(n)
            strand = '+'
        start = int(tab[1]) + 1
        end = tab[2]
        des = 'gene_id "{}"; gene_name "{}";'.format(name, name)
        gtf = '\t'.join([
            tab[0],
            'bed',
            'gene',
            str(start),
            end,
            '.',
            strand,
            '.',
            des])
        print(gtf)