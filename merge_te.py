#!/usr/bin/env python

"""
Merge Insertions by window

extract, merge insertions
full version
output of anno_te.py
16-column format, tab-delimited
read_id start end chr start end ins_start ins_end ins_pct te_id te_start te_end te_length te_pct strand tag

criteria
merge insertions within 100bp window
sorted by chr,start: sort -k4,4 -k5,5n 
"""

import os
import sys
import argparse
from collections import Counter


class MergeINS(object):
    """
    Merge Insertion by window
    """
    def __init__(self, x, **kwargs):
        self.x = x
        args = {
            'window': 100,
            'out_dir': None
        }
        args.update(kwargs)
        for k,v in args.items():
            if not hasattr(self, k):
                setattr(self, k, v)


    def is_sorted(self, a, b):
        """
        sorted by chr,start
        """
        if len(a) == len(b):
            f1 = a[3] <= b[3]
            f2 = int(a[4]) <= int(b[4]) if a[3] == b[3] else True
            out = f1 and f2
            # out = a[3] <= b[3] and int(a[4]) <= int(b[4])
        else:
            out = False
        return out


    def in_window(self, a, b, w=100):
        """
        check if a,b within the same window
        chr   : column-4
        start : column-5
        """
        if len(a) == len(b):
            out = a[3] <= b[3] and abs(int(a[4]) - int(b[4])) <= w
        else:
            out = len(a) == 0 or len(b) == 0
        return out


    def merge_chunk(self, x):
        """
        x : list of records (see above message)
        
        chunk of INS records:
        merge ins, add read count
        output:
        chr start end name 255 strand count
        """
        # one TE per chunk
        te = Counter([i[9] for i in x]).most_common(1) # top ranked TE 
        top_te = te[0][0] # most common ids
        x2 = [i for i in x if i[9] == top_te] # filter chunk by TEs
        if len(x2) > 0:
            s = min([i[4] for i in x2]) # start
            e = max([i[5] for i in x2]) # end
            t_pct = max([i[13] for i in x2]) # max te_pct
            t_tag = [i[-1] for i in x2] #
            if 'full' in t_tag:
                tag = 'full'
            else:
                gmax = Counter(t_tag).most_common(1)
                tag = gmax[0][0]        
            out = [
                x2[0][3],
                s,
                e,
                '{}:{}-{}'.format(x2[0][3], s, e),
                255,
                '+',
                top_te,
                x2[0][12],
                t_pct,
                len(x2),
                tag,            
            ]
            # convert to str
            out = list(map(str, out))
            return '\t'.join(out)
        
    
    def merge_ins(self):
        i = 0
        j = 0
        out = []
        with open(self.x) as r:
            chunk = []
            p_pre = []
            for line in r:
                i += 1
                p = line.strip().split('\t')
                if not p[-1] == 'full':
                    continue # only full
                if len(p_pre) == 0:
                    p_pre = p # init
                    chunk.append(p) # 1st item
                else:
                    # check
                    if not self.is_sorted(p_pre, p):
                        sys.exit('file is not sorted by: sort -k4,4 -k5,5n in line-{}'.format(i))
                    # add to chunk
                    if self.in_window(p_pre, p, self.window):
                        p_pre = p # update p_pre
                        chunk.append(p) # add to chunk
                    else:
                        # do-something in the chunk
                        ins = self.merge_chunk(chunk)
                        j += 1
                        if isinstance(ins, str):
                            out.append(ins)
                        p_pre = p # init
                        chunk = [p] # init
            # last one
            if len(chunk) > 0:
                ins = self.merge_chunk(chunk)
                if isinstance(ins, str):
                    out.append(ins)
        return out
    
    
    def run(self):
        return self.merge_ins()
        

def get_args():
    """ 
    Parsing arguments
    """
    example = '\n'.join([
        'Examples:',
    ])  
    parser = argparse.ArgumentParser(
        prog='merge_te_ins.py',
        description='merge te insertions by window',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('x', help='insertion file, output of anno_te.py')
    parser.add_argument('-w', '--window', dest='window', type=int, default=100,
        help='window size to merge TEs, default [100]')
    return parser


def main():
    args = get_args().parse_args()
    out = MergeINS(args.x, window=args.window).run()
    print('\n'.join(out))
            
            
if __name__ == '__main__':
    main()

