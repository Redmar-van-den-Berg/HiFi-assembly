#!/usr/bin/env python3

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as fin, open(outfile, 'wt') as fout:
    for line in fin:
        if line.startswith('S'):
            spline=line.strip().split('\t')
            print(f'>{spline[1]}\n{spline[2]}', file=fout)
