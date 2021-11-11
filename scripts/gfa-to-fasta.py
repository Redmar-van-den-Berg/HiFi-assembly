#!/usr/bin/env python3

import sys

infiles = sys.argv[1:]

for filename in infiles:
    with open(filename) as fin:
        for line in fin:
            if line.startswith('S'):
                spline=line.strip().split('\t')
                print(f'>{spline[1]}\n{spline[2]}')
