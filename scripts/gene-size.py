#!/usr/bin/env python3

import argparse

from pysam import FastaFile


def main(args):
    record = FastaFile(args.reference)

    for gene, region in zip(args.genes, args.regions):
        size = len(record.fetch(region=region))
        print(f'{gene}\t{size}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', required=True,
                        help='Reference file in FASTA format')
    parser.add_argument('--genes', required=True, nargs='+',
                        help='One or more genes')
    parser.add_argument('--regions', required=True, nargs='+',
                        help='One or more regions, same order as --genes')

    args = parser.parse_args()

    # Check to see if we get a region for every gene
    msg =  'Unequal number of genes and regions, aborting..'
    assert len(args.genes) == len(args.regions), msg

    main(args)
