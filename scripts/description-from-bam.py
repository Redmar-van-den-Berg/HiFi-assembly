#!/usr/bin/env python3

import argparse
import pysam

from normalizer import description_extractor

extractor = description_extractor.description_extractor

def get_reference_region(record, region):
    """ Extract region from the FastaFile record"""
    return record.fetch(region=region)

def get_contigs_region(bamfile, region):
    """ Extract contigs that overlap region from the bamfile """
    bam = pysam.AlignmentFile(bamfile)
    for contig in bam.fetch(region=region):
        yield contig

def main(args):
    # Get the sequence of the reference
    reference = pysam.FastaFile(args.reference)
    ref_seq = get_reference_region(reference, args.region)

    # Get the contigs that overlap region from the bamfile
    for contig in get_contigs_region(args.bam, args.region):
        name = contig.qname
        seq = contig.seq
        # seq can be empty, and supplementary alignments are usually badly
        # mapped
        if seq and not contig.is_supplementary:
            description = extractor(ref_seq, contig.seq)
            print(f'{name}\t{description}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, help='Input bam file')
    parser.add_argument('--reference', required=True, help='Reference fasta file')
    parser.add_argument('--region', required=True,
                        help='Region to extract from reference and bam file')

    args = parser.parse_args()
    main(args)
