#!/usr/bin/env python3

import argparse
import json
import sys

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO

from pyBlast import pyBlastFlat

def jsonify(record):
    """ Fetch some useful values from the alignment record """
    # Blast records are pretty deeply nested
    alignment = record.alignment
    hsp = alignment.hsp

    # Information about the hit
    data = {
            'score': hsp.score,
            'expect': hsp.expect,
            'positives': hsp.positives,
            'gaps': hsp.gaps,
            'aln_length': hsp.align_length,
            'mismatch': hsp.align_length - hsp.identities
    }

    # Information about the query
    query = {
            'name': record.query,
            'start': hsp.query_start,
            'end': hsp.query_end,
            'sequence': hsp.query
    }

    # Information about the subject
    sbjct = {
            'name': alignment.hit_def,
            'start': hsp.sbjct_start,
            'end': hsp.sbjct_end,
            'sequence': hsp.sbjct
    }

    data['query'] = query
    data['sbjct'] = sbjct

    return data


def main(args):
    cmd = NcbiblastnCommandline(query=args.query, db=args.database)

    with pyBlastFlat(cmd) as pb:
        # Just yank all records into memory
        records = list(pb)

        # If json output was requested
        if args.json:
            data = [jsonify(record) for record in records]
            with open(args.json, 'w') as fout:
                print(json.dumps(data, indent=True), file=fout)

        # If fasta output was requested
        if args.fasta:
            with open(args.fasta, 'w') as fout:
                for record in records:
                    SeqIO.write(pyBlastFlat.fasta(record), fout, 'fasta')
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', required=True,
        help='Fasta file to blast against the database')
    parser.add_argument('--database', required=True,
        help='Fasta file to use as a database')
    parser.add_argument('--json', required=False,
        help='Output blast results in JSON format')
    parser.add_argument('--fasta', required=False,
        help='Output blast results in FASTA format')

    args = parser.parse_args()
    main(args)

