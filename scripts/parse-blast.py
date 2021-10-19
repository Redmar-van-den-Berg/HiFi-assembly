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


def better(record1, record2):
    """ Is record1 better than record2 """
    return record1.alignment.hsp.score > record2.alignment.hsp.score


def get_best_hits(pb):
    """ Return the best hit for each contig in --database """
    best = dict()

    for record in pb:
        # The subject name is the name of the contig
        name = record.alignment.hit_def
        if name not in best:
            best[name] = record
        elif better(record, best[name]):
            best[name] = record

    # Once we are done, we can just return it as a list of records
    return list(best.values())


def record_to_fasta(record):
    """ Return a record in fasta format

    Includes some messing about with the blast data to put information about
    the hit into the fasta header

    A hit for CYP2D6 on contig utg1 for the first 1000 bases would get the
    following fasta header:
    >utg1:1-1000 (CYP2D6)
    """
    alignment = record.alignment
    hsp = alignment.hsp
    name = f'{alignment.hit_def}:{hsp.sbjct_start}-{hsp.sbjct_end}'
    seq = hsp.sbjct
    return f'>{name} ({record.query})\n{seq}'

def main(args):
    cmd = NcbiblastnCommandline(query=args.query, db=args.database)

    with pyBlastFlat(cmd) as pb:
        # Get the best hit for each contig (each contig is only present once)
        records = get_best_hits(pb)

        # If json output was requested
        if args.json:
            data = [jsonify(record) for record in records]
            with open(args.json, 'w') as fout:
                print(json.dumps(data, indent=True), file=fout)

        # If fasta output was requested
        if args.fasta:
            with open(args.fasta, 'w') as fout:
                for record in records:
                    print(record_to_fasta(record), file=fout)

        # If we need to write the genes
        if args.genes:
            for record in records:
                # Determine the name of the query, up to the first space
                query_name = record.query.split(' ')[0]
                fname = f'{args.genes}/{query_name}.fasta'
                with open(fname, 'a') as fout:
                    print(record_to_fasta(record), file=fout)


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
    parser.add_argument('--genes', required=False,
        help=(
            'Output the best hit for each sequence in database to this folder.'
            'Use the names of the sequences in query as file name.')
    )

    args = parser.parse_args()
    main(args)
