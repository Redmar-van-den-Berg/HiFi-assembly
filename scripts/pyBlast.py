#!/usr/bin/env python3

# Author: Redmar van den Berg
# Copyright: 2017-2018 NVWA
#
# NVWA: Nederlandse Voedsel- en Warenautoriteit
# Netherlands Food and Consumer Product Safety Authority
# https://english.nvwa.nl/
#
# Licence: GPLv3

import os
import sys
import subprocess
import shutil
import tempfile
import unittest
import copy
import itertools
from functools import partial
from contextlib import contextmanager
from pprint import pprint

import Bio
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class pyBlast():
    """ A full python wrapper around makeblastdb and ncbi-blast+

    Uses temporary files and context manager to make sure that nothing
    remains on the file system after it has ran.

    Please note, this should not be used for very large datasets since
    both the database creation and the blast search itself are done on the fly
    for every invocation

    """
    def __init__(self, cmd, rm_tmp=True, verbose=False):
        self.cmd = cmd
        # We copy the fasta target over, since we will override the db variable
        # with the actual path to the makeblastdb database
        self.target = cmd.db
        self.rm_tmp = rm_tmp
        self.verbose = verbose

        if cmd.program_name == 'blastn':
            self.dbtype = 'nucl'
        elif cmd.program_name in ['blastp', 'blastx', 'deltablast', 'psiblast'
            'rpsblast', 'rpstblastn', 'tblastn']:
            self.dbtype = 'prot'
        else:
            raise RuntimeError(
                'Unknown blast command "{}"'.format(cmd.program_name)
            )

        # Query file name
        file_name = os.path.basename(self.target)
        # Temporary files
        self.tempdir = tempfile.mkdtemp()
        self.blastout = os.path.join(self.tempdir, file_name + '.xml')
        self.blastdb = os.path.join(self.tempdir, file_name + '.blastdb')

        # Update the cmd
        self.cmd.outfmt = 5
        self.cmd.db = self.blastdb
        self.cmd.out = self.blastout


    def __enter__(self):
        self.makeblastdb(self.target, self.dbtype, self.blastdb)
        self.run_blast(self.cmd)#, self.query, self.target, self.tempdir)
        self.file_in = open(self.blastout, 'r')
        return NCBIXML.parse(self.file_in)

    def __exit__(self, exc_type, exc_value, traceback):
        if self.rm_tmp:
            shutil.rmtree(self.tempdir)
        self.file_in.close()

    def makeblastdb(self, target, dbtype, blastdb, makeblastdb='makeblastdb'):

        command = '{} -in {} -input_type fasta -out {} -dbtype {}'
        command = command.format(makeblastdb, target, blastdb, dbtype)
        command = command.split()

        if dbtype not in ['nucl', 'prot']:
            raise SyntaxError('dbtype must be "nucl" or "prot"')

        if not os.path.isfile(target):
            msg = '{} target "{}" does not exist or is not a file'
            raise FileNotFoundError(msg.format(makeblastdb, target))

        try:
            if self.verbose:
                subprocess.run(command)
            else:
                subprocess.run(command, stdout=subprocess.PIPE)
        except FileNotFoundError as e:
            msg = 'Unknown command to create blast database: {}'
            e.args=(msg.format(makeblastdb),)
            raise e

    def run_blast(self, cmd):
        """ Run the ncbi blast command and return the location of the file """
        if self.verbose:
            print(cmd)
        cmd()


class pyBlastFlat(pyBlast):
    """ Yield flattened Blast Records

    Each record is guaranteed to have:
    - One Alignment per Record
    - One Hsp per Alignment

    Records without alignments are silently dropped
    Multiple Alignments or Hsps are returned in individual Records
    """

    def __init__(self, cmd, rm_tmp=True, add_mismatch=False, min_id=0,
        min_cov=0, verbose=False
    ):
        super().__init__(cmd, rm_tmp, verbose)
        self.add_mismatch = add_mismatch
        self.min_id = min_id #minimum identity to report a hit
        self.min_cov = min_cov # minimum coverage to report a hit

    def __enter__(self):
        return (record for record in self.flatten(super().__enter__()))

    def flatten(self, pb):
        for record in pb:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    flat_record = copy.deepcopy(record)
                    flat_record.alignment = copy.deepcopy(alignment)
                    flat_record.alignment.hsp = copy.deepcopy(hsp)

                    del flat_record.alignments
                    del flat_record.alignment.hsps

                    if self.add_mismatch:
                        self.mismatch(flat_record)

                    # Skip records that overlap less then min_cov
                    align_length = flat_record.alignment.hsp.align_length
                    query_length = flat_record.query_length
                    if align_length/query_length < self.min_cov:
                        continue

                    # Skip record that have an identity less then min_id
                    identities = flat_record.alignment.hsp.identities
                    query_length = flat_record.query_length
                    if identities/query_length < self.min_id:
                        continue

                    yield flat_record

    def mismatch(self, record):
        """ Return the number of mismatches in a blast hit """
        # This function only works with flat blast Records
        if hasattr(record, 'alignments'):
            raise RuntimeError('only flat Records are supported')

        alignment = record.alignment
        hsp = record.alignment.hsp

        # I don't know what the difference between these is
        assert hsp.identities == hsp.positives

        # The alignment and the query are the same length
        # We only count the mismatches
        if record.query_length == hsp.align_length:
            record.mismatch = record.query_length - hsp.identities
        # When the hit is longer then the query, we also count the gaps as
        # mismatches
        elif record.query_length < hsp.align_length:
            record.mismatch = hsp.align_length - hsp.identities
        elif record.query_length > hsp.align_length:
            record.mismatch = record.query_length - hsp.identities
        else:
            raise NotImplementedError('unknown composition of blast hit')

    def records_overlap(record1, record2):
        """ Records overlap when their targets overlap """
        if hasattr(record1, 'alignments') or hasattr(record2, 'alignments'):
            raise RuntimeError('only flat Records are supported')

        # Records cannot overlap if they are to different sequences
        if record1.alignment.hit_def != record2.alignment.hit_def:
            return False

        hsp1 = record1.alignment.hsp
        hsp2 = record2.alignment.hsp
        return _hit_overlap(hsp1, hsp2)

    @staticmethod
    def best_hits_pos(records):
        """ Find the best hits for each position on the genome """

        def overlap(record1, record2):
            return pyBlastFlat.records_overlap(record1, record2)

        def list_overlap(record1, lijst):
            """ Do any of the hits in lijst overlap hit1 """
            for record in lijst:
                if overlap(record1, record):
                    return True
            else:
                return False

        def better_hit(record1, lijst):
            """ Is hit better then every hits in lijst it overlaps with """
            # If we find record1 is better then another record,
            # we can't return True yet, because it might be worse then another
            # record it *also* overlaps with. That is why we store the boolean
            # and only return it once we have exhausted lijst
            better=False
            for record in lijst:
                if overlap(record1, record):
                    if record1.mismatch >= record.mismatch:
                        # record1 is worse or equal to an existing
                        # overlapping record
                        return False
                    elif record1.mismatch < record.mismatch:
                        # record1 is better than an existing
                        # overlapping record
                        better=True
            else:
                return better

        def remove_worse(record1, lijst):
            """ Remove all hits that overlap hit1 but have more mismatches """
            # Be sure to iterate over a copy of lijst
            for record in lijst[:]:
                is_better = record1.mismatch < record.mismatch
                if overlap(record1, record) and is_better:
                    lijst.remove(record)

        def equal_hit(hit1, lijst):
            """ Is hit equal to one or more hits in lijst """
            for hit in lijst:
                if overlap(hit1, hit) and hit1.mismatch == hit.mismatch:
                    return True
            else:
                return False

        def test_invariant(lijst):
            """ Invariant: no hit in lijst is better then an overlapping hit """
            for hit1, hit2 in itertools.combinations(lijst, 2):
                # If hit1 and 2 overlap, and they have different mismatch scores
                msg='{} and {} overlap and have different mismatch scores ({} and {})'
                if overlap(hit1, hit2) and hit1.mismatch != hit2.mismatch:
                    raise AssertionError(
                        (msg.format(hit1.query_id,
                            hit2.query_id,
                            hit1.mismatch,
                            hit2.mismatch)
                        )
                    )

        best = list()

        for record in records:
            test_invariant(best)

            # If hit is new
            if not list_overlap(record, best):
                best.append(record)
            # If hit is better then an existing hit
            elif better_hit(record, best):
                remove_worse(record, best)
                best.append(record)

        return sorted(best, key=lambda x: x.mismatch)

    @staticmethod
    def best_hits_gene(records, name=str):
        """ Return only the best hit for each gene

            The Optional function 'name' will be called on record.query to
            create the name for the gene that will be used to find the best hit
            for each gene. This comes in handy when your database contains many
            alleles of the same gene, and you are only interested in the best
            hit.

            For example, you might want to treat CMY-2_1 and CMY-2_23 the same
            if the underscore only denotes a different allele.
            In that case, specify name=lambda x: x.split('_')[0]
        """

        best = dict()
        for record in records:
            gene_name = name(record.query)
            if gene_name not in best:
                best[gene_name] = record
            elif record.mismatch < best[gene_name].mismatch:
                best[record.query] = record
        return best.values()

    @staticmethod
    def fasta(record):
        """ Return the Blast object as SeqRecord """
        # This function only works with flat blast Records
        if hasattr(record, 'alignments'):
            raise RuntimeError('only flat Records are supported')
        return SeqRecord(
            Seq(record.alignment.hsp.sbjct),
            id = record.query,
            description='{}:{}-{}'.format(
                record.alignment.hit_def,
                record.alignment.hsp.sbjct_start,
                record.alignment.hsp.sbjct_end
            )
        )


def _minmax(*args):
    """ Return the min and max of the input arguments """
    min_ = min(*args)
    max_ = max(*args)
    return(min_, max_)

def _srange(begin, end):
    """ Return a set based on range """
    return set(range(begin, end))

def _hit_overlap(hsp1, hsp2):
    """ Determine whether the hits of two hsps overlap """
    hit1_begin, hit1_end = _minmax(hsp1.sbjct_start, hsp1.sbjct_end)
    hit2_begin, hit2_end = _minmax(hsp2.sbjct_start, hsp2.sbjct_end)

    hit1_range = _srange(hit1_begin, hit1_end)
    hit2_range = _srange(hit2_begin, hit2_end)

    return not hit1_range.isdisjoint(hit2_range)

if __name__ == '__main__':
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio import SeqIO
    cmd=NcbiblastnCommandline(
        query = 'test/sul2_1_AF542061.fasta',
        db = 'test/102637-001-018_k64-contigs.fa',
        evalue = 0.001
    )
    with pyBlastFlat(cmd, rm_tmp=False, min_cov=0.5, verbose=True) as pb:
        for record in pb:
            fasta = pyBlastFlat.fasta(record)
            print(fasta)
            print(SeqIO.write(fasta, sys.stdout, 'fasta'))
