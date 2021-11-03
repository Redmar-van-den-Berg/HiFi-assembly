#!/usr/bin/env python3

def add_padding(seq1, start1, end1, seq2, start2, end2):
    """ 'Align' seq1 to seq2 by adding padding. This makes use of the sequence
        itself, as wel as  the 'start' and 'end' of the blast hits

        Note that the start and end are in the format from blast and NOT python, namely:
        - 1 based
        - closed
    """

    # If both sequences are the same length
    if len(seq1) == len(seq2):
        return seq1

    # If seq1 is longer
    if len(seq1) > len(seq2):
        return seq1[start1-1:end1]
