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
        # If the hit is not the full length of seq2
        if end1-start1+1 < len(seq2):
            # How much of seq2 is missing before the hit?
            missing_before = end2 - start2 - 1
            seq_start = start1 - missing_before
            # How much of seq2 is missing after the hit
            missing_after = len(seq2) - end2
            seq_end = end1 + missing_after
            return seq1[seq_start-1:seq_end]
        else:
            return seq1[start1-1:end1]
