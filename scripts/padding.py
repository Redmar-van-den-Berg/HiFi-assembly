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

    # If seq1 is longer, we need to truncate it
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
    # If seq1 is shorter, we have to padd it
    if len(seq1) < len(seq2):
        # If the full length of seq1 matches
        if len(seq1) == end1 - start1 + 1:
            missing_before = start2 - 1
            return '-'*missing_before + seq1

def extract_hit_region(seq1, start1, end1, seq2, start2, end2):
    """ Extract the hit region for seq2 from seq1.

        Briefly, this means that if part of seq1 is a blast hit for seq2, we
        extend the region around this hit for seq1 untill it is the same size
        as seq2.

        This way, we can fin insertions or deletions in the assembly where the
        target gene should be, which are otherwise missed by BLAST.
    """
    # How much of the beginning of seq2 is missing from the blast hit
    extended_begin = start1 - min(start1, start2)

    # How much of the end of seq2 is missing from the blast hit
    extended_end = end1 + len(seq2) - end2

    # Return the slice of seq1 that contains both the original hit and the
    # extension
    return seq1[extended_begin:extended_end]
