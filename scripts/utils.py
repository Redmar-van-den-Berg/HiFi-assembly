#!/usr/bin/env python3

def extract_hit_region(seq1, start1, end1, seq2_length, start2, end2):
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
    extended_end = end1 + seq2_length - end2

    # Return the slice of seq1 that contains both the original hit and the
    # extension
    return seq1[extended_begin:extended_end]
