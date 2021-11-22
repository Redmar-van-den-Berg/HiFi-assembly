#!/usr/bin/env python3

def extract_hit_region(seq1, start1, end1, seq2_length, start2, end2):
    """ Extract the hit region for seq2 from seq1.

        Briefly, this means that if part of seq1 is a blast hit for seq2, we
        extend the region around this hit for seq1 until it is the same size
        as seq2.

        This way, we can find insertions or deletions in the assembly where the
        target gene should be, which are otherwise missed by BLAST.
    """
    # How much of the beginning of seq2 is missing from the blast hit
    extended_begin = start1 - min(start1, start2)

    # How much of the end of seq2 is missing from the blast hit
    extended_end = end1 + seq2_length - end2

    # Return the slice of seq1 that contains both the original hit and the
    # extension
    return seq1[extended_begin:extended_end]

def extend_hit_reference(seq1, start1, end1, seq2, start2, end2):
    """ Extend the hit on seq1, assuming missing content is identical to seq2


        If the hit covers all of seq1, but seq1 is shorter than seq2, assume
        the missing sequences from seq1 are identical, and add those sequences.
    """
    # If the hit covers all of seq2, we simply return the region from seq1 that
    # matches
    if len(seq2) == end1 - start1 + 1:
        return seq1[start1 -1: end1]

    # If the full length of seq1 matches
    if len(seq1) == end1 - start1 + 1:
        # If seq2 extend beyond the blast hit
        if len(seq2) > end2:
            return seq1[start1 - 1: end1] + seq2[end2:]

    # If the hit is smaller than seq1
    return seq1[start1 - 1: end1]
