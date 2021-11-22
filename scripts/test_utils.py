import pytest

from utils import extract_hit_region, extend_hit_reference

# These are the test cases for extracting regions of seq2 from seq1. This
# assumes that seq2 has already been blasted against seq1, which is the source
# of start1, end1 and start2, end2
# The format is:
## - seq1
## - start1
## - end1
## - length seq2
## - start2
## - end2
## - expected
extract_cases = [
        # Test case where both seq1 and seq2 are identical, and fully match
        ('AAAAA', 1, 5, 5, 1, 5, 'AAAAA'),
        # Test case where seq1 has a suffix that does not match seq2
        ('AAAAAGGG', 1, 5, 5, 1, 5, 'AAAAA'),
        # Test case where seq1 has a prefix that does not match seq2
        ('GGGAAAAA', 4, 8, 5, 1, 5, 'AAAAA'),
        # Test case where seq1 should be extended towards the beginning
        ('GAAAAA', 2, 5, 5, 2, 5, 'GAAAA'),
        # Test case where seq1 should be extended towards the beginning a lot
        ('GGGGGAAAAA', 6, 10, 10, 6, 10, 'GGGGGAAAAA'),
        # Test case where seq2 should be extended towards the end
        ('AAAAAGGG', 1, 5, 8, 1, 5, 'AAAAAGGG'),
        # Test case where content from seq2 is missing from the start of seq1
        ('AAAAA', 1, 5, 8, 4, 8, 'AAAAA'),
        # Test case where content from seq2 is missing from the end of seq1
        ('AAAAA', 1, 5, 8, 1, 5, 'AAAAA'),
        # Test case where seq2 should be extended towards the end, but not to
        # full lenght of seq1
        ('AAAAAGGG', 1, 5, 11, 1, 5, 'AAAAAGGG'),
        # Test case where seq2 should be extended towards the beginning, but not to
        # full lenght of seq1
        ('GGGGGAAAAA', 6, 10, 14, 10, 14, 'GGGGGAAAAA'),
]

@pytest.mark.parametrize('seq1,start1,end1,seq2_len,start2,end2,expected', extract_cases)
def test_extract_hit_region(seq1, start1, end1, seq2_len, start2, end2, expected):
    assert extract_hit_region(seq1, start1, end1, seq2_len, start2, end2) == expected


# These are the test cases for extending seq1 with sequences from seq2. This
# assumes that seq2 has already been blasted against seq1, which is the source
# of start1, end1 and start2, end2
# The format is:
## - seq1
## - start1
## - end1
## - seq2
## - start2
## - end2
## - expected
extend_cases = [
        # Test case where both seq1 and seq2 are identical, and fully match
        ('AAAAA', 1, 5, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where the hit is smaller than seq1. In this case,
        # we should only return the sequence content of the blast hit
        ('AAAAA', 2, 5, 'TAAAA', 2, 5, 'AAAA'),
        # Test case where the seq1 matches, but is larger, than the reference
        ('TCGAAAAATCG', 4, 8, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 matches, but is larger (at the beginning) than
        # the reference
        ('TCGAAAAA', 4, 8, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 matches, but is larger (at the end) than
        # the reference
        ('AAAAATCG', 1, 5, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where the end of seq1 should be extended with content from seq2
        ('AAAAA', 1, 5, 'AAAAATCG', 1, 5, 'AAAAATCG'),
        # Test case where the start of seq1 should be extended with content from seq2
        ('AAAAA', 1, 5, 'CTGAAAAA', 4, 8, 'CTGAAAAA'),
]

@pytest.mark.parametrize('seq1,start1,end1,seq2,start2,end2,expected', extend_cases)
def test_extend_hit_reference(seq1, start1, end1, seq2, start2, end2, expected):
    assert extend_hit_reference(seq1, start1, end1, seq2, start2, end2) == expected
