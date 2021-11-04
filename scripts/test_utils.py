import pytest

from utils import extract_hit_region

# These are the test cases for extracting regions of seq2 from seq1. This
# assumes that seq2 has already been blasted against seq1, which is the source
# of start1, end1 and start2, end2
extract_cases = [
        # Test case where both seq1 and seq2 are identical, and fully match
        ('AAAAA', 1, 5, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 has a suffix that does not match seq2
        ('AAAAAGGG', 1, 5, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 has a prefix that does not match seq2
        ('GGGAAAAA', 4, 8, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 should be extended towards the beginning
        ('GAAAAA', 2, 5, 'AAAAA', 2, 5, 'GAAAA'),
        # Test case where seq1 should be extended towards the beginning a lot
        ('GGGGGAAAAA', 6, 10, 'AAAAAAAAAA', 6, 10, 'GGGGGAAAAA'),
        # Test case where seq2 should be extended towards the end
        ('AAAAAGGG', 1, 5, 'AAAAATTT', 1, 5, 'AAAAAGGG'),
        # Test case where content from seq2 is missing from the start of seq1
        ('AAAAA', 1, 5, 'GGGAAAAA', 4, 8, 'AAAAA'),
        # Test case where content from seq2 is missing from the end of seq1
        ('AAAAA', 1, 5, 'AAAAAGGG', 1, 5, 'AAAAA'),
        # Test case where seq2 should be extended towards the end, but not to
        # full lenght of seq1
        ('AAAAAGGG', 1, 5, 'AAAAATTTCCC', 1, 5, 'AAAAAGGG'),
        # Test case where seq2 should be extended towards the beginning, but not to
        # full lenght of seq1
        ('GGGGGAAAAA', 6, 10, 'CCCCAAAAAAAAAA', 10, 14, 'GGGGGAAAAA'),
]

@pytest.mark.parametrize('seq1,start1,end1,seq2,start2,end2,expected', extract_cases)
def test_extract_hit_region(seq1, start1, end1, seq2, start2, end2, expected):
    assert extract_hit_region(seq1, start1, end1, seq2, start2, end2) == expected
