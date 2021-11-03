import pytest

from padding import add_padding

# These are the test cases for aligning seq1 to seq2 by using padding. This
# assumes that seq2 has already been blasted against seq1, which is the source
# of start1, end1 and start2, end2
cases = [ 
        # Test case where both seq1 and seq2 are identical, and fully match
        ('AAAAA', 1, 5, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where both seq1 and seq2 are the same length, but seq1
        # overhangs the beginning
        ('GGAAA', 3, 5, 'AAAAA', 3, 5, 'GGAAA'),
        # Test case where both seq1 and seq2 are the same length, but seq1
        # overhangs the end
        ('AAAGG', 1, 3, 'AAAAA', 1, 3, 'AAAGG'),
        # Test case where seq1 is longer, and the beginning of the sequence
        # matches
        ('AAAAAGG', 1, 5, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 is longer, and the end of the sequence matches
        ('GGAAAAA', 3, 7, 'AAAAA', 1, 5, 'AAAAA')
        ]

@pytest.mark.parametrize('seq1,start1,end1,seq2,start2,end2,expected', cases)
def test_add_padding(seq1, start1, end1, seq2, start2, end2, expected):
    assert add_padding(seq1, start1, end1, seq2, start2, end2) == expected

