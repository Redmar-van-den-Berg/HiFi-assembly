import pytest

from padding import add_padding

# These are the test cases for aligning seq1 to seq2 by using padding. This
# assumes that seq2 has already been blasted against seq1, which is the source
# of start1, end1 and start2, end2
cases = [ 
        # Test case where both seq1 and seq2 are identical, and fully match
        ('AAAAAAAAAA', 1, 10, 'AAAAAAAAAA', 1, 10, 'AAAAAAAAAA')
        ]

@pytest.mark.parametrize('seq1,start1,end1,seq2,start2,end2,expected', cases)
def test_add_padding(seq1, start1, end1, seq2, start2, end2, expected):
    assert add_padding(seq1, start1, end1, seq2, start2, end2) == expected

