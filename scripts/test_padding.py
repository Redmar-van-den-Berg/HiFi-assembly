import pytest

from padding import add_padding, extract_hit_region

# These are the test cases for aligning seq1 to seq2 by using padding. This
# assumes that seq2 has already been blasted against seq1, which is the source
# of start1, end1 and start2, end2
padding_cases = [
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
        ('GGAAAAA', 3, 7, 'AAAAA', 1, 5, 'AAAAA'),
        # Test case where seq1 is longer, but only the center of the sequence
        # matches
        ('GGGAAAGGG', 4, 6, 'AAAAA', 2, 4, 'GAAAG'),
        # Test case where seq1 is shorter, and must be padded at the beginning
        ('AAA', 1, 3, 'GGAAA', 3, 5, '--AAA'),
        # Test case where seq1 is shorter, and must be padded at the end
        #('AAA', 1, 3, 'AAAGG', 1, 3, 'AAA--'),
]

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

@pytest.mark.parametrize('seq1,start1,end1,seq2,start2,end2,expected', padding_cases)
def test_add_padding(seq1, start1, end1, seq2, start2, end2, expected):
    assert add_padding(seq1, start1, end1, seq2, start2, end2) == expected

@pytest.mark.parametrize('seq1,start1,end1,seq2,start2,end2,expected', extract_cases)
def test_extract_hit_region(seq1, start1, end1, seq2, start2, end2, expected):
    assert extract_hit_region(seq1, start1, end1, seq2, start2, end2) == expected
