#!/usr/bin/env python3

import argparse
import re


def trim_single_description(description, ref_size):
    """ Trim a single description, if safe
    >>> identical = '='

    trim_single_description(identical, 20)
    '='
    """
    # If the description is that the sequence is identical to the reference, we
    # have nothing to trim
    if description == '=':
        return '='

    # If a single variant is a deletion that spans the end of the referene
    if re.match(f'\d+_{ref_size}', description):
        return '='

def trim_multiple_descriptions(description, ref_size):
    """ Trim the first and last description, if safe
    """
    # If there is only a single variant
    # Get the list of variants
    variants = description[1:-1].split(';')
    # If the first variant is a deletion starting at 1
    if re.match('1_\d+del', variants[0]):
        variants.pop(0)

    # If the second variant is a deletion till the end of the reference
    if re.match(f'\d+_{ref_size}del', variants[-1]):
        variants.pop(-1)

    # If there are now no variants left
    if not variants:
        return '='

def trim_description(description, ref_size):
    """ Trim insertions / deletions at the edges

    >>> identical = '='
    >>> contained = '[1_3del;13_20del]'
    >>> end_missing = '17_20del'

    >>> trim_description(identical, 20)
    '='

    >>> trim_description(contained, 20)
    '='

    >>> trim_description(end_missing, 20)
    '='
    """
    # If there are multiple descriptions
    if '[' in description:
        return trim_multiple_descriptions(description, ref_size)
    else:
        return trim_single_description(description, ref_size)


def main(args):
    with open(args.descriptions) as fin:
        for line in fin:
            name, description = line.strip('\n').split('\t')
            trimmed = trim_description(description)
            print(f'{name}\t{trimmed}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--descriptions', required=True,
                        help='TSV file of descriptions')

    args = parser.parse_args()

    main(args)
