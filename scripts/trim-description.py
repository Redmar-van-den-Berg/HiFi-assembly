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

    # If a single variant is a deletion that spans the end of the reference
    if re.match(f'\d+_{ref_size}', description):
        return '='

    # If a single variant is a deletion at the beginning of the reference
    if re.match('1_\d+del', description):
        return '='

def trim_first_variant(variant):
    """ Determine wether the first variant can safely be trimmed
    """
    # If the first variant is a deletion starting at 1
    if re.match('1_\d+del', variant):
        return True
    # If the first variant is an extention of the reference starting at 0
    elif re.match('0_1ins\w+', variant):
        return True

    return False


def trim_last_variant(variant, ref_size):
    """ Determine wether the last variant can safely be trimmed
    """
    # If the last variant is a deletion till the end of the reference
    if re.match(f'\d+_{ref_size}del', variant):
        return True
    elif re.match(f'{ref_size}_{ref_size+1}ins\w+', variant):
        return True
    return False


def trim_multiple_descriptions(description, ref_size):
    """ Trim the first and last description, if safe
    """
    # Get the list of variants
    variants = description[1:-1].split(';')

    # Get the first and last variants
    first_variant = variants[0]
    last_variant = variants[-1]

    if trim_first_variant(first_variant):
        variants.pop(0)
    if trim_last_variant(last_variant, ref_size):
        variants.pop(-1)

    # If there are now no variants left
    if not variants:
        return '='
    # If there is only a single variant left
    if len(variants) == 1:
        return variants[0]


def trim_description(description, ref_size):
    """ Trim insertions / deletions at the edges

    >>> identical = '='
    >>> contained = '[1_3del;13_20del]'
    >>> end_missing = '17_20del'
    >>> begin_missing = '1_3del'
    >>> ref_contained = '[0_1insWWW;20_21insWWW]'
    >>> has_prefix = '[0_1insWWW;13_20del]'

    >>> trim_description(identical, 20)
    '='

    >>> trim_description(contained, 20)
    '='

    >>> trim_description(end_missing, 20)
    '='

    >>> trim_description(begin_missing, 20)
    '='

    >>> trim_description(ref_contained, 20)
    '='

    >>> trim_description(has_prefix, 20)
    '='

    >>> trim_description(has_prefix, 21)
    '13_20del'
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
