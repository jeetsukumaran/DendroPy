#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Various bitwise utilities.
"""

def bit_length(n):
    """
    Return the number of bits necessary to represent an integer in binary,
    excluding the sign and leading zeros. Also can be used to determine the
    index of the highest set bit, or the width of the bitstring
    representing the integer.
    """
    try:
        s = bin(n)          # binary representation:  bin(-37) --> '-0b100101'
        s = s.lstrip('-0b') # remove leading zeros and minus sign
        return len(s)       # len('100101') --> 6
    except TypeError:  # if n is None
        return 0

def int_as_bitstring(n, length=None, symbol0=None, symbol1=None, reverse=False):
    if length is None:
        length = bit_length(n)
    s = bin(n)[2:].rjust(length, "0")
    if symbol0 is not None:
        s = s.replace("0", symbol0)
    if symbol1 is not None:
        s = s.replace("1", symbol1)
    if reverse:
        return s[::-1]
    else:
        return s

def num_set_bits(n):
    return bin(n).count("1")

def least_significant_set_bit(n):
    """
    Returns least-significant bit in integer 'n' that is set.
    """
    m = n & (n - 1)
    return m ^ n

def indexes_of_set_bits(s, fill_bitmask=-1, one_based=False, ordination_in_mask=False):
    return [i for i in set_bit_index_iter(s, fill_bitmask, one_based, ordination_in_mask)]

def set_bit_index_iter(s, fill_bitmask=-1, one_based=False, ordination_in_mask=False):
    """
    Returns the index of each bit that is on in ``s`` and the ``fill_bitmask``

        If 'one_based` is True then the 0x01 bit is returned as 1 instead of 0.
        If ``ordination_in_mask`` is True then the indices returned will be the
            count of the 1's in the fill_bitmask that are to the right of the bit rather
            than the total number of digits to the right of the bit. Thus, the
            index will be the index in a taxon block that is the subset of the
            full set of taxa).
    """
    currBitIndex = one_based and 1 or 0
    test_bit = 1
    maskedSplitRep = s & fill_bitmask
    standard_ordination = not ordination_in_mask
    while test_bit <= maskedSplitRep:
        if maskedSplitRep & test_bit:
            yield currBitIndex
        if standard_ordination or (fill_bitmask & test_bit):
            currBitIndex += 1
        test_bit <<= 1
