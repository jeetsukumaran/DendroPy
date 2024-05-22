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
Combinatoric and related functionCombinatoric and related functions.
"""

from dendropy.utility import deprecate

def factorial(num):
    """factorial(n): return the factorial of the integer num.
    factorial(0) = 1
    factorial(n) with n<0 is -factorial(abs(n))
    """
    deprecate.dendropy_deprecation_warning(
        preamble="Deprecated since Dendropy 5:",
        old_construct="combinatorics.factorial",
        new_construct="math.factorial"
    )
    result = 1
    for i in range(1, abs(num)+1):
        result *= i
    return result

def choose(population, sample):
    """
    Returns  ``population`` choose ``sample``, given
    by: n! / k!(n-k)!, where n == ``population`` and
    k == ``sample``.
    """
    if sample > population:
        return 0
    s = max(sample, population - sample)
    assert s <= population
    assert population > -1
    if s == population:
        return 1
    numerator = 1
    denominator = 1
    for i in range(s+1, population + 1):
        numerator *= i
        denominator *= (i - s)
    return numerator/denominator

def num_edges_on_tree(num_leaves, is_rooted=True):
    if is_rooted:
        return (2 * num_leaves) - 2
    else:
        return (2 * num_leaves) - 3

def num_internal_nodes_on_tree(num_leaves, is_rooted=True):
    if is_rooted:
        return num_leaves - 1
    else:
        return num_leaves - 2

def num_internal_edges_on_tree(num_leaves, is_rooted=True):
    if is_rooted:
        return num_leaves - 2
    else:
        return num_leaves - 3
