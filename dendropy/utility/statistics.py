#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Functions to calculate some general statistics.
"""

def _calc_mean_and_variance_pop_n(values):
    n = 0
    s = 0.0
    ss = 0.0
    for v in values:
        n += 1
        s += v
        ss += v*v
    if n == 0:
        raise IndexError("values in calc_mean_and_variance cannot be empty")
    mean = float(s)/n
    var = (ss - mean*s)/n
    return mean, var, n

def calc_mean_and_population_variance(values):
    """Returns the mean and population variance while only passing over the
    elements in values once."""
    return _calc_mean_and_variance_pop_n(values)[:2]

def calc_mean_and_sample_variance(values):
    """Returns the mean and sample variance while only passing over the
    elements in values once."""
    mean, pop_var, n = _calc_mean_and_variance_pop_n(values)
    if n == 1:
        return mean, float('inf')
    samp_var = n*pop_var/(n-1)
    return mean, samp_var

