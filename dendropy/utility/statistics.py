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

import math

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

def calc_unimodal_emp_hpd(values, conf=0.05):
    """
    Assuming a **unimodal** distribution, returns the 0.95 highest posterior
    density (HPD) interval for a set of samples from a posterior distribution.
    Adapted from `emp.hpd` in the "TeachingDemos" R package (Copyright Greg
    Snow; licensed under the Artistic License).
    """
    conf = min([conf, 1.0 - conf])
    n = len(values)
    nn = int(round(n * conf))
    x = sorted(values)
    xx = []
    for i in range(nn):
        Z1 = x[n-(nn-i)]
        Z2 = x[i]
        #print "==>", Z1, Z2
        xx.append(Z1 - Z2)
    m = min(xx)
    nnn = xx.index(m)
    #print "n=", n
    #print "nn=", nn
    #print "xx=", xx
    #print "m=", m
    #print "nnn=", n
    return (x[nnn], x[n-nn+nnn])

def median(pool):
    """
    Returns median of sample. From: http://wiki.python.org/moin/SimplePrograms
    """
    copy = sorted(pool)
    size = len(copy)
    if size % 2 == 1:
        return copy[(size - 1) / 2]
    else:
        return (copy[size/2 - 1] + copy[size/2]) / 2

def summarize_sample(values):
    """
    Summarizes a sample of values, returing a dictionary with:

        - `range`       : tuple pair representing minimum and maximum values
        - `mean`        : mean of sample
        - `median`      : median of sample
        - `quart_5_95`  : tuple pair representing 5% and 95% quartile
        - `hpd95`       : tuple pair representing 5% and 95% HPD

    """
    pass
