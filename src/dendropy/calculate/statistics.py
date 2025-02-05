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
Functions to calculate some general statistics.
"""

from dendropy.calculate import probability
from operator import itemgetter

def _mean_and_variance_pop_n(values):
    n = 0
    s = 0.0
    ss = 0.0
    for v in values:
        n += 1
        s += v
        ss += v*v
    if n == 0:
        raise IndexError("values in mean_and_variance cannot be empty")
    mean = float(s)/n
    var = (ss - mean*s)/n
    return mean, var, n

def mean_and_population_variance(values):
    """Returns the mean and population variance while only passing over the
    elements in values once."""
    return _mean_and_variance_pop_n(values)[:2]

def mean_and_sample_variance(values):
    """Returns the mean and sample variance while only passing over the
    elements in values once."""
    mean, pop_var, n = _mean_and_variance_pop_n(values)
    if n == 1:
        return mean, float('inf')
    samp_var = n*pop_var/(n-1)
    return mean, samp_var

def mode(values, bin_size=0.1):
    """
    Returns the mode of a set of values.
    """
    bins = {}
    for v in values:
        if bin_size is not None:
            idx = int(round(float(v)/bin_size))
        else:
            idx = v
        if idx in bins:
            bins[idx] += 1
        else:
            bins[idx] = 1
    sorted_bins = sorted(bins.items(), key=itemgetter(1), reverse=True)
    max_count = sorted_bins[0][1]
    results = [(sorted_bins[i][0] * bin_size) for i in range(len(sorted_bins)) if sorted_bins[i][1] >= max_count]
    return results

def median(pool):
    """
    Returns median of sample. From: http://wiki.python.org/moin/SimplePrograms
    """
    copy = sorted(pool)
    size = len(copy)
    if size % 2 == 1:
        idx = int((size - 1) / 2)
        return copy[idx]
    else:
        idx1 = int(size/2) - 1
        idx2 = int(size/2)
        return (copy[idx1] + copy[idx2]) / 2

def empirical_hpd(values, conf=0.05):
    """
    Assuming a **unimodal** distribution, returns the 0.95 highest posterior
    density (HPD) interval for a set of samples from a posterior distribution.
    Adapted from ``emp.hpd`` in the "TeachingDemos" R package (Copyright Greg
    Snow; licensed under the Artistic License).
    """
    conf = min([conf, 1.0 - conf])
    n = len(values)
    nn = int(round(n * conf))
    x = sorted(values)
    xx = []
    if nn == 0:
        raise ValueError("Sample size too small: %s" % len(values))
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

def empirical_cdf(values, v):
    """
    Returns the proportion of values in ``values`` <= ``v``.
    """
    count = 0.0
    for idx, v0 in enumerate(values):
        if v0 < v:
            count += 1
    return count / len(values)

def quantile(values, q):
    """
    Returns q-th quantile.
    """
    values = sorted(values)
    size = len(values)
    idx = int(round(size * q)) - 1
    if idx == 0:
        raise ValueError("Sample size too small: %s" % len(values))
    return values[idx]



# http://adorio-research.org/wordpress/?p=125
# File    quantile.py
# Desc    computes sample quantiles
# Author  Ernesto P. Adorio, PhD.
#         UPDEPP (U.P. at Clarkfield)
# Version 0.0.1 August 7. 2009
def quantile(x, q,  qtype = 7, issorted = False):
    from math import modf, floor
    """
    Args:
       x - input data
       q - quantile
       qtype - algorithm
       issorted- True if x already sorted.

    Compute quantiles from input array x given q.For median,
    specify q=0.5.

    References:
       http://reference.wolfram.com/mathematica/ref/Quantile.html
       http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile

    Author:
	Ernesto P.Adorio Ph.D.
	UP Extension Program in Pampanga, Clark Field.
    """
    if not issorted:
        y = sorted(x)
    else:
        y = x
    if not (1 <= qtype <= 9):
       return None  # error!

    # Parameters for the Hyndman and Fan algorithm
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3

            (0,   0, 0, 1), # California linear interpolation, R type 4
            (0.5, 0, 0, 1), # hydrologists method, R type 5
            (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
            (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
            (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
           ]

    a, b, c, d = abcd[qtype-1]
    n = len(x)
    g, j = modf( a + (n+b) * q -1)
    if j < 0:
        return y[0]
    elif j >= n:
        return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!

    j = int(floor(j))
    if g ==  0:
       return y[j]
    else:
       return y[j] + (y[j+1]- y[j])* (c + d * g)

def quantile_5_95(values):
    """
    Returns 5% and 95% quantiles.
    """
    values = sorted(values)
    size = len(values)
    idx5 = int(round(size * 0.05)) - 1
    idx95 = int(round(size * 0.95)) - 1
    if idx5 == 0:
        raise ValueError("Sample size too small: %s" % len(values))
    return values[idx5], values[idx95]

def variance_covariance(data, population_variance=False):
    """
    Returns the Variance-Covariance matrix for ``data``.
    From: http://www.python-forum.org/pythonforum/viewtopic.php?f=3&t=17441
    """
    N = len(data) # number of vectors
    D = len(data[0]) # dimensions per vector
    if population_variance:
        denom = N
    else:
        denom = N-1.0

    means = [0.0 for i in range(D)] # intialize 1xD mean vector
    for i in range(N):
        for j in range(D):
            means[j] += data[i][j]
    means = [i/N for i in means]
    # print "Means:"," ".join(map(str,means)),"\n"

    covar = [[0.0 for i in range(D)] for j in range(D)] # initialize DxD covariance matrix

    for i in range(D):
        for j in range(i+1): #  covariance symmetric, only do lower triangle of matrix
            sum = 0.0
            for k in range(N):
                sum += data[k][i]*data[k][j]
            covar[i][j] = sum/denom - means[i]*means[j]*N/denom

    for j in range(D):
            for k in range(j+1):
                covar[k][j] = covar[j][k]

    # print "covariance:"
    # for i in range(D):
    #     print " ".join(map(str,covar[i]))
    # print ""
    return covar

def rank(value_to_be_ranked, value_providing_rank):
    """
    Returns the rank of ``value_to_be_ranked`` in set of values, ``values``.
    Works even if ``values`` is a non-orderable collection (e.g., a set).
    A binary search would be an optimized way of doing this if we can constrain
    ``values`` to be an ordered collection.
    """
    num_lesser = [v for v in value_providing_rank if v < value_to_be_ranked]
    return len(num_lesser)

class FishersExactTest(object):
    """
    Given a 2x2 table:

        +---+---+
        | a | b |
        +---+---+
        | c | d |
        +---+---+

    represented by a list of lists::

        [[a,b],[c,d]]

    this calculates the sum of the probability of this table and all others
    more extreme under the null hypothesis that there is no association between
    the categories represented by the vertical and horizontal axes.
    """

    @staticmethod
    def probability_of_table(table):
        """
        Given a 2x2 table:

            +---+---+
            | a | b |
            +---+---+
            | c | d |
            +---+---+

        represented by a list of lists::

            [[a,b],[c,d]]

        this returns the probability of this table under the null hypothesis of
        no association between rows and columns, which was shown by Fisher to be
        a hypergeometric distribution:

            p = ( choose(a+b, a) * choose(c+d, c) ) / choose(a+b+c+d, a+c)

        """
        a = table[0][0]
        b = table[0][1]
        c = table[1][0]
        d = table[1][1]
        return probability.hypergeometric_pmf(a, a+b, c+d, a+c)

    def __init__(self, table):
        self.table = table
        self.flat_table = [table[0][0], table[0][1], table[1][0], table[1][1]]
        self.min_value = min(self.flat_table)
        self.max_value = max(self.flat_table)

    def _rotate_cw(self, table):
        """
        Returns a copy of table such that all the values
        are rotated clockwise once.
        """
        return [ [ table[1][0], table[0][0] ],
                [table[1][1], table[0][1] ] ]

    def _min_rotation(self):
        """
        Returns copy of self.table such that the smallest value is in the first
        (upper left) cell.
        """
        table = [list(self.table[0]), list(self.table[1])]
        while table[0][0] != self.min_value:
            table = self._rotate_cw(table)
        return table

    def _max_rotation(self):
        """
        Returns copy of self.table such that the largest value is in the first
        (upper left) cell.
        """
        table = [list(self.table[0]), list(self.table[1])]
        while table[0][0] != self.max_value:
            table = self._rotate_cw(table)
        return table

    def _sum_left_tail(self):
        """
        Returns the sum of probabilities of tables that are *more* extreme than
        the current table.
        """
        # left_tail_tables = self._get_left_tail_tables()
        # p_vals = [ self.probability_of_table(t) for t in left_tail_tables ]
        p_vals = self._get_left_tail_probs()
        return sum(p_vals)

    def _sum_right_tail(self):
        """
        Returns the sum of probabilities of tables that are *less* extreme than
        the current table.
        """
        # right_tail_tables = self._get_right_tail_tables()
        # p_vals = [ self.probability_of_table(t) for t in right_tail_tables ]
        p_vals = self._get_right_tail_probs()
        return sum(p_vals)

    def _get_left_tail_probs(self):
        """
        Returns list of probabilities of all tables *more* extreme than the
        current table.
        """
        table = self._min_rotation()
        row_totals = [sum(table[0]), sum(table[1])]
        col_totals = [table[0][0] + table[1][0], table[0][1] + table[1][1]]
        p_vals = []
        while True:
            table[0][0] -= 1
            if table[0][0] < 0:
                break
            table[0][1] = row_totals[0] - table[0][0]
            table[1][0] = col_totals[0] - table[0][0]
            table[1][1] = row_totals[1] - table[1][0]
            p_vals.append(self.probability_of_table(table))
        return p_vals

    def _get_right_tail_probs(self):
        """
        Returns list of probabilities of all tables *less* extreme than the
        current table.
        """
        table = self._min_rotation()
        row_totals = [sum(table[0]), sum(table[1])]
        col_totals = [table[0][0] + table[1][0], table[0][1] + table[1][1]]
        p_vals = []
        while True:
            table[0][0] += 1
            table[0][1] = row_totals[0] - table[0][0]
            if table[0][1] < 0:
                break
            table[1][0] = col_totals[0] - table[0][0]
            if table[1][0] < 0:
                break
            table[1][1] = row_totals[1] - table[1][0]
            if table[1][1] < 0:
                break
            p_vals.append(self.probability_of_table(table))
        return p_vals

    def _get_left_tail_tables(self):
        """
        Returns all tables that are *more* extreme than the current table.
        """
        table = self._min_rotation()
        row_totals = [sum(table[0]), sum(table[1])]
        col_totals = [table[0][0] + table[1][0], table[0][1] + table[1][1]]
        left_tail_tables = []
        while True:
            table[0][0] -= 1
            if table[0][0] < 0:
                break
            table[0][1] = row_totals[0] - table[0][0]
            table[1][0] = col_totals[0] - table[0][0]
            table[1][1] = row_totals[1] - table[1][0]
            left_tail_tables.append([list(table[0]), list(table[1])])
        return left_tail_tables

    def _get_right_tail_tables(self):
        """
        Returns all tables that are *less* extreme than the current table.
        """
        table = self._min_rotation()
        row_totals = [sum(table[0]), sum(table[1])]
        col_totals = [table[0][0] + table[1][0], table[0][1] + table[1][1]]
        right_tail_tables = []
        while True:
            table[0][0] += 1
            table[0][1] = row_totals[0] - table[0][0]
            if table[0][1] < 0:
                break
            table[1][0] = col_totals[0] - table[0][0]
            if table[1][0] < 0:
                break
            table[1][1] = row_totals[1] - table[1][0]
            if table[1][1] < 0:
                break
            right_tail_tables.append([list(table[0]), list(table[1])])
        return right_tail_tables

    def left_tail_p(self):
        """
        Returns the sum of probabilities of this table and all others more
        extreme.
        """
        return self.probability_of_table(self.table) + self._sum_left_tail()

    def right_tail_p(self):
        """
        Returns the sum of probabilities of this table and all others more
        extreme.
        """
        return self.probability_of_table(self.table) + self._sum_right_tail()

    def two_tail_p(self):
        """
        Returns the sum of probabilities of this table and all others more
        extreme.
        """
        p0 = self.probability_of_table(self.table)
        all_p_vals = self._get_left_tail_probs() + self._get_right_tail_probs()
        p_vals = []
        for p in all_p_vals:
            if p <= p0:
                p_vals.append(p)
        return sum(p_vals) + p0

def summarize(values):
    """
    Summarizes a sample of values:

        - ``range``       : tuple pair representing minimum and maximum values
        - ``mean``        : mean of sample
        - ``median``      : median of sample
        - ``var``         : (sample) variance
        - ``sd``          : (sample) standard deviation
        - ``hpd95``       : tuple pair representing 5% and 95% HPD
        - ``quant_5_95``  : tuple pair representing 5% and 95% quantile

    """
    summary = {}
    if len(values) == 0:
        raise ValueError("No values in data")
    try:
        summary['range'] = (min(values), max(values))
    except (ValueError, OverflowError):
        summary['range'] = None
    try:
        summary['mean'], summary['var'] = mean_and_sample_variance(values)
        try:
            #summary['sd'] = math.sqrt(summary['var'])
            summary['sd'] = summary['var'] ** 0.5
        except ValueError:
            summary['sd'] = 0.0
        except OverflowError:
            summary['sd'] = None
    except (ValueError, OverflowError, IndexError):
        summary['mean'], summary['var'], summary['sd'] = None, None, None
    try:
        summary['median'] = median(values)
    except (ValueError, OverflowError):
        summary['median'] = None
    try:
        summary['hpd95'] = empirical_hpd(values, conf=0.95)
    except (ValueError, OverflowError):
        summary['hpd95'] = None
    try:
        summary['quant_5_95'] = quantile_5_95(values)
    except (ValueError, OverflowError):
        summary['quant_5_95'] = None
    return summary
