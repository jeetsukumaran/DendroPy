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
Calculate distances between points, vectors, and matrices.
"""

import math
from dendropy.mathlib import linearalg

def squared_mahalanobis(u, v, cov=None, population_variance=False):
    """
    Returns the *squared* Mahalanobis distance between matrices `u` and `v`.
    `u` and `v` must be 2-dimensional, and have the same number of columns
    (though they can have different number of rows).
    That is, they must be a list of lists, and the lengths of the inner lists
    must be equal, though the length of the outer lists can be differnt). If
    there simple vectors of values, they should be represented as lists of
    multiple single-element lists, or use the convenience function
    `squared_mahalanobis_1d(u, v)`.
    `cov` is the covariance matrix. If not given, the pooled covariance of `u`
    and `v` are used.
    If the pooled covariance is calculated (i.e., `cov` is not given), then
    if `population_variance` is False (default), the data are treated as samples
    rather than a population.
    If only relative distances are needed (as they are in most cases), then
    this function should be preferred over `mahalanobis(u, v)` which has the
    added computational expense of taking the square root.

    The following examples calculate the Mahalanobis distances between
    matrices, using the pooled covariances::

        >>> u = [ [2, 3.14, 1.3], [1, 4, 5] ]
        >>> v = [ [4, 1, 1], [5, 3, 2], [1, 3, 4], [1, 4, 4] ]
        >>> print squared_mahalanobis(u, v)
        2.15009570304
        >>> u = [ [2, 5], [5, 7] ]
        >>> v = [ [1, 5], [4, 8] ]
        >>> print squared_mahalanobis(u, v)
        28.8888888889
        >>> # for single column, you can use
        ... `squared_mahalanobis_vec(u, v)` or:
        >>> u = [ [1], [3], [5] ]
        >>> v = [ [2], [3] ]
        >>> print squared_mahalanobis(u, v)
        0.425

    The followed extended example shows a different approach. Here, we have
    a known vector of means, and we are interested in calculating the distances
    of different datasets to these means. Instead of using the pooled covariances
    of the two matrices, only the covariance of the datasets are used (the means
    are taken to be the truth)::

        #! /usr/bin/env python

        import random
        from dendropy.mathlib import linearalg
        from dendropy.mathlib import distance

        nrows = 10
        ncols = 4
        v1 = []
        for i in range(nrows):
            v1.append([random.gauss(0, 10) for j in range(ncols)])
        v2 = []
        for i in range(nrows):
            v2.append([random.gauss(10, 10) for j in range(ncols)])
        v3 = []
        for i in range(nrows):
            v3.append([random.gauss(-10, 10) for j in range(ncols)])

        c1 = [ [0] * 4 ]
        c2 = [ [10] * 4 ]
        c3 = [ [-10] * 4 ]
        v1 = linearalg.new_matrix(v1)
        s1 = v1.covariance_by_cols()
        v2 = linearalg.new_matrix(v2)
        s2 = v2.covariance_by_cols()
        v3 = linearalg.new_matrix(v3)
        s3 = v3.covariance_by_cols()

        print
        print "-- v1 --"
        print "d(c1, v1) = {}".format(distance.squared_mahalanobis(c1, v1, cov=s1))
        print "d(c2, v1) = {}".format(distance.squared_mahalanobis(c2, v1, cov=s1))
        print "d(c3, v1) = {}".format(distance.squared_mahalanobis(c3, v1, cov=s1))

        print
        print "-- v2 --"
        print "d(c1, v2) = {}".format(distance.squared_mahalanobis(c1, v2, cov=s2))
        print "d(c2, v2) = {}".format(distance.squared_mahalanobis(c2, v2, cov=s2))
        print "d(c3, v2) = {}".format(distance.squared_mahalanobis(c3, v2, cov=s2))

        print
        print "-- v3 --"
        print "d(c1, v3) = {}".format(distance.squared_mahalanobis(c1, v3, cov=s3))
        print "d(c2, v3) = {}".format(distance.squared_mahalanobis(c2, v3, cov=s3))
        print "d(c3, v3) = {}".format(distance.squared_mahalanobis(c3, v3, cov=s3))

    This results in::

        -- v1 --
        d(c1, v1) = 0.170092378552
        d(c2, v1) = 8.59447583779
        d(c3, v1) = 8.898973355

        -- v2 --
        d(c1, v2) = 30.7463150693
        d(c2, v2) = 0.434906936737
        d(c3, v2) = 121.212523478

        -- v3 --
        d(c1, v3) = 4.19154157624
        d(c2, v3) = 18.4689462227
        d(c3, v3) = 0.128840873046

    The `mahal` function of MATLAB calculates the Mahalanobis distance as well.
    Its implementation and usage are a little different:

        d = mahal(Y,X) computes the Mahalanobis distance (in squared units) of
        each observation in Y from the reference sample in matrix X. If Y is
        n-by-m, where n is the number of observations and m is the dimension of
        the data, d is n-by-1. X and Y must have the same number of columns,
        but can have different numbers of rows. X must have more rows than
        columns.

        For observation I, the Mahalanobis distance is defined by d(I) =
        (Y(I,:)-mu)*inv(SIGMA)*(Y(I,:)-mu)', where mu and SIGMA are the sample
        mean and covariance of the data in X.

    Thus, the following MATLAB code::

        >> g1 = [ 2 2; 2 5; 6 5; 7 3; 4 7; 6 4; 5 3; 4 6; 2 5; 1 3; ];
        >> g2 = [ 6 5; 7 4; 8 7; 5 6; 5 4; ];
        >> mahal(g1, g2)

    can be replicated by::

        #! /usr/bin/env python

        from dendropy.mathlib import distance
        from dendropy.mathlib import linearalg

        g1 = [ [2, 2], [2, 5], [6, 5], [7, 3], [4, 7], [6, 4], [5, 3], [4, 6], [2, 5], [1, 3], ]
        g2 = [ [6, 5], [7, 4], [8, 7], [5, 6], [5, 4], ]
        s = linearalg.new_matrix(g2).covariance_by_cols(population_variance=False)
        for g in g1:
            print distance.squared_mahalanobis([g], g2, cov=s)

    """
    if not isinstance(u, linearalg.Matrix):
        u = linearalg.new_matrix(u)
    if not isinstance(v, linearalg.Matrix):
        v = linearalg.new_matrix(v)
    assert len(u[0]) == len(v[0])
    if cov is None:
        cov = linearalg.pooled_covariance(u, v, population_variance=population_variance)
    if len(cov) == 1:
        cov_inv = cov
    else:
        cov_inv = cov.inverse()

    # TODO: column means have already been calculated to
    # get the covariances: avoid duplicating calcs
    mean_diffs = linearalg.new_matrix([u.col_means() - v.col_means()])

    # `mean_diffs` is now a row vector
    # for consistency with the standard formulation, we should:
    #
    #     mean_diffs = mean_diffs.tr()
    #     d = mean_diffs.tr().mmul(pooled_cov_inv).mmul(mean_diffs)
    #
    # Instead ...
    d = mean_diffs.mmul(cov_inv).mmul(mean_diffs.tr())
    return d[0][0]

def squared_mahalanobis_1d(u, v, cov=None, population_variance=False):
    """
    Returns the *squared* Mahalanobis distance between vectors `u` and `v`.
    `u` and `v` are assumed to be single-dimension, i.e., a simple list or
    Vector of values.
    `cov` is the covariance matrix. If not given, the pooled covariance of `u`
    and `v` are used.
    If the pooled covariance is calculated (i.e., `cov` is not given), then
    if `population_variance` is False (default), the data are treated as samples
    rather than a population.
    If only relative distances are needed (as they are in most cases), then
    this function should be preferred over `mahalanobis_1d(u, v)` which has the
    added computational expense of taking the square root.
    For more details and examples, see `dendropy.mathlib.distance.squared_mahalanobis()`.

       >>> u = [1, 2, 2, 4, 1, 4]
       >>> v = [2, 1, 1, 0, 2, 1]
       >>> print squared_mahalanobis_1d(u, v)
       1.3800154321

    """
    u = linearalg.new_matrix([[i] for i in u])
    v = linearalg.new_matrix([[i] for i in v])
    return squared_mahalanobis(u, v, cov=cov, population_variance=population_variance)

def mahalanobis(u, v, cov=None, population_variance=False):
    """
    Returns the Mahalanobis distance between matrices `u` and `v`.
    `u` and `v` must be 2-dimensional, and have the same number of columns
    (though they can have different number of rows).
    That is, they must be a list of lists, and the lengths of the inner lists
    must be equal, though the length of the outer lists can be differnt). If
    there simple vectors of values, they should be represented as lists of
    multiple single-element lists, or use the convenience function
    `mahalanobis_1d(u, v)`.
    `cov` is the covariance matrix. If not given, the pooled covariance of `u`
    and `v` are used.
    If the pooled covariance is calculated (i.e., `cov` is not given), then
    if `population_variance` is False (default), the data are treated as samples
    rather than a population.
    If only relative distances are needed (as they are in most cases), then
    `squared_mahalanobis(u, v) should be preferred over this one as this has
    the added computational expense of taking the square root.
    For more details and examples, see `dendropy.mathlib.distance.squared_mahalanobis()`.

        >>> u = [ [2, 3.14, 1.3], [1, 4, 5] ]
        >>> v = [ [4, 1, 1], [5, 3, 2], [1, 3, 4], [1, 4, 4] ]
        >>> print mahalanobis(u, v)
        1.46632046397
        >>> u = [ [2, 5], [5, 7] ]
        >>> v = [ [1, 5], [4, 8] ]
        >>> print mahalanobis(u, v)
        5.37483849887
        >>> # for single column, you can use
        ... `mahalanobis_vec(u, v)` or:
        >>> u = [ [1], [3], [5] ]
        >>> v = [ [2], [3] ]
        >>> print mahalanobis(u, v)
        0.65192024052

    """
    return math.sqrt(squared_mahalanobis(u, v, cov=cov, population_variance=population_variance))


def mahalanobis_1d(u, v, cov=None, population_variance=False):
    """
    Returns the  Mahalanobis distance between vectors `u` and `v`.
    `u` and `v` are assumed to be single-dimension, i.e., a simple list or
    Vector of values.
    `cov` is the covariance matrix. If not given, the pooled covariance of `u`
    and `v` are used.
    If the pooled covariance is calculated (i.e., `cov` is not given), then
    if `population_variance` is False (default), the data are treated as samples
    rather than a population.
    If only relative distances are needed (as they are in most cases), then
    `squared_mahalanobis_1d(u, v) should be preferred over this one as this has
    the added computational expense of taking the square root.
    For more details and examples, see `dendropy.mathlib.distance.squared_mahalanobis()`.

       >>> u = [1, 2, 2, 4, 1, 4]
       >>> v = [2, 1, 1, 0, 2, 1]
       >>> print squared_mahalanobis_1d(u, v)
       1.3800154321

    """
    return math.sqrt(squared_mahalanobis_1d(u, v, cov=cov, population_variance=population_variance))
