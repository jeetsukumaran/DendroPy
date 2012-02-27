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

def _pooled_covariance(u, v):
    """
    Returns pooled covariance matrix of u, v.
    """
    assert len(u[0]) == len(v[0]), "Number of columns in matrices not equal"
    nrow1 = len(u)
    nrow2 = len(v)
    total_rows = nrow1 + nrow2
    f1 = float(nrow1) / total_rows
    f2 = float(nrow2) / total_rows
    s1 = u.covariance_by_cols()
    s2 = v.covariance_by_cols()
    pooled_cov = []
    for i, r1 in enumerate(s1):
        pooled_cov.append([])
        for j, c1 in enumerate(r1):
            pooled_cov[-1].append(f1 * s1[i][j] + f2 *s2[i][j])
    pooled_cov = linearalg.new_matrix(pooled_cov)
    return pooled_cov

def squared_mahalanobis(u, v):
    """
    Returns the *squared* Mahalanobis distance between matrices `u` and `v`.
    `u` and `v` must be 2-dimensional, and have the same number of columns.
    That is, they must be a list of lists, and the lengths of the inner lists
    must be equal. If there simple vectors of values, they should be
    represented as lists of multiple single-element lists, or use the
    convenience function `squared_mahalanobis_vec(u, v)`.  If only relative
    distances are needed (as they are in most cases), then this function should
    be preferred over `mahalanobis(u, v)` which has the added computational
    expense of taking the square root.

        >>> u = [ [2, 3.14, 1.3], [1, 4, 5] ]
        >>> v = [ [4, 1, 1], [5, 3, 2], [1, 3, 4], [1, 4, 4] ]
        >>> print squared_mahalanobis(u, v)
        1.46632046397
        >>> u = [ [2, 5], [5, 7] ]
        >>> v = [ [1, 5], [4, 8] ]
        >>> print squared_mahalanobis(u, v)
        5.37483849887
        >>> # for single column, you can:
        >>> u = [ [1], [3], [5] ]
        >>> v = [ [2], [3] ]
        >>> print squared_mahalanobis(u, v)
        0.65192024052
        >>> # or:
        >>> u = [ 1, 3, 5 ]
        >>> v = [ 2, 3 ]
        >>> print squared_mahalanobis(u, v)
        0.65192024052
    """
    if not isinstance(u, linearalg.Matrix):
        u = linearalg.new_matrix(u)
    if not isinstance(v, linearalg.Matrix):
        v = linearalg.new_matrix(v)
    assert len(u[0]) == len(v[0])
    pooled_cov = _pooled_covariance(u, v)
    if len(pooled_cov) == 1:
        pooled_cov_inv = pooled_cov
    else:
        pooled_cov_inv = pooled_cov.inverse()

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
    d = mean_diffs.mmul(pooled_cov_inv).mmul(mean_diffs.tr())
    return d[0][0]

def mahalanobis(u, v):
    """
    Returns Mahalanobis distance between matrices `u` and `v`.
    `u` and `v` must be 2-dimensional, and have the same number of columns.
    That is, they must be a list of lists, and the lengths of the inner lists
    must be equal. If there simple vectors of values, they should be
    represented as lists of multiple single-element lists, or use the
    convenience function `squared_mahalanobis_vec(u, v)`.  If only relative
    distances are needed (as they are in most cases), then
    `squared_mahalanobis(u, v) should be preferred over this one as this has
    the added computational expense of taking the square root.
    """
    return math.sqrt(squared_mahalanobis(u,v))


