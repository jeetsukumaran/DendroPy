#! /usr/bin/env python

############################################################################
##  util_for_testing.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################
"Functions that are useful for writing tests."
import itertools
import unittest
import inspect

def approx_equal(x, y, tol=1e-5):
    "Returns True if x and y differ by less than tol"
    return (abs(x - y) < tol)

def vec_approx_equal(x, y, tol=1e-5):
    """Returns True if each element in the iterable `x` differs by less than
    `tol` from the corresponding element in `y`
    """
    if len(x) != len(y):
        return False
    for i, j in itertools.izip(x, y):
        if abs(i - j) >= tol:
            return False
    return True

def mat_approx_equal(x, y, tol=1e-5):
    """Returns True if each cell in 2D iterable matrix `x` differs by less than
    `tol` from the corresponding element in `y`
    """
    if len(x) != len(y):
        return False
    for row_x, row_y in itertools.izip(x, y):
        if len(row_x) != len(row_y):
            return False
        for i, j in itertools.izip(row_x, row_y):
            if abs(i - j) >= tol:
                return False
    return True
    

def _failure(tester, msg):
    calling_frame = inspect.currentframe().f_back.f_back
    co = calling_frame.f_code
    emsg = "%s\nCalled from file %s, line %d, in %s" % (msg, co.co_filename, calling_frame.f_lineno, co.co_name)
    if isinstance(tester, unittest.TestCase):
        tester.assertTrue(False, emsg)
    else:
        raise AssertionError(emsg)

def assert_approx_equal(x, y, tester=None, tol=1e-5):
    """Asserts that x and y are approximately equal.  

    If `tester` is a unittest.TestCase then assertTrue is used; otherwise 
    AssertionErrors are raised.
    """
    if abs(x - y) >= tol:
        _failure(tester, "%f != %f" % (x, y))

def assert_vec_approx_equal(x, y, tester=None, tol=1e-5):
    """Returns True if each element in the iterable `x` differs by less than
    `tol` from the corresponding element in `y`

    If `tester` is a unittest.TestCase then assertTrue is used; otherwise 
    AssertionErrors are raised.
    """
    if len(x) != len(y):
        _failure(tester, "vectors of numbers differ in length (%d vs %d)" % (len(x), len(y)))
    for n, itup in enumerate(itertools.izip(x, y)):
        i, j = itup
        if abs(i - j) >= tol:
            _failure(tester, "%f != %f at element %d" % (i, j, n))
    

def assert_mat_approx_equal(x, y, tester=None, tol=1e-5):
    """Returns True if each cell in 2D iterable matrix `x` differs by less than
    `tol` from the corresponding element in `y`
    If `tester` is a unittest.TestCase then assertTrue is used; otherwise 
    AssertionErrors are raised.
    """
    if len(x) != len(y):
        _failure(tester, "Matrices differs in length (%d vs %d)" % (len(x), len(y)))
    for n, row_tup in enumerate(itertools.izip(x, y)):
        row_x, row_y = row_tup
        if len(row_x) != len(row_y):
            _failure(tester, "row %d of matrix differs in length (%d vs %d)" % (n, len(row_x), len(row_y)))
        for col, cell_tup in enumerate(itertools.izip(row_x, row_y)):
            i, j = cell_tup
            if abs(i - j) >= tol:
                _failure(tester, "%f != %f for column %d of row %d" % (i, j, col, n))

