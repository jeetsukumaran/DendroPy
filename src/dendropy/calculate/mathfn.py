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
Some common mathematical functions.
"""

import functools
from dendropy.utility import deprecate
import math

def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    deprecate.dendropy_deprecation_warning(
        preamble="Deprecated since Dendropy 5:",
        old_construct="mathfn.gcd",
        new_construct="math.gcd"
    )
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Return lowest common multiple."""
    return a * b // math.gcd(a, b)

def LCM(*args):
    """Return lcm of args."""
    return functools.reduce(lcm, args)
