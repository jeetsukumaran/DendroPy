#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
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
