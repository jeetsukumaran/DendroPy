#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
DEPRECATED IN DENDROPY 4: USE `dendropy.model.continuous` instead.
"""

from dendropy.utility import error
from dendropy.model import continuous

def simulate_continuous(node, rng, **kwargs):
    error.dendropy_construct_migration_warning(
            "dendropy.continuous.simulate_continuous",
            "dendropy.model.continuous.evolve_continuous_char")
    return continuous.evolve_continuous_char(node, rng, **kwargs)
