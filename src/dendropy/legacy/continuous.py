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
DEPRECATED IN DENDROPY 4: USE `dendropy.model.continuous` instead.
"""

from dendropy.model import continuous
from dendropy.utility import deprecate

def simulate_continuous(node, rng=None, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.continuous' module has moved to 'dendropy.model.continuous', and this function has been renamed 'evolve_continuous_char()'.",
            old_construct="from dendropy import continuous\ncontinuous.simulate_continuous(...)",
            new_construct="from dendropy.model import continuous\ncontinuous.evolve_continuous_char(...)",
            epilog="Note that this function is also available through 'dendropy.simulate.charsim.evolve_continuous_char(...)'.")
    return continuous.evolve_continuous_char(node, rng, **kwargs)

class PhylogeneticIndependentContrasts(continuous.PhylogeneticIndependentContrasts):

    def __init__(self,
            tree,
            char_matrix,
            polytomy_strategy=None):
        deprecate.dendropy_deprecation_warning(
                preamble="The 'dendropy.continuous' module has moved to 'dendropy.model.continuous'.",
                old_construct="from dendropy import continuous\ncontinuous.PhylogeneticIndependentContrasts(...)",
                new_construct="from dendropy.model import continuous\ncontinuous.PhylogeneticIndependentContrasts(...)",
                )
        continuous.PhylogeneticIndependentContrasts.__init__(self,
                tree=tree,
                char_matrix=char_matrix,
                polytomy_strategy=polytomy_strategy)

