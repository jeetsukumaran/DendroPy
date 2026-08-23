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

