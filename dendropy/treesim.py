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
DEPRECATED IN DENDROPY 4: USE `dendropy.simulate.treesim` instead.
"""

from dendropy.utility import error
error.dendropy_module_migration_warning("dendropy.treesim", "dendropy.simulate.treesim")

# legacy
from dendropy.simulate.popgensim import pop_gen_tree
from dendropy.model.birthdeath import birth_death_tree as birth_death
from dendropy.model.birthdeath import discrete_birth_death_tree as discrete_birth_death
from dendropy.model.birthdeath import uniform_pure_birth_tree as uniform_pure_birth
from dendropy.model.coalescent import contained_coalescent_tree as contained_coalescent
from dendropy.model.coalescent import pure_kingman_tree as pure_kingman
from dendropy.model.coalescent import mean_kingman_tree as mean_kingman
from dendropy.model.coalescent import constrained_kingman_tree as constrained_kingman
from dendropy.model.treeshape import star_tree
