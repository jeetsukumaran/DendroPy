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
This module provides a convenient interface that aggregates, wraps, and/or
implements functions and classes that simulate trees under various
models and processes.
"""

import dendropy

###############################################################################
## Import tree generation functions

from dendropy.model.birthdeath import birth_death_tree
from dendropy.model.birthdeath import discrete_birth_death_tree
from dendropy.model.birthdeath import uniform_pure_birth_tree
from dendropy.model.coalescent import contained_coalescent_tree
from dendropy.model.coalescent import pure_kingman_tree
from dendropy.model.coalescent import mean_kingman_tree
from dendropy.model.coalescent import constrained_kingman_tree
from dendropy.model.treeshape import star_tree

