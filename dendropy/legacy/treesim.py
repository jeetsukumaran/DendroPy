#! /usr/bin/env python

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
DEPRECATED IN DENDROPY 4: USE `dendropy.simulate.treesim` instead.
"""

from dendropy.simulate import treesim
from dendropy.simulate import popgensim
from dendropy.utility import deprecate

def pop_gen_tree(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.pop_gen_tree()' function has moved to 'dendropy.simulate.popgensim.pop_gen_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.pop_gen_tree(...)",
            new_construct="from dendropy.simulate import popgensim\ntree = popgensim.pop_gen_tree(...)")
    return popgensim.pop_gen_tree(*args, **kwargs)

def birth_death(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.birth_death()' function has moved to 'dendropy.simulate.treesim.birth_death_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.birth_death(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.birth_death_tree(...)")
    return treesim.birth_death_tree(*args, **kwargs)

# from dendropy.model.birthdeath import discrete_birth_death_tree as discrete_birth_death
def discrete_birth_death(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.discrete_birth_death()' function has moved to 'dendropy.simulate.treesim.discrete_birth_death_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.discrete_birth_death(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.discrete_birth_death_tree(...)")
    return treesim.discrete_birth_death_tree(*args, **kwargs)

# from dendropy.model.birthdeath import uniform_pure_birth_tree as uniform_pure_birth
def uniform_pure_birth(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.uniform_pure_birth()' function has moved to 'dendropy.simulate.treesim.uniform_pure_birth_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.uniform_pure_birth(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.uniform_pure_birth_tree(...)")
    return treesim.uniform_pure_birth_tree(*args, **kwargs)

# from dendropy.model.coalescent import contained_coalescent_tree as contained_coalescent
def contained_coalescent(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.contained_coalescent()' function has moved to 'dendropy.simulate.treesim.contained_coalescent_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.contained_coalescent(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.contained_coalescent_tree(...)")
    return treesim.contained_coalescent_tree(*args, **kwargs)

# from dendropy.model.coalescent import pure_kingman_tree as pure_kingman
def pure_kingman(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.pure_kingman()' function has moved to 'dendropy.simulate.treesim.pure_kingman_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.pure_kingman(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.pure_kingman_tree(...)")
    return treesim.pure_kingman_tree(*args, **kwargs)

# from dendropy.model.coalescent import mean_kingman_tree as mean_kingman
def mean_kingman(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.mean_kingman()' function has moved to 'dendropy.simulate.treesim.mean_kingman_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.mean_kingman(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.mean_kingman_tree(...)")
    return treesim.mean_kingman_tree(*args, **kwargs)

# from dendropy.model.coalescent import constrained_kingman_tree as constrained_kingman
def constrained_kingman(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.constrained_kingman()' function has moved to 'dendropy.simulate.treesim.constrained_kingman_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.constrained_kingman(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.constrained_kingman_tree(...)")
    return treesim.constrained_kingman_tree(*args, **kwargs)

# from dendropy.model.treeshape import star_tree
def star_tree(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesim.star()' function has moved to 'dendropy.simulate.treesim.star_tree()'.",
            old_construct="from dendropy import treesim\ntree = treesim.star(...)",
            new_construct="from dendropy.simulate import treesim\ntree = treesim.star_tree(...)")
    return treesim.star_tree(*args, **kwargs)
