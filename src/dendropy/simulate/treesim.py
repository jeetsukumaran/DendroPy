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
This module provides a convenient interface that aggregates, wraps, and/or
implements functions and classes that simulate trees under various
models and processes. This module just exposes these function and classes under
the ``dendropy.simulate.treesim`` namespace. The actual functions and classes
are defined under the the appropriate model namespace in the ``dendropy.model``
sub-package.
"""

###############################################################################
## Import tree generation functions

import random
import collections
import dendropy
from dendropy.model.birthdeath import birth_death_tree
from dendropy.model.birthdeath import discrete_birth_death_tree
from dendropy.model.birthdeath import uniform_pure_birth_tree
from dendropy.model.coalescent import contained_coalescent_tree
from dendropy.model.coalescent import pure_kingman_tree
from dendropy.model.coalescent import mean_kingman_tree
from dendropy.model.coalescent import constrained_kingman_tree
from dendropy.model.treeshape import star_tree
from dendropy.simulate import treesim
from dendropy.calculate import treemeasure

## Required for Sphix auto-documentation of this module
__all__ = [
    "birth_death_tree",
    "discrete_birth_death_tree",
    "contained_coalescent_tree",
    "pure_kingman_tree",
    "mean_kingman_tree",
    "constrained_kingman_tree",
    "star_tree",
    "rand_trees",
    ]


## Following are wrappers to for more efficient vectorized sampling, all using
## the general interface, `f([rng], <model/sim parameters>, <size>)`
#
#   The model parameters may be specified as:
#   -   A single dict or map, in which case it will be repeated
#       through the replicates,
#   -   A function, in which case it will be called for each replicate
#       with two positional arguments: the 0-based replicate index and
#       the random numberg generator object to use; the function should
#       return a dict or mapping of keyword-value pairs for the model
#       simulation call.
#   -   An iterable of dicts or maps, for *each* of which
#       `n_replicates` simulations will be generated, in order.
#
#  Example::
#
#    import itertools
#    from dendropy.simulate import treesim
#
#    fixed_params = { "birth_rate": 1.0, "death_rate": 0.0, "num_extant_tips": 20, }
#
#    fixed_params_list_data = [
#        { "birth_rate": 1.0, "death_rate": 0.0, "num_extant_tips": 20, }
#        for birth_rate in (0.1, 0.2, 0.5, 1.0)
#    ]
#
#    def kwargs_generator_fn(rep_idx, rng):
#        return { "birth_rate": rng.uniform(0.1, 1.0), "death_rate": 0.0, "num_extant_tips": rng.randint(10, 20), }
#
#    def run(n_replicates):
#        print(treesim.birthdeath_coalescence_ages(None, fixed_params, n_replicates))
#        print(treesim.birthdeath_coalescence_ages(None, kwargs_generator_fn, n_replicates))
#        print(treesim.birthdeath_coalescence_ages(None, fixed_params_list_data, n_replicates))
#
#    run(10)

def rand_trees(
    rng,
    model_fn,
    model_kwargs,
    n_replicates,
):
    """
    The model parameters may be specified as:
    -   A single dict or map, in which case it will be repeated
        through the replicates,
    -   A function, in which case it will be called for each replicate
        with two positional arguments: the 0-based replicate index and
        the random numberg generator object to use; the function should
        return a dict or mapping of keyword-value pairs for the model
        simulation call.
    -   An iterable of dicts or maps, for *each* of which
        `n_replicates` simulations will be generated, in order.
    """
    for rep_idx in range(n_replicates):
        if rng is None:
            rng = random.Random()
        if isinstance(model_kwargs, collections.abc.Mapping):
            model_kwargs_data = model_kwargs
        elif callable(model_kwargs):
            model_kwargs_data = model_kwargs(rep_idx, rng)
        elif not isinstance(model_kwargs, str):
            # raise ValueError(f"Assuming `model_kwargs` as iterable of keyword maps: expecting same number of model keyword maps ({len(model_kwargs)}) as number of replicates ({n_replicates})")
            for model_kwargs_data in model_kwargs:
                tree = birth_death_tree(**model_kwargs_data)
                yield tree
            continue
        else:
            model_kwargs_data = {}
        tree = model_fn(**model_kwargs_data)
        yield tree

def coalescence_ages(rng, model_fn, model_kwargs, n_replicates):
    result = [*map(
                treemeasure.coalescence_ages,
                rand_trees(rng, model_fn, model_kwargs, n_replicates,),
                )]
    return result

def birthdeath_coalescence_ages(rng, model_kwargs, n_replicates):
    result = [*map(
                treemeasure.coalescence_ages,
                rand_trees(rng, birth_death_tree, model_kwargs, n_replicates,),
                )]
    return result


