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

## Required for Sphix auto-documentation of this module
__all__ = [
    "birth_death_tree",
    "discrete_birth_death_tree",
    "contained_coalescent_tree",
    "pure_kingman_tree",
    "mean_kingman_tree",
    "constrained_kingman_tree",
    "star_tree",
    "birthdeath_trees_iter",
    ]

def _normalize_kwargs(kwargs, model_args_and_defaults):
    args = []
    model_arg_values = {}
    for model_arg_name, default in model_args_and_defaults:
        value = kwargs.pop(model_arg_name, default)
        model_arg_values[model_arg_name] = value
        args.append(value)
    return args, kwargs, model_arg_values


def birthdeath_trees_iter(
    rng,
    model_kwargs,
    n_replicates,
):
    args, kwargs, model_arg_values = _normalize_kwargs(
                model_kwargs,
                (
                    ("birth_rate", 1.0),
                    ("death_rate", 1.0),
                    ("birth_rate_sd", 0.0),
                    ("death_rate_sd", 0.0),
                ))
    kwargs["rng"] = rng
    kwargs["taxon_namespace"] = dendropy.TaxonNamespace()
    return_value = []
    for rep_idx in range(n_replicates):
        tree = birth_death_tree(*args, **kwargs)
        yield tree

