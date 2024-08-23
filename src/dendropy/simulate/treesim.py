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

import dendropy
import random
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
    ]


# helps each common functional interface implementation wrap call to wrapped
# library implementation
# def _normalize_args_and_kwargs(kwargs, model_args_and_defaults):
#     args = []
#     model_arg_values = {}
#     for model_arg_name, default in model_args_and_defaults:
#         value = kwargs.pop(model_arg_name, default)
#         model_arg_values[model_arg_name] = value
#         args.append(value)
#     return args, kwargs, model_arg_values

# def _setup_shared_context(rng, kwargs):
#     context = {}
#     context["rng"] = rng
#     if "taxon_namespace" in kwargs:
#         context["taxon_namespace"] = kwargs["taxon_namespace"]
#     else:
#         context["taxon_namespace"] = dendropy.TaxonNamespace()
#     return context

def iter_birthdeath_trees(
    rng,
    model_kwargs_fn,
    n_replicates,
):
    return_value = []
    for rep_idx in range(n_replicates):
        if rng is None:
            rng = random.Random()
        model_kwargs = model_kwargs_fn(rep_idx, rng)
        tree = birth_death_tree(**model_kwargs)
        yield tree

def map_birthdeath_trees(
    fn,
    rng,
    model_kwargs_fn,
    n_replicates,
):
    for tree in iter_birthdeath_trees(
        rng=rng,
        model_kwargs_fn=model_kwargs_fn,
        n_replicates=n_replicates,
    ):
        yield fn(tree)

def mapped_birthdeath_trees(
    fn,
    rng,
    model_kwargs_fn,
    n_replicates,
):
    """
    Usage:

        #! /usr/bin/env python
        # -*- coding: utf-8 -*-

        import random
        from dendropy.simulate import treesim
        from dendropy.calculate import treemeasure

        def run(n_replicates, min_leaves=10, max_leaves=1000, ):
            rng = random.Random()
            n_replicates = 10
            result = treesim.mapped_birthdeath_trees(
                treemeasure.coalescence_ages,
                rng=rng,
                model_kwargs_fn=lambda rep_idx, rrng: {
                    "birth_rate": 1.0,
                    "death_rate": 0.0,
                    "num_extant_tips": rrng.randint(min_leaves, max_leaves),
                },
                n_replicates=n_replicates,
            )
            return result

        print(run(1000))

    """
    return [*map_birthdeath_trees(fn, rng, model_kwargs_fn, n_replicates)]

def birthdeath_coalescence_ages(rng, model_kwargs_fn, n_replicates):
    """
    Usage:

        def run(n_replicates):
            print(treesim.birthdeath_coalescence_ages(None, lambda rep_idx, rng: {
                "birth_rate": 1.0,
                "death_rate": 0.0,
                "num_extant_tips": rng.randint(10, 100),
                # "gsa_ntax": 5000,
                }, n_replicates))

        run(100)
    """
    result = mapped_birthdeath_trees(
        treemeasure.coalescence_ages,
        rng=rng,
        model_kwargs_fn=model_kwargs_fn,
        n_replicates=n_replicates,
    )
    return result


