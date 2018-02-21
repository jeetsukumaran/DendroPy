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
DEPRECATED IN DENDROPY 4: USE `dendropy.model.coalescent` instead.
"""

from dendropy.model import coalescent
from dendropy.utility import deprecate
from dendropy.utility import constants

def discrete_time_to_coalescence(n_genes,
                                 pop_size=None,
                                 n_to_coalesce=2,
                                 rng=None):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.discrete_time_to_coalescence(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.discrete_time_to_coalescence(...)")
    return coalescent.discrete_time_to_coalescence(
            n_genes=n_genes,
            pop_size=pop_size,
            n_to_coalesce=n_to_coalesce,
            rng=rng)

def time_to_coalescence(n_genes,
        pop_size=None,
        n_to_coalesce=2,
        rng=None):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.time_to_coalescence(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.time_to_coalescence(...)")
    return coalescent.time_to_coalescence(
        n_genes=n_genes,
        pop_size=pop_size,
        n_to_coalesce=n_to_coalesce,
        rng=rng)

def expected_tmrca(n_genes, pop_size=None, n_to_coalesce=2):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.expected_tmrca(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.expected_tmrca(...)")
    return coalescent.expected_tmrca(n_genes, pop_size=pop_size, n_to_coalesce=n_to_coalesce)

def coalesce(nodes,
             pop_size=None,
             period=None,
             rng=None,
             use_expected_tmrca=False):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent', and this function has been renamed 'coalesce_nodes'.",
            old_construct="from dendropy import coalescent\ncoalescent.coalesce(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.coalesce_nodes(...)")
    return coalescent.coalesce_nodes(
            nodes=nodes,
            pop_size=pop_size,
            period=period,
            rng=rng,
            use_expected_tmrca=use_expected_tmrca)

def node_waiting_time_pairs(tree, ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.node_waiting_time_pairs(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.node_waiting_time_pairs(...)")
    return coalescent.node_waiting_time_pairs(tree=tree,
            ultrametricity_precision=ultrametricity_precision)

def extract_coalescent_frames(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.extract_coalescent_frames(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.extract_coalescent_frames(...)")
    return coalescent.extract_coalescent_frames(*args, **kwargs)

def log_probability_of_coalescent_frames(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.log_probability_of_coalescent_frames(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.log_probability_of_coalescent_frames(...)")
    return log_probability_of_coalescent_frames(*args, **kwargs)

def log_probability_of_coalescent_tree(*args, **kwargs):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.coalescent' module has moved to 'dendropy.model.coalescent'.",
            old_construct="from dendropy import coalescent\ncoalescent.log_probability_of_coalescent_tree(...)",
            new_construct="from dendropy.model import coalescent\ncoalescent.log_probability_of_coalescent_tree(...)")
    return log_probability_of_coalescent_tree(*args, **kwargs)

