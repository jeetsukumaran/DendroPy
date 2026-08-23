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
DEPRECATED IN DENDROPY 4: USE `dendropy.model.reconcile` instead.
"""

from dendropy.utility import deprecate
from dendropy.model import reconcile

def reconciliation_discordance(gene_tree, species_tree):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.reconcile' module has moved to 'dendropy.model.reconcile'.",
            old_construct="from dendropy import reconcile\nreconcile.reconciliation_discordance(...)",
            new_construct="from dendropy.model import reconcile\nreconcile.reconciliation_discordance(...)",
            )
    return reconcile.reconciliation_discordance(gene_tree, species_tree)

def monophyletic_partition_discordance(tree, taxon_namespace_partition):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.reconcile' module has moved to 'dendropy.model.reconcile'.",
            old_construct="from dendropy import reconcile\nreconcile.monophyletic_partition_discordance(...)",
            new_construct="from dendropy.model import reconcile\nreconcile.monophyletic_partition_discordance(...)",
            )
    return reconcile.monophyletic_partition_discordance(tree, taxon_namespace_partition)

class ContainingTree(reconcile.ContainingTree):
    def __init__(self, *args, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="The 'dendropy.reconcile' module has moved to 'dendropy.model.reconcile'.",
                old_construct="from dendropy import reconcile\nreconcile.ContainingTree(...)",
                new_construct="from dendropy.model import reconcile\nreconcile.ContainingTree(...)",
                )
        reconcile.ContainingTree.__init__(self, *args, **kwargs)
