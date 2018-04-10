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
