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
Models and operations with tree shapes.
"""

import dendropy

##############################################################################
### Treeshape Generation

def star_tree(taxon_namespace, **kwargs):
    "Builds and returns a star tree from the given taxa block."
    star_tree = dendropy.Tree(taxon_namespace=taxon_namespace, **kwargs)
    for taxon in taxon_namespace:
        star_tree.seed_node.new_child(taxon=taxon)
    return star_tree

