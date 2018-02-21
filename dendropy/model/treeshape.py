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

