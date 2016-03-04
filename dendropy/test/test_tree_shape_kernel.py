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

import math
import unittest
import dendropy
from dendropy.calculate.treecompare import TreeShapeKernel

class PhyloNode(dendropy.Node):

    def __init__(self, *args, **kwargs):
        dendropy.Node.__init__(self, *args, **kwargs)

    def _get_clades(self):
        return self._child_nodes
    clades = property(_get_clades)

    def _get_branch_length(self):
        return self.edge.length
    branch_length = property(_get_branch_length)

class PhyloTree(dendropy.Tree):

    def node_factory(cls, *args, **kwargs):
        return PhyloNode(*args, **kwargs)

    def get_nonterminals(self, order):
        if order == "postorder":
            return list(self.postorder_internal_node_iter())
        else:
            raise ValueError(order)

    def get_terminals(self):
        return list(self.leaf_node_iter())

class PhyloTreeList(dendropy.TreeList):

    DEFAULT_TREE_TYPE = PhyloTree

    def tree_factory(*args, **kwargs):
        return PhyloTree(*args, **kwargs)

class TreeShapeKernelTests(unittest.TestCase):

    def test_small_single(self):
        ## Test based on:
        ##      KAMPHIR
        ##      By Art F.Y. Poon and Rosemary McCloskey
        ##      https://github.com/ArtPoon/kamphir.git
        T1_str = "( ( A:0.5, B:0.25 )E:0.5, ( C:0.25, D:0.25 )F:0.5 )G;"
        T2_str = "( ( ( A:0.25, B:0.25 )E:0.5, C:0.25 )F:0.5, D:0.25 )G;"
        taxon_namespace = dendropy.TaxonNamespace()
        T1 = PhyloTree.get_from_string(T1_str, "newick", taxon_namespace=taxon_namespace)
        T2 = PhyloTree.get_from_string(T2_str, "newick", taxon_namespace=taxon_namespace)
        tree_shape_kernel = TreeShapeKernel(decay_factor=0.5, gauss_factor=1)
        assert tree_shape_kernel(T1, T2) == 1.125 * (1+math.exp(-0.0625))

if __name__ == "__main__":
    unittest.main()

