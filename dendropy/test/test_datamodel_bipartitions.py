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
Testing of calculation of and operations with bipartitions.
"""

import warnings
import unittest
import re
import sys
import json
from dendropy.test.support import pathmap
from dendropy.test.support import paupsplitsreference
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.utility import messaging
from dendropy.interop import paup
from dendropy.calculate import treecompare
import dendropy

# @unittest.skip('BipartitionEncodingTestCase skipped. Test in development')
class BipartitionEncodingTestCase(ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        ref_path = pathmap.tree_source_path("bipartition_encoding_fixture.json")
        with open(ref_path, "r") as src:
            cls.reference = json.load(src)

    def test_encoding(self):
        for source_name in self.reference:
            # if "multifurcating" in source_name:
            #     continue
            tree_filepath = pathmap.tree_source_path(source_name)
            for rooting in self.reference[source_name]:
                for collapse_unrooted_basal_bifurcation_desc in self.reference[source_name][rooting]:
                    if "collapse_unrooted_basal_bifurcation=True" in collapse_unrooted_basal_bifurcation_desc:
                        collapse_unrooted_basal_bifurcation = True
                    elif "collapse_unrooted_basal_bifurcation=False" in collapse_unrooted_basal_bifurcation_desc:
                        collapse_unrooted_basal_bifurcation = False
                    else:
                        raise ValueError(collapse_unrooted_basal_bifurcation_desc)
                    for suppress_unifurcations_desc in self.reference[source_name][rooting][collapse_unrooted_basal_bifurcation_desc]:
                        if "suppress_unifurcations=True" in suppress_unifurcations_desc:
                            suppress_unifurcations = True
                        elif "suppress_unifurcations=False" in suppress_unifurcations_desc:
                            suppress_unifurcations = False
                        else:
                            raise ValueError(suppress_unifurcations_desc)
                        trees_bipartitions_ref = self.reference[source_name][rooting][collapse_unrooted_basal_bifurcation_desc][suppress_unifurcations_desc]
                        source_path = pathmap.tree_source_path(source_name)
                        trees = dendropy.TreeList.get_from_path(
                                source_path,
                                "nexus",
                                rooting=rooting,
                                suppress_leaf_node_taxa=False,
                                suppress_internal_node_taxa=False,
                                )
                        for tree_idx, tree in enumerate(trees):
                            tree_bipartitions_ref = trees_bipartitions_ref[str(tree_idx)]
                            bipartition_encoding = tree.encode_bipartitions(
                                    suppress_unifurcations=suppress_unifurcations,
                                    collapse_unrooted_basal_bifurcation=collapse_unrooted_basal_bifurcation,
                                    )
                            seen = set()
                            for edge in tree.postorder_edge_iter():
                                bipartition = edge.bipartition
                                assert edge.head_node.taxon is not None
                                assert edge.head_node.taxon.label is not None
                                label = edge.head_node.taxon.label
                                # print("{}: {}: {}: {}".format(source_name, tree_idx, rooting, label, ))
                                # print("    {}".format(tree_bipartitions_ref[label]))
                                # print("    {} ({}), {}({})".format(
                                #     bipartition.split_bitmask,
                                #     bipartition.as_bitstring(),
                                #     bipartition.leafset_bitmask,
                                #     bipartition.leafset_as_bitstring(),
                                #     ))
                                expected_leafset_bitmask = int(tree_bipartitions_ref[label]["leafset_bitmask"])
                                self.assertEqual(bipartition.leafset_bitmask, expected_leafset_bitmask)
                                expected_split_bitmask = int(tree_bipartitions_ref[label]["split_bitmask"])
                                self.assertEqual(bipartition.split_bitmask, expected_split_bitmask)

if __name__ == "__main__":
    unittest.main()

