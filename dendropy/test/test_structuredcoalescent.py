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

import unittest
import dendropy
from dendropy.model import structuredcoalescent

# class SimpleStructuredCoalescentTestCase(unittest.TestCase):

#     def test_simple1(self):
        # sp_tree = dendropy.Tree.get(
        #         data="((a:2.0000,b:2.0000):3.0000,c:5.0000);",
        #         schema="newick",)
        # gene_tree = dendropy.Tree.get(
        #         data="((a:3.0000,b:3.0000):3.0000,c:6.0000);",
        #         schema="newick",)
        # msc = structuredcoalescent.StructuredCoalescent(sp_tree)

        # gs_map = {}
        # for nd in gene_tree.leaf_node_iter():
        #     gs_map[nd.taxon] = sp_tree.taxon_namespace.get_taxon(label=nd.taxon.label)
        #     assert gs_map[nd.taxon] is not None
        # gfn = lambda x: gs_map[x]

        # print(msc.score_coalescent_tree(
        #     coalescent_tree=gene_tree,
        #     coalescent_to_structure_map_fn=gfn))

class SimpleStructuredTreeFittingTestCase(unittest.TestCase):

    def test_simple1(self):
        sp_taxa = dendropy.TaxonNamespace(["a","b","c"])
        sp_tree = dendropy.Tree.get(
                data="((a:2.0000,b:2.0000):3.0000,c:5.0000);",
                schema="newick",
                rooting="force-rooted",
                taxon_namespace=sp_taxa)
        sp_tree.encode_bipartitions()
        gene_taxa = dendropy.TaxonNamespace(["a","b","c"])
        gene_tree = dendropy.Tree.get(
                data="((a:3.0000,b:3.0000):3.0000,c:6.0000);",
                schema="newick",
                rooting="force-rooted",
                taxon_namespace=gene_taxa)
        gene_tree.encode_bipartitions()
        msc = structuredcoalescent.StructuredCoalescent(sp_tree)

        gs_map = {}
        for nd in gene_tree.leaf_node_iter():
            gs_map[nd.taxon] = sp_tree.taxon_namespace.get_taxon(label=nd.taxon.label)
            assert gs_map[nd.taxon] is not None
        gfn = lambda x: gs_map[x]

        edge_head_coalescent_edges, edge_tail_coalescent_edges = msc._fit_coalescent_tree(
                coalescent_tree=gene_tree,
                coalescent_to_structure_map_fn=gfn)

        expected_head_coalescent_edges = {
            "001": set(["001"]),
            "010": set(["010"]),
            "100": set(["100"]),
            "011": set(["001", "010"]),
            "111": set(["011", "100"]),
        }
        expected_tail_coalescent_edges = {
            "001": set(["001"]),
            "010": set(["010"]),
            "100": set(["100"]),
            "011": set(["010"]),
            "111": set(["011", "100"]),
        }
        expected_coalescing_edges = {
            "001": set([]),
            "010": set([]),
            "100": set([]),
            "011": set(["010", "001"]),
            "111": set(["100", "011"]),
        }
        for structure_tree_edge in edge_head_coalescent_edges:
            ss = structure_tree_edge.bipartition.leafset_as_bitstring()
            # print("-- {}".format(ss))

            cs_head = set([ce.bipartition.leafset_as_bitstring() for ce in edge_head_coalescent_edges[structure_tree_edge]])
            # print(cs_head, expected_head_coalescent_edges[ss])
            self.assertEqual(cs_head, expected_head_coalescent_edges[ss])

            cs_head = set([ce.bipartition.leafset_as_bitstring() for ce in edge_head_coalescent_edges[structure_tree_edge]])
            # print(cs_head, expected_head_coalescent_edges[ss])
            self.assertEqual(cs_head, expected_head_coalescent_edges[ss])

            coalescing_edges = edge_head_coalescent_edges[structure_tree_edge] - edge_tail_coalescent_edges[structure_tree_edge]
            cs = set([ce.bipartition.leafset_as_bitstring() for ce in coalescing_edges])
            self.assertEqual(cs, expected_coalescing_edges[ss])

if __name__ == "__main__":
    unittest.main()

