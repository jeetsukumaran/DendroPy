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

# class SimpleStructuredTreeFittingTestCase(unittest.TestCase):

#     def test_simple1(self):
#         sp_taxa = dendropy.TaxonNamespace(["a","b","c"])
#         sp_tree = dendropy.Tree.get(
#                 data="((a:2.0000,b:2.0000):3.0000,c:5.0000);",
#                 schema="newick",
#                 rooting="force-rooted",
#                 taxon_namespace=sp_taxa)
#         sp_tree.encode_bipartitions()
#         gene_taxa = dendropy.TaxonNamespace(["a","b","c"])
#         gene_tree = dendropy.Tree.get(
#                 data="((a:3.0000,b:3.0000):3.0000,c:6.0000);",
#                 schema="newick",
#                 rooting="force-rooted",
#                 taxon_namespace=gene_taxa)
#         gene_tree.encode_bipartitions()
#         msc = structuredcoalescent.StructuredCoalescent(sp_tree)

#         gs_map = {}
#         for nd in gene_tree.leaf_node_iter():
#             gs_map[nd.taxon] = sp_tree.taxon_namespace.get_taxon(label=nd.taxon.label)
#             assert gs_map[nd.taxon] is not None
#         gfn = lambda x: gs_map[x]

#         edge_head_coalescent_edges, edge_tail_coalescent_edges = msc._fit_coalescent_tree(
#                 coalescent_tree=gene_tree,
#                 coalescent_to_structure_map_fn=gfn)

#         expected_head_coalescent_edges = {
#             "001": set(["001"]),
#             "010": set(["010"]),
#             "100": set(["100"]),
#             "011": set(["001", "010"]),
#             "111": set(["011", "100"]),
#         }
#         expected_tail_coalescent_edges = {
#             "001": set(["001"]),
#             "010": set(["010"]),
#             "100": set(["100"]),
#             "011": set(["010"]),
#             "111": set(["011", ]),
#         }
#         expected_coalescing_edges = {
#             "001": set([]),
#             "010": set([]),
#             "100": set([]),
#             "011": set(["010", "001"]),
#             "111": set(["100", "011"]),
#         }
#         for structure_tree_edge in edge_head_coalescent_edges:
#             ss = structure_tree_edge.bipartition.leafset_as_bitstring()
#             # print("-- {}".format(ss))

#             cs_head = set([ce.bipartition.leafset_as_bitstring() for ce in edge_head_coalescent_edges[structure_tree_edge]])
#             # print(cs_head, expected_head_coalescent_edges[ss])
#             self.assertEqual(cs_head, expected_head_coalescent_edges[ss])

#             cs_tail = set([ce.bipartition.leafset_as_bitstring() for ce in edge_tail_coalescent_edges[structure_tree_edge]])
#             print(cs_tail, expected_tail_coalescent_edges[ss])
#             self.assertEqual(cs_tail, expected_tail_coalescent_edges[ss])

#             coalescing_edges = edge_head_coalescent_edges[structure_tree_edge] - edge_tail_coalescent_edges[structure_tree_edge]
#             cs = set([ce.bipartition.leafset_as_bitstring() for ce in coalescing_edges])
#             self.assertEqual(cs, expected_coalescing_edges[ss])

class StructuredCoalescentBasicTestCase(unittest.TestCase):

    def generate_system(self,
            speciation_ages,
            coalescent_ages,
            ):
        """
        Generates a species tree and a coalescent tree based on Figure 1 of:

            Rannala B, Yang Z. 2003. Bayesian estimation of species divergence
            ages and ancestral population sizes using DNA sequences from
            multiple loci. Genetics 164L 1645-1656.

        """
        assert len(speciation_ages) == 3
        assert len(coalescent_ages) == 6
        speciation_ages = sorted(speciation_ages)
        coalescent_ages = sorted(coalescent_ages)

        species_taxa = dendropy.TaxonNamespace(["H","C","G","O"])
        species_tree_str = "(((H,C)HC,G)HCG,O)HCGO;"
        species_tree = dendropy.Tree.get(
                data=species_tree_str,
                schema="newick",
                taxon_namespace=species_taxa,
                )
        species_taxa.is_mutable = False
        for nd in species_tree.leaf_node_iter():
            nd.age = 0.0
            nd.label = nd.taxon.label
        for nd_idx, nd in enumerate(species_tree.postorder_internal_node_iter()):
            nd.age = speciation_ages[nd_idx]
        species_tree.set_edge_lengths_from_node_ages()

        gene_taxa = dendropy.TaxonNamespace(["H1", "H2", "H3", "C1", "C2", "G", "O"])
        coalescent_tree_str = "(((H1, ((H2, H3)a,(C1, C2)b)c)d,G)e,O)f;"
        coalescent_tree = dendropy.Tree.get(
                data=coalescent_tree_str,
                schema="newick",
                taxon_namespace=gene_taxa,
                )
        gene_taxa.is_mutable = False
        for nd in coalescent_tree.leaf_node_iter():
            nd.age = 0.0
            nd.label = nd.taxon.label
        for nd_idx, nd in enumerate(coalescent_tree.postorder_internal_node_iter()):
            nd.age = coalescent_ages[nd_idx]
        coalescent_tree.set_edge_lengths_from_node_ages()

        coalescent_to_species_taxon_map = {
            gene_taxa.require_taxon("H1"): species_taxa.require_taxon("H"),
            gene_taxa.require_taxon("H2"): species_taxa.require_taxon("H"),
            gene_taxa.require_taxon("H3"): species_taxa.require_taxon("H"),
            gene_taxa.require_taxon("C1"): species_taxa.require_taxon("C"),
            gene_taxa.require_taxon("C2"): species_taxa.require_taxon("C"),
            gene_taxa.require_taxon("G"): species_taxa.require_taxon("G"),
            gene_taxa.require_taxon("O"): species_taxa.require_taxon("O"),
        }

        return species_tree, coalescent_tree, coalescent_to_species_taxon_map

    def test_fixed_species_tree_fitting(self):
        species_tree, coalescent_tree, coalescent_to_species_taxon_map = self.generate_system(
                speciation_ages=[10, 20, 30],
                coalescent_ages=[5, 6, 15, 16, 35, 36]
                )
        # print(species_tree.as_string("newick"))
        # print(coalescent_tree.as_string("newick"))
        msc = structuredcoalescent.StructuredCoalescent(species_tree)
        edge_head_coalescent_edges, edge_tail_coalescent_edges = msc._fit_coalescent_tree(
                coalescent_tree=coalescent_tree,
                coalescent_to_structure_map_fn=lambda x: coalescent_to_species_taxon_map[x])

        _get_edge = lambda t,s: t.find_node(filter_fn=lambda n: n.label==s).edge

        expected_head_coalescent_edges = {
            _get_edge(species_tree, "H"): set([
                                               _get_edge(coalescent_tree, "H1"),
                                               _get_edge(coalescent_tree, "H2"),
                                               _get_edge(coalescent_tree, "H3"),
                                               ]),
            _get_edge(species_tree, "C"): set([
                                               _get_edge(coalescent_tree, "C1"),
                                               _get_edge(coalescent_tree, "C2"),
                                               ]),
            _get_edge(species_tree, "G"): set([
                                               _get_edge(coalescent_tree, "G"),
                                              ]),
            _get_edge(species_tree, "O"): set([
                                               _get_edge(coalescent_tree, "O"),
                                              ]),
            _get_edge(species_tree, "HC"): set([
                                               _get_edge(coalescent_tree, "H1"),
                                               _get_edge(coalescent_tree, "a"),
                                               _get_edge(coalescent_tree, "b"),
                                               ]),
            _get_edge(species_tree, "HCG"): set([
                                               _get_edge(coalescent_tree, "d"),
                                               _get_edge(coalescent_tree, "G"),
                                               ]),
            _get_edge(species_tree, "HCGO"): set([
                                               _get_edge(coalescent_tree, "d"),
                                               _get_edge(coalescent_tree, "G"),
                                               _get_edge(coalescent_tree, "O"),
                                               ]),
        }

        for structure_tree_edge in edge_head_coalescent_edges:
            # print("{}: {} vs. {}".format(
            #     structure_tree_edge.head_node.label,
            #     [ce.head_node.label for ce in edge_head_coalescent_edges[structure_tree_edge]],
            #     [ce.head_node.label for ce in expected_head_coalescent_edges[structure_tree_edge]]))
            self.assertEqual(
                    set(edge_head_coalescent_edges[structure_tree_edge]),
                    expected_head_coalescent_edges[structure_tree_edge]
                    )

            # cs_tail = set([ce.bipartition.leafset_as_bitstring() for ce in edge_tail_coalescent_edges[structure_tree_edge]])
            # # print(cs_tail, expected_tail_coalescent_edges[ss])
            # self.assertEqual(cs_tail, expected_tail_coalescent_edges[ss])

            # coalescing_edges = edge_head_coalescent_edges[structure_tree_edge] - edge_tail_coalescent_edges[structure_tree_edge]
            # cs = set([ce.bipartition.leafset_as_bitstring() for ce in coalescing_edges])
            # self.assertEqual(cs, expected_coalescing_edges[ss])

if __name__ == "__main__":
    unittest.main()

