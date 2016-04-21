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
import json
import dendropy
from dendropy.model import multispeciescoalescent
from dendropy.test.support import pathmap

def generate_multispecies_coalescent_system(
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
    speciation_ages = sorted(float(i) for i in speciation_ages)
    coalescent_ages = sorted(float(i) for i in coalescent_ages)

    species_taxa = dendropy.TaxonNamespace(["H","C","G","O"])
    species_tree_str = "(((H,C)HC,G)HCG,O)HCGO;"
    species_tree = dendropy.Tree.get(
            data=species_tree_str,
            schema="newick",
            taxon_namespace=species_taxa,
            rooting="force-rooted",
            )
    species_taxa.is_mutable = False
    for nd in species_tree.leaf_node_iter():
        nd.age = 0.0
        nd.label = nd.taxon.label
    species_tree.find_node_with_label("HC").age = speciation_ages[0]
    species_tree.find_node_with_label("HCG").age = speciation_ages[1]
    species_tree.find_node_with_label("HCGO").age = speciation_ages[2]
    species_tree.set_edge_lengths_from_node_ages()

    gene_taxa = dendropy.TaxonNamespace(["H1", "H2", "H3", "C1", "C2", "G", "O"])
    coalescent_tree_str = "(((H1, ((H2, H3)a,(C1, C2)b)c)d,G)e,O)f;"
    coalescent_tree = dendropy.Tree.get(
            data=coalescent_tree_str,
            schema="newick",
            taxon_namespace=gene_taxa,
            rooting="force-rooted",
            )
    gene_taxa.is_mutable = False
    for nd in coalescent_tree.leaf_node_iter():
        nd.age = 0.0
        nd.label = nd.taxon.label
    coalescent_tree.find_node_with_label("a").age = coalescent_ages[0]
    coalescent_tree.find_node_with_label("b").age = coalescent_ages[1]
    coalescent_tree.find_node_with_label("c").age = coalescent_ages[2]
    coalescent_tree.find_node_with_label("d").age = coalescent_ages[3]
    coalescent_tree.find_node_with_label("e").age = coalescent_ages[4]
    coalescent_tree.find_node_with_label("f").age = coalescent_ages[5]
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

class MultispeciesCoalescentFixedSingleTreesCalculationTestCase(unittest.TestCase):

    def test1(self):
        with open(pathmap.other_source_path("multispecies_coalescent_test_data.json")) as src:
            test_regimes = json.load(src)
        for test_regime in test_regimes:
            species_tree = dendropy.Tree.get(
                    data=test_regime["species_tree"],
                    schema="newick",
                    rooting="force-rooted",
                    )
            species_tree.taxon_namespace.is_mutable = False
            msc = multispeciescoalescent.MultispeciesCoalescent(species_tree=species_tree)
            coalescent_species_lineage_label_map = test_regime["coalescent_species_lineage_label_map"]
            coalescent_species_lineage_map_fn = lambda x: species_tree.taxon_namespace.require_taxon(coalescent_species_lineage_label_map[x.label])
            coalescent_taxa = dendropy.TaxonNamespace(sorted(coalescent_species_lineage_label_map.keys()))
            coalescent_taxa.is_mutable = False
            for sub_regime in test_regime["coalescent_trees"]:
                coalescent_tree = dendropy.Tree.get(
                        data=sub_regime["coalescent_tree"],
                        schema="newick",
                        rooting="force-rooted",
                        taxon_namespace=coalescent_taxa,
                        )
                obs_ln_likelihood = msc.score_coalescent_tree(
                        coalescent_tree=coalescent_tree,
                        coalescent_species_lineage_map_fn=coalescent_species_lineage_map_fn,
                        )
                exp_ln_likelihood = sub_regime["log_likelihood"]
                self.assertAlmostEqual(obs_ln_likelihood, exp_ln_likelihood, 2)

class MultispeciesCoalescentBasicTestCase(unittest.TestCase):

    def calc_log_likelihood(self,
            species_tree,
            coalescent_tree,
            thetas=None,
            default_theta=1.0,
            ):
        tau_HC = species_tree.find_node_with_label("HC").age
        tau_HCG = species_tree.find_node_with_label("HCG").age
        tau_HCGO = species_tree.find_node_with_label("HCGO").age
        t3_H = coalescent_tree.find_node_with_label("a").age
        t2_C = coalescent_tree.find_node_with_label("b").age
        t3_HC = coalescent_tree.find_node_with_label("c").age - species_tree.find_node_with_label("HC").age
        t2_HC = coalescent_tree.find_node_with_label("d").age - coalescent_tree.find_node_with_label("c").age
        t3_HCGO = coalescent_tree.find_node_with_label("e").age - species_tree.find_node_with_label("HCGO").age
        t2_HCGO = coalescent_tree.find_node_with_label("f").age - coalescent_tree.find_node_with_label("e").age
        if thetas is None:
            thetas = {}
        theta_H = thetas.get("H", default_theta)
        theta_C = thetas.get("C", default_theta)
        theta_HC = thetas.get("HC", default_theta)
        theta_HCG = thetas.get("HCG", default_theta)
        theta_HCGO = thetas.get("HCGO", default_theta)
        p1 = 2.0/theta_H * math.exp(-6 * t3_H/theta_H) * math.exp(-2 * (tau_HC-t3_H)/theta_H)
        p2 = 2.0/theta_C * math.exp(-2 * t2_C/theta_C)
        p3 = 2.0/theta_HC * math.exp(-6 * t3_HC/theta_HC) * 2.0/theta_HC * math.exp(-2 * t2_HC/theta_HC)

        # below is as originally given in the Rannala and Yang paper,
        # but is incorrect
        # p4 = math.exp(-2 * (tau_HCG - tau_HC - (t3_HC + t2_HC)) / theta_HCG)
        p4 = math.exp(-2 * (tau_HCGO - tau_HCG) / theta_HCG)

        p5 = 2.0/theta_HCGO * math.exp(-6 * t3_HCGO/theta_HCGO)
        p6 = 2.0/theta_HCGO * math.exp(-2 * t2_HCGO/theta_HCGO)
        p = p1 * p2 * p3 * p4 * p5 * p6

        q1 = math.log(2.0/theta_H) + (-6 * t3_H/theta_H) + (-2 * (tau_HC-t3_H)/theta_H)
        q2 = math.log(2.0/theta_C) + (-2 * t2_C/theta_C)
        q3 = math.log(2.0/theta_HC) + (-6 * t3_HC/theta_HC) + math.log(2.0/theta_HC) + (-2 * t2_HC/theta_HC)
        # below is as originally given in the Rannala and Yang paper,
        # but is incorrect
        # q4 = (-2 * (tau_HCG - tau_HC - (t3_HC + t2_HC)) / theta_HCG)
        q4 = (-2 * (tau_HCGO - tau_HCG) / theta_HCG)
        q5 = math.log(2.0/theta_HCGO) + (-6 * t3_HCGO/theta_HCGO)
        q6 = math.log(2.0/theta_HCGO) + (-2 * t2_HCGO/theta_HCGO)
        q = q1 + q2 + q3 + q4 + q5 + q6

        self.assertAlmostEqual(math.log(p), q, 7)

        return q

    def get_node(self, tree, label):
        return tree.find_node(filter_fn=lambda n: n.label==label)

    def get_edge(self, tree, label):
        return tree.find_node(filter_fn=lambda n: n.label==label).edge

    def test_fixed_species_tree_fitting(self):
        species_tree, coalescent_tree, coalescent_to_species_taxon_map = generate_multispecies_coalescent_system(
                speciation_ages=[10, 20, 30],
                coalescent_ages=[5, 6, 15, 16, 35, 36]
                )
        msc = multispeciescoalescent.MultispeciesCoalescent(species_tree)
        edge_head_coalescent_edges, edge_tail_coalescent_edges, edge_coalescent_nodes = msc._fit_coalescent_tree(
                coalescent_tree=coalescent_tree,
                coalescent_species_lineage_map_fn=lambda x: coalescent_to_species_taxon_map[x])

        expected_head_coalescent_edges = {
            self.get_edge(species_tree, "H"): set([
                                               self.get_edge(coalescent_tree, "H1"),
                                               self.get_edge(coalescent_tree, "H2"),
                                               self.get_edge(coalescent_tree, "H3"),
                                               ]),
            self.get_edge(species_tree, "C"): set([
                                               self.get_edge(coalescent_tree, "C1"),
                                               self.get_edge(coalescent_tree, "C2"),
                                               ]),
            self.get_edge(species_tree, "G"): set([
                                               self.get_edge(coalescent_tree, "G"),
                                              ]),
            self.get_edge(species_tree, "O"): set([
                                               self.get_edge(coalescent_tree, "O"),
                                              ]),
            self.get_edge(species_tree, "HC"): set([
                                               self.get_edge(coalescent_tree, "H1"),
                                               self.get_edge(coalescent_tree, "a"),
                                               self.get_edge(coalescent_tree, "b"),
                                               ]),
            self.get_edge(species_tree, "HCG"): set([
                                               self.get_edge(coalescent_tree, "d"),
                                               self.get_edge(coalescent_tree, "G"),
                                               ]),
            self.get_edge(species_tree, "HCGO"): set([
                                               self.get_edge(coalescent_tree, "d"),
                                               self.get_edge(coalescent_tree, "G"),
                                               self.get_edge(coalescent_tree, "O"),
                                               ]),
        }
        expected_tail_coalescent_edges = {
            self.get_edge(species_tree, "H"): set([
                                               self.get_edge(coalescent_tree, "H1"),
                                               self.get_edge(coalescent_tree, "a"),
                                               ]),
            self.get_edge(species_tree, "C"): set([
                                               self.get_edge(coalescent_tree, "b"),
                                               ]),
            self.get_edge(species_tree, "G"): set([
                                               self.get_edge(coalescent_tree, "G"),
                                              ]),
            self.get_edge(species_tree, "O"): set([
                                               self.get_edge(coalescent_tree, "O"),
                                              ]),
            self.get_edge(species_tree, "HC"): set([
                                               self.get_edge(coalescent_tree, "d"),
                                               ]),
            self.get_edge(species_tree, "HCG"): set([
                                               self.get_edge(coalescent_tree, "d"),
                                               self.get_edge(coalescent_tree, "G"),
                                               ]),
            self.get_edge(species_tree, "HCGO"): set([
                                               ]),
        }
        expected_coalescing_nodes = {
            self.get_edge(species_tree, "H"): set([
                                               self.get_node(coalescent_tree, "a"),
                                               ]),
            self.get_edge(species_tree, "C"): set([
                                               self.get_node(coalescent_tree, "b"),
                                               ]),
            self.get_edge(species_tree, "G"): set([
                                              ]),
            self.get_edge(species_tree, "O"): set([
                                              ]),
            self.get_edge(species_tree, "HC"): set([
                                               self.get_node(coalescent_tree, "d"),
                                               self.get_node(coalescent_tree, "c"),
                                               ]),
            self.get_edge(species_tree, "HCG"): set([
                                               ]),
            self.get_edge(species_tree, "HCGO"): set([
                                               self.get_node(coalescent_tree, "f"),
                                               self.get_node(coalescent_tree, "e"),
                                               ]),
        }

        for species_tree_edge in edge_head_coalescent_edges:
            # print("-- {} --".format(species_tree_edge.head_node.label if species_tree_edge.head_node else "<root>"))
            # print("{}: {} vs. {}".format(
            #     species_tree_edge.head_node.label,
            #     [ce.head_node.label for ce in edge_head_coalescent_edges[species_tree_edge]],
            #     [ce.head_node.label for ce in expected_head_coalescent_edges[species_tree_edge]]))
            self.assertEqual(
                    set(edge_head_coalescent_edges[species_tree_edge]),
                    expected_head_coalescent_edges[species_tree_edge]
                    )
            # print("{}: {} vs. {}".format(
            #     species_tree_edge.head_node.label if species_tree_edge.head_node else "<root>",
            #     [ce.head_node.label for ce in edge_tail_coalescent_edges[species_tree_edge]],
            #     [ce.head_node.label for ce in expected_tail_coalescent_edges[species_tree_edge]]))
            self.assertEqual(
                    set(edge_tail_coalescent_edges[species_tree_edge]),
                    expected_tail_coalescent_edges[species_tree_edge]
                    )
            # print("{}: {} vs. {}".format(
            #     species_tree_edge.head_node.label if species_tree_edge.head_node else "<root>",
            #     [nd.label for nd in edge_coalescent_nodes[species_tree_edge]],
            #     [nd.label for nd in expected_coalescing_nodes[species_tree_edge]]))
            self.assertEqual(
                    set(edge_coalescent_nodes[species_tree_edge]),
                    expected_coalescing_nodes[species_tree_edge]
                    )

        expected_lnL = self.calc_log_likelihood(
                species_tree=species_tree,
                coalescent_tree=coalescent_tree,)
        s = msc.score_coalescent_tree(
                coalescent_tree=coalescent_tree,
                coalescent_species_lineage_map_fn=lambda x: coalescent_to_species_taxon_map[x],
                )
        self.assertAlmostEqual(s, expected_lnL)

if __name__ == "__main__":
    unittest.main()

