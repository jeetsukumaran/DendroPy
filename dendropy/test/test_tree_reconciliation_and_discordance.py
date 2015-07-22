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
Tests reconciliation calculations.
"""

import os
import unittest
import dendropy
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

from dendropy.model import reconcile

class ContainingTreeDeepCoalescenceSmall(unittest.TestCase):

    def setUp(self):
        dataset = dendropy.DataSet.get_from_path(pathmap.tree_source_path(filename="deepcoal1.nex"), "nexus")
        self.species_tree = dataset.get_tree_list(label="ContainingTree")[0]
        self.gene_trees = dataset.get_tree_list(label="EmbeddedTrees")
        self.species_tree.taxon_namespace.is_mutable = False
        self.gene_taxon_to_population_taxon_map = dendropy.TaxonNamespaceMapping(
                domain_taxon_namespace=self.gene_trees.taxon_namespace,
                range_taxon_namespace=self.species_tree.taxon_namespace,
                mapping_fn=lambda t: self.species_tree.taxon_namespace.require_taxon(label=t.label[0].upper()))
        self.expected_under_original_brlens = [4, 6, 4, 2, 4, 3, 3, 4, 5, 4]

    def testFixedEdgesDeepCoalCount(self):
        results = []
        for idx, gt in enumerate(self.gene_trees):
            ct = reconcile.ContainingTree(containing_tree=self.species_tree,
                    contained_taxon_namespace=self.gene_trees.taxon_namespace,
                    contained_to_containing_taxon_map=self.gene_taxon_to_population_taxon_map,
                    contained_trees=[gt],
                    fit_containing_edge_lengths=False,
                    )
            dc = ct.num_deep_coalescences()
            results.append(dc)

            ## FOR DEBUGGING
            # mesqf = pathmap.named_output_stream("ContainingTreeDeepCoalescence_Small_FixedEdges_t%02d_dc%02d.nex" % (idx+1, dc), False)
            # with mesqf:
            #     ct.write_as_mesquite(mesqf)

        self.assertEqual(results, self.expected_under_original_brlens)

    def testFittedEdgesDeepCoalCount(self):
        for idx, gt in enumerate(self.gene_trees):
            gt.encode_bipartitions()
            ct = reconcile.ContainingTree(containing_tree=self.species_tree,
                    contained_taxon_namespace=self.gene_trees.taxon_namespace,
                    contained_to_containing_taxon_map=self.gene_taxon_to_population_taxon_map,
                    contained_trees=[gt],
                    fit_containing_edge_lengths=True,
                    )
            dc = ct.num_deep_coalescences()

            ## FOR DEBUGGING
            # mesqf = pathmap.named_output_stream("ContainingTreeDeepCoalescence_Small_FittedEdges_t%02d_dc%02d.nex" % (idx+1, dc), False)
            # with mesqf:
            #     ct.write_as_mesquite(mesqf)

class DeepCoalTest(unittest.TestCase):

    def testFittedDeepCoalCounting(self):

        taxa = dendropy.TaxonNamespace()

        gene_trees = dendropy.TreeList.get_from_string("""
            [&R] (A,(B,(C,D))); [&R] ((A,C),(B,D)); [&R] (C,(A,(B,D)));
            """, "newick", taxon_namespace=taxa)

        species_trees = dendropy.TreeList.get_from_string("""
            [&R] (A,(B,(C,D)));
            [&R] (A,(C,(B,D)));
            [&R] (A,(D,(C,B)));
            [&R] (B,(A,(C,D)));
            [&R] (B,(C,(A,D)));
            [&R] (B,(D,(C,A)));
            [&R] (C,(A,(B,D)));
            [&R] (C,(B,(A,D)));
            [&R] (C,(D,(B,A)));
            [&R] (D,(A,(B,C)));
            [&R] (D,(B,(A,C)));
            [&R] (D,(C,(B,A)));
            [&R] ((A,B),(C,D));
            [&R] ((A,C),(B,D));
            [&R] ((A,D),(C,B));
            """, "NEWICK", taxon_namespace=taxa)

        # expected results, for each gene tree / species tree pairing, with
        # cycling through species trees for each gene tree
        expected_deep_coalescences = [ 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 1, 2, 2,
                                            2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 2, 0, 2,
                                            2, 1, 2, 3, 3, 3, 0, 1, 1, 3, 3, 3, 2, 1, 2 ]
        assert len(expected_deep_coalescences) == len(gene_trees) * len(species_trees)

        for t in gene_trees + species_trees:
            t.update_bipartitions()
        idx = 0
        _LOG.info("Species\t\tGene\t\tDC\t\tExp.DC\t\tDiff")
        for gt in gene_trees:
            gt.update_bipartitions()
            for st in species_trees:
                st.update_bipartitions()
                dc = reconcile.reconciliation_discordance(gt, st)
                _LOG.info("%s\t\t%s\t\t%s\t\t%s\t\t%s"
                    % (st._as_newick_string(),
                       gt._as_newick_string(),
                       dc,
                       expected_deep_coalescences[idx],
                       dc - expected_deep_coalescences[idx]))
                assert dc == expected_deep_coalescences[idx]
                idx += 1

    def testGroupedDeepCoalCounting(self):
            src_trees = [
                        "((a1,a2)x,b1)y;",
                        "((a1, (a2, a3), b1), (b2,(b3,b4)));",
                        "(((((a1, a2),a3), b1), b2), (b3, ((b4,b5),b6)));",
                        "((b1, (b2, b3), a1), (a2,(a3, a4)));",
                        "(((((b1, b2),b3), a1), a2), (a3, ((a4,a5),a6)));",
                        "((a1,a2),(b1,b2),(c1,c2));",
                        "((a1,a2),(b1,b2,c3),(c1,c2));",
                        "(((a1,a2),(b1,b2),c1),c2);",
                        ]
            scores = [ 0, 1, 2, 1, 2, 0, 1, 1 ]
            for src_tree, expected in zip(src_trees, scores):
                tree = dendropy.Tree.get_from_string(src_tree, "NEWICK")
                groups = dendropy.TaxonNamespacePartition(tree.taxon_namespace,
                    membership_fn=lambda x: x.label[0])
                dc = reconcile.monophyletic_partition_discordance(tree, groups)
            assert dc == expected, \
                "deep coalescences by groups: expecting %d, but found %d" % (expected, dc)

if __name__ == "__main__":
    unittest.main()

