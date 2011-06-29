#! /usr/bin/env python

############################################################################
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Tests reconciliation calculations.
"""

import os
import unittest
import dendropy
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

from dendropy import reconcile

class ContainingTreeDeepCoalescenceSmall(unittest.TestCase):

    def setUp(self):
        dataset = dendropy.DataSet.get_from_path(pathmap.tree_source_path(filename="deepcoal1.nex"), "nexus")
        self.species_tree = dataset.get_tree_list(label="ContainingTree")[0]
        self.gene_trees = dataset.get_tree_list(label="EmbeddedTrees")
        self.species_tree.taxon_set.lock()
        self.gene_taxon_to_population_taxon_map = dendropy.TaxonSetMapping(
                domain_taxon_set=self.gene_trees.taxon_set,
                range_taxon_set=self.species_tree.taxon_set,
                mapping_func=lambda t: self.species_tree.taxon_set.require_taxon(label=t.label[0].upper()))
        self.expected_under_original_brlens = [4, 6, 4, 2, 4, 3, 3, 4, 5, 4]

    def testFixedEdgesDeepCoalCount(self):
        results = []
        for idx, gt in enumerate(self.gene_trees):
            ct = reconcile.ContainingTree(containing_tree=self.species_tree,
                    contained_taxon_set=self.gene_trees.taxon_set,
                    contained_to_containing_taxon_map=self.gene_taxon_to_population_taxon_map,
                    contained_trees=[gt],
                    fit_containing_edge_lengths=False,
                    )
            dc = ct.num_deep_coalescences()
            results.append(dc)
            mesqf = pathmap.named_output_stream("ContainingTreeDeepCoalescence_Small_FixedEdges_t%02d_dc%02d.nex" % (idx+1, dc), False)
            ct.write_as_mesquite(mesqf)
        self.assertEqual(results, self.expected_under_original_brlens)

    def testFittedEdgesDeepCoalCount(self):
        for idx, gt in enumerate(self.gene_trees):
            ct = reconcile.ContainingTree(containing_tree=self.species_tree,
                    contained_taxon_set=self.gene_trees.taxon_set,
                    contained_to_containing_taxon_map=self.gene_taxon_to_population_taxon_map,
                    contained_trees=[gt],
                    fit_containing_edge_lengths=True,
                    )
            dc = ct.num_deep_coalescences()
            mesqf = pathmap.named_output_stream("ContainingTreeDeepCoalescence_Small_FittedEdges_t%02d_dc%02d.nex" % (idx+1, dc), False)
            ct.write_as_mesquite(mesqf)

class DeepCoalTest(unittest.TestCase):

    def testFittedDeepCoalCounting(self):

        taxa = dendropy.TaxonSet()

        gene_trees = dendropy.TreeList.get_from_string("""
            [&R] (A,(B,(C,D))); [&R] ((A,C),(B,D)); [&R] (C,(A,(B,D)));
            """, "newick", taxon_set=taxa)

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
            """, "NEWICK", taxon_set=taxa)

        # expected results, for each gene tree / species tree pairing, with
        # cycling through species trees for each gene tree
        expected_deep_coalescences = [ 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 1, 2, 2,
                                            2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 2, 0, 2,
                                            2, 1, 2, 3, 3, 3, 0, 1, 1, 3, 3, 3, 2, 1, 2 ]
        assert len(expected_deep_coalescences) == len(gene_trees) * len(species_trees)

        for t in gene_trees + species_trees:
            t.update_splits()
        idx = 0
        _LOG.info("Species\t\tGene\t\tDC\t\tExp.DC\t\tDiff")
        for gt in gene_trees:
            for st in species_trees:
                dc = reconcile.reconciliation_discordance(gt, st)
                _LOG.info("%s\t\t%s\t\t%s\t\t%s\t\t%s"
                    % (st.compose_newick(),
                       gt.compose_newick(),
                       dc,
                       expected_deep_coalescences[idx],
                       dc - expected_deep_coalescences[idx]))
                assert dc == expected_deep_coalescences[idx]
                idx += 1

    def testGroupedDeepCoalCounting(self):
        src_trees = { "((a1,a2)x,b1)y;" : 0,
                      "((a1, (a2, a3), b1), (b2,(b3,b4)))" : 1,
                      "(((((a1, a2),a3), b1), b2), (b3, ((b4,b5),b6)))" : 2,
                      "((b1, (b2, b3), a1), (a2,(a3, a4)))" : 1,
                      "(((((b1, b2),b3), a1), a2), (a3, ((a4,a5),a6)))" : 2,
                      "((a1,a2),(b1,b2),(c1,c2))" : 0,
                      "((a1,a2),(b1,b2,c3),(c1,c2))" : 1,
                      "(((a1,a2),(b1,b2),c1),c2)" : 1
                    }
        for src_tree, expected in src_trees.items():
            tree = dendropy.Tree.get_from_string(src_tree, "NEWICK")
            groups = dendropy.TaxonSetPartition(tree.taxon_set,
                membership_func=lambda x: x.label[0])
#            groups = [[],[]]
#            for taxon in tree.taxon_set:
#                if taxon.label.startswith('a'):
#                    groups[0].append(taxon)
#                elif taxon.label.startswith('b'):
#                    groups[1].append(taxon)
#                elif taxon.label.startswith('c'):
#                    if len(groups) < 3:
#                        groups.append([])
#                    groups[2].append(taxon)
            dc = reconcile.monophyletic_partition_discordance(tree, groups)
            assert dc == expected, \
                "deep coalescences by groups: expecting %d, but found %d" % (expected, dc)

if __name__ == "__main__":
    unittest.main()

