############################################################################
##  test_coal.py
##
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
Tests coalescence calculations.
"""

import unittest
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
import dendropy
from dendropy import coalescent

class CalcIntervalsTest(unittest.TestCase):

    def testSimple1(self):

        t = dendropy.Tree.get_from_string("((((a:1, b:1):1, c:2):1, d:3, e:3):2, (f:4, g:4):1)", "newick")
        i1 = coalescent.node_waiting_time_pairs(t)
        assert [x[1] for x in i1] == [1.0, 1.0, 1.0, 1.0, 1.0]
        i2 = coalescent.extract_coalescent_frames(t)
        assert i2 == {7: 1.0, 6:1.0, 5:1.0, 3:1.0, 2:1.0}
        check = coalescent.log_probability_of_coalescent_tree(t, 10)

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
                dc = coalescent.num_deep_coalescences_with_fitted_tree(gt, st)
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
            groups = [[],[]]
            for taxon in tree.taxon_set:
                if taxon.label.startswith('a'):
                    groups[0].append(taxon)
                elif taxon.label.startswith('b'):
                    groups[1].append(taxon)
                elif taxon.label.startswith('c'):
                    if len(groups) < 3:
                        groups.append([])
                    groups[2].append(taxon)
            dc = coalescent.num_deep_coalescences_with_grouping(tree, groups)
            assert dc == expected, \
                "deep coalescences by groups: expecting %d, but found %d" % (expected, dc)

#    if coalescent.de_hoon_statistics:
#        def testKLDiv(self):
#            from dendropy import treegen
#            taxa_block = treegen.random_taxa_block(100)
#            ctrees = [treegen.pure_kingman(taxa_block, 20000) for i in range(10)]
#            _LOG.info("KL divergence from Kingman coalescent of trees generated under pure Kingman process: %s" % coalescent.kl_divergence_coalescent_trees(ctrees, 20000))
#            ytrees = [treegen.uniform_pure_birth(taxa_block) for i in range(10)]
#            _LOG.info("KL divergence from Kingman coalescent of trees generated under Yule process: %s" % coalescent.kl_divergence_coalescent_trees(ytrees, 20000))

if __name__ == "__main__":
    unittest.main()

