#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests of parsimony scoring.
"""

import unittest
import dendropy
from dendropy.calculate import treescore
from dendropy.test.support import pathmap

class ParsimonyScoringTest(unittest.TestCase):

    def test_pscores_with_gaps_as_new_state(self):
        # #NEXUS
        # begin paup;
        #     set warnroot = no;
        #     exe apternodus.chars.nexus;
        #     gett file = apternodus.tre;
        #     set criterion = parsimony;
        #     pset gap = newstate;
        #     pscore;
        # end;
        expected_scores = [396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 713, 715, 723, 733, 672, 719, 734, 709, 695, 686]
        self.verify_pscores(
                trees_fname="apternodus.tre",
                chars_fname="apternodus.chars.nexus",
                matrix_type=dendropy.StandardCharacterMatrix,
                gaps_as_missing=False,
                expected_scores=expected_scores)

    def test_pscores_with_gaps_as_missing(self):
        # #NEXUS
        # begin paup;
        #     set warnroot = no;
        #     exe apternodus.chars.nexus;
        #     gett file = apternodus.tre;
        #     set criterion = parsimony;
        #     pset gap = missing;
        #     pscore;
        # end;
        expected_scores = [ 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 671, 670, 678, 687, 633, 675, 689, 668, 652, 644]
        self.verify_pscores(
                trees_fname="apternodus.tre",
                chars_fname="apternodus.chars.nexus",
                matrix_type=dendropy.StandardCharacterMatrix,
                gaps_as_missing=True,
                expected_scores=expected_scores)


    def verify_pscores(self,
            trees_fname,
            chars_fname,
            matrix_type,
            gaps_as_missing,
            expected_scores):
        taxon_namespace = dendropy.TaxonNamespace()
        chars = matrix_type.get(
                path=pathmap.char_source_path(chars_fname),
                schema="nexus",
                taxon_namespace=taxon_namespace)
        trees = dendropy.TreeList.get(
                path=pathmap.tree_source_path(trees_fname),
                schema="nexus",
                taxon_namespace=taxon_namespace)
        self.assertEqual(len(expected_scores), len(trees))
        for n, tree in enumerate(trees):
            score_by_character_list = []
            pscore = treescore.parsimony_score(
                    tree,
                    chars,
                    gaps_as_missing=gaps_as_missing,
                    score_by_character_list=score_by_character_list)
            self.assertEqual(expected_scores[n], pscore)

if __name__ == "__main__":
    unittest.main()


