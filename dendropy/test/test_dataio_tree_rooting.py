#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Tests rooting interpretation.
"""

import unittest
import dendropy
from cStringIO import StringIO
from dendropy.dataio.nexustokenizer import RootingInterpreter
from dendropy.dataio import tree_source_iter

class RootingIntepreterInstantiationTest(unittest.TestCase):

    def testDefault(self):
        ri = RootingInterpreter()
        self.assertEqual(ri._as_rooted, None)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsRootedTrue(self):
        ri = RootingInterpreter(as_rooted=True)
        self.assertEqual(ri._as_rooted, True)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsRootedFalse(self):
        ri = RootingInterpreter(as_rooted=False)
        self.assertEqual(ri._as_rooted, False)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsUnrootedTrue(self):
        ri = RootingInterpreter(as_unrooted=True)
        self.assertEqual(ri._as_rooted, False)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsUnrootedFalse(self):
        ri = RootingInterpreter(as_unrooted=False)
        self.assertEqual(ri._as_rooted, True)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsDefaultRootedTrue(self):
        ri = RootingInterpreter(default_as_rooted=True)
        self.assertEqual(ri._as_rooted, None)
        self.assertEqual(ri._default_as_rooted, True)

    def testAsDefaultRootedFalse(self):
        ri = RootingInterpreter(default_as_rooted=False)
        self.assertEqual(ri._as_rooted, None)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsDefaultUnrootedTrue(self):
        ri = RootingInterpreter(default_as_unrooted=True)
        self.assertEqual(ri._as_rooted, None)
        self.assertEqual(ri._default_as_rooted, False)

    def testAsDefaultUnrootedFalse(self):
        ri = RootingInterpreter(default_as_unrooted=False)
        self.assertEqual(ri._as_rooted, None)
        self.assertEqual(ri._default_as_rooted, True)

class TreeRootingIntepretationTest(unittest.TestCase):

    def setUp(self):
        self.trees_all_rooted = [
            "[&R] ((A,B),(C,(D,E)))",
            "[&R] ((A,B),(C,(D,E)))",
            "[&R] ((A,B),(C,(D,E)))",
            "[&R] ((A,B),(C,(D,E)))",
            "[&R] ((A,B),(C,(D,E)))",
        ]
        self.trees_all_unrooted = [
            "[&U] ((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
        ]
        self.trees_all_unspecified = [
            "((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
        ]
        self.trees_mixed = [
            "[&R] ((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "[&R] ((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
            "[&R] ((A,B),(C,(D,E)))",
            "[&U] ((A,B),(C,(D,E)))",
            "((A,B),(C,(D,E)))",
        ]
        self.trees_mixed_is_rooted = [True, None, True, False, None, False, None, None, True, False, None]

    def get_as_nexus(self, trees):
        nx = []
        nx.append("#NEXUS")
        nx.append("begin taxa;")
        nx.append("dimensions ntax=5;")
        nx.append("end;")
        nx.append("begin trees;")
        for i, t in enumerate(trees):
            nx.append("tree %d = %s;" % (i, t))
        nx.append("end;")
        return "\n".join(nx)

    def get_as_newick(self, trees):
        return "\n".join([(t + ";") for t in trees])

    def check_all_rooted(self, trees):
        v1 = [tree.is_rooted for tree in trees]
        self.assertEqual(v1, [True] * len(trees))

    def check_all_unrooted(self, trees):
        v1 = [tree.is_rooted for tree in trees]
        self.assertEqual(v1, [False] * len(trees))

    def check_mixed_with_default(self, trees, default):
        v1 = [tree.is_rooted for tree in trees]
        for i, r in enumerate(v1):
            if self.trees_mixed_is_rooted[i] is None:
                self.assertEqual(v1[i], default)
            else:
                self.assertEqual(v1[i], self.trees_mixed_is_rooted[i])

    def iterated_trees(self, tree_source_string, schema, **kwargs):
        ti = tree_source_iter(StringIO(tree_source_string), schema, **kwargs)
        return [t for t in ti]

    ###########################################################################
    ############################ NEXUS Reader  ################################
    ###########################################################################

    ###########################################################################
    ## NEXUS / as_rooted = True

    def testNexusReaderAllRootedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnrootedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnspecifiedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNexusReaderMixedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_rooted(trees)

    ###########################################################################
    ## NEXUS / as_rooted = False

    def testNexusReaderAllRootedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnrootedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnspecifiedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNexusReaderMixedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_unrooted(trees)

    ###########################################################################
    ## NEXUS / default_as_rooted = True

    def testNexusReaderAllRootedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnrootedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnspecifiedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNexusReaderMixedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, True)

    ###########################################################################
    ## NEXUS / default_as_rooted = False

    def testNexusReaderAllRootedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnrootedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnspecifiedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNexusReaderMixedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, False)

    ###########################################################################
    ############################ NEXUS Iterator ###############################
    ###########################################################################

    ###########################################################################
    ## NEXUS / as_rooted = True

    def testNexusReaderAllRootedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnrootedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnspecifiedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNexusReaderMixedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_rooted(trees)

    ###########################################################################
    ## NEXUS / as_rooted = False

    def testNexusReaderAllRootedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnrootedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnspecifiedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNexusReaderMixedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_unrooted(trees)

    ###########################################################################
    ## NEXUS / default_as_rooted = True

    def testNexusReaderAllRootedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnrootedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnspecifiedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNexusReaderMixedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, True)

    ###########################################################################
    ## NEXUS / default_as_rooted = False

    def testNexusReaderAllRootedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_rooted),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNexusReaderAllUnrootedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unrooted),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNexusReaderAllUnspecifiedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_all_unspecified),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNexusReaderMixedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_nexus(self.trees_mixed),
                "nexus",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, False)

    ###########################################################################
    ############################ NEWICK Reader  ################################
    ###########################################################################

    ###########################################################################
    ## NEWICK / as_rooted = True

    def testNewickReaderAllRootedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnrootedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnspecifiedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNewickReaderMixedAsRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_mixed),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_rooted(trees)

    ###########################################################################
    ## NEWICK / as_rooted = False

    def testNewickReaderAllRootedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnrootedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnspecifiedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNewickReaderMixedAsUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_mixed),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_unrooted(trees)

    ###########################################################################
    ## NEWICK / default_as_rooted = True

    def testNewickReaderAllRootedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnrootedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnspecifiedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNewickReaderMixedAsDefaultRooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_mixed),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, True)

    ###########################################################################
    ## NEWICK / default_as_rooted = False

    def testNewickReaderAllRootedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnrootedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnspecifiedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNewickReaderMixedAsDefaultUnrooted(self):
        trees = dendropy.TreeList.get_from_string(
                self.get_as_newick(self.trees_mixed),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, False)

    ###########################################################################
    ############################ NEWICK Iterator ###############################
    ###########################################################################

    ###########################################################################
    ## NEWICK / as_rooted = True

    def testNewickReaderAllRootedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnrootedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnspecifiedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNewickReaderMixedAsRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_mixed),
                "newick",
                as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_rooted(trees)

    ###########################################################################
    ## NEWICK / as_rooted = False

    def testNewickReaderAllRootedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnrootedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnspecifiedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNewickReaderMixedAsUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_mixed),
                "newick",
                as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_all_unrooted(trees)

    ###########################################################################
    ## NEWICK / default_as_rooted = True

    def testNewickReaderAllRootedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnrootedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnspecifiedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_rooted(trees)

    def testNewickReaderMixedAsDefaultRooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_mixed),
                "newick",
                default_as_rooted=True)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, True)

    ###########################################################################
    ## NEWICK / default_as_rooted = False

    def testNewickReaderAllRootedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_rooted),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_rooted))
        self.check_all_rooted(trees)

    def testNewickReaderAllUnrootedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unrooted),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unrooted))
        self.check_all_unrooted(trees)

    def testNewickReaderAllUnspecifiedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_all_unspecified),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_all_unspecified))
        self.check_all_unrooted(trees)

    def testNewickReaderMixedAsDefaultUnrooted(self):
        trees = self.iterated_trees(
                self.get_as_newick(self.trees_mixed),
                "newick",
                default_as_rooted=False)
        self.assertEqual(len(trees), len(self.trees_mixed))
        self.check_mixed_with_default(trees, False)

if __name__ == "__main__":
    unittest.main()
