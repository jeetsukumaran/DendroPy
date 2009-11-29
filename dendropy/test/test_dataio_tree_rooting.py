#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Tests rooting interpretation.
"""

import unittest
import dendropy
from dendropy.dataio.nexustokenizer import RootingInterpreter

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





if __name__ == "__main__":
    unittest.main()
