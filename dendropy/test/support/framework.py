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
Extensions of the unittest framework.
"""

import unittest
import dendropy
from dendropy.utility import messaging

# defining this here means that unittest will exclude all lines from this
# module in the traceback report when an assertion is raised, allowing
# for the point of failure to be the point where the assertion statement
# was made, rather than the point where an exception was raised because
# the assertion as false
#__unittest = True

class ExtendedTestCase(unittest.TestCase):
    """
    Extends unittest.TestCase with various new assertion tests.
    """

    def _get_logger(self):
        if not hasattr(self, "_logger") or self._logger is None:
            self._logger = messaging.get_logger(self.__class__.__name__)
        return self._logger

    def _set_logger(self, logger):
        self._logger = logger

    logger = property(_get_logger, _set_logger)

    def assertIsSame(self, obj1, obj2, message=None):
        if message is None:
            message = "%s (%d) is not %s (%d)" % (obj1, id(obj1), obj2, id(obj2))
        self.assertTrue(obj1 is obj2, message)

    def assertIsNotSame(self, obj1, obj2, message=None):
        if message is None:
            message = "%s (%d) is %s (%d)" % (obj1, id(obj1), obj2, id(obj2))
        self.assertTrue(obj1 is not obj2, message)

    def assertIsContainedIn(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is not in: %s" % (obj1, obj2)
        self.assertTrue(obj1 in obj2, message)

    def assertIsNotContainedIn(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is in: %s" % (obj1, obj2)
        self.assertTrue(obj1 not in obj2, message)

class DataObjectVerificationTestCase(ExtendedTestCase):
    """
    Extends ExtendedTestCase with tests for data object comparisons.
    """

    def assertDistinctButEqual(self, data_object1, data_object2, distinct_taxa=True, equal_oids=None):
        """
        Verifies that two DendroPy phylogenetic data objects (Tree, TreeList,
        CharArray, Dataset etc.) are independent objects, but equal. That is,
        if `data_object1` and `data_object` are the same object, then the
        distinction criterion is failed. If `distinct_taxa` is False, then the
        distinct/independence is criterion is enforced down to Taxon and
        TaxonSet objects. If `distinct_taxa` is True, then the test fails if
        the Taxon and TaxonSet objects are NOT distinct. Same holds for
        `equal_oids`, except if `equal_oids` is None then oid's are not
        checked.

        Equality is assessed on the basis of identical representation of the
        objects in domain space.
        """
        if type(data_object1) != type(data_object2):
            raise TypeError("Objects to be compared must be of the same type, but was given %s and %s objects" \
                % (type(data_object1), type(data_object2)))
        if isinstance(data_object1, dendropy.Tree):
            self.assertDistinctButEqualTrees(data_object1, data_object2, distinct_taxa=distinct_taxa, equal_oids=equal_oids)
        else:
            raise TypeError("Unsupported type for comparison: %s" % type(data_object1))

    def assertDistinctButEqualTrees(self, tree1, tree2, distinct_taxa=True, equal_oids=None):
        """
        `tree1` and `tree2` must be different objects, but with the same
        structure, returning nodes with the same labels and associated with
        taxa with the same labels in the same order when traversed in the same
        way. If `distinct_taxa==False` then corresponding nodes on the two
        trees must have the same `Taxon` object (or None). If
        `distinct_taxa==True` then corresponding nodes must have different
        `Taxon` objects or None, but, if they have `Taxon` objects, then the
        labels of the `Taxon` objects must be the same even though the objects
        themselves are different.
        """
        self.logger.debug(tree1.as_newick_str())
        tree1.debug_check_tree(logger=self.logger)
        self.logger.debug(tree2.as_newick_str())
        tree2.debug_check_tree(logger=self.logger)

        self.assertIsNotSame(tree1, tree2)
        if distinct_taxa:
            self.assertIsNotSame(tree1.taxon_set, tree2.taxon_set)
        if equal_oids is True:
            self.assertEqual(tree1.oid, tree2.oid)
        elif equal_oids is False:
            self.assertNotEqual(tree1.oid, tree2.oid)

        tree1_nodes = [nd for nd in tree1.postorder_node_iter()]
        tree2_nodes = [nd for nd in tree2.postorder_node_iter()]
        self.assertEqual(len(tree1_nodes), len(tree2_nodes))
        for nd_idx, node1 in enumerate(tree1_nodes):
            node2 = tree2_nodes[nd_idx]
            if node1.taxon is not None:
                self.assert_(node2.taxon is not None)
                if distinct_taxa:
                    self.assertIsNotSame(node1.taxon, node2.taxon)
                else:
                    self.assertIsSame(node1.taxon, node2.taxon)
                if equal_oids is True:
                    self.assertEqual(node1.oid, node2.oid)
                elif equal_oids is False:
                    self.assertNotEqual(node1.oid, node2.oid)
                self.assertEqual(node1.taxon.label, node2.taxon.label)
                self.assertIsContainedIn(node1.taxon, tree1.taxon_set)
                self.assertIsContainedIn(node2.taxon, tree2.taxon_set)
            else:
                self.assertIsSame(node2.taxon, None)
            if node1.edge.length is not None:
                self.assertIsNotSame(node2.edge.length, None)
                self.assertAlmostEqual(node1.edge.length, node2.edge.length, 3)
            else:
                self.assertIsSame(node2.edge.length, None)
            self.assertEqual(len(node1.child_nodes()), len(node2.child_nodes()))
        tree1_edges = [edge for edge in tree1.postorder_edge_iter()]
        tree2_edges = [edge for edge in tree2.postorder_edge_iter()]
        self.assertEqual(len(tree1_edges), len(tree2_edges))
        for edge_idx, edge1 in enumerate(tree1_edges):
            edge2 = tree2_edges[edge_idx]
            self.assertIsNotSame(edge1, edge2)
            if edge1.length is None:
                self.assertIsSame(edge1.length, None)
            else:
                self.assertEqual(edge1.length, edge2.length)
            if equal_oids is True:
                self.assertEqual(edge1.oid, edge2.oid)
            elif equal_oids is False:
                self.assertNotEqual(edge1.oid, edge2.oid)

#    def compare_datasets(ds1, ds2, tester, distinct_taxa=True, equal_oids=False):
#        self.logger.info("Comparing dataset taxon sets ...")
#        compare_dataset_taxon_sets(ds1, ds2, tester, distinct_taxa, equal_oids)
#        self.logger.info("Comparing dataset tree lists ...")
#        compare_dataset_tree_lists(ds1, ds2, tester, distinct_taxa, equal_oids)
#        self.logger.info("Comparing dataset character arrays ...")
#        compare_dataset_char_arrays(ds1, ds2, tester, distinct_taxa, equal_oids)
#
#    def compare_dataset_taxon_sets(ds1, ds2, tester, distinct_taxa=True, equal_oids=False):
#        self.assertEqual(len(ds1.taxon_sets), len(ds2.taxon_sets))
#        if distinct_taxa:
#            self.assertTrue(ds1.taxon_sets is not ds2.taxon_sets)
#        for ts_idx, ts1 in enumerate(ds1.taxon_sets):
#            ts2 = ds2.taxon_sets[ts_idx]
#            self.logger.info("Comparing taxa of taxon set %d: %d taxa vs. %d taxa" \
#                % (ts_idx, len(ts1), len(ts2)))
#            compare_individual_taxon_sets(ts1, ts2, tester, distinct_taxa, equal_oids)
#
#    def compare_individual_taxon_sets(ts1, ts2, tester, distinct_taxa=True, equal_oids=False):
#        if distinct_taxa:
#            self.assertTrue(ts1 is not ts2)
#        self.assertEqual(len(ts1), len(ts2))
#        if equal_oids:
#            self.assertNotEqual(ts1.oid, ts2.oid)
#        self.assertEqual(ts1.label, ts2.label)
#        for taxon_idx, taxon1 in enumerate(ts1):
#            self.logger.debug("Taxon %d: '%s' == '%s'" % (taxon_idx, taxon1.label, ts2[taxon_idx].label))
#            taxon2 = ts2[taxon_idx]
#            if distinct_taxa:
#                self.assertTrue(taxon1 is not taxon2)
#            self.assertEqual(taxon1.label, taxon2.label)
#            if equal_oids:
#                self.assertNotEqual(taxon1.oid, taxon2.oid)
#
#    def compare_dataset_tree_lists(ds1, ds2, tester, distinct_taxa=True, equal_oids=False):
#        self.assertTrue(ds1.tree_lists is not ds2.tree_lists)
#        self.assertEqual(len(ds1.tree_lists), len(ds2.tree_lists))
#        for tree_list_idx, tree_list1 in enumerate(ds1.tree_lists):
#            tree_list2 = ds2.tree_lists[tree_list_idx]
#            if distinct_taxa:
#                self.assertTrue(tree_list1.taxon_set is not tree_list2.taxon_set)
#                self.assertTrue(tree_list1.taxon_set in ds1.taxon_sets)
#                self.assertTrue(tree_list2.taxon_set in ds2.taxon_sets)
#                self.assertTrue(tree_list1.taxon_set not in ds2.taxon_sets)
#                self.assertTrue(tree_list2.taxon_set not in ds1.taxon_sets)
#            compare_individual_tree_lists(tree_list1, tree_list2, tester, distinct_taxa, equal_oids)
#
#    def compare_individual_tree_lists(tree_list1, tree_list2, tester, distinct_taxa=True, equal_oids=False):
#        self.assertTrue(tree_list1 is not tree_list2)
#        self.assertEqual(len(tree_list1), len(tree_list2))
#        if distinct_taxa:
#            self.assertTrue(tree_list1.taxon_set is not tree_list2.taxon_set)
#        compare_individual_taxon_sets(tree_list1.taxon_set, tree_list2.taxon_set, tester, distinct_taxa, equal_oids)
#        self.assertEqual(tree_list1.label, tree_list2.label)
#        if equal_oids:
#            self.assertNotEqual(tree_list1.oid, tree_list2.oid)
#        for tree_idx, tree1 in enumerate(tree_list1):
#            tree2 = tree_list2[tree_idx]
#            compare_individual_trees(tree1, tree2, tester, distinct_taxa, equal_oids)
#
#    def compare_individual_trees(tree1, tree2, tester, distinct_taxa=True, equal_oids=False):
#            self.logger.debug(tree1.as_newick_str())
#            tree1.debug_check_tree(logger=self.logger)
#            self.logger.debug(tree2.as_newick_str())
#            tree2.debug_check_tree(logger=self.logger)
#
#            self.assertTrue(tree1 is not tree2)
#            if distinct_taxa:
#                self.assertTrue(tree1.taxon_set is not tree2.taxon_set)
#            self.assertTrue(tree1.taxon_set is tree_list1.taxon_set)
#            self.assertTrue(tree2.taxon_set is tree_list2.taxon_set)
#
#            tree1_nodes = [nd for nd in tree1.postorder_node_iter()]
#            tree2_nodes = [nd for nd in tree2.postorder_node_iter()]
#            self.assertEqual(len(tree1_nodes), len(tree2_nodes))
#            for nd_idx, node1 in enumerate(tree1_nodes):
#                node2 = tree2_nodes[nd_idx]
#                if node1.taxon is not None:
#                    self.assert_(node2.taxon is not None)
#                    if distinct_taxa:
#                        self.assertTrue(node1.taxon is not node2.taxon)
#                    else:
#                        self.assertTrue(node1.taxon is node2.taxon)
#                    if equal_oids:
#                        self.assertNotEqual(node1.oid, node2.oid)
#                    self.assertEqual(node1.taxon.label, node2.taxon.label)
#                    self.assertTrue(node1.taxon in tree1.taxon_set)
#                    self.assertTrue(node2.taxon in tree2.taxon_set)
#                else:
#                    self.assert_(node2.taxon is None)
#                if node1.edge.length is not None:
#                    self.assert_(node2.edge.length is not None)
#                    self.assertAlmostEqual(node1.edge.length, node2.edge.length, 3)
#                else:
#                    self.assert_(node2.edge.length is None)
#                self.assertEqual(len(node1.child_nodes()), len(node2.child_nodes()))
#
#    def compare_dataset_char_arrays(ds1, ds2, tester, distinct_taxa=True, equal_oids=False):
#        self.assertEqual(len(ds1.char_arrays), len(ds2.char_arrays))
#        for char_array_idx, char_array1 in enumerate(ds1.char_arrays):
#            char_array2 = ds2.char_arrays[char_array_idx]
#            compare_individual_char_arrays(char_array1, char_array2, tester, distinct_taxa, equal_oids)
#
#    def compare_individual_char_arrays(char_array1, char_array2, tester, distinct_taxa=True, equal_oids=False):
#        self.assertEqual(len(char_array1), len(char_array2))
#        self.assertEqual(len(char_array1.taxon_set), len(char_array2.taxon_set))
#        for taxon_idx, taxon1 in enumerate(char_array1.taxon_set):
#            self.assertEqual(char_array1.taxon_set[taxon_idx].label,
#                    char_array2.taxon_set[taxon_idx].label)
#            seq1 = char_array1[taxon_idx]
#            seq2 = char_array2[taxon_idx]
#            self.assertEqual(len(seq1), len(seq2))
#            for cell_idx, cell1 in enumerate(seq1):
#                cell2 = seq2[cell_idx]
#                state1 = cell1.value
#                state2 = cell2.value
#                self.assertEqual(state1.symbol, state2.symbol)
#                self.assertEqual(state1.token, state2.token)
#                self.assertEqual(state1.multistate, state2.multistate)
#                self.assertEqual(state1.fundamental_symbols, state2.fundamental_symbols)

