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
Extensions of the unittest framework for data comparison and testing.
"""

import unittest
import dendropy
from dendropy.test.support import extendedtest

class DataObjectVerificationTestCase(extendedtest.ExtendedTestCase):
    """
    Extends ExtendedTestCase with tests for data object comparisons.
    """

    def assertDistinctButEqual(self, data_object1, data_object2, **kwargs):
        """
        Verifies that two DendroPy phylogenetic data objects (Tree, TreeList,
        CharArray, DataSet etc.) are independent objects, but equal. That is,
        if `data_object1` and `data_object` are the same object, then the
        distinction criterion is failed. If `distinct_taxa` is True, then the
        distinct/independence is criterion is enforced down to Taxon and
        TaxonSet objects. If `distinct_taxa` is False, then the test fails if
        the Taxon and TaxonSet objects are NOT distinct. Same holds for
        `equal_oids`, except if `equal_oids` is None then oid's are not
        checked.

        Equality is assessed on the basis of identical representation of the
        objects in domain space.
        """
        distinct_taxa = kwargs.get("distinct_taxa", True)
        equal_oids = kwargs.get("equal_oids", None)
        if type(data_object1) != type(data_object2):
            raise ValueError("Objects to be compared must be of the same type, but was given %s and %s objects" \
                % (type(data_object1), type(data_object2)))
        if isinstance(data_object1, dendropy.Taxon):
            self.assertDistinctButEqualTaxon(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.TaxonSet):
            self.assertDistinctButEqualTaxonSet(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.Tree):
            self.assertDistinctButEqualTree(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.TreeList):
            self.assertDistinctButEqualTreeList(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.DiscreteCharacterArray):
            self.assertDistinctButEqualDiscreteCharArray(data_object1, data_object2, **kwargs)
        else:
            raise ValueError("Unsupported type for comparison: %s" % type(data_object1))

    def assertDistinctButEqualTaxon(self, taxon1, taxon2, **kwargs):
        equal_oids = kwargs.get("equal_oids", None)
        ignore_underscore_substitution = kwargs.get("ignore_underscore_substitution", False)
        self.assertNotSame(taxon1, taxon2)
        if taxon1.label is None:
            self.assertSame(taxon2.label, None)
        elif ignore_underscore_substitution:
            self.assertEqual(taxon1.label.replace(" ", "_"), taxon2.label.replace(" ", "_"))
        else:
            self.assertEqual(taxon1.label, taxon2.label)
        if equal_oids is True:
            self.assertEqual(taxon1.oid, taxon2.oid)
        elif equal_oids is False:
            self.assertNotEqual(taxon1.oid, taxon2.oid)

    def assertDistinctButEqualTaxonSet(self, taxon_set1, taxon_set2, **kwargs):
        equal_oids = kwargs.get("equal_oids", None)
        ignore_underscore_substitution = kwargs.get("ignore_underscore_substitution", False)
        self.assertNotSame(taxon_set1, taxon_set2)
        self.assertEqual(len(taxon_set1), len(taxon_set2))
        for tidx, taxon1 in enumerate(taxon_set1):
            taxon2 = taxon_set2[tidx]
            self.assertDistinctButEqualTaxon(taxon1, taxon2, **kwargs)

    def assertDistinctButEqualTreeList(self, tree_list1, tree_list2, **kwargs):
        """
        `tree_list1` and `tree_list2` must be distinct but equivalent objects.
        """
        distinct_taxa = kwargs.get("distinct_taxa", True)
        equal_oids = kwargs.get("equal_oids", None)
        self.assertTrue(tree_list1 is not tree_list2)
        self.assertEqual(len(tree_list1), len(tree_list2))
        if distinct_taxa:
            self.assertNotSame(tree_list1.taxon_set, tree_list2.taxon_set)
            self.assertDistinctButEqualTaxonSet(tree_list1.taxon_set, tree_list2.taxon_set, **kwargs)
        else:
            self.assertSame(tree_list1.taxon_set, tree_list2.taxon_set)
        if distinct_taxa:
            self.assertDistinctButEqualTaxonSet(tree_list1.taxon_set, tree_list2.taxon_set, **kwargs)
        else:
            self.assertSame(tree_list1.taxon_set, tree_list2.taxon_set)
        self.assertEqual(tree_list1.label, tree_list2.label)
        if equal_oids is True:
            self.assertEqual(tree_list1.oid, tree_list1.oid)
        elif equal_oids is False:
            self.assertNotEqual(tree_list1.oid, tree_list2.oid)
        distinct_trees = kwargs.get("distinct_trees", True)
        for tree_idx, tree1 in enumerate(tree_list1):
            tree2 = tree_list2[tree_idx]
            if distinct_trees:
                self.assertSame(tree1.taxon_set, tree_list1.taxon_set)
                self.assertSame(tree2.taxon_set, tree_list2.taxon_set)
                self.assertDistinctButEqualTree(tree1, tree2, **kwargs)
            else:
                self.assertSame(tree1, tree2)

    def assertDistinctButEqualTree(self, tree1, tree2, **kwargs):
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
        distinct_taxa = kwargs.get("distinct_taxa", True)
        equal_oids = kwargs.get("equal_oids", None)
        self.logger.debug(tree1.as_newick_str())
        tree1.debug_check_tree(logger=self.logger)
        self.logger.debug(tree2.as_newick_str())
        tree2.debug_check_tree(logger=self.logger)

        self.assertNotSame(tree1, tree2)
        if distinct_taxa:
            self.assertNotSame(tree1.taxon_set, tree2.taxon_set)
            self.assertDistinctButEqualTaxonSet(tree1.taxon_set, tree2.taxon_set, **kwargs)
        else:
            self.assertSame(tree1.taxon_set, tree2.taxon_set)
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
                    self.assertNotSame(node1.taxon, node2.taxon)
                else:
                    self.assertSame(node1.taxon, node2.taxon)
                if equal_oids is True:
                    self.assertEqual(node1.oid, node2.oid)
                elif equal_oids is False:
                    self.assertNotEqual(node1.oid, node2.oid)
                self.assertEqual(node1.taxon.label, node2.taxon.label)
                self.assertContained(node1.taxon, tree1.taxon_set)
                self.assertContained(node2.taxon, tree2.taxon_set)
            else:
                self.assertSame(node2.taxon, None)
            if node1.edge.length is not None:
                self.assertNotSame(node2.edge.length, None)
                self.assertAlmostEqual(node1.edge.length, node2.edge.length, 3)
            else:
                self.assertSame(node2.edge.length, None)
            self.assertEqual(len(node1.child_nodes()), len(node2.child_nodes()))
        tree1_edges = [edge for edge in tree1.postorder_edge_iter()]
        tree2_edges = [edge for edge in tree2.postorder_edge_iter()]
        self.assertEqual(len(tree1_edges), len(tree2_edges))
        for edge_idx, edge1 in enumerate(tree1_edges):
            edge2 = tree2_edges[edge_idx]
            self.assertNotSame(edge1, edge2)
            if edge1.length is None:
                self.assertSame(edge1.length, None)
            else:
                self.assertAlmostEqual(edge1.length, edge2.length, 2)
            if equal_oids is True:
                self.assertEqual(edge1.oid, edge2.oid)
            elif equal_oids is False:
                self.assertNotEqual(edge1.oid, edge2.oid)

    def assertDistinctButEqualStateAlphabetElement(self, sae1, sae2, **kwargs):
        pass

    def assertDistinctButEqualStateAlphabet(self, state_alphabet1, state_alphabet2, **kwargs):
        equal_oids = kwargs.get("equal_oids", None)
        if equal_oids is True:
            self.assertEqual(state_alphabet1.oid, state_alphabet2.oid)
        elif equal_oids is False:
            self.assertNotEqual(state_alphabet1.oid, state_alphabet2.oid)
        if state_alphabet1.missing is not None:
            self.assertDistinctButEqualStateAlphabetElement(state_alphabet1.missing, state_alphabet2.missing, **kwargs)
        else:
            self.assertSame(state_alphabet2.missing, None)
        self.assertEqual(len(state_alphabet1), len(state_alphabet2))
        for state_idx, state1 in enumerate(state_alphabet1):
            state2 = state_alphabet2[state_idx]
            self.assertDistinctButEqualStateAlphabetElement(state1, state2)

    def assertDistinctButEqualDiscreteCharArray(self, char_array1, char_array2, **kwargs):
        """
        `char_array1` and `char_array2` must be distinct but equivalent objects
         (if `distinct_state_alphabets` is False, as it should be for fixed state
         alphabets such as ).
        """
        distinct_state_alphabets = kwargs.get("distinct_state_alphabets", None)
        distinct_taxa = kwargs.get("distinct_taxa", True)
        equal_oids = kwargs.get("equal_oids", None)
        if distinct_taxa:
            self.assertNotSame(char_array1.taxon_set, char_array2.taxon_set)
            self.assertDistinctButEqualTaxonSet(char_array1.taxon_set, tree2.taxon_set, **kwargs)
        else:
            self.assertSame(char_array1.taxon_set, char_array2.taxon_set)
        if equal_oids is True:
            self.assertEqual(char_array1.oid, char_array2.oid)
        elif equal_oids is False:
            self.assertNotEqual(char_array1.oid, char_array2.oid)
        if distinct_state_alphabets is True:
            self.assertEqual(len(char_array1.state_alphabets), len(char_array2.state_alphabets))
            for sai, sa1 in enumerate(char_array1.state_alphabets):
                sa2 = char_array2.state_alphabets[sai]
                self.assertDistinctButEqualStateAlphabet(sa1, sa2)
        elif distinct_state_alphabets is False:
            self.assertSame(char_array1.state_alphabets, char_array2.state_alphabets)
            self.assertEqual(char_array1.state_alphabets, char_array2.state_alphabets)
            self.assertSame(char_array1.default_state_alphabet, char_array2.default_state_alphabet)
        self.assertEqual(len(char_array1.column_types), len(char_array2.column_types))
        for coli, col1 in enumerate(char_array1.column_types):
            if distinct_state_alphabets is True:
                col2 = char_array2.column_types[coli]
                self.assertDistinctButEqualStateAlphabet(col1.state_alphabet, col2.state_alphabet)
            elif distinct_state_alphabets is False:
                self.assertSame(col1.state_alphabet, col2.state_alphabet)
        for vi, vec1 in enumerate(char_array1):
            vec2 = char_array2[vi]
            taxon1 = vec1.taxon
            taxon2 = vec2.taxon
            if distinct_taxa:
                self.assertDistinctButEqualTaxon(taxon1, taxon2, **kwargs)
            else:
                self.assertSame(taxon1, taxon2)
            for i, c1 in enumerate(vec1):
                c2 = vec2[i]
                self.assertNotSame(c1, c2)
                if c1.column_type is not None:
                    self.assertNotSame(c1.column_type, c2.column_type)
                    if distinct_state_alphabets is True:
                        self.assertDistinctButEqualStateAlphabet(c1.column_type.state_alphabet, c2.column_type.state_alphabet)
                    elif distinct_state_alphabets is False:
                        self.assertSame(c1.column_type.state_alphabet, c2.column_type.state_alphabet)
                    self.assertContained(c1.column_type.state_alphabet, char_array1.state_alphabets)
                    self.assertContained(c2.column_type.state_alphabet, char_array2.state_alphabets)
                    self.assertContained(c1.value, c1.column_type.state_alphabet)
                    self.assertContained(c2.value, c2.column_type.state_alphabet)
                else:
                    self.assertSame(c2.column_type, None)
                    self.assertContained(c1.value, char_array1.default_state_alphabet)
                    self.assertContained(c2.value, char_array2.default_state_alphabet)
                if distinct_state_alphabets is True:
                    self.assertNotSame(c1.value, c2.value)
                elif distinct_state_alphabets is False:
                    self.assertSame(c1.value, c2.value)

