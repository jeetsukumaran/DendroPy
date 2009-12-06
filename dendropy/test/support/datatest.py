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
from dendropy.test.support import pathmap
from dendropy.test.support import extendedtest

class DataObjectVerificationTestCase(extendedtest.ExtendedTestCase):
    """
    Extends ExtendedTestCase with tests for data object comparisons.
    """

    def roundTripData(self, dataset, format):
        output_path = pathmap.named_output_path(filename="roundtrip_test.%s" % format, suffix_timestamp=True)
        self.logger.info("Writing DataSet object to file '%s' and back" % output_path)
        dataset.write(stream=open(output_path, "w"), format=format)
        self.logger.info("Reading file '%s'" % output_path)
        return dendropy.DataSet(stream=open(output_path, "rU"), format=format)

    def roundTripDataSetTest(self,
                      dataset,
                      format,
                      **kwargs):
        rt_dataset = self.roundTripData(dataset, format)
        self.logger.info("Comparing original and round-tripped DataSet objects")
        self.assertDistinctButEqual(dataset, rt_dataset, **kwargs)

    def assertDistinctButEqual(self, data_object1, data_object2, **kwargs):
        """
        Verifies that two DendroPy phylogenetic data objects (Tree, TreeList,
        CharMatrix, DataSet etc.) are independent objects, but equal. That is,
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
        if isinstance(data_object1, dendropy.DataSet):
            self.assertDistinctButEqualDataSet(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.Taxon):
            self.assertDistinctButEqualTaxon(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.TaxonSet):
            self.assertDistinctButEqualTaxonSet(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.Tree):
            self.assertDistinctButEqualTree(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.TreeList):
            self.assertDistinctButEqualTreeList(data_object1, data_object2, **kwargs)
        elif isinstance(data_object1, dendropy.CharacterMatrix):
            self.assertDistinctButEqualCharMatrix(data_object1, data_object2, **kwargs)
        else:
            raise ValueError("Unsupported type for comparison: %s" % type(data_object1))

    def assertDistinctButEqualDataSet(self, dataset1, dataset2, **kwargs):
        if "distinct_taxa" in kwargs and not kwargs["distinct_taxa"]:
            raise Exception("Distinct TaxonSet objects criterion must be enforced when comparing DataSet objects")
        self.logger.info("Comparing DataSet objects %d and %d" % (id(dataset1), id(dataset2)))

        self.assertEqual(len(dataset1.taxon_sets), len(dataset2.taxon_sets))
        for idx, taxon_set1 in enumerate(dataset1.taxon_sets):
            taxon_set2 = dataset2.taxon_sets[idx]
            self.assertDistinctButEqualTaxonSet(taxon_set1, taxon_set2, **kwargs)

        self.assertEqual(len(dataset1.char_matrices), len(dataset2.char_matrices))
        for tsi, char_matrix1 in enumerate(dataset1.char_matrices):
            char_matrix2 = dataset2.char_matrices[tsi]
            self.assertDistinctButEqualCharMatrix(char_matrix1, char_matrix2, **kwargs)

        self.assertEqual(len(dataset1.tree_lists), len(dataset2.tree_lists))
        for tsi, tree_list1 in enumerate(dataset1.tree_lists):
            tree_list2 = dataset2.tree_lists[tsi]
            self.assertDistinctButEqualTreeList(tree_list1, tree_list2, **kwargs)

    def assertDistinctButEqualTaxon(self, taxon1, taxon2, **kwargs):
        equal_oids = kwargs.get("equal_oids", None)
        ignore_underscore_substitution = kwargs.get("ignore_underscore_substitution", False)
        self.logger.info("Comparing Taxon objects %d ('%s') and %d ('%s')" % (id(taxon1), taxon1.label, id(taxon2), taxon2.label))
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
        ignore_taxon_order = kwargs.get("ignore_taxon_order", False)
        self.logger.info("Comparing TaxonSet objects %d and %d" % (id(taxon_set1), id(taxon_set2)))
        self.assertNotSame(taxon_set1, taxon_set2)
        self.assertEqual(len(taxon_set1), len(taxon_set2))
        if not ignore_taxon_order:
            for tidx, taxon1 in enumerate(taxon_set1):
                taxon2 = taxon_set2[tidx]
                self.assertDistinctButEqualTaxon(taxon1, taxon2, **kwargs)
        else:
            if ignore_underscore_substitution:
                get_label = lambda t: str(t.label).replace("_", " ")
            else:
                get_label = lambda t: str(t.label)
            labels1 = set([get_label(t) for t in taxon_set1])
            labels2 = set([get_label(t) for t in taxon_set2])
            self.assertEqual(labels1, labels2)
            if equal_oids is not None:
                oids1 = set([t.oid for t in taxon_set1])
                oids2 = set([t.oid for t in taxon_set2])
                if equal_oids:
                    self.assertEqual(oids1, oids2)
                else:
                    self.assertNotEqual(oids1, oids2)

    def assertDistinctButEqualTreeList(self, tree_list1, tree_list2, **kwargs):
        """
        `tree_list1` and `tree_list2` must be distinct but equivalent objects.
        """
        distinct_taxa = kwargs.get("distinct_taxa", True)
        equal_oids = kwargs.get("equal_oids", None)
        self.logger.info("Comparing TreeList objects %d and %d" % (id(tree_list1), id(tree_list2)))
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
        self.logger.info("Comparing Tree objects %d and %d" % (id(tree1), id(tree2)))
        self.logger.info(tree1.as_newick_str())
        tree1.debug_check_tree(logger=self.logger)
        self.logger.info(tree2.as_newick_str())
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
                    if equal_oids is True:
                        self.assertEqual(node1.oid, node2.oid)
                    elif equal_oids is False:
                        self.assertNotEqual(node1.oid, node2.oid)
                    self.assertEqual(node1.taxon.label, node2.taxon.label)
                    self.assertContained(node1.taxon, tree1.taxon_set)
                    self.assertContained(node2.taxon, tree2.taxon_set)
                else:
                    self.assertSame(node1.taxon, node2.taxon)
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
        equal_oids = kwargs.get("equal_oids", None)
        equal_symbols = kwargs.get("equal_symbols", True)
        self.assertNotSame(sae1, sae2)
        if equal_oids is True:
            self.assertEqual(sae1.oid, sae2.oid)
        elif equal_oids is False:
            self.assertNotEqual(sae1.oid, sae2.oid)
        if equal_symbols:
            self.assertEqual(sae1.symbol, sae2.symbol)
        self.assertEqual(sae1.multistate, sae2.multistate)
        if sae1.member_states is not None:
            self.assertNotSame(sae2.member_states, None)
            for mi, ms1 in enumerate(sae1.member_states):
                ms2 = sae2.member_states[mi]
                self.assertDistinctButEqualStateAlphabetElement(ms1, ms2, **kwargs)
        else:
            self.assertSame(sae2.member_states, None)

    def assertDistinctButEqualStateAlphabet(self, state_alphabet1, state_alphabet2, **kwargs):
        equal_oids = kwargs.get("equal_oids", None)
        self.assertNotSame(state_alphabet1, state_alphabet2)
        self.logger.info("Comparing StateAlphabet objects %d and %d" % (id(state_alphabet1), id(state_alphabet2)))
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

    def assertDistinctButEqualDiscreteCharMatrix(self, char_matrix1, char_matrix2, **kwargs):
        """
        `char_matrix1` and `char_matrix2` must be distinct but equivalent objects
         (if `distinct_state_alphabets` is False, as it should be for fixed state
         alphabets such as ).
        """
        distinct_state_alphabets = kwargs.get("distinct_state_alphabets", None)
        distinct_taxa = kwargs.get("distinct_taxa", True)
        equal_oids = kwargs.get("equal_oids", None)
        ignore_chartypes = kwargs.get("ignore_chartypes", True)
        self.logger.info("Comparing DiscreteCharacterMatrix objects %d and %d" % (id(char_matrix1), id(char_matrix2)))
        self.assertNotSame(char_matrix1, char_matrix2)
        if distinct_taxa:
            self.assertNotSame(char_matrix1.taxon_set, char_matrix2.taxon_set)
            self.assertDistinctButEqualTaxonSet(char_matrix1.taxon_set, char_matrix2.taxon_set, **kwargs)
        else:
            self.assertSame(char_matrix1.taxon_set, char_matrix2.taxon_set)
        if equal_oids is True:
            self.assertEqual(char_matrix1.oid, char_matrix2.oid)
        elif equal_oids is False:
            self.assertNotEqual(char_matrix1.oid, char_matrix2.oid)
        if distinct_state_alphabets is True:
            self.assertEqual(len(char_matrix1.state_alphabets), len(char_matrix2.state_alphabets))
            for sai, sa1 in enumerate(char_matrix1.state_alphabets):
                sa2 = char_matrix2.state_alphabets[sai]
                self.assertDistinctButEqualStateAlphabet(sa1, sa2)
        elif distinct_state_alphabets is False:
            for sai, sa1 in enumerate(char_matrix1.state_alphabets):
                sa2 = char_matrix2.state_alphabets[sai]
                self.assertSame(sa1, sa2)
                self.assertSame(char_matrix1.default_state_alphabet, char_matrix2.default_state_alphabet)
        if not ignore_chartypes:
            self.assertEqual(len(char_matrix1.character_types), len(char_matrix2.character_types))
        for coli, col1 in enumerate(char_matrix1.character_types):
            if distinct_state_alphabets is True:
                col2 = char_matrix2.character_types[coli]
                self.assertDistinctButEqualStateAlphabet(col1.state_alphabet, col2.state_alphabet)
            elif distinct_state_alphabets is False:
                self.assertSame(col1.state_alphabet, col2.state_alphabet)
        self.assertEqual(len(char_matrix1), len(char_matrix2))

        for ti, taxon1 in enumerate(char_matrix1):
            vec1 = char_matrix1[taxon1]
            taxon2 = char_matrix2.taxon_set[ti]
            vec2 = char_matrix2[taxon2]
            self.logger.info("Comparing CharacterDataVector objects %d and %d" % (id(vec2), id(vec2)))
            if distinct_taxa:
                self.assertDistinctButEqualTaxon(taxon1, taxon2, **kwargs)
            else:
                self.assertSame(taxon1, taxon2)
            self.logger.info("%s: %s" % (taxon1.label, vec1.symbols_as_string()))
            self.logger.info("%s: %s" % (taxon2.label, vec2.symbols_as_string()))
            self.assertEqual(len(vec1), len(vec2))
            for i, c1 in enumerate(vec1):
                c2 = vec2[i]
                self.assertNotSame(c1, c2)
                if len(char_matrix1.state_alphabets) == 1:
                    self.assertContained(c1.value, char_matrix1.state_alphabets[0])
                    self.assertContained(c2.value, char_matrix2.state_alphabets[0])
                else:
                    # assume mapped by columns, and will be checked there
                    pass
                if distinct_state_alphabets is True:
                    self.assertDistinctButEqualStateAlphabetElement(c1.value, c2.value)
                elif distinct_state_alphabets is False:
                    self.assertSame(c1.value, c2.value)
                if not ignore_chartypes:
                    if c1.character_type is not None:
                        self.assertNotSame(c1.character_type, c2.character_type)
                        if distinct_state_alphabets is True:
                            self.assertDistinctButEqualStateAlphabet(c1.character_type.state_alphabet, c2.character_type.state_alphabet)
                        elif distinct_state_alphabets is False:
                            self.assertSame(c1.character_type.state_alphabet, c2.character_type.state_alphabet)
                        self.assertContained(c1.character_type.state_alphabet, char_matrix1.state_alphabets)
                        self.assertContained(c2.character_type.state_alphabet, char_matrix2.state_alphabets)
                        self.assertContained(c1.value, c1.character_type.state_alphabet)
                        self.assertContained(c2.value, c2.character_type.state_alphabet)
                    else:
                        self.assertSame(c2.character_type, None)

    def assertDistinctButEqualCharMatrix(self, char_matrix1, char_matrix2, **kwargs):
        if isinstance(char_matrix1, dendropy.DiscreteCharacterMatrix):
            self.assertDistinctButEqualDiscreteCharMatrix(char_matrix1, char_matrix2, **kwargs)
        else:
            raise NotImplementedError()

    def text_to_label_symbol_tuples(self, text):
        """
        Takes a tab-delimited string in the form of:
            <TAXON_NAME>\t<CHARACTERS>
        and returns a list of pairs with first element the
        taxon label and the second a string of state symbols.
        """
        data = []
        for i in text.split("\n"):
            if i:
                j = i.split("\t")
                assert len(j) == 2
                if j:
                    if j[0] and j[1]:
                        values = [j[0], j[1]]
                        data.append(values)
        return data

    def text_to_label_value_list(self, text, val_type=float):
        """
        Takes a tab-delimited string in the form of:
            <TAXON_NAME>\t<CHAR1> <CHAR2> <CHAR3> ...
        and returns a list of pairs with first element the
        taxon label and the second a list of values.
        """
        data = []
        for i in text.split("\n"):
            if i:
                j = i.split("\t")
                assert len(j) == 2
                if j:
                    if j[0] and j[1]:
                        val_list = [val_type(x) for x in j[1].split(' ')]
                        values = [j[0], val_list]
                        data.append(values)
        return data

    def char_matrix_to_label_symbol_tuples(self, char_matrix):
        """
        Takes a `char_matrix` and returns a list of pairs with first element the
        taxon label and the second a string of state symbols.
        """
        data = []
        for t, s in char_matrix.items():
            data.append((t.label, s.symbols_as_string()))
        return data

    def assertEqualCharMatrixLabelSymbols(self, char_matrix, **kwargs):
        """
        Takes a CharacterMatrix object, extracts tuples in the form of
        (<Taxon label string>, <sequence symbols string>), and compares it with
        *one* of the following passed as keyword arguments:

            - `expected_label_symbol_tuples`, a list of tuples of strings
               corresponding to the tuples expected from the CharacterMatrix object.

            - `expected_label_symbol_text`, a tab-delimited string with rows
               describing the tuples expected from the CharacterMatrix object
               in the form of::

                    <Taxon label string>\t<sequence symbols string>

            - `expected_label_symbol_stream`, a file-like object with the file
               containing rows describing the tuples expected from the
               CharacterMatrix object in the form of::

                    <Taxon label string>\t<sequence symbols string>

        """
        ignore_underscore_substitution = kwargs.get("ignore_underscore_substitution", False)
        if "expected_label_symbol_tuples" in kwargs:
            expected_label_symbol_tuples = kwargs["expected_label_symbol_tuples"]
        elif "expected_label_symbol_text" in kwargs:
            expected_label_symbol_tuples = self.text_to_label_symbol_tuples(kwargs["expected_label_symbol_text"])
        elif "expected_label_symbol_stream" in kwargs:
            expected_label_symbol_text = kwargs["expected_label_symbol_stream"].read()
            expected_label_symbol_tuples = self.text_to_label_symbol_tuples(expected_label_symbol_text)

        obs_label_symbol_tuples = self.char_matrix_to_label_symbol_tuples(char_matrix)

        self.assertEqual(len(obs_label_symbol_tuples), len(expected_label_symbol_tuples))
        for i, x1 in enumerate(expected_label_symbol_tuples):
            tax_label1 = x1[0]
            seq_symbols1 = x1[1]
            tax_label2 = obs_label_symbol_tuples[i][0]
            seq_symbols2 = obs_label_symbol_tuples[i][1]
            if ignore_underscore_substitution:
                tax_label1 = tax_label1.replace("_", " ")
                tax_label2 = tax_label2.replace("_", " ")
            self.assertEqual(tax_label1, tax_label2)
            self.assertEqual(seq_symbols1, seq_symbols2)

    def assertEqualCharMatrixLabelContinuousValues(self, char_matrix, **kwargs):
        ignore_underscore_substitution = kwargs.get("ignore_underscore_substitution", False)
        if "expected_label_symbol_tuples" in kwargs:
            expected_label_value_list = kwargs["expected_label_symbol_tuples"]
        elif "expected_label_symbol_text" in kwargs:
            expected_label_value_list = self.text_to_label_value_list(kwargs["expected_label_symbol_text"])
        elif "expected_label_symbol_stream" in kwargs:
            expected_label_value_text = kwargs["expected_label_symbol_stream"].read()
            expected_label_value_list = self.text_to_label_value_list(expected_label_value_text)
        val_type = kwargs.get("val_type", float)

        obs_label_value_list = []
        for t in char_matrix.taxon_set:
            obs_label_value_list.append([t.label, [val_type(v.value) for v in char_matrix[t]]])

        self.assertEqual(len(obs_label_value_list), len(expected_label_value_list))
        for i, x1 in enumerate(expected_label_value_list):
            tax_label1 = x1[0]
            seq_values1 = x1[1]
            tax_label2 = obs_label_value_list[i][0]
            seq_values2 = obs_label_value_list[i][1]
            if ignore_underscore_substitution:
                tax_label1 = tax_label1.replace("_", " ")
                tax_label2 = tax_label2.replace("_", " ")
            self.assertEqual(tax_label1, tax_label2)
            for j, v1 in enumerate(seq_values1):
                v2 = seq_values2[j]
                self.assertAlmostEqual(v1, v2, 4)

class ComplexMultiTaxonSetDataVerificationTest(DataObjectVerificationTestCase):

    def setUp(self):
        self.taxon_set_names = [
            "Pythonidae",
            "Caenophidia",
            "Primates"
        ]
        self.taxon_set_taxon_labels = [
            ["Python molurus", "Python sebae", "Morelia tracyae", "Morelia amethistina",
             "Morelia nauta", "Morelia kinghorni", "Morelia clastolepis", "Morelia boeleni",
             "Python reticulatus", "Python timoriensis", "Morelia oenpelliensis",
             "Morelia viridis", "Morelia carinata", "Morelia spilota", "Morelia bredli",
             "Antaresia maculosa", "Antaresia childreni", "Antaresia stimsoni",
             "Antaresia perthensis", "Leiopython albertisii", "Bothrochilus boa",
             "Liasis olivaceus", "Liasis mackloti", "Liasis fuscus", "Apodora papuana",
             "Aspidites ramsayi", "Aspidites melanocephalus", "Python brongersmai", "Python regius"],
            ["Lystrophis dorbignyi", "Waglerophis merremi", "Lystrophis histricus",
              "Xenoxybelis argenteus", "Liophis jaegeri", "Liophis elegantissimus",
              "Liophis meridionalis", "Liophis typhlus", "Erythrolamprus aesculapii",
              "Atractus albuquerquei", "Atractus trihedrurus", "Heterodon nasicus",
              "Sibynomorphus mikanii", "Ninia atrata", "Dipsas indica", "Dipsas neivai",
              "Sibynomorphus garmani", "Hydrops triangularis", "Helicops infrataeniatus",
              "Pseudoeryx plicatilis", "Helicops pictiventris", "Helicops gomesi",
              "Helicops angulatus", "Oxyrhopus rhombifer", "Oxyrhopus clathratus",
              "Lycognathophis seychellensis", "Taeniophallus brevirostris",
              "Hydrodynastes gigas", "Hydrodynastes bicinctus", "Pseudoboa nigra",
              "Pseudoboa coronata", "Boiruna maculata", "Drepanoides anomalus",
              "Clelia bicolor", "Phimophis guerini", "Psomophis joberti",
              "Psomophis genimaculatus", "Leptodeira annulata", "Thamnodynastes rutilus",
              "Pseudotomodon trigonatus", "Tachymenis peruviana", "Tomodon dorsatus",
              "Ptychophis flavovirgatus", "Calamodontophis paucidens", "Pseudablabes agassizi",
              "Philodryas patagoniensis", "Philodryas aestiva", "Philodryas mattogrossensis",
              "Taeniophallus affinis", "Liophis amarali", "Imantodes cenchoa",
              "Apostolepis dimidiata", "Phalotris nasutus", "Phalotris lemniscatus",
              "Apostolepis assimilis", "Elapomorphus quinquelineatus", "Siphlophis pulcher",
              "Siphlophis compressus", "Dendroaspis polylepis", "Bothrolycus ater",
              "Bothrophthalmus lineatus", "Bothrophthalmus brunneus", "Lamprophis virgatus",
              "Lamprophis fuliginosus mentalis", "Lamprophis lineatus",
              "Lamprophis fuliginosus", "Lamprophis olivaceus", "Lamprophis capensis",
              "Lamprophis inornatus", "Lycodonomorphus whytii", "Lamprophis fiskii",
              "Lycodonomorphus rufulus", "Lamprophis guttatus", "Homoroselaps lacteus",
              "Atractaspis bibronii", "Atractaspis corpulenta", "Atractaspis boulengeri",
              "Atractaspis micropholis", "Aparallactus modestus", "Polemon collaris",
              "Polemon acanthias", "Polemon notatum", "Xenocalamus transvaalensis",
              "Amblyodipsas polylepis", "Macrelaps microlepidotus", "Mehelya nyassae",
              "Psammophis phillipsi", "Psammophis praeornatus", "Psammophis sp.",
              "Psammophis schokari", "Prosymna janii", "Dromicodryas bernieri",
              "Buhoma depressiceps", "Liophidium chabaudi", "Heteroliodon occipitalis",
              "Lycophidion capense", "Rhamphiophis oxyrhynchus", "Hemirhagerrhis hildebrandtii",
              "Amplorhinus multimaculatus", "Lycophidion laterale", "Buhoma procterae",
              "Mehelya stenophthalmus", "Psammodynastes sp.", "Rhamphiophis rostratus",
              "Pseudaspis cana", "Gonionotophis brussauxi", "Psammophylax rhombeatus",
              "Psammophylax variabilis", "Mehelya poensis", "Duberria variegata",
              "Duberria lutrix", "Stenophis citrinus", "Lycodryas sanctijohannis",
              "Ithycyphus oursi"],
            ["Lemur catta", "Homo sapiens", "Pan", "Gorilla", "Pongo",
              "Hylobates", "Macaca fuscata", "Macaca mulatta", "Macaca fascicularis",
              "Macaca sylvanus", "Saimiri sciureus", "Tarsius syrichta"]
        ]
        self.tree_list_labels = [
            "Pythonidae Random Trees",
            "Caenophidia Random Yule 1",
            "Caenophidia Random Yule 2",
            "Primates Random Yule",
            "Pythonidae CytB MLE",
        ]
        self.tree_list_taxon_set_labels = {
            "Pythonidae Random Trees": "Pythonidae",
            "Caenophidia Random Yule 1": "Caenophidia",
            "Caenophidia Random Yule 2": "Caenophidia",
            "Primates Random Yule": "Primates",
            "Pythonidae CytB MLE": "Pythonidae",
        }
        self.tree_list_strings = [
            """
            #NEXUS
            BEGIN TREES;
                Title Pythonidae_Random_Trees;
                TRANSLATE
                    1 Python_molurus,
                    2 Python_sebae,
                    3 Morelia_tracyae,
                    4 Morelia_amethistina,
                    5 Morelia_nauta,
                    6 Morelia_kinghorni,
                    7 Morelia_clastolepis,
                    8 Morelia_boeleni,
                    9 Python_reticulatus,
                    10 Python_timoriensis,
                    11 Morelia_oenpelliensis,
                    12 Morelia_viridis,
                    13 Morelia_carinata,
                    14 Morelia_spilota,
                    15 Morelia_bredli,
                    16 Antaresia_maculosa,
                    17 Antaresia_childreni,
                    18 Antaresia_stimsoni,
                    19 Antaresia_perthensis,
                    20 Leiopython_albertisii,
                    21 Bothrochilus_boa,
                    22 Liasis_olivaceus,
                    23 Liasis_mackloti,
                    24 Liasis_fuscus,
                    25 Apodora_papuana,
                    26 Aspidites_ramsayi,
                    27 Aspidites_melanocephalus,
                    28 Python_brongersmai,
                    29 Python_regius;
                TREE Tree01 = ((11:5.5010411471,14:5.5010411471)Node36995248:61.991593,((((((25:0.7074902224,(27:0.2369844919,13:0.2369844919)Node36995920:0.470506)Node36995600:2.097596,(21:1.0695698693,4:1.0695698693)Node36995952:1.735517)Node36995536:0.35162,5:3.1567060575)Node36995472:1.093757,(3:0.585328938,29:0.585328938)Node36996400:3.665134)Node36995440:33.019561,((((9:1.3252196473,(28:0.6261312708,22:0.6261312708)Node36996880:0.699088)Node36996464:3.583644,10:4.9088635169)Node36996560:4.903818,((7:0.5202931818,16:0.5202931818)Node36997072:3.896642,(24:0.1932288198,23:0.1932288198)Node36997360:4.223707)Node36997168:5.395746)Node36996176:3.492068,8:13.3047497112)Node36996528:23.965274)Node36995344:0.377328,((1:12.8622487575,((18:1.2013217474,2:1.2013217474)Node36997712:3.830584,6:5.03190605)Node36997840:7.830343)Node36997616:21.168528,((26:2.7873255922,((17:0.1054281369,20:0.1054281369)Node36998032:1.300192,19:1.4056203892)Node36998320:1.381705)Node36997968:4.806974,(12:2.4211779654,15:2.4211779654)Node36998704:5.173121)Node36998064:26.436477)Node36997392:3.616575)Node36962224:29.845282)Node36995184;
                TREE Tree02 = (8:26.0384706379,((((20:15.2611605145,21:15.2611605145)Node36999024:2.956938,(13:16.0625616183,(((14:7.0505997998,15:7.0505997998)Node37061072:7.474769,(16:13.0854581588,(19:10.2267471766,(17:4.8707020453,18:4.8707020453)Node37061584:5.356045)Node37061456:2.858711)Node37061392:1.439911)Node37060944:0.166805,12:14.6921741266)Node37061040:1.370387)Node37060720:2.155537)Node36962160:1.648835,((22:10.850645624,(24:2.8304830038,23:2.8304830038)Node37062064:8.020163)Node37061872:1.940949,25:12.7915950606)Node37061744:7.075338)Node36995152:0.952266,((11:19.051793087,(((4:8.1204292308,(7:2.0657194471,(6:0.9996517243,5:0.9996517243)Node37062896:1.066068)Node37062768:6.05471)Node37062544:1.954125,3:10.0745540733)Node37062448:7.198644,(9:6.6158514927,10:6.6158514927)Node37063184:10.657346)Node37062512:1.778595)Node37062288:1.426579,((2:14.3734153538,(29:11.7636662419,(28:10.5318200836,1:10.5318200836)Node37063664:1.231846)Node37063536:2.609749)Node37063248:3.933241,(26:9.727198476,27:9.727198476)Node37063760:8.579458)Node37063088:2.171715)Node37062192:0.340828)Node36999056:5.219271)Node36995120;
                TREE Tree03 = (((((2:16.9224723285,((1:13.9521024428,28:13.9521024428)Node37064432:0.218629,29:14.170731254)Node37064496:2.751741)Node37064048:2.973954,(26:7.827543828,27:7.827543828)Node37654640:12.068883)Node37063696:3.261106,(((13:11.7896562791,(14:4.0544526578,15:4.0544526578)Node37655120:7.735204)Node37064656:2.978776,(16:13.5106259111,(19:11.2360265701,(17:3.4884792675,18:3.4884792675)Node37655536:7.747547)Node37655408:2.274599)Node37655152:1.257806)Node37654896:2.903479,12:17.6719111069)Node37654704:5.485621)Node36998896:1.069297,(20:13.3580917481,21:13.3580917481)Node37655728:10.868738)Node37064208:0.036473,(((22:6.2148795185,(24:2.6053610856,23:2.6053610856)Node37656272:3.609518)Node37655632:1.914716,25:8.1295951417)Node37655984:14.87802,(8:20.4730397427,((((4:4.4332464164,((5:0.8115194045,6:0.8115194045)Node37656912:0.622029,7:1.4335488208)Node37656976:2.999698)Node37656304:3.286251,3:7.719497883)Node37656688:9.885018,(9:5.7810239951,10:5.7810239951)Node37657392:11.823492)Node37656528:0.296425,11:17.9009406942)Node37656656:2.572099)Node37656560:2.534576)Node37656048:1.255687)Node37064144;
                TREE Tree04 = (((2:16.5921325784,(29:14.3799352093,(1:13.3821288803,28:13.3821288803)Node37658128:0.997806)Node37658000:2.212197)Node36962032:3.430395,(26:11.703431691,27:11.703431691)Node37658224:8.319096)Node37657776:3.092724,((((9:6.4632684179,10:6.4632684179)Node37658160:13.34671,(((4:4.6709208052,(5:2.0888008587,(6:1.4812304636,7:1.4812304636)Node37712496:0.60757)Node37712368:2.58212)Node37711984:1.717541,3:6.3884621114)Node37712016:3.030316,8:9.4187780678)Node37658576:10.391201)Node37658288:0.121668,11:19.9316466464)Node37657648:2.717608,((((12:15.9488284198,((16:14.0002983201,(19:11.4477602235,(17:5.7284429434,18:5.7284429434)Node37713488:5.719317)Node37713360:2.552538)Node37713072:1.553888,(14:4.4953163706,15:4.4953163706)Node37713648:11.05887)Node37713200:0.394642)Node37712400:2.234726,13:18.1835548451)Node37712816:1.848444,(20:12.7456406695,21:12.7456406695)Node37714000:7.286358)Node37712848:1.892237,((22:13.5304034382,(24:1.7445458007,23:1.7445458007)Node37714352:11.785858)Node37714128:1.710768,25:15.2411713811)Node37714192:6.683064)Node37712912:0.725019)Node37658384:0.465997)Node37657712;
                TREE Tree05 = ((9:10.8300265774,10:10.8300265774)Node37714768:17.634323,((((12:14.0077690635,((13:11.7785955814,(14:4.9757977016,15:4.9757977016)Node37715504:6.802798)Node37715184:1.825003,(16:12.3274764924,(19:8.0033214502,(17:4.0152243781,18:4.0152243781)Node37715760:3.988097)Node37715792:4.324155)Node37715536:1.276123)Node37715312:0.40417)Node37714992:4.885081,(20:11.1450426742,21:11.1450426742)Node37715600:7.747808)Node37714960:1.317722,((22:10.9742846954,(24:3.1666996368,23:3.1666996368)Node37773840:7.807585)Node37773616:3.787122,25:14.7614067829)Node37773680:5.449165)Node37714512:4.529061,(((((4:5.5031961516,(7:1.5725765811,(5:1.4357285473,6:1.4357285473)Node37774608:0.136848)Node37774480:3.93062)Node37774160:4.509387,3:10.0125827763)Node37773904:6.858868,8:16.8714507361)Node37774032:5.485833,11:22.3572839238)Node37774128:1.80695,((2:16.535465064,(29:14.962419097,(28:13.1354320851,1:13.1354320851)Node37775312:1.826987)Node37775184:1.573046)Node37774928:4.520647,(26:14.6640030997,27:14.6640030997)Node37775408:6.392109)Node37775024:3.108122)Node37773872:0.575398)Node37715024:3.724717)Node37714704;
                TREE Tree06 = (((((((16:12.6196623606,(19:10.1068643474,(17:4.3931890637,18:4.3931890637)Node37776432:5.713675)Node37776304:2.512798)Node37776016:1.834271,(12:13.41617982,(14:4.0398271364,15:4.0398271364)Node37776784:9.376353)Node37776144:1.037753)Node37775952:0.835746,13:15.2896785188)Node37775472:1.651414,(20:12.7629580607,21:12.7629580607)Node37777072:4.178134)Node37775568:0.754225,((22:10.6502403013,(24:3.6716556879,23:3.6716556879)Node37777360:6.978585)Node37777232:4.064978,25:14.7152178676)Node37777264:2.980099)Node37064112:4.099713,((((4:4.6447123226,(7:2.1717490576,(6:0.9040957981,5:0.9040957981)Node37831344:1.267653)Node37831216:2.472963)Node37830928:4.270569,3:8.9152813163)Node37777040:3.165779,8:12.0810604926)Node37830704:6.584465,11:18.6655255697)Node37830864:3.129505)Node37775856:1.011175,((9:11.5585932535,10:11.5585932535)Node37831504:11.17184,((2:15.9676320675,((28:13.8599706405,1:13.8599706405)Node37832048:1.044872,29:14.9048428404)Node37832176:1.062789)Node37831888:5.314147,(26:15.5786974573,27:15.5786974573)Node37832400:5.703082)Node37832016:1.448654)Node37831600:0.075772)Node37775792;
                TREE Tree07 = ((2:29.0711850036,(1:22.4124826989,(29:21.9264827074,28:21.9264827074)Node37833136:0.486)Node37833008:6.658702)Node37832848:2.122614,((((9:14.5492404256,10:14.5492404256)Node37833168:5.892517,((((12:14.834011977,((14:6.2738399067,15:6.2738399067)Node37833904:7.331733,13:13.6055732883)Node37834032:1.228439)Node37833520:0.322649,((19:10.1091415386,(17:3.4832007905,18:3.4832007905)Node37834512:6.625941)Node37834288:4.225992,16:14.3351337411)Node37834256:0.821527)Node37833680:2.938028,(20:9.2842258749,21:9.2842258749)Node37834736:8.810464)Node37833456:0.826871,((22:11.0151449224,(24:4.0880111593,23:4.0880111593)Node37101968:6.927134)Node37101808:2.70301,25:13.7181551735)Node37101712:5.203405)Node37833744:1.520198)Node37833328:1.10888,(11:17.6009445592,((3:5.8027959732,(4:4.3810721744,((5:1.8029155441,6:1.8029155441)Node37102608:0.077449,7:1.8803643554)Node37102672:2.500708)Node37102544:1.421724)Node37102288:6.371263,8:12.174058564)Node37102352:5.426886)Node37101936:3.949694)Node37833360:1.246336,(26:9.8470454886,27:9.8470454886)Node37103024:12.949929)Node37833200:8.396824)Node37832784;
                TREE Tree08 = (((((((3:8.3143208605,(((5:0.9126763167,7:0.9126763167)Node37103920:0.714408,6:1.6270842916)Node37103824:2.881324,4:4.5084082245)Node37103888:3.805913)Node37103568:4.138097,8:12.452417468)Node37103504:9.507796,((25:17.2449081924,(22:12.2505582513,(24:6.2921986908,23:6.2921986908)Node37104720:5.95836)Node37104592:4.99435)Node37104336:1.509383,((20:13.0633246673,21:13.0633246673)Node37104752:5.193371,(((19:10.3499331345,(17:3.3857328938,18:3.3857328938)Node37105392:6.9642)Node37105200:3.437917,16:13.7878501853)Node37104880:2.063714,(12:15.2096070352,(13:13.6526925966,(14:4.9216311908,15:4.9216311908)Node37163248:8.731061)Node37105232:1.556914)Node37105584:0.641957)Node37104976:2.405131)Node37104784:0.497596)Node37104432:3.205922)Node37103120:0.983173,((9:14.7952666012,10:14.7952666012)Node37163312:6.613431,11:21.4086976551)Node37103280:1.534689)Node37102896:4.226216,(26:6.6352256469,27:6.6352256469)Node37163504:20.534377)Node37714672:6.419829,(1:23.4003588995,(28:23.3469263821,29:23.3469263821)Node37164080:0.053433)Node37163696:10.189073)Node37103408:8.436564,2:42.025995763)Node37103344;
                TREE Tree09 = (((((12:14.5888315231,(((14:4.6124916508,15:4.6124916508)Node37164880:5.065344,13:9.6778360646)Node37164720:4.896004,(16:14.05750906,(19:7.6084936983,(17:1.8635827152,18:1.8635827152)Node37165456:5.744911)Node37165328:6.449015)Node37165264:0.516331)Node37164848:0.014991)Node37164304:2.826242,((9:8.600122812,10:8.600122812)Node37165616:8.738044,(((25:13.0669390597,(20:9.3384061584,21:9.3384061584)Node37166192:3.728533)Node37165808:1.387248,(22:10.9703192915,(24:7.1343513196,23:7.1343513196)Node37166512:3.835968)Node37166384:3.483868)Node37165840:1.989938,(((3:8.3209994023,(((6:1.0739843066,5:1.0739843066)Node37166992:0.774827,7:1.8488112596)Node37166896:2.025233,4:3.8740441117)Node37166960:4.446955)Node37166704:1.741944,8:10.0629436951)Node37166544:3.527281,11:13.5902249452)Node37166576:2.853901)Node37165968:0.894041)Node37165712:0.076906)Node37775760:4.091583,(26:9.1538829568,27:9.1538829568)Node37220432:12.352774)Node37164368:3.290072,((29:22.6656294361,28:22.6656294361)Node37220944:0.519654,1:23.1852834479)Node37220976:1.611445)Node37164496:4.512234,2:29.3089629268)Node37164432;
                TREE Tree10 = (2:41.037475788,((29:32.1860160849,(1:28.8078489217,28:28.8078489217)Node37221776:3.378167)Node37221616:7.880428,(((11:20.2736124855,((3:8.9277450894,(((6:3.1533524892,5:3.1533524892)Node37222448:0.318849,7:3.4722017773)Node37222288:4.468139,4:7.9403403274)Node37222416:0.987405)Node37222096:3.537623,8:12.4653684465)Node37222224:7.808244)Node37221872:3.045662,(((22:14.4158757568,(24:4.302477921,23:4.302477921)Node37223184:10.113398)Node37222864:6.094015,((9:11.1261504133,10:11.1261504133)Node37223376:9.215066,25:20.3412168975)Node37223248:0.168673)Node37222960:1.171445,((20:16.554642571,21:16.554642571)Node37223408:4.801492,(((13:13.6491727777,(14:5.9537557687,15:5.9537557687)Node37224208:7.695417)Node37223888:0.593107,12:14.242279557)Node37223728:1.515694,(16:13.2829808571,(19:9.8084470504,(17:5.8114705139,18:5.8114705139)Node37277904:3.996977)Node37224368:3.474534)Node37224304:2.474992)Node37223760:5.598162)Node37223632:0.3252)Node37222736:1.63794)Node37221968:5.612331,(26:9.3156335341,27:9.3156335341)Node37224432:19.615972)Node37221936:11.134838)Node37221584:0.971032)Node37221424;
            END;
        """,
        """
            #NEXUS
            BEGIN TREES;
                Title Caenophidia_Random_Yule_1;
            [!Parameters: Tree simulator: Uniform speciation (Yule). [seed: 1259785426822]]
                TRANSLATE
                    1 Lystrophis_dorbignyi,
                    2 Waglerophis_merremi,
                    3 Lystrophis_histricus,
                    4 Xenoxybelis_argenteus,
                    5 Liophis_jaegeri,
                    6 Liophis_elegantissimus,
                    7 Liophis_meridionalis,
                    8 Liophis_typhlus,
                    9 Erythrolamprus_aesculapii,
                    10 Atractus_albuquerquei,
                    11 Atractus_trihedrurus,
                    12 Heterodon_nasicus,
                    13 Sibynomorphus_mikanii,
                    14 Ninia_atrata,
                    15 Dipsas_indica,
                    16 Dipsas_neivai,
                    17 Sibynomorphus_garmani,
                    18 Hydrops_triangularis,
                    19 Helicops_infrataeniatus,
                    20 Pseudoeryx_plicatilis,
                    21 Helicops_pictiventris,
                    22 Helicops_gomesi,
                    23 Helicops_angulatus,
                    24 Oxyrhopus_rhombifer,
                    25 Oxyrhopus_clathratus,
                    26 Lycognathophis_seychellensis,
                    27 Taeniophallus_brevirostris,
                    28 Hydrodynastes_gigas,
                    29 Hydrodynastes_bicinctus,
                    30 Pseudoboa_nigra,
                    31 Pseudoboa_coronata,
                    32 Boiruna_maculata,
                    33 Drepanoides_anomalus,
                    34 Clelia_bicolor,
                    35 Phimophis_guerini,
                    36 Psomophis_joberti,
                    37 Psomophis_genimaculatus,
                    38 Leptodeira_annulata,
                    39 Thamnodynastes_rutilus,
                    40 Pseudotomodon_trigonatus,
                    41 Tachymenis_peruviana,
                    42 Tomodon_dorsatus,
                    43 Ptychophis_flavovirgatus,
                    44 Calamodontophis_paucidens,
                    45 Pseudablabes_agassizi,
                    46 Philodryas_patagoniensis,
                    47 Philodryas_aestiva,
                    48 Philodryas_mattogrossensis,
                    49 Taeniophallus_affinis,
                    50 Liophis_amarali,
                    51 Imantodes_cenchoa,
                    52 Apostolepis_dimidiata,
                    53 Phalotris_nasutus,
                    54 Phalotris_lemniscatus,
                    55 Apostolepis_assimilis,
                    56 Elapomorphus_quinquelineatus,
                    57 Siphlophis_pulcher,
                    58 Siphlophis_compressus,
                    59 Dendroaspis_polylepis,
                    60 Bothrolycus_ater,
                    61 Bothrophthalmus_lineatus,
                    62 Bothrophthalmus_brunneus,
                    63 Lamprophis_virgatus,
                    64 Lamprophis_fuliginosus_mentalis,
                    65 Lamprophis_lineatus,
                    66 Lamprophis_fuliginosus,
                    67 Lamprophis_olivaceus,
                    68 Lamprophis_capensis,
                    69 Lamprophis_inornatus,
                    70 Lycodonomorphus_whytii,
                    71 Lamprophis_fiskii,
                    72 Lycodonomorphus_rufulus,
                    73 Lamprophis_guttatus,
                    74 Homoroselaps_lacteus,
                    75 Atractaspis_bibronii,
                    76 Atractaspis_corpulenta,
                    77 Atractaspis_boulengeri,
                    78 Atractaspis_micropholis,
                    79 Aparallactus_modestus,
                    80 Polemon_collaris,
                    81 Polemon_acanthias,
                    82 Polemon_notatum,
                    83 Xenocalamus_transvaalensis,
                    84 Amblyodipsas_polylepis,
                    85 Macrelaps_microlepidotus,
                    86 Mehelya_nyassae,
                    87 Psammophis_phillipsi,
                    88 Psammophis_praeornatus,
                    89 Psammophis_sp.,
                    90 Psammophis_schokari,
                    91 Prosymna_janii,
                    92 Dromicodryas_bernieri,
                    93 Buhoma_depressiceps,
                    94 Liophidium_chabaudi,
                    95 Heteroliodon_occipitalis,
                    96 Lycophidion_capense,
                    97 Rhamphiophis_oxyrhynchus,
                    98 Hemirhagerrhis_hildebrandtii,
                    99 Amplorhinus_multimaculatus,
                    100 Lycophidion_laterale,
                    101 Buhoma_procterae,
                    102 Mehelya_stenophthalmus,
                    103 Psammodynastes_sp.,
                    104 Rhamphiophis_rostratus,
                    105 Pseudaspis_cana,
                    106 Gonionotophis_brussauxi,
                    107 Psammophylax_rhombeatus,
                    108 Psammophylax_variabilis,
                    109 Mehelya_poensis,
                    110 Duberria_variegata,
                    111 Duberria_lutrix,
                    112 Stenophis_citrinus,
                    113 Lycodryas_sanctijohannis,
                    114 Ithycyphus_oursi;
                TREE 'Tree # 1 simulated by Uniform speciation (Yule)' = (((((39:2.5803913536185727,(((52:0.2357864240123999,(102:0.08282509658987092,93:0.08282509658987092):0.152961327422529):0.2214802459365109,24:0.4572666699489108):1.945893456645554,84:2.4031601265944644):0.17723122702410854):0.4908444328116085,11:3.0712357864301825):2.5923132893962073,(63:1.2307914403400522,75:1.2307914403400522):4.432757635486336):1.0200745430701508,((((12:1.0579747675353615,(86:0.15378232386940346,91:0.15378232386940346):0.9041924436659581):0.5743314663735045,((88:0.006294520924597713,72:0.006294520924597713):0.12468858138623497,13:0.13098310231083268):1.5013231315980335):2.1317992308430798,(14:2.3758552182409325,(104:1.293245787183993,5:1.293245787183993):1.0826094310569394):1.3882502465110136):2.88700689532952,((50:4.636765296575451,(47:2.837093884632398,(97:0.8651570723310044,(28:0.19139889193677923,55:0.19139889193677923):0.6737581803942251):1.9719368123013938):1.7996714119430528):2.009001819244041,(101:4.835318430733998,((((68:2.0848955335180084,(31:1.5716404803657122,8:1.5716404803657122):0.5132550531522961):0.7955041381082514,(4:1.3077858310399064,(100:0.8873985846072808,70:0.8873985846072808):0.4203872464326258):1.5726138405863532):0.6586501189330831,(((26:1.891067590912069,92:1.891067590912069):1.1149435576870839,(82:0.25746327691214477,87:0.25746327691214477):2.74854787168701):0.4356397934181371,(112:2.495190552072538,((99:0.17106572727158684,64:0.17106572727158684):1.3867287730143514,((10:0.2824094492118545,6:0.2824094492118545):0.28972985119081,(53:0.5196727389007325,73:0.5196727389007325):0.05246656150193206):0.9856551998832738):0.9373960517866001):0.9464603899447522):0.0973988485420526):1.2867452337992613,((((29:1.7215589086054652,45:1.7215589086054652):0.11544685579407793,((27:0.11556314744938696,32:0.11556314744938696):1.5140987386161997,(89:0.5602461951446934,60:0.5602461951446934):1.0694156909208932):0.2073438783339571):1.1821218778042855,56:3.0191276422038293):1.8013909015341303,((((80:0.4159012620467458,38:0.4159012620467458):0.9455586524540294,98:1.3614599145007749):0.21760713188022085,106:1.5790670463809957):0.8297501576538641,(109:1.1835562005995945,(16:0.1966764138290763,54:0.1966764138290763):0.9868797867705182):1.2252610034352658):2.4117013397030984):0.005276480620645585):0.009523406375393796):1.8104486850854935):0.005345244261973421):0.03251125881507414):3.3163763811034594,((((((((61:0.29440964113392615,69:0.29440964113392615):0.3615622076311924,22:0.6559718487651186):0.28053922627863004,90:0.9365110750437485):3.989121139718452,(((18:0.7715314621196668,57:0.7715314621196668):1.003961873499669,((67:1.418732654788288,(51:0.1607722077267308,46:0.1607722077267308):1.2579604470615573):0.11556644564416302,19:1.5342991004324509):0.24119423518688496):3.140439366954856,((1:1.1323042044192035,(21:0.8629264757424177,85:0.8629264757424177):0.26937772867678555):3.4633382314176044,(17:0.8225231007965416,78:0.8225231007965416):3.7731193350402665):0.3202902667373836):0.00969951218800986):1.2914460263151295,((((((3:0.020875692978952264,7:0.020875692978952264):2.1217587452008093,((95:0.3987219008076808,37:0.3987219008076808):1.1413282820743367,103:1.5400501828820172):0.602584255297744):1.3067037265101777,(81:3.4483043230149955,74:3.4483043230149955):0.0010338416749450605):0.431766610328737,(((114:0.08486282510842975,71:0.08486282510842975):0.11900197549055996,36:0.2038648005989897):2.6070280983562655,25:2.8108928989552546):1.070211876063421):0.27680725078377577,(41:0.5825650296382654,43:0.5825650296382654):3.5753469961641855):0.5836966091466412,(110:0.4788832122520656,76:0.4788832122520656):4.262725422697028):1.475469606128238):2.030370458504919,96:8.247448699582245):0.8644549455025734,(((77:3.59013030665239,((65:0.9047126423132336,2:0.9047126423132336):0.4248151716462242,49:1.3295278139594575):2.260602492692932):0.37693314788474425,(44:2.0510372399763175,(62:1.5500434045105569,(94:0.07016826986024709,42:0.07016826986024709):1.47987513465031):0.5009938354657605):1.9160262145608158):0.2776031998200735,(105:2.459096378917186,111:2.459096378917186):1.7855702754400207):4.867236990727616):0.6510634659297476,((((30:1.9179447691748817,((((20:0.5451740985817892,23:0.5451740985817892):0.3462820991359227,107:0.8914561977177119):0.32762008624409256,(108:0.1452226939665266,40:0.1452226939665266):1.073853589995278):0.189861452181888,58:1.408937736143692):0.5090070330311898):0.3045649989769858,(15:1.4256139877331226,66:1.4256139877331226):0.7968957804187454):3.1573702908983234,((83:2.510008469276216,(35:1.2005356875244955,9:1.2005356875244955):1.3094727817517215):0.42604935479417244,(113:0.5206752711837962,33:0.5206752711837962):2.415382552886594):2.443822234979802):3.490716408791597,((59:1.4376801579218637,79:1.4376801579218637):2.990193278827107,(34:0.6175676472858066,48:0.6175676472858066):3.8103057894631642):4.442723031092817):0.8923706431727824):0.2370328889854287);
                TREE 'Tree # 2 simulated by Uniform speciation (Yule)' = ((((((((37:0.8888284529991368,90:0.8888284529991368):0.6179845460761253,((14:0.15636412648927278,(39:0.1466713801724733,(49:0.11690699973443669,52:0.11690699973443669):0.029764380438036597):0.00969274631679949):1.0415233834134567,24:1.1978875099027304):0.3089254891725323):1.948383113025092,(5:1.0879808217341231,(92:0.0028865657703142497,79:0.0028865657703142497):1.085094255963809):2.367215290366231):0.41957316158538777,((((68:1.3862889621309114,((114:0.02657053988394304,94:0.02657053988394304):1.1428861759253193,33:1.1694567158092624):0.21683224632164863):0.3442793379034314,47:1.730568300034343):0.15830357302397058,((25:0.632207415184395,61:0.632207415184395):0.5766974391716106,((113:0.2242826258945202,(70:0.02348381951083763,93:0.02348381951083763):0.20079880638368255):0.14764210480175918,87:0.37192473069627924):0.8369801236597264):0.679967018702307):0.5752037933357892,(29:0.316652669104458,80:0.316652669104458):2.147422997289643):1.4106936072916396):0.8648168516922232,(((91:0.06429624910226073,8:0.06429624910226073):0.7047100274204471,(111:0.21276772343696312,82:0.21276772343696312):0.5562385530857447):1.0504388047915683,(42:1.253941407708918,(35:0.4304844929893645,54:0.4304844929893645):0.8234569147195527):0.5655036736053586):2.920141044063689):1.5874793368447784,(13:6.301055161934009,((99:3.2078324932195503,(107:1.0425346163444287,(22:0.20396588907793925,(19:0.10503480550833101,20:0.10503480550833101):0.09893108356960821):0.8385687272664891):2.165297876875123):2.7208331691433987,(21:2.621083718976404,((31:0.8591250177291859,50:0.8591250177291859):0.36538516211842287,58:1.2245101798476092):1.3965735391287966):3.3075819433865448):0.37238949957105977):0.026010300288734155):0.9700832668604472,((62:2.0941754643952475,(((36:1.3660529709566918,((76:0.417426133312241,112:0.417426133312241):0.6628608688280649,9:1.0802870021403064):0.285765968816385):0.09329755326555632,(40:1.3575595578244244,(109:0.7140729948922789,(81:0.14711502885691005,(59:0.10683852914506871,(15:0.04146173578105493,95:0.04146173578105493):0.0653767933640138):0.0402764997118413):0.566957966035369):0.6434865629321449):0.10179096639782385):0.3324479480976391,((106:0.697373205717216,53:0.697373205717216):0.7792481657967958,((28:0.048117913828986766,34:0.048117913828986766):1.3134406001119145,75:1.3615585139409014):0.11506285757311052):0.31517710080587513):0.30237699207536123):3.9726564853505484,((64:0.25603881286514946,4:0.25603881286514946):5.800182892249088,(((108:1.2833013830661986,(65:0.24376614099400176,(98:0.15459044737166538,67:0.15459044737166538):0.08917569362233642):1.0395352420721962):0.6660526017318172,(((96:0.5894998846264646,44:0.5894998846264646):0.020763083992489363,38:0.610262968618954):0.657476763143336,103:1.2677397317622905):0.6816142530357252):0.5743725102783195,(((56:0.17787678711967408,6:0.17787678711967408):0.83705786945194,((26:0.10760242801545113,1:0.10760242801545113):0.4765936197285708,((46:0.042083660311608206,(83:0.02298266851554873,30:0.02298266851554873):0.019100991796059476):0.1780326088932885,27:0.22011626920489674):0.3640797785391252):0.4307386088275923):1.067313782357906,((101:0.6520560738191765,(104:0.23845976966589805,57:0.23845976966589805):0.41359630415327847):0.5243495592435838,43:1.176405633062761):0.90584280586676):0.4414780561468143):3.532495210037903):0.010610244631558368):1.230316779337394):2.7028512709168075,((((41:2.2510998624245264,((2:1.4882408922452242,10:1.4882408922452242):0.7607188746412992,((66:1.8052984197318778,(105:1.4818912845562495,(60:0.5068000366194512,(45:0.17396805505502094,18:0.17396805505502094):0.3328319815644304):0.9750912479367977):0.32340713517562863):0.2702753099096074,((73:0.675141932654651,(((12:0.31168638400657395,3:0.31168638400657395):0.31908013460542983,7:0.6307665186120038):0.02074994742856444,74:0.6515164660405683):0.023625466614082442):0.7694869264002483,(89:0.5706943906632249,11:0.5706943906632249):0.8739344683916743):0.6309448705865859):0.17338603724503776):0.002140095538004554):0.47533898125935287,86:2.7264388436838782):0.2608566305493007,(17:1.1340964316697055,100:1.1340964316697055):1.853199042563476):1.352341075424337,(((51:0.52226272333101,(16:0.2353706237113944,84:0.2353706237113944):0.28689209961961576):2.8967876264492256,78:3.419050349780235):0.43341356126819075,((((88:0.7232099347469628,110:0.7232099347469628):2.22176282193872,((77:0.11718331436936599,32:0.11718331436936599):0.7077594908142963,(97:0.013841635782994697,63:0.013841635782994697):0.8111011694006677):2.120029951502021):0.41190810979214215,((((55:0.3823494768836099,102:0.3823494768836099):0.0887039417601812,85:0.471053418643791):1.5686585956254357,48:2.0397120142692264):0.6641349855763554,((71:0.3430299369698703,23:0.3430299369698703):0.09833316627643436,69:0.44136310324630473):2.262483896599276):0.6530338666322432):0.034751948961034514,72:3.3916328154388586):0.46083109560956764):0.48717263860909044):5.660363450342481);
                TREE 'Tree # 3 simulated by Uniform speciation (Yule)' = ((((((((76:2.348287103354214,(107:0.24319850581253818,41:0.24319850581253818):2.1050885975416764):1.054477032777944,(((100:0.19425842996623838,67:0.19425842996623838):1.3547592927277783,84:1.5490177226940165):1.200660144721675,(((10:0.04380740269560295,19:0.04380740269560295):0.2858267627712927,97:0.3296341654668956):1.4868926180617423,((71:1.2175071010416778,((15:0.11307695834004278,81:0.11307695834004278):0.11287377613775015,62:0.22595073447779299):0.9915563665638849):0.23349763234156878,(24:0.3961705750577298,80:0.3961705750577298):1.0548341583255176):0.36552205014539113):0.9331510838870531):0.6530862687164684):1.599611982516063,(103:1.8332901508503316,((27:0.5757063846616234,109:0.5757063846616234):0.7934696370935191,(59:1.011020629836568,39:1.011020629836568):0.35815539191857393):0.4641141290951894):3.1690859677978906):0.896690396519931,((((94:1.0471300840857205,64:1.0471300840857205):1.0073164084470059,((8:0.09436824549805703,16:0.09436824549805703):0.1819711478060492,98:0.27633939330410623):1.77810709922862):1.7398766105751582,(35:1.6805061028632629,((46:0.5614913827210902,(50:0.22264075336059558,113:0.22264075336059558):0.3388506293604946):0.778959500190291,112:1.340450882911381):0.34005521995188215):2.1138170002446213):0.5718962640892663,(((((((17:0.4341772548641502,53:0.4341772548641502):0.15772181261413196,49:0.5918990674782824):0.04770585347771558,86:0.639604920955998):0.943507258575378,(87:0.7510621525256367,34:0.7510621525256367):0.8320500270057394):0.04869113697983617,48:1.6318033165112118):0.5030227487363705,((63:0.07530402663572705,104:0.07530402663572705):1.3613194663009278,23:1.4366234929366548):0.6982025723109276):1.7556251756160361,(((1:0.5000247099741794,7:0.5000247099741794):2.2187248631623886,(54:1.3482094775754916,(13:1.0424804827554106,101:1.0424804827554106):0.305728994820081):1.3705400955610767):0.4778919463626469,((73:0.4078730209381425,37:0.4078730209381425):0.821628071574718,(11:0.8009122398760583,52:0.8009122398760583):0.42858885263680224):1.9671404269863553):0.6938097213644032):0.475768126333532):1.5328471479710029):0.6125227345970107,((110:1.1501692688011096,((92:0.06560513795986168,96:0.06560513795986168):1.084529880655818,105:1.1501350186156793):3.425018543003032E-5):0.10760818966311793,(85:0.3528827711956035,79:0.3528827711956035):0.9048946872686242):5.253811791300936):1.0223805578417418,((74:2.7720959177924143,((72:0.02139111994988565,9:0.02139111994988565):2.1645803457452883,((26:0.002358493876065758,61:0.002358493876065758):0.3076104138472471,36:0.30996890772331287):1.8760025579718616):0.5861244520972411):0.13584173613229655,(21:1.0169260434657317,((78:0.23539193195712593,83:0.23539193195712593):0.15399630367368242,6:0.3893882356308084):0.6275378078349234):1.8910116104589807):4.626032153682193):1.0659856198524016,((((57:3.045088467726223,((20:0.7249013499668998,4:0.7249013499668998):0.1265552413486197,3:0.8514565913155195):2.1936318764107043):1.4042870993405587,((77:0.678957019126178,(95:0.18976147871035706,18:0.18976147871035706):0.48919554041582086):1.322982425762536,(93:2.0009469805354176,29:2.0009469805354176):9.92464353296332E-4):2.4474361221780683):0.2438202921657327,(60:3.2180494269667053,((42:0.3624829431578298,75:0.3624829431578298):0.10998407379020719,90:0.472467016948037):2.745582410018669):1.475146432265809):2.5694173643647407,((((2:1.7166648335588965,99:1.7166648335588965):0.1867910719941611,((((33:0.763361685950789,88:0.763361685950789):0.04285180535598959,45:0.8062134913067786):0.18963647303336356,69:0.995849964340142):0.8174911479733239,((47:0.8784824233927232,25:0.8784824233927232):0.4522513436452115,((32:0.2998688294453695,58:0.2998688294453695):0.33158450730493083,28:0.6314533367503006):0.6992804302876344):0.4826073452755312):0.09011479323959182):0.3336384774035872,12:2.237094382956644):1.6767613321347081,((111:2.1962217172481253,(51:0.7806356606134565,31:0.7806356606134565):1.4155860566346699):1.716757510632607,43:3.9129792278807325):8.764872106199858E-4):3.3487575085059027):1.337342203862052):1.4000445725406905,((((91:0.12305275796633904,106:0.12305275796633904):0.7686518055296836,89:0.8917045634960228):3.160106709349426,(((40:1.108668737567735,114:1.108668737567735):0.34988141067631257,5:1.4585501482440477):0.04413372325066822,(68:1.5025363992472711,66:1.5025363992472711):1.4747224744508685E-4):2.549127401350733):1.5228864568162637,(102:4.958835091165425,(((108:1.3790356467842662,(38:0.2647014122538049,30:0.2647014122538049):1.1143342345304612):1.8946242950296461,((44:0.5331347280795637,14:0.5331347280795637):0.28875960109064336,(((82:0.7201014640635851,56:0.7201014640635851):0.002412651099174804,70:0.72251411516276):0.05710954642758855,(22:0.29507318194558857,65:0.29507318194558857):0.48455047964475967):0.04227066757985867):2.451765612643705):0.06314790525591132,55:3.336807847069822):1.6220272440956005):0.6158626384962894):4.425302270338285);
                TREE 'Tree # 4 simulated by Uniform speciation (Yule)' = (((((((((11:1.211588771200587,104:1.211588771200587):0.8405984569478285,(40:1.5360923006897762,(19:1.1004340829605928,44:1.1004340829605928):0.4356582177291836):0.5160949274586394):0.2437520535366676,36:2.2959392816850834):0.7631241869309465,(102:0.07668230278388447,20:0.07668230278388447):2.9823811658321446):0.926164639986822,(18:3.909805581828459,(((38:0.24099911230362706,72:0.24099911230362706):0.8981260924470523,103:1.139125204750679):2.4438060912547104,(((51:1.5095323919880959,(89:0.9870127582821163,(52:0.9440474169700938,(69:0.48472913189355,(12:0.3192593768166802,43:0.3192593768166802):0.16546975507686978):0.4593182850765437):0.042965341312022474):0.5225196337059796):0.15072216868269847,45:1.6602545606707944):1.206305510628824,(85:0.32126906424875656,26:0.32126906424875656):2.545291007050861):0.716371224705771):0.3268742858230702):0.0754225267743921):1.4332820019289487,(15:2.511228581962446,(53:0.8255631410175753,27:0.8255631410175753):1.6856654409448708):2.9072815285693525):2.0540819616793637,((100:1.7541275972626111,30:1.7541275972626111):4.44668551925303,(((54:0.2720585924634813,64:0.2720585924634813):0.2622783950893201,4:0.5343369875528013):1.1509555118583596,(92:0.3438813245087821,108:0.3438813245087821):1.3414111749023796):4.515520617104479):1.2717789556955235):1.186148714820006,((((((57:0.104391056799186,90:0.104391056799186):0.8950967063357733,39:0.9994877631349594):3.3611439176012072,(((99:0.17747534523603226,35:0.17747534523603226):0.7674163791218436,59:0.9448917243578758):0.6749586337288728,(97:0.7191211844518219,91:0.7191211844518219):0.9007291736349265):2.7407813226494175):1.5100486304084466,(((24:0.776106227741738,50:0.776106227741738):0.6741963552689836,(23:0.4931010580490366,14:0.4931010580490366):0.9572015249616849):2.3980274114435898,(94:2.7729787349688566,47:2.7729787349688566):1.0753512594854544):2.022350316690302):0.46581388208344804,(1:6.1603471811292065,((5:2.6712984847453574,(106:1.1070173221071977,80:1.1070173221071977):1.5642811626381605):1.2443429014014484,(((98:0.6552856061423131,21:0.6552856061423131):2.0854001540761486,((65:0.4275420679039235,74:0.4275420679039235):0.7692361866407648,79:1.1967782545446883):1.543907505673774):0.20110487130074028,((28:0.028276831495780162,71:0.028276831495780162):1.401043529815148,(((34:0.44912492779434654,58:0.44912492779434654):0.5449514804100969,(77:0.9149246211486786,55:0.9149246211486786):0.07915178705576495):0.1421427673640798,(37:0.984178465031758,33:0.984178465031758):0.1520407105367654):0.2931011857424051):1.5124702702082744):0.9738507546276037):2.244705794982402):0.1761470120988528):2.1614822399832447,((((((16:1.8577998782437375,48:1.8577998782437375):1.8145237787132371,(32:0.8413619971588142,(17:0.8412788637238218,101:0.8412788637238218):8.313343499252645E-5):2.8309616597981604):0.8976689545356605,((2:0.7366545791938748,((75:0.49025315569616906,66:0.49025315569616906):0.22822824978759143,107:0.7184814054837604):0.018173173710114395):2.861072703167212,((63:0.6894077689683678,(56:0.2355587012495178,78:0.2355587012495178):0.45384906771884986):2.1591706919079865,96:2.8485784608763534):0.7491488214847332):0.9722653291315475):0.34557742140288733,(((109:0.8726342594534149,49:0.8726342594534149):1.9270239791531198,86:2.799658238606534):0.16776606773745026,(((111:0.004871775392833259,6:0.004871775392833259):0.4552590419128098,113:0.46013081730564304):0.7135272880483688,(67:0.5862802054745229,22:0.5862802054745229):0.587377899879489):1.793766200989973):1.9481457265515367):0.1847957075049406,((29:2.893619551950081,83:2.893619551950081):0.11878953897765755,(42:0.21109357163672565,61:0.21109357163672565):2.801315519291012):2.0879566494727237):2.2149477428028397,((((93:2.57206673722916,((76:1.1219453640143149,(70:0.7307531084327519,82:0.7307531084327519):0.39119225558156306):1.449182422688505,73:2.5711277867028195):9.389505263406936E-4):2.202571047146109,7:4.774637784375269):0.5214465864538722,(8:3.6265009237151604,((110:0.5043656432005792,60:0.5043656432005792):0.8946913186976839,9:1.3990569618982631):2.2274439618168977):1.669583447113982):1.3291980073379228,((13:1.5418067693415118,(10:0.41429749491078316,105:0.41429749491078316):1.1275092744307287):0.4295877394621273,62:1.971394508803639):4.653887869363427):0.6900311050362367):1.1826629500080035):0.16076435381986487):1.341259212968828,((81:3.3659624747447996,(((112:0.0510243681698741,68:0.0510243681698741):0.8659728075714744,114:0.9169971757413485):1.9653209291756122,87:2.8823181049169597):0.4836443698278396):3.3491505673099016,((46:0.35827697171273576,41:0.35827697171273576):0.23650395891039916,(95:0.44707533762270507,(88:0.34293376790159813,((31:0.10265584393319595,25:0.10265584393319595):0.13736271577427397,(3:0.10177429190362695,84:0.10177429190362695):0.13824426780384297):0.1029152081941282):0.10414156972110698):0.14770559300042976):6.120332111431566):3.2848869579452975);
                TREE 'Tree # 5 simulated by Uniform speciation (Yule)' = ((((22:1.0194700467071671,49:1.0194700467071671):0.3079902339115988,55:1.3274602806187659):4.557710901824335,((((((114:1.8716972058303185,48:1.8716972058303185):1.5206382193213936,(1:0.6887168472148287,83:0.6887168472148287):2.7036185779368846):0.6885051193126961,((40:0.9852728127020866,27:0.9852728127020866):0.21100585982000547,(36:0.43670135091753326,99:0.43670135091753326):0.7595773216045587):2.884561871942317):0.26662337474691244,((51:0.8109544595615031,17:0.8109544595615031):2.668489890560083,(82:0.786025957770023,31:0.786025957770023):2.693418392351563):0.8680195690897355):0.08750708731270895,((94:0.8718693957655052,93:0.8718693957655052):2.9242327639673937,92:3.7961021597328997):0.6388688467911314):1.4380488586002005,((((43:0.7132813069965125,((68:0.4460789471610896,(15:0.10534508200255688,53:0.10534508200255688):0.34073386515853277):0.12322030287247977,46:0.5692992500335695):0.14398205696294306):3.1099882369185092,(65:2.5134602640283528,(84:0.6048021019455533,86:0.6048021019455533):1.9086581620827996):1.3098092798866687):1.7983202527946491,((69:3.23475511763795,104:3.23475511763795):1.6959601729201552,((32:2.833249264213957,60:2.833249264213957):1.3869438820148574,(((77:1.2911259410128317,8:1.2911259410128317):2.074532612195817,54:3.36565855320865):0.5370588590693343,(62:2.7429343328363,112:2.7429343328363):1.1597830794416821):0.3174757339508321):0.7105221443292897):0.6908745061515666):0.23929097733623553,105:5.860880774045907):0.01213909107832451):0.012151317318870406):4.114828817556897,(((((((13:0.27236833516755554,58:0.27236833516755554):1.7355827543248235,((3:0.11055579736277124,12:0.11055579736277124):0.5656934809477352,89:0.6762492783105065):1.3317018111818726):0.5273566227725699,78:2.5353077122649488):0.34869657827450673,38:2.884004290539456):4.365581317459739,(56:4.569768180009242,((14:2.4691262230290323,37:2.4691262230290323):0.8581874898134838,(25:0.9131617326654251,44:0.9131617326654251):2.4141519801770914):1.2424544671667255):2.6798174279899545):1.8548498490622323,((20:6.92764552408139,34:6.92764552408139):1.0480569911776998,(((80:8.718693512701108E-5,26:8.718693512701108E-5):6.692425640406069,(71:6.691606369509737,((((28:3.985103125572778,23:3.985103125572778):0.2664201201852206,(66:0.21549526738953065,10:0.21549526738953065):4.036027978368468):0.5872095168213411,102:4.838732762579338):0.12759484097429513,(50:2.7657055696974604,((109:0.45054555068432683,74:0.45054555068432683):0.5202606225427523,21:0.9708061732270794):1.794899396470381):2.2006220338561735):1.7252787659561035):9.06457831460134E-4):0.13128567866025484,((45:2.3645268056243185,((79:0.5543045823247924,100:0.5543045823247924):0.4011376412001895,41:0.9554422235249821):1.4090845820993363):2.357037571657857,(96:0.17335017214591172,90:0.17335017214591172):4.548214205136265):2.102234128719276):1.1519040092576383):1.128732941802338):0.8339684088076916,(((((85:2.199719892111033,73:2.199719892111033):1.3797219320518774,29:3.5794418241629113):1.7145019078341526,((97:0.5020178075754378,39:0.5020178075754378):3.666198828467737,((9:3.302392781420432,((((59:0.8502397575963401,24:0.8502397575963401):0.08443491565906196,(111:0.6725832278390843,76:0.6725832278390843):0.26209144541631757):0.7042592521890801,6:1.6389339254444812):0.8184147984858895,33:2.4573487239303717):0.8450440574900597):0.7745759978218745,(((91:0.4035658491578398,5:0.4035658491578398):0.16798546502497552,18:0.5715513141828155):2.9945267162405647,(((35:1.7700724241616534,(72:0.7767513582350571,95:0.7767513582350571):0.993321065926597):0.07921185699253852,64:1.8492842811541919):0.9696513443289988,(110:0.25210594071183934,67:0.25210594071183934):2.5668296847713523):0.7471424049401874):0.5108907488189268):0.09124785680086872):1.1257270959538888):2.93430446307928,((16:1.6559215229181359,61:1.6559215229181359):3.99073387540266,(((19:1.5753950893255524,(7:0.9664746195114635,107:0.9664746195114635):0.6089204698140891):2.342921399624504,((108:1.231159102141442,4:1.231159102141442):0.030160262813296517,((2:0.46966326188339663,((63:0.0750977211790196,81:0.0750977211790196):0.03280538187109167,47:0.10790310305011129):0.36176015883328533):0.009520736696578382,103:0.479183998579975):0.7821353663747631):2.656997123995318):0.5959980806788564,87:4.5143145696289135):1.132340828691884):2.581592796755546):0.4070058360927919,((42:3.081365679973879,(30:3.0810430343027444,52:3.0810430343027444):3.2264567113447047E-4):1.970436572588453,(((((75:0.6207992022752251,70:0.6207992022752251):1.1715055114533186,106:1.7923047137285433):0.15621027765805867,(101:0.48838325973734087,57:0.48838325973734087):1.4601317316492612):1.5587194982010875,(98:0.07848123137900034,113:0.07848123137900034):3.4287532582086904):1.536737458056622,(88:0.12665073875670196,11:0.12665073875670196):4.917321208887608):0.00783030491802021):3.583451778606802):1.3031498346999861):0.06159613413087751);
                TREE 'Tree # 6 simulated by Uniform speciation (Yule)' = ((((((99:3.20400723118286,(26:2.347615012710926,(28:1.8029415299400249,((104:0.2163425410905689,52:0.2163425410905689):1.3612841663090212,42:1.5776267073995904):0.22531482254043464):0.544673482770901):0.8563922184719335):1.3908013074727192,((35:1.1495527490517885,((27:0.3467481776704405,(61:0.23278323809006682,5:0.23278323809006682):0.1139649395803738):0.2677491378206175,110:0.614497315491058):0.5350554335607305):1.5579157477368017,(46:0.00493274409653687,21:0.00493274409653687):2.7025357526920537):1.8873400418669883):1.1302690100771693,23:5.725077548732747):0.6354817520380314,((50:0.5468725547283521,65:0.5468725547283521):0.8266606267769194,11:1.3735331815052716):4.987026119265508):0.7919042564370797,114:7.1524635572078585):2.8475364427921392,((((106:1.228350060916384,25:1.228350060916384):0.6620343440369881,(96:1.5705784429031773,((15:0.9207005087731518,111:0.9207005087731518):0.09510034311622627,24:1.0158008518893782):0.5547775910137985):0.3198059620501953):1.7968234992324104,(2:0.45899847153946655,(85:0.18686189148229787,91:0.18686189148229787):0.2721365800571688):3.228209432646317):6.016230966033873,((((((((58:0.355906831225317,100:0.355906831225317):1.7870190369084846,67:2.142925868133802):0.7422899133435755,(((81:0.4944412808383662,(77:0.23388990802486356,60:0.23388990802486356):0.2605513728135027):0.5275027725946112,12:1.0219440534329773):1.7336776539237713,71:2.755621707356748):0.1295940741206289):2.372051065049385,((((((8:1.582108300168188,87:1.582108300168188):0.8473642500999874,39:2.429472550268175):0.18290913255898267,((32:0.8984056065247499,83:0.8984056065247499):0.23041547162888906,(38:0.9285110192739813,(56:0.81846423998205,13:0.81846423998205):0.1100467792919313):0.20031005887965778):1.4835606046735188):0.21701963380196013,((97:1.4969558994759475,(33:0.6046208702175546,55:0.6046208702175546):0.8923350292583927):0.1412840179274488,(48:0.6218162125034119,(44:0.2356124893566791,34:0.2356124893566791):0.38620372314673285):1.0164237048999845):1.1911613992257215):0.27684284261618886,((86:1.1986005373123048,1:1.1986005373123048):0.5491497518980509,72:1.7477502892103556):1.358493870034951):0.04987147495352853,((36:1.4662503190662177,82:1.4662503190662177):0.3081668631887433,7:1.774417182254961):1.3816984519438742):2.101151212327927):2.51944050513608,109:7.776707351662842):0.3592321897531146,(92:1.6613172491601036,108:1.6613172491601036):6.474622292255855):0.8665562543256309,(((78:4.252890842253196,(113:0.6344795207600264,(6:0.634138516390291,(57:0.5294831054819263,68:0.5294831054819263):0.10465541090836458):3.410043697355444E-4):3.6184113214931695):2.7380673255198653,(((((76:0.07763534031455922,70:0.07763534031455922):4.098652465042318,((((49:0.6770788577888316,88:0.6770788577888316):1.1710533239972707,(20:0.13122092469774896,(41:0.004161503766527513,10:0.004161503766527513):0.12705942093122144):1.7169112570883533):0.4595043900683393,(9:2.1407653189011766,75:2.1407653189011766):0.16687125295326544):0.36288209472553495,54:2.670518666579977):1.505769138776899):0.3583208158148901,((4:1.521865866202031,74:1.521865866202031):1.2643315310630852,(69:2.371981311249843,(84:0.7177086480402061,31:0.7177086480402061):1.6542726632096365):0.4142160860152728):1.7484112239066496):0.6667756652429073,(((19:0.9311948019697905,112:0.9311948019697905):1.1552756361283298,18:2.086470438098121):2.0431604268566512,((63:1.811004580521817,((53:0.18951940688584876,94:0.18951940688584876):0.3416466645429463,89:0.531166071428795):1.2798385090930215):1.2726838563591045,(95:1.438646957601102,62:1.438646957601102):1.64504147927982):1.0459424280738505):1.071753421459901):0.4865016818089487,(101:0.8203593429578636,98:0.8203593429578636):4.86752662526576):1.303072199549438):1.7223831522723851,(((90:0.7373753825875102,107:0.7373753825875102):3.2548883921164853,(17:2.7621303573802365,30:2.7621303573802365):1.2301334173237586):1.5641299414980239,(((102:0.030489971967820977,105:0.030489971967820977):3.552615136361267,(64:0.9814184390063282,(47:0.5231425224087788,45:0.5231425224087788):0.45827591659754935):2.6016866693227594):0.053583596270966385,73:3.6366887046000533):1.9197050116019658):3.156947603843427):0.28915447569614317):0.5961035166050964,((3:0.25751530777667053,40:0.25751530777667053):9.244171506672995,((14:8.509049687455281,((29:1.0497368521204749,(51:0.20679474056272576,103:0.20679474056272576):0.8429421115577492):3.374069807922634,((66:2.385144497242506,(80:2.079013664440029,(43:0.23927519209465475,79:0.23927519209465475):1.8397384723453736):0.30613083280247744):1.776931885980389,(22:0.14915489516509192,93:0.14915489516509192):4.012921488057805):0.2617302768202133):4.085243027412173):0.3100052091366133,(37:2.4871057961446112,(16:0.8771158757427711,59:0.8771158757427711):1.6099899204018397):6.331949100447282):0.6826319178577688):0.09691249789702135):0.10483955787297079):0.29656112978034255);
                TREE 'Tree # 7 simulated by Uniform speciation (Yule)' = (((((((3:1.7773733499703177,(62:1.539740826923612,(29:0.18975587040928835,93:0.18975587040928835):1.3499849565143236):0.23763252304670612):1.8578006105242237,(((((12:0.2214646504946761,63:0.2214646504946761):1.3415486405254033,90:1.5630132910200794):0.6175050755392009,((10:0.6594674818818371,106:0.6594674818818371):0.11771978774145056,8:0.7771872696232879):1.4033310969359922):0.7739502605165907,(22:0.09462232563167076,110:0.09462232563167076):2.8598463014441977):0.6794646241904925,(92:3.234249574034215,1:3.234249574034215):0.3996836772321485):0.0012407092281784338):1.419471442296315,((26:0.553410322846109,52:0.553410322846109):2.438030473236479,((((98:0.0014353093394337,23:0.0014353093394337):1.3227535400031243,(76:0.44856466646808374,24:0.44856466646808374):0.8756241828744741):0.9002342834702994,38:2.224423132812857):0.7649710689682401,((20:1.2617084938727838,55:1.2617084938727838):1.2851885170617066,(101:2.546090721748062,91:2.546090721748062):8.062891864279749E-4):0.44249719084660705):0.0020465943014913055):2.063204606708268):1.5339562737135735,(68:0.86256137110357,104:0.86256137110357):5.726040305400862):1.1130099704142,(36:2.2470747183890367,42:2.2470747183890367):5.454536928529592):0.9337075536761591,((((44:0.674078920973472,80:0.674078920973472):3.6570619544166747,((81:0.3611493653703629,(27:0.20171014607397425,14:0.20171014607397425):0.15943921929638866):0.14464534542731686,100:0.5057947107976797):3.825346164592467):1.6982873238864977,(((84:2.405469893661435,((37:1.1601698311835118,43:1.1601698311835118):0.1908424212237685,65:1.3510122524072803):1.0544576412541549):0.3922439668363911,51:2.797713860497826):0.25929830235674994,((71:1.7503118800842759,96:1.7503118800842759):0.9942656867721443,97:2.74457756685642):0.31243459599815604):2.9724160364220675):2.4533590448775042,((((105:4.553623443150942,(((74:0.26513379147981314,82:0.26513379147981314):0.05361394684022307,(107:0.0446790340554512,19:0.0446790340554512):0.2740687042645851):2.3751749949283942,83:2.6939227332484306):1.8597007099025098):1.219914068096936,((((21:0.45697748468792054,50:0.45697748468792054):2.0255786839963594,((69:0.1612897474912156,18:0.1612897474912156):0.772911547196707,(112:0.2597106456042082,113:0.2597106456042082):0.6744906490837145):1.5483548739963577):0.07903088972891159,((58:0.14011189084038084,6:0.14011189084038084):2.0308252686231207,(28:0.7585651269349798,87:0.7585651269349798):1.412372032528522):0.39064989894969054):0.6924052039422974,((16:2.1329320799076745,(79:0.02844165610418828,53:0.02844165610418828):2.1044904238034863):0.23697680188239423,(54:0.2032992706037837,9:0.2032992706037837):2.166609611186285):0.8840833805654207):2.519545248892386):0.4731604437463032,(((((46:0.6010196856050586,111:0.6010196856050586):0.4935384115139135,33:1.0945580971189721):0.520105209414153,(59:0.12639808896208396,13:0.12639808896208396):1.4882652175710414):4.347816031979979,((85:1.688416757495176,(56:1.3643108330284313,66:1.3643108330284313):0.3241059244667446):1.4775578319281355,(77:3.050381259246752,102:3.050381259246752):0.11559333017655916):2.7965047490897916):0.036091735782565935,((17:2.623761302515578,72:2.623761302515578):0.25572438582285634,(40:2.7947353782894933,((31:1.1753930100988474,(11:0.547742563518091,75:0.547742563518091):0.6276504465807564):1.1808339084860797,((7:1.1981457216251739,(25:0.7725283046298715,4:0.7725283046298715):0.4256174169953025):0.4079531550885891,89:1.6060988767137634):0.7501280418711641):0.4385084597045668):0.08475031004894047):3.1190853859572343):0.2481268806985106):0.167158799683618,49:6.4138567546777985):2.06893048947635):0.1525319564406414):1.3646807994052095,(((114:3.3081173311971432,((47:2.7602996209100406,((45:1.0726421096619212,70:1.0726421096619212):0.4714158525962284,(39:1.0603910289240819,(108:0.08682985503063752,94:0.08682985503063752):0.9735611738934442):0.48366693333406807):1.2162416586518907):0.48714632860074253,((48:1.2874635833506116,67:1.2874635833506116):1.941035726125457,88:3.2284993094760686):0.018946640034714073):0.060671381686360534):1.4227337710362689,(((99:0.40328607733351096,73:0.40328607733351096):0.808416706037634,(30:0.18058782092283687,(2:0.08541909158517104,78:0.08541909158517104):0.09516872933766587):1.031114962448308):0.15711771658141738,(86:1.0673091027757189,35:1.0673091027757189):0.3015113971768435):3.3620306022808504):5.087079193128547,(((((60:0.9515012344807722,(57:0.9336059473994116,61:0.9336059473994116):0.017895287081360876):0.46976462731112295,(34:0.30804691710706755,64:0.30804691710706755):1.1132189446848277):0.5671213521382718,103:1.9883872139301664):2.0914833189625064,((41:0.04800560084319876,32:0.04800560084319876):3.0871867652724014,95:3.135192366115601):0.9446781667770714):3.3043893815054477,(15:0.7910473285947146,(5:0.12037544235371771,109:0.12037544235371771):0.6706718862409967):6.593212585803408):2.4336703809638394):0.18206970463803834);
                TREE 'Tree # 8 simulated by Uniform speciation (Yule)' = (((72:5.2204150525713615,(((((75:0.41024724910081306,2:0.41024724910081306):0.21121128336290088,12:0.6214585324637139):1.373859818495205,((42:0.2395139073677459,43:0.2395139073677459):1.3029150384539345,87:1.5424289458216804):0.4528894051372389):1.2172074665655173,((((6:0.19708734085428462,98:0.19708734085428462):0.4823515106530021,65:0.6794388515072868):1.2888918746514173,((77:0.02071060880036789,103:0.02071060880036789):0.23455078632440346,111:0.25526139512477136):1.7130693310339324):0.5169310111770848,(114:2.409157553785707,(31:2.3659351071092725,94:2.3659351071092725):0.04322244667643495):0.07610418355008125):0.7272640801886475):1.3622973534098706,(((40:0.07757407952634972,17:0.07757407952634972):1.5747704784230845,(71:1.4404758665279596,53:1.4404758665279596):0.21186869142147485):2.539600048192359,(((55:1.552773199336903,59:1.552773199336903):1.7458492257586218,((7:0.7097009637413992,52:0.7097009637413992):0.39148662996396993,((112:0.06790986369838425,69:0.06790986369838425):0.3983996076636445,44:0.4663094713620287):0.6348781223433406):2.1974348313901553):0.28501110211890135,((68:1.781398979748029,((28:1.4089055178046381,39:1.4089055178046381):0.14067759543741568,(99:0.6121198947266245,1:0.6121198947266245):0.9374632185154292):0.23181586650597527):0.033820795271214205,(101:1.400052196178604,9:1.400052196178604):0.4151675788406393):1.7684137521951833):0.6083110789273662):0.38287856479251436):0.6455918816370526):2.0360068625212278,((((5:2.790246276566559,((((96:0.04149593865160155,91:0.04149593865160155):0.44669461373171765,(78:0.03671654862515874,11:0.03671654862515874):0.45147400375816044):0.9949730264593303,(23:0.009356417391800056,48:0.009356417391800056):1.4738071614508494):1.2116428820826373,((102:1.4591724889978221,(97:1.416573403617036,(80:0.30825650159443585,89:0.30825650159443585):1.1083169020225998):0.04259908538078606):0.06059479456656252,107:1.519767283564385):1.175039177360902):0.09543981564127302):0.5916749498830455,46:3.3819212264496064):0.8412333291174724,(30:2.539133046979453,83:2.539133046979453):1.6840215085876236):1.9948585727624735,(54:0.5154445942170504,(66:0.003989380811657699,56:0.003989380811657699):0.5114552134053926):5.702568534112497):1.0384087867630363):2.7435780849074107,(((((92:0.9720623177999801,14:0.9720623177999801):0.16071463074956455,(10:0.12021884727807847,35:0.12021884727807847):1.012558101271466):2.534758194972247,(26:0.3747914028516058,84:0.3747914028516058):3.292743740670187):2.795804233265915,((((16:1.250170857351656,93:1.250170857351656):3.3934557800638556,(((70:1.5960932392631475,(27:0.29612274494583113,(73:0.07568066001825076,45:0.07568066001825076):0.22044208492758038):1.299970494317316):1.7299300718969814,(85:0.7408561735916854,58:0.7408561735916854):2.5851671375684435):1.238169647493784,((95:2.2514242222458862,((41:1.9620589115555918,(33:0.3161482141345429,8:0.3161482141345429):1.6459106974210485):0.1099348664926264,34:2.071993778048218):0.1794304441976686):0.7392787740652391,(88:2.3464784763098128,21:2.3464784763098128):0.6442245200013127):1.5734899623427865):0.07943367876159958):0.2644880984478119,(((((81:2.075664545394056,76:2.075664545394056):0.10834203828423654,61:2.1840065836782925):0.22802741002573648,109:2.412033993704029):1.1935489560416055,(64:0.12789657063726015,82:0.12789657063726015):3.4776863791083756):0.5223841302238662,(((((90:0.16100408449389206,105:0.16100408449389206):0.16950051438947994,100:0.33050459888337197):0.6020698569559035,36:0.9325744558392753):1.372341382546282,(79:1.4456077727402985,25:1.4456077727402985):0.8593080656452594):0.06353743743275123,(86:0.9140060031447576,15:0.9140060031447576):1.4544472726735511):1.759513804151192):0.7801476558938232):0.893726653684549,((29:1.7065927112806227,(50:0.22607879475310505,113:0.22607879475310505):1.4805139165275174):0.9917134922015391,((63:0.27043831094723675,62:0.27043831094723675):0.3154004967020743,(57:0.48736524824157657,104:0.48736524824157657):0.09847355940773458):2.1124673958328506):3.1035351860657103):0.661497987239834):0.3990951718594348,(((((110:0.8258144649394806,49:0.8258144649394806):1.0240066423740142,74:1.849821107313495):1.1912182267976612,60:3.0410393341111566):0.9034048487442257,((22:2.1444909165500703,((67:0.4709429005237257,51:0.4709429005237257):1.0383285558607818,38:1.5092714563845078):0.6352194601655633):0.358456330723591,106:2.502947247273661):1.4414969355817202):2.018371147111596,((24:5.11021152360978,(19:1.3382368907945168,(47:0.3683135133091489,3:0.3683135133091489):0.9699233774853678):3.771974632815261):0.8481431768418217,((108:4.45568111260699,(13:2.648048357269368,4:2.648048357269368):1.8076327553376206):0.36721939779295426,((18:0.7733866209660952,20:0.7733866209660952):1.810125465202255,(37:1.4185034712495082,32:1.4185034712495082):1.1650086149188423):2.239388424231593):1.1354541900516573):0.004460629515376992):0.899619218680164):3.137565451352856);
                TREE 'Tree # 9 simulated by Uniform speciation (Yule)' = (((((23:1.017130945067286,75:1.017130945067286):0.2763239962627771,(97:0.7383434218589645,78:0.7383434218589645):0.5551115194710984):2.6570800130917482,6:3.950534954421812):3.794191164218941,((40:2.6953763743189194,(64:1.6919394933036336,(((56:0.3519009971441078,32:0.3519009971441078):1.1668614989009876,8:1.5187624960450947):0.17315583528962003,(20:1.6912941038509621,(95:0.7847160568197388,61:0.7847160568197388):0.9065780470312247):6.242274837523716E-4):2.116196891912303E-5):1.003436881015285):4.75074939932664,((((99:0.22558375570882383,25:0.22558375570882383):2.007245514581259,(100:1.7529318910781768,77:1.7529318910781768):0.4798973792119054):3.30314797010882,(((47:2.412127585090697,42:2.412127585090697):0.06381801624833994,39:2.4759456013390375):1.9988347406856974,(((60:0.7432856899978842,44:0.7432856899978842):0.852847153640885,67:1.5961328436387678):0.11369229695192654,(((4:0.33763351452836454,41:0.33763351452836454):0.05512013726642647,85:0.39275365179479105):0.3304382550944581,79:0.7231919068892493):0.9866332337014465):2.7649552014340406):1.0611968983741689):1.883861048480488,((((((90:2.9750039079058963,(96:0.36482810447284125,7:0.36482810447284125):2.610175803433054):0.21268166611868503,45:3.1876855740245804):0.3467140311102542,(112:3.3346806090961354,(((((22:0.08868784754909187,48:0.08868784754909187):0.8673697522426342,((54:0.11795945726491863,103:0.11795945726491863):0.7868152554871702,74:0.904774712752089):0.05128288703963703):0.3070359299032765,84:1.2630935296950023):0.29624579222789255,(((106:0.8724202476269016,63:0.8724202476269016):0.6742388177831953,52:1.5466590654100956):0.0015972387672537928,18:1.5482563041773496):0.01108301774554423):1.5770486887314876,((50:0.662943154560071,21:0.662943154560071):0.21563299244357964,26:0.8785761470036507):2.2578118636507316):0.19829259844175298):0.1997189960386991):0.2412523606447941,(24:0.19119650297432514,5:0.19119650297432514):3.584455462805304):2.1696253617528183,(114:4.518259057526128,(((62:1.9117039054651033,((28:1.5883367314074492,58:1.5883367314074492):0.17294192790304666,35:1.7612786593104959):0.15042524615460715):1.432340637940378,((46:0.9860684908337861,(55:0.734345261409901,38:0.734345261409901):0.251723229423885):0.9721556984092379,37:1.9582241892430228):1.3858203541624587):0.03359492017413697,((76:0.3067260791051034,(12:0.3066649680278747,104:0.3066649680278747):6.111107722876423E-5):1.5742906714311293,65:1.8810167505362319):1.4966227130433862):1.1406195939465078):1.4270182700063205):1.3681545388255332,(70:0.12605856389810155,16:0.12605856389810155):7.187373302459879):0.10640642252141169):0.02628748476616751):0.2986003449951932):2.255273881359247,(((((15:0.5099041291360453,68:0.5099041291360453):4.308918959729866,(((110:1.9854652390255247,(113:1.5491825445188476,(88:0.4058435961978725,2:0.4058435961978725):1.1433389483209762):0.43628269450667667):0.05655391266137123,87:2.042019151686896):0.056255197095051855,(1:0.9058245321965291,9:0.9058245321965291):1.1924498165854198):2.720548740083961):0.3655832584541272,(17:2.9150987516232134,(91:1.5241705299365993,86:1.5241705299365993):1.3909282216866121):2.2693075956968243):1.033858243644117,(34:4.664372841854619,94:4.664372841854619):1.5538917491095365):0.5156336727474565,(((((73:1.7884984320569532,(98:0.7937515188270203,69:0.7937515188270203):0.9947469132299342):1.618762718665666,((107:2.3791769809361103,10:2.3791769809361103):0.7400536709099146,89:3.119230651846026):0.2880304988765942):0.9255504207303228,((101:0.4718010395819499,105:0.4718010395819499):2.4080017868065022,72:2.8798028263884525):1.4530087450644908):1.4475488343737912,(((36:0.610951962357272,14:0.610951962357272):4.0002929850592785,(81:3.1450192736149387,((27:1.7640580531708168,33:1.7640580531708168):1.3801252316941013,((31:0.44464105025371575,(53:0.30411460291765374,19:0.30411460291765374):0.14052644733606207):0.25135601193312385,3:0.6959970621868397):2.448186222678079):8.359887500193317E-4):1.4662256738016115):0.17534499299832046,(80:1.6393344257852065,(108:0.4335715018231018,111:0.4335715018231018):1.2057629239621057):3.1472555146296632):0.9937704654118635):0.3551904589282929,((59:0.8058774121717827,66:0.8058774121717827):5.0899477840211045,(((30:2.1543496987016644,(109:0.41540451522121635,57:0.41540451522121635):1.7389451834804492):0.38618053920082596,(13:2.3735492860672047,49:2.3735492860672047):0.16698095183528538):3.3539576458551,(((82:1.4569937634052825,(11:0.6898110408713436,(102:0.3902073644546026,(93:0.35124225136378895,71:0.35124225136378895):0.03896511309081358):0.29960367641674107):0.7671827225339397):1.1622193345884577,((29:0.06589883986766422,83:0.06589883986766422):2.3410845861359464,92:2.4069834260036105):0.21222967199012968):1.0779431963125545,(43:0.9203944544243726,51:0.9203944544243726):2.776761839881923):2.1973315894512955):0.0013373124352959022):0.2397256685621401):0.5983473989565837):3.2661017362883893);
                TREE 'Tree # 10 simulated by Uniform speciation (Yule)' = ((((84:3.9758265194174327,(92:2.811799977483899,68:2.811799977483899):1.1640265419335347):2.6099321095076884,(((103:0.7408270530709783,(18:0.7408015123807465,((77:0.6023781542577221,49:0.6023781542577221):0.13252076751955966,2:0.7348989217772817):0.005902590603464842):2.5540690231769523E-5):2.9461589191179485,(((6:0.32255228914200396,16:0.32255228914200396):0.8166959546838475,11:1.1392482438258513):1.7046717010507755,(65:0.8389437579849384,81:0.8389437579849384):2.0049761868916884):0.8430660273123011):0.8903226382251646,(((((44:0.9734632937379658,104:0.9734632937379658):0.36436412832991594,((96:0.03809384733416016,23:0.03809384733416016):0.697216994049078,7:0.7353108413832381):0.6025165806846439):0.5251437341339154,43:1.8629711562017972):0.5647756259970461,86:2.4277467821988425):1.096932410440006,(21:3.443311694246106,((79:0.4887078713988029,48:0.4887078713988029):1.8175704540612858,5:2.306278325460089):1.1370333687860181):0.08136749839274277):1.052629417775243):2.0084500185110303):1.7568264730036167,((80:3.079949599126683,(42:0.1551045213778532,94:0.1551045213778532):2.9248450777488304):3.161559456980693,(((((((62:0.04773013293792959,3:0.04773013293792959):0.4687059394631754,66:0.5164360724011049):0.8644350769354897,107:1.380871149336595):1.5643370712561024,(46:2.924887335299998,(58:2.1657445548141894,8:2.1657445548141894):0.7591427804858089):0.020320885292699098):0.3059390391394366,90:3.251147259732133):2.7228204842829125,29:5.973967744015043):0.24265956044133447,(((((9:4.377882559393299,(((59:1.6302525790282003,67:1.6302525790282003):0.040010737046462845,25:1.6702633160746632):0.903041357870075,((((30:0.13453900525794873,14:0.13453900525794873):0.9739626584431993,24:1.108501663701148):0.7248310210725202,47:1.8333326847736684):0.014685155797610594,87:1.848017840571279):0.725286833373459):1.8045778854485628):0.5999874422359001,((75:2.3285823445551466,(109:1.17381686659414,101:1.17381686659414):1.1547654779610066):2.643299856208181,(((111:2.0587228276033405,(99:0.617047396181564,27:0.617047396181564):1.4416754314217761):0.053612541871283045,31:2.1123353694746236):0.7581224642663948,63:2.870457833741017):2.1014243670223105):0.00598780086587163):0.3313441984400522,((((93:0.41932786822894863,91:0.41932786822894863):1.3204267245145784,72:1.7397545927435272):0.208009437621344,(69:1.1273234268168264,28:1.1273234268168264):0.8204406035480444):0.2578240187512625,45:2.2055880491161335):3.1036261509531182):0.5631693673879757,((60:3.881801993387112,((35:1.0153598152322105,54:1.0153598152322105):2.540305953868706,((97:2.278022798301517,((32:0.2586843773057431,108:0.2586843773057431):0.14401466189026685,56:0.4026990391960099):1.8753237591055067):0.2711648385015559,(70:0.46886699609455756,89:0.46886699609455756):2.0803206407085146):1.006478132297844):0.32613622428619704):1.3718195429580395,(((((53:0.27551473097404605,12:0.27551473097404605):0.35492780265590884,100:0.6304425336299548):3.6927442104610257,(4:1.9084665014509226,((41:0.8290905483899873,74:0.8290905483899873):0.35903781736264523,20:1.188128365752633):0.7203381356982894):2.4147202426400596):0.2029533407408023,((57:0.3541346155039802,61:0.3541346155039802):2.0313397994941207,17:2.385474414998101):2.1406656698336826):0.40488487425619446,37:4.931024959087975):0.32259657725717494):0.6187620311120753):0.27763927150374645,(((55:3.5630191525453174,((82:1.8003742135818652,((114:0.20878969445532472,50:0.20878969445532472):0.1828498288465328,110:0.3916395233018575):1.4087346902800075):1.1264846440623035,88:2.9268588576441683):0.6361602949011501):1.2398616686607513,(((102:0.22038812284149883,95:0.22038812284149883):1.7906269596143836,((26:0.5895375526779978,78:0.5895375526779978):0.014867503946891097,76:0.604405056624889):1.4066100258309935):1.3903218137377917,((19:1.604039244460153,38:1.604039244460153):0.5780395928004002,105:2.182078837260553):1.2192580589331212):1.4015439250123956):0.34043044507621484,73:5.143311266282281):1.00671157267869):0.06660446549540541):0.024881751650995848):2.1010760458213635):1.6574148980712593,(((((36:0.8606484952025208,112:0.8606484952025208):0.3511441470539232,39:1.2117926422564445):0.25846121563685165,(40:1.1361320822354615,13:1.1361320822354615):0.3341217756578342):1.1413415182594915,(34:2.471638085634848,33:2.471638085634848):0.1399572905179385):2.8611667912586607,(83:5.453632817653242,((((10:0.6800351846447861,(51:0.6299925584886462,71:0.6299925584886462):0.05004262615613981):0.37027317453549136,85:1.0503083591802775):1.4258807137370746,(((98:0.11483581960963396,(1:0.03758104624691612,22:0.03758104624691612):0.07725477336271784):1.0805376787248526,64:1.195373498334487):1.0764640884166898,15:2.271837586751177):0.2043514861661765):0.9705804052190817,((52:0.1807575816493472,113:0.1807575816493472):1.5300856991465288,106:1.7108432807958762):1.7359261973405586):2.006863339516812):0.01912934975820202):4.5272378325885505);
            END;
        """,
        """
            #NEXUS
            BEGIN TREES;
                Title Caenophidia_Random_Yule_2;
            [!Parameters: Tree simulator: Coalescent Trees. [seed: 1259785463881]]
                TRANSLATE
                    1 Lystrophis_dorbignyi,
                    2 Waglerophis_merremi,
                    3 Lystrophis_histricus,
                    4 Xenoxybelis_argenteus,
                    5 Liophis_jaegeri,
                    6 Liophis_elegantissimus,
                    7 Liophis_meridionalis,
                    8 Liophis_typhlus,
                    9 Erythrolamprus_aesculapii,
                    10 Atractus_albuquerquei,
                    11 Atractus_trihedrurus,
                    12 Heterodon_nasicus,
                    13 Sibynomorphus_mikanii,
                    14 Ninia_atrata,
                    15 Dipsas_indica,
                    16 Dipsas_neivai,
                    17 Sibynomorphus_garmani,
                    18 Hydrops_triangularis,
                    19 Helicops_infrataeniatus,
                    20 Pseudoeryx_plicatilis,
                    21 Helicops_pictiventris,
                    22 Helicops_gomesi,
                    23 Helicops_angulatus,
                    24 Oxyrhopus_rhombifer,
                    25 Oxyrhopus_clathratus,
                    26 Lycognathophis_seychellensis,
                    27 Taeniophallus_brevirostris,
                    28 Hydrodynastes_gigas,
                    29 Hydrodynastes_bicinctus,
                    30 Pseudoboa_nigra,
                    31 Pseudoboa_coronata,
                    32 Boiruna_maculata,
                    33 Drepanoides_anomalus,
                    34 Clelia_bicolor,
                    35 Phimophis_guerini,
                    36 Psomophis_joberti,
                    37 Psomophis_genimaculatus,
                    38 Leptodeira_annulata,
                    39 Thamnodynastes_rutilus,
                    40 Pseudotomodon_trigonatus,
                    41 Tachymenis_peruviana,
                    42 Tomodon_dorsatus,
                    43 Ptychophis_flavovirgatus,
                    44 Calamodontophis_paucidens,
                    45 Pseudablabes_agassizi,
                    46 Philodryas_patagoniensis,
                    47 Philodryas_aestiva,
                    48 Philodryas_mattogrossensis,
                    49 Taeniophallus_affinis,
                    50 Liophis_amarali,
                    51 Imantodes_cenchoa,
                    52 Apostolepis_dimidiata,
                    53 Phalotris_nasutus,
                    54 Phalotris_lemniscatus,
                    55 Apostolepis_assimilis,
                    56 Elapomorphus_quinquelineatus,
                    57 Siphlophis_pulcher,
                    58 Siphlophis_compressus,
                    59 Dendroaspis_polylepis,
                    60 Bothrolycus_ater,
                    61 Bothrophthalmus_lineatus,
                    62 Bothrophthalmus_brunneus,
                    63 Lamprophis_virgatus,
                    64 Lamprophis_fuliginosus_mentalis,
                    65 Lamprophis_lineatus,
                    66 Lamprophis_fuliginosus,
                    67 Lamprophis_olivaceus,
                    68 Lamprophis_capensis,
                    69 Lamprophis_inornatus,
                    70 Lycodonomorphus_whytii,
                    71 Lamprophis_fiskii,
                    72 Lycodonomorphus_rufulus,
                    73 Lamprophis_guttatus,
                    74 Homoroselaps_lacteus,
                    75 Atractaspis_bibronii,
                    76 Atractaspis_corpulenta,
                    77 Atractaspis_boulengeri,
                    78 Atractaspis_micropholis,
                    79 Aparallactus_modestus,
                    80 Polemon_collaris,
                    81 Polemon_acanthias,
                    82 Polemon_notatum,
                    83 Xenocalamus_transvaalensis,
                    84 Amblyodipsas_polylepis,
                    85 Macrelaps_microlepidotus,
                    86 Mehelya_nyassae,
                    87 Psammophis_phillipsi,
                    88 Psammophis_praeornatus,
                    89 Psammophis_sp.,
                    90 Psammophis_schokari,
                    91 Prosymna_janii,
                    92 Dromicodryas_bernieri,
                    93 Buhoma_depressiceps,
                    94 Liophidium_chabaudi,
                    95 Heteroliodon_occipitalis,
                    96 Lycophidion_capense,
                    97 Rhamphiophis_oxyrhynchus,
                    98 Hemirhagerrhis_hildebrandtii,
                    99 Amplorhinus_multimaculatus,
                    100 Lycophidion_laterale,
                    101 Buhoma_procterae,
                    102 Mehelya_stenophthalmus,
                    103 Psammodynastes_sp.,
                    104 Rhamphiophis_rostratus,
                    105 Pseudaspis_cana,
                    106 Gonionotophis_brussauxi,
                    107 Psammophylax_rhombeatus,
                    108 Psammophylax_variabilis,
                    109 Mehelya_poensis,
                    110 Duberria_variegata,
                    111 Duberria_lutrix,
                    112 Stenophis_citrinus,
                    113 Lycodryas_sanctijohannis,
                    114 Ithycyphus_oursi;
                TREE Tree_#_1_simulated_by_Coalescent_Trees = (7:10994.0,((((92:93.0,43:93.0):264.0,(36:1.0,62:1.0):356.0):6815.0,((((87:0.0,50:0.0):1754.0,((53:301.0,(((89:76.0,(3:14.0,79:14.0):62.0):56.0,(81:129.0,(46:20.0,109:20.0):109.0):3.0):3.0,88:135.0):166.0):395.0,65:696.0):1058.0):987.0,((((84:45.0,(15:3.0,38:3.0):42.0):194.0,(100:36.0,74:36.0):203.0):197.0,55:436.0):47.0,(((69:14.0,49:14.0):5.0,51:19.0):200.0,(91:8.0,73:8.0):211.0):264.0):2258.0):3065.0,((90:519.0,(((104:48.0,24:48.0):6.0,26:54.0):7.0,18:61.0):458.0):3313.0,((12:174.0,(11:158.0,60:158.0):16.0):2931.0,(78:2823.0,(((10:93.0,(((6:15.0,112:15.0):33.0,103:48.0):21.0,5:69.0):24.0):41.0,((14:60.0,(111:57.0,(35:47.0,66:47.0):10.0):3.0):58.0,108:118.0):16.0):279.0,(80:174.0,114:174.0):239.0):2410.0):282.0):727.0):1974.0):1366.0):1375.0,(((58:45.0,(110:28.0,85:28.0):17.0):1992.0,((((22:60.0,105:60.0):17.0,(32:76.0,(17:70.0,((95:6.0,67:6.0):7.0,33:13.0):57.0):6.0):1.0):415.0,(8:86.0,102:86.0):406.0):827.0,((((113:17.0,63:17.0):43.0,16:60.0):446.0,((54:283.0,(((23:14.0,39:14.0):74.0,2:88.0):182.0,30:270.0):13.0):200.0,(((((96:38.0,29:38.0):55.0,25:93.0):46.0,(64:125.0,34:125.0):14.0):55.0,93:194.0):207.0,(47:32.0,(13:2.0,20:2.0):30.0):369.0):82.0):23.0):161.0,(42:244.0,40:244.0):423.0):652.0):718.0):63.0,(((61:137.0,((101:30.0,70:30.0):49.0,(71:13.0,82:13.0):66.0):58.0):1498.0,(((((9:395.0,(37:89.0,((59:67.0,45:67.0):5.0,98:72.0):17.0):306.0):6.0,(57:11.0,28:11.0):390.0):80.0,(52:244.0,97:244.0):237.0):696.0,(((41:158.0,99:158.0):36.0,(27:20.0,48:20.0):174.0):526.0,(((1:45.0,(21:29.0,75:29.0):16.0):35.0,(86:48.0,106:48.0):32.0):222.0,76:302.0):418.0):457.0):152.0,((19:196.0,(31:155.0,77:155.0):41.0):232.0,(72:151.0,44:151.0):277.0):901.0):306.0):17.0,((4:481.0,(68:17.0,107:17.0):464.0):31.0,((83:50.0,94:50.0):121.0,56:171.0):341.0):1140.0):448.0):6447.0):2447.0):0.0;
                TREE Tree_#_2_simulated_by_Coalescent_Trees = (((((95:1050.0,((((10:54.0,32:54.0):149.0,18:203.0):143.0,((13:42.0,42:42.0):291.0,56:333.0):13.0):103.0,(105:377.0,(36:218.0,113:218.0):159.0):72.0):601.0):402.0,(((55:66.0,(92:55.0,74:55.0):11.0):80.0,20:146.0):223.0,(84:167.0,44:167.0):202.0):1083.0):2474.0,((((25:154.0,104:154.0):295.0,114:449.0):624.0,(((79:53.0,53:53.0):10.0,52:63.0):312.0,(((40:137.0,34:137.0):52.0,(35:116.0,((87:13.0,41:13.0):29.0,67:42.0):74.0):73.0):112.0,93:301.0):74.0):698.0):510.0,(((((11:40.0,2:40.0):202.0,38:242.0):725.0,(((76:178.0,(12:144.0,(((66:28.0,72:28.0):30.0,68:58.0):9.0,39:67.0):77.0):34.0):535.0,28:713.0):141.0,6:854.0):113.0):133.0,1:1100.0):288.0,((49:118.0,((46:26.0,5:26.0):78.0,(78:67.0,14:67.0):37.0):14.0):71.0,54:189.0):1199.0):195.0):2343.0):1347.0,(((((((85:100.0,58:100.0):236.0,110:336.0):28.0,((22:34.0,43:34.0):144.0,70:178.0):186.0):1741.0,(98:528.0,((45:89.0,(7:44.0,(83:30.0,(30:2.0,107:2.0):28.0):14.0):45.0):93.0,77:182.0):346.0):1577.0):117.0,(51:368.0,(27:89.0,(82:28.0,75:28.0):61.0):279.0):1854.0):451.0,(89:1142.0,80:1142.0):1531.0):912.0,(50:411.0,(106:213.0,((17:57.0,23:57.0):127.0,(108:56.0,63:56.0):128.0):29.0):198.0):3174.0):1688.0):16115.0,(((94:588.0,(((19:203.0,(((81:6.0,57:6.0):36.0,(103:2.0,9:2.0):40.0):49.0,(((8:11.0,3:11.0):15.0,(4:2.0,24:2.0):24.0):63.0,73:89.0):2.0):112.0):131.0,26:334.0):146.0,88:480.0):108.0):137.0,((21:19.0,48:19.0):317.0,71:336.0):389.0):3806.0,((((99:303.0,(101:13.0,47:13.0):290.0):2.0,((15:14.0,31:14.0):136.0,109:150.0):155.0):1779.0,(((91:415.0,((29:149.0,59:149.0):47.0,111:196.0):219.0):64.0,33:479.0):1392.0,(((90:36.0,96:36.0):53.0,(37:75.0,16:75.0):14.0):166.0,86:255.0):1616.0):213.0):709.0,((65:281.0,(69:203.0,97:203.0):78.0):437.0,((100:307.0,64:307.0):9.0,((62:122.0,(112:116.0,(102:58.0,61:58.0):58.0):6.0):183.0,60:305.0):11.0):402.0):2075.0):1738.0):16857.0):0.0;
                TREE Tree_#_3_simulated_by_Coalescent_Trees = (((((94:260.0,((112:3.0,5:3.0):174.0,(106:34.0,76:34.0):143.0):83.0):1167.0,(((((107:8.0,21:8.0):10.0,64:18.0):7.0,90:25.0):290.0,(27:158.0,(19:122.0,(55:2.0,78:2.0):120.0):36.0):157.0):605.0,((((1:28.0,(26:2.0,99:2.0):26.0):125.0,(17:81.0,72:81.0):72.0):235.0,(74:274.0,(75:222.0,39:222.0):52.0):114.0):339.0,(62:36.0,41:36.0):691.0):193.0):507.0):386.0,(113:637.0,(18:596.0,(61:128.0,96:128.0):468.0):41.0):1176.0):2437.0,(20:3088.0,((((84:89.0,(46:49.0,111:49.0):40.0):287.0,(36:136.0,(34:37.0,66:37.0):99.0):240.0):180.0,(97:389.0,((44:16.0,104:16.0):223.0,((53:58.0,38:58.0):49.0,69:107.0):132.0):150.0):167.0):1518.0,((16:475.0,((11:108.0,(77:3.0,88:3.0):105.0):255.0,91:363.0):112.0):371.0,((((15:32.0,7:32.0):413.0,(114:8.0,85:8.0):437.0):3.0,((82:196.0,(47:75.0,9:75.0):121.0):113.0,((((43:14.0,100:14.0):34.0,108:48.0):107.0,98:155.0):30.0,(42:110.0,52:110.0):75.0):124.0):139.0):303.0,((23:58.0,35:58.0):353.0,(40:207.0,83:207.0):204.0):340.0):95.0):1228.0):1014.0):1162.0):26422.0,(((((((51:104.0,((71:34.0,54:34.0):5.0,45:39.0):65.0):21.0,(28:3.0,13:3.0):122.0):763.0,((63:36.0,(50:22.0,10:22.0):14.0):143.0,14:179.0):709.0):1275.0,((95:133.0,(29:8.0,24:8.0):125.0):1262.0,(70:1237.0,((6:332.0,(((((56:22.0,3:22.0):125.0,(2:136.0,(101:76.0,49:76.0):60.0):11.0):6.0,87:153.0):16.0,(12:36.0,67:36.0):133.0):74.0,57:243.0):89.0):892.0,(80:28.0,33:28.0):1196.0):13.0):158.0):768.0):7.0,((8:143.0,30:143.0):277.0,(86:147.0,81:147.0):273.0):1750.0):801.0,(60:51.0,89:51.0):2920.0):5079.0,(((59:810.0,(93:450.0,73:450.0):360.0):2302.0,((48:102.0,103:102.0):2310.0,((109:235.0,(79:3.0,65:3.0):232.0):494.0,(92:49.0,25:49.0):680.0):1683.0):700.0):1867.0,(((102:75.0,32:75.0):469.0,((37:189.0,58:189.0):47.0,22:236.0):308.0):2218.0,(((68:22.0,110:22.0):37.0,(4:3.0,31:3.0):56.0):649.0,105:708.0):2054.0):2217.0):3071.0):22622.0):0.0;
                TREE Tree_#_4_simulated_by_Coalescent_Trees = ((48:56.0,27:56.0):31205.0,((((52:187.0,(54:180.0,(32:168.0,41:168.0):12.0):7.0):62.0,((73:44.0,(68:29.0,85:29.0):15.0):4.0,84:48.0):201.0):183.0,(24:269.0,((92:37.0,53:37.0):207.0,((13:29.0,(63:27.0,(45:14.0,100:14.0):13.0):2.0):46.0,49:75.0):169.0):25.0):163.0):11966.0,(((((7:37.0,6:37.0):18.0,26:55.0):712.0,(104:272.0,(22:42.0,71:42.0):230.0):495.0):2176.0,((14:73.0,35:73.0):1429.0,(61:364.0,96:364.0):1138.0):1441.0):2126.0,(((((((((((((102:17.0,75:17.0):25.0,101:42.0):45.0,((10:3.0,23:3.0):81.0,65:84.0):3.0):62.0,(88:14.0,4:14.0):135.0):38.0,105:187.0):125.0,(57:197.0,113:197.0):115.0):283.0,((1:37.0,99:37.0):481.0,25:518.0):77.0):181.0,89:776.0):973.0,36:1749.0):202.0,(((94:777.0,(21:86.0,64:86.0):691.0):189.0,((39:289.0,110:289.0):238.0,(109:407.0,(76:69.0,((81:29.0,(107:13.0,19:13.0):16.0):11.0,37:40.0):29.0):338.0):120.0):439.0):278.0,(31:10.0,58:10.0):1234.0):707.0):252.0,(((((79:161.0,(42:30.0,(5:5.0,62:5.0):25.0):131.0):41.0,(16:31.0,70:31.0):171.0):381.0,(111:55.0,106:55.0):528.0):43.0,((56:196.0,((93:14.0,87:14.0):60.0,(15:56.0,47:56.0):18.0):122.0):234.0,44:430.0):196.0):1521.0,(((98:11.0,86:11.0):30.0,43:41.0):362.0,55:403.0):1744.0):56.0):1123.0,((((17:65.0,(82:20.0,59:20.0):45.0):88.0,34:153.0):382.0,50:535.0):258.0,((103:30.0,80:30.0):86.0,33:116.0):677.0):2533.0):986.0,((((78:44.0,2:44.0):200.0,(114:21.0,38:21.0):223.0):1580.0,(12:1755.0,(77:1363.0,((((((8:29.0,40:29.0):14.0,74:43.0):30.0,(9:21.0,83:21.0):52.0):111.0,28:184.0):358.0,((95:160.0,(51:147.0,29:147.0):13.0):326.0,((69:130.0,11:130.0):146.0,((112:29.0,30:29.0):156.0,20:185.0):91.0):210.0):56.0):364.0,(60:534.0,(46:191.0,108:191.0):343.0):372.0):457.0):392.0):69.0):1544.0,(3:729.0,(((18:14.0,72:14.0):236.0,(97:153.0,90:153.0):97.0):270.0,((91:348.0,67:348.0):70.0,66:418.0):102.0):209.0):2639.0):944.0):757.0):7329.0):18863.0):0.0;
                TREE Tree_#_5_simulated_by_Coalescent_Trees = ((((92:28.0,72:28.0):922.0,((48:21.0,98:21.0):383.0,(((18:34.0,110:34.0):144.0,(90:128.0,33:128.0):50.0):133.0,((10:3.0,107:3.0):68.0,84:71.0):240.0):93.0):546.0):187.0,((((86:57.0,4:57.0):2.0,94:59.0):221.0,((((55:19.0,64:19.0):20.0,6:39.0):86.0,62:125.0):57.0,43:182.0):98.0):3.0,59:283.0):854.0):9130.0,((((95:2508.0,(((((9:87.0,37:87.0):1.0,(113:70.0,112:70.0):18.0):468.0,(((87:34.0,(54:18.0,104:18.0):16.0):370.0,(56:128.0,41:128.0):276.0):39.0,(24:214.0,42:214.0):229.0):113.0):496.0,((20:321.0,(1:159.0,46:159.0):162.0):84.0,82:405.0):647.0):500.0,(111:128.0,40:128.0):1424.0):956.0):939.0,(((23:200.0,32:200.0):76.0,((105:73.0,(80:9.0,114:9.0):64.0):29.0,((50:64.0,11:64.0):0.0,63:64.0):38.0):174.0):637.0,(3:319.0,79:319.0):594.0):2534.0):224.0,(((((100:204.0,12:204.0):66.0,(109:42.0,8:42.0):228.0):463.0,21:733.0):1934.0,(((((15:17.0,47:17.0):234.0,(38:7.0,49:7.0):244.0):301.0,88:552.0):12.0,45:564.0):1454.0,(75:1867.0,(((((97:64.0,89:64.0):76.0,30:140.0):35.0,(53:111.0,((81:40.0,36:40.0):11.0,(74:49.0,(16:17.0,93:17.0):32.0):2.0):60.0):64.0):0.0,60:175.0):296.0,((39:16.0,27:16.0):411.0,((((101:14.0,69:14.0):227.0,(13:18.0,68:18.0):223.0):22.0,((66:89.0,83:89.0):39.0,(96:80.0,((58:18.0,29:18.0):37.0,34:55.0):25.0):48.0):135.0):111.0,(((52:18.0,7:18.0):237.0,51:255.0):4.0,(25:2.0,61:2.0):257.0):115.0):53.0):44.0):1396.0):151.0):649.0):205.0,(((((2:104.0,71:104.0):6.0,70:110.0):27.0,67:137.0):104.0,(((76:7.0,31:7.0):0.0,5:7.0):190.0,17:197.0):44.0):2324.0,((57:57.0,35:57.0):71.0,(108:9.0,65:9.0):119.0):2437.0):307.0):799.0):3787.0,(((44:776.0,((102:24.0,19:24.0):55.0,28:79.0):697.0):0.0,77:776.0):591.0,((106:450.0,(((26:30.0,103:30.0):59.0,78:89.0):76.0,73:165.0):285.0):479.0,(((14:89.0,(85:61.0,22:61.0):28.0):470.0,99:559.0):274.0,91:833.0):96.0):438.0):6091.0):2809.0):0.0;
                TREE Tree_#_6_simulated_by_Coalescent_Trees = (((113:828.0,((((2:6.0,40:6.0):285.0,(53:103.0,14:103.0):188.0):172.0,114:463.0):220.0,(74:13.0,9:13.0):670.0):145.0):8858.0,(((((((90:82.0,(((52:7.0,100:7.0):7.0,28:14.0):47.0,83:61.0):21.0):135.0,(84:183.0,82:183.0):34.0):184.0,((37:76.0,(11:61.0,96:61.0):15.0):10.0,(64:31.0,16:31.0):55.0):315.0):1189.0,103:1590.0):110.0,(6:22.0,43:22.0):1678.0):689.0,70:2389.0):1486.0,(((((((65:13.0,31:13.0):6.0,34:19.0):112.0,110:131.0):487.0,(13:27.0,63:27.0):591.0):194.0,((15:27.0,49:27.0):778.0,(5:295.0,((108:7.0,89:7.0):61.0,25:68.0):227.0):510.0):7.0):1517.0,((12:141.0,(30:48.0,107:48.0):93.0):487.0,39:628.0):1701.0):613.0,(((41:113.0,46:113.0):283.0,(((102:44.0,(86:29.0,24:29.0):15.0):71.0,29:115.0):203.0,(((85:29.0,56:29.0):248.0,((48:99.0,58:99.0):78.0,(77:141.0,((66:44.0,45:44.0):39.0,19:83.0):58.0):36.0):100.0):30.0,(106:24.0,78:24.0):283.0):11.0):78.0):2182.0,(((79:46.0,81:46.0):113.0,71:159.0):192.0,(93:39.0,(62:24.0,61:24.0):15.0):312.0):2227.0):364.0):933.0):5811.0):20101.0,(((8:340.0,(1:169.0,((21:27.0,105:27.0):29.0,76:56.0):113.0):171.0):1002.0,(((((91:412.0,(27:352.0,(97:131.0,((104:27.0,38:27.0):44.0,(94:14.0,88:14.0):57.0):60.0):221.0):60.0):93.0,(67:220.0,69:220.0):285.0):42.0,((((7:44.0,(33:22.0,98:22.0):22.0):103.0,59:147.0):208.0,((3:71.0,57:71.0):27.0,42:98.0):257.0):41.0,(32:144.0,60:144.0):252.0):151.0):619.0,(112:238.0,35:238.0):928.0):143.0,((23:45.0,((80:14.0,44:14.0):6.0,50:20.0):25.0):420.0,((68:75.0,72:75.0):30.0,(36:54.0,47:54.0):51.0):360.0):844.0):33.0):7284.0,((((4:7.0,92:7.0):1027.0,(((111:89.0,20:89.0):19.0,10:108.0):231.0,(73:157.0,17:157.0):182.0):695.0):1100.0,(55:823.0,(22:211.0,99:211.0):612.0):1311.0):1510.0,(((109:115.0,18:115.0):174.0,((75:99.0,(54:41.0,87:41.0):58.0):110.0,(26:205.0,(101:178.0,95:178.0):27.0):4.0):80.0):296.0,51:585.0):3059.0):4982.0):21161.0):0.0;
                TREE Tree_#_7_simulated_by_Coalescent_Trees = (((88:1326.0,((38:22.0,84:22.0):308.0,((6:43.0,(59:25.0,34:25.0):18.0):259.0,14:302.0):28.0):996.0):2691.0,((((((3:13.0,81:13.0):862.0,114:875.0):178.0,(91:888.0,(64:46.0,31:46.0):842.0):165.0):937.0,((99:239.0,77:239.0):1281.0,((((19:2.0,22:2.0):2.0,97:4.0):126.0,(((75:33.0,(20:6.0,80:6.0):27.0):48.0,((104:38.0,(11:21.0,65:21.0):17.0):16.0,(56:16.0,109:16.0):38.0):27.0):24.0,(66:4.0,94:4.0):101.0):25.0):754.0,(42:285.0,(74:92.0,(51:88.0,(32:4.0,24:4.0):84.0):4.0):193.0):599.0):636.0):470.0):1288.0,48:3278.0):463.0,(21:224.0,(106:165.0,(102:43.0,96:43.0):122.0):59.0):3517.0):276.0):2887.0,(((((1:42.0,26:42.0):151.0,(7:54.0,112:54.0):139.0):4160.0,(((35:61.0,17:61.0):575.0,((107:8.0,87:8.0):8.0,18:16.0):620.0):3707.0,(101:277.0,(10:97.0,25:97.0):180.0):4066.0):10.0):604.0,(((((8:19.0,113:19.0):163.0,((12:81.0,((78:10.0,60:10.0):6.0,72:16.0):65.0):67.0,(37:79.0,4:79.0):69.0):34.0):164.0,62:346.0):3426.0,(((((((39:11.0,103:11.0):182.0,(58:171.0,13:171.0):22.0):1216.0,54:1409.0):111.0,((((98:25.0,83:25.0):169.0,((27:16.0,89:16.0):7.0,(29:9.0,41:9.0):14.0):171.0):128.0,40:322.0):1122.0,(15:658.0,(((63:68.0,100:68.0):3.0,108:71.0):81.0,30:152.0):506.0):786.0):76.0):410.0,(52:1009.0,(((61:97.0,(111:46.0,49:46.0):51.0):124.0,(71:83.0,5:83.0):138.0):455.0,23:676.0):333.0):921.0):215.0,((36:717.0,((16:150.0,(57:4.0,46:4.0):146.0):143.0,(((95:76.0,2:76.0):72.0,((69:44.0,43:44.0):23.0,93:67.0):81.0):111.0,(85:68.0,105:68.0):191.0):34.0):424.0):383.0,47:1100.0):1045.0):39.0,(((28:63.0,76:63.0):5.0,68:68.0):89.0,44:157.0):2027.0):1588.0):757.0,(((82:286.0,73:286.0):70.0,(55:34.0,53:34.0):322.0):639.0,(92:368.0,((67:308.0,110:308.0):38.0,((90:34.0,45:34.0):150.0,(9:96.0,(50:4.0,33:4.0):92.0):88.0):162.0):22.0):627.0):3534.0):428.0):425.0,(79:330.0,(70:77.0,86:77.0):253.0):5052.0):1522.0):0.0;
                TREE Tree_#_8_simulated_by_Coalescent_Trees = (((((((111:84.0,112:84.0):279.0,(82:200.0,(104:17.0,52:17.0):183.0):163.0):473.0,(5:15.0,8:15.0):821.0):1410.0,((((30:23.0,23:23.0):89.0,29:112.0):1268.0,((109:78.0,2:78.0):278.0,(((16:38.0,75:38.0):108.0,88:146.0):50.0,(72:62.0,92:62.0):134.0):160.0):1024.0):502.0,((70:1071.0,108:1071.0):43.0,((67:132.0,(107:126.0,(12:33.0,38:33.0):93.0):6.0):571.0,94:703.0):411.0):768.0):364.0):316.0,(((((45:100.0,14:100.0):7.0,(35:14.0,56:14.0):93.0):712.0,((((22:191.0,71:191.0):118.0,(20:284.0,(17:273.0,((36:1.0,69:1.0):9.0,43:10.0):263.0):11.0):25.0):167.0,64:476.0):22.0,(106:227.0,41:227.0):271.0):321.0):393.0,((((33:220.0,((102:8.0,63:8.0):153.0,((101:40.0,(81:15.0,59:15.0):25.0):3.0,34:43.0):118.0):59.0):50.0,(((113:30.0,(58:23.0,51:23.0):7.0):123.0,(44:32.0,114:32.0):121.0):74.0,(9:41.0,73:41.0):186.0):43.0):6.0,76:276.0):173.0,15:449.0):763.0):170.0,((((74:427.0,(100:318.0,(90:78.0,((7:52.0,93:52.0):7.0,(49:1.0,66:1.0):58.0):19.0):240.0):109.0):12.0,(105:67.0,(39:27.0,4:27.0):40.0):372.0):380.0,((10:20.0,(87:13.0,25:13.0):7.0):106.0,((18:20.0,50:20.0):67.0,98:87.0):39.0):693.0):362.0,(78:199.0,(19:18.0,77:18.0):181.0):982.0):201.0):1180.0):761.0,((79:90.0,11:90.0):1131.0,(((97:3.0,60:3.0):175.0,65:178.0):262.0,(62:12.0,47:12.0):428.0):781.0):2102.0):9721.0,((((40:270.0,55:270.0):383.0,((85:13.0,21:13.0):563.0,((91:69.0,13:69.0):233.0,(6:2.0,24:2.0):300.0):274.0):77.0):2326.0,((37:15.0,54:15.0):2373.0,(89:1952.0,((95:849.0,((80:133.0,53:133.0):162.0,(57:240.0,(26:125.0,103:125.0):115.0):55.0):554.0):704.0,(((32:35.0,31:35.0):230.0,99:265.0):683.0,(1:462.0,28:462.0):486.0):605.0):399.0):436.0):591.0):3573.0,(((((83:14.0,46:14.0):294.0,((68:94.0,48:94.0):75.0,96:169.0):139.0):150.0,3:458.0):702.0,((61:2.0,110:2.0):57.0,86:59.0):1101.0):2449.0,((27:3.0,42:3.0):605.0,84:608.0):3001.0):2943.0):6492.0):0.0;
                TREE Tree_#_9_simulated_by_Coalescent_Trees = (((((54:336.0,92:336.0):396.0,76:732.0):491.0,((((106:116.0,3:116.0):41.0,40:157.0):91.0,((4:20.0,7:20.0):222.0,111:242.0):6.0):34.0,6:282.0):941.0):14711.0,((((1:1196.0,(49:1170.0,((((38:16.0,108:16.0):443.0,(71:24.0,(83:24.0,19:24.0):0.0):435.0):88.0,((42:62.0,110:62.0):31.0,(107:88.0,20:88.0):5.0):454.0):16.0,(((((112:15.0,78:15.0):122.0,34:137.0):72.0,((109:20.0,31:20.0):50.0,105:70.0):139.0):140.0,((9:151.0,48:151.0):3.0,35:154.0):195.0):136.0,(102:210.0,(47:29.0,98:29.0):181.0):275.0):78.0):607.0):26.0):1075.0,(((81:8.0,(45:5.0,97:5.0):3.0):156.0,(100:9.0,89:9.0):155.0):328.0,(56:222.0,46:222.0):270.0):1779.0):97.0,(((5:138.0,36:138.0):387.0,113:525.0):550.0,(((33:243.0,((99:114.0,96:114.0):60.0,60:174.0):69.0):115.0,73:358.0):101.0,((52:72.0,61:72.0):39.0,2:111.0):348.0):616.0):1293.0):628.0,((((39:106.0,80:106.0):72.0,(51:101.0,70:101.0):77.0):649.0,((57:181.0,72:181.0):113.0,(((84:26.0,(8:21.0,10:21.0):5.0):86.0,(24:9.0,22:9.0):103.0):111.0,41:223.0):71.0):533.0):80.0,(((69:70.0,(29:16.0,55:16.0):54.0):37.0,(12:90.0,(11:80.0,77:80.0):10.0):17.0):416.0,44:523.0):384.0):2089.0):12938.0):2397.0,((((((59:32.0,(82:5.0,17:5.0):27.0):47.0,(21:2.0,13:2.0):77.0):14.0,((94:34.0,90:34.0):38.0,53:72.0):21.0):464.0,(50:16.0,103:16.0):541.0):3114.0,((((((88:50.0,93:50.0):22.0,37:72.0):26.0,16:98.0):79.0,25:177.0):359.0,(((30:16.0,91:16.0):34.0,23:50.0):138.0,114:188.0):348.0):2036.0,(75:551.0,62:551.0):2021.0):1099.0):976.0,(((67:465.0,(((32:43.0,26:43.0):11.0,58:54.0):146.0,((43:29.0,63:29.0):63.0,64:92.0):108.0):265.0):79.0,(27:404.0,(104:157.0,(((101:26.0,(18:21.0,((65:12.0,95:12.0):9.0,68:21.0):0.0):5.0):64.0,79:90.0):36.0,87:126.0):31.0):247.0):140.0):715.0,(((15:79.0,74:79.0):13.0,((66:9.0,86:9.0):12.0,(14:9.0,85:9.0):12.0):71.0):98.0,28:190.0):1069.0):3388.0):13684.0):0.0;
                TREE Tree_#_10_simulated_by_Coalescent_Trees = (((14:1646.0,((62:193.0,(92:43.0,104:43.0):150.0):1180.0,((95:40.0,76:40.0):978.0,(109:852.0,(107:490.0,(((16:102.0,41:102.0):250.0,(74:146.0,((108:53.0,112:53.0):41.0,(96:6.0,25:6.0):88.0):52.0):206.0):39.0,48:391.0):99.0):362.0):166.0):355.0):273.0):1480.0,(33:1427.0,(((85:43.0,52:43.0):96.0,111:139.0):251.0,(87:383.0,80:383.0):7.0):1037.0):1699.0):7819.0,(((110:128.0,((81:104.0,(82:34.0,37:34.0):70.0):6.0,31:110.0):18.0):1245.0,(90:584.0,((36:1.0,91:1.0):143.0,53:144.0):440.0):789.0):7288.0,((((102:23.0,((105:12.0,29:12.0):5.0,23:17.0):6.0):170.0,(((83:38.0,79:38.0):24.0,(42:10.0,26:10.0):52.0):5.0,(97:20.0,(8:3.0,93:3.0):17.0):47.0):126.0):1363.0,(((((24:1.0,11:1.0):11.0,38:12.0):122.0,(44:132.0,51:132.0):2.0):169.0,((19:175.0,(27:104.0,(61:92.0,64:92.0):12.0):71.0):82.0,13:257.0):46.0):330.0,((((77:95.0,35:95.0):80.0,(68:20.0,4:20.0):155.0):156.0,((86:142.0,(10:34.0,65:34.0):108.0):44.0,(((50:15.0,5:15.0):85.0,(89:80.0,(6:23.0,21:23.0):57.0):20.0):64.0,((71:43.0,((45:20.0,72:20.0):19.0,59:39.0):4.0):39.0,100:82.0):82.0):22.0):145.0):244.0,((55:338.0,(34:330.0,54:330.0):8.0):210.0,(75:191.0,(3:25.0,56:25.0):166.0):357.0):27.0):58.0):923.0):3833.0,(((39:1164.0,((((9:160.0,58:160.0):35.0,((30:44.0,73:44.0):145.0,(88:183.0,43:183.0):6.0):6.0):100.0,((66:20.0,78:20.0):87.0,84:107.0):188.0):683.0,(106:10.0,22:10.0):968.0):186.0):1111.0,(((17:598.0,(((15:78.0,40:78.0):239.0,99:317.0):230.0,46:547.0):51.0):219.0,((114:10.0,94:10.0):286.0,60:296.0):521.0):1069.0,((69:53.0,32:53.0):0.0,47:53.0):1833.0):389.0):254.0,(((((70:27.0,113:27.0):164.0,(49:82.0,12:82.0):109.0):50.0,(((2:42.0,57:42.0):100.0,101:142.0):78.0,20:220.0):21.0):613.0,18:854.0):458.0,((63:452.0,(67:243.0,((1:156.0,(7:108.0,98:108.0):48.0):79.0,103:235.0):8.0):209.0):566.0,28:1018.0):294.0):1217.0):2860.0):3272.0):2284.0):0.0;

            END;
        """,
        """
            #NEXUS
            BEGIN TREES;
                Title Primates_Random_Yule;
            [!Parameters: Tree simulator: Uniform speciation (Yule). [seed: 1259785889153]]
                TRANSLATE
                    1 Lemur_catta,
                    2 Homo_sapiens,
                    3 Pan,
                    4 Gorilla,
                    5 Pongo,
                    6 Hylobates,
                    7 Macaca_fuscata,
                    8 Macaca_mulatta,
                    9 Macaca_fascicularis,
                    10 Macaca_sylvanus,
                    11 Saimiri_sciureus,
                    12 Tarsius_syrichta;
                TREE 'Tree # 1 simulated by Uniform speciation (Yule)' = (((2:2.9422567834231814,((10:0.444699078418728,(3:0.41779614238296503,(12:0.061659151274651755,11:0.061659151274651755):0.35613699110831337):0.02690293603576299):0.5091364998727755,1:0.9538355782915036):1.9884212051316779):3.5402216677894742,9:6.482478451212656):3.517521548787344,((4:1.4347879453096342,(8:1.3314808989756342,(5:0.1645667139624229,7:0.1645667139624229):1.166914185013211):0.10330704633400022):3.9444946173265514,6:5.3792825626361855):4.620717437363813);
                TREE 'Tree # 2 simulated by Uniform speciation (Yule)' = ((4:1.5613553157054088,(11:0.23886218406593138,5:0.23886218406593138):1.3224931316394772):8.438644684294589,((6:7.384533879140683,(3:3.8519695895747725,12:3.8519695895747725):3.532564289565911):2.582439221207633,(((2:0.38149390569100095,8:0.38149390569100095):4.61913048557425,9:5.000624391265251):1.0813448191432584,(10:4.531357186405765,(1:1.3423253659972916,7:1.3423253659972916):3.1890318204084727):1.5506120240027446):3.885003889939807):0.03302689965168175);
                TREE 'Tree # 3 simulated by Uniform speciation (Yule)' = ((8:0.4870727355699301,11:0.4870727355699301):9.512927264430065,((6:5.387377829292574,(((9:3.707387876969255,12:3.707387876969255):0.5638841724812745,3:4.27127204945053):0.38758006202517276,(2:4.125361263356914,1:4.125361263356914):0.5334908481187889):0.7285257178168713):3.3720343278831013,(5:8.370074656396419,(10:7.8222724989804115,(7:4.90057284806337,4:4.90057284806337):2.9216996509170414):0.547802157416008):0.3893375007792559):1.2405878428243207);
                TREE 'Tree # 4 simulated by Uniform speciation (Yule)' = (((((9:0.8485168762756372,7:0.8485168762756372):1.494656625827443,((10:0.28587778527112556,(11:0.27618714785933973,8:0.27618714785933973):0.009690637411785838):0.11131771026531428,6:0.39719549553643985):1.9459780065666403):1.751917441724805,5:4.095090943827885):2.815857030080398,((4:0.48206811914970554,3:0.48206811914970554):3.211419610149455,12:3.6934877292991604):3.2174602446091227):3.089052026091716,(2:6.024291320545672,1:6.024291320545672):3.975708679454328);
                TREE 'Tree # 5 simulated by Uniform speciation (Yule)' = (((6:2.182225802535574,((2:1.705821071987443,5:1.705821071987443):0.306140079204767,((((11:0.6784194845929409,12:0.6784194845929409):0.25911267241651365,9:0.9375321570094546):0.2927799954461013,3:1.2303121524555558):0.25235293773994905,(7:0.5320378489177201,8:0.5320378489177201):0.9506272412777849):0.529296060996705):0.1702646513433639):2.300237256268224,1:4.482463058803798):5.517536941196201,(10:9.629626416739569,4:9.629626416739569):0.3703735832604297);
                TREE 'Tree # 6 simulated by Uniform speciation (Yule)' = ((11:5.992148561545486,3:5.992148561545486):4.007851438454514,(((8:1.2768199873499768,5:1.2768199873499768):2.4927472622850493,((1:0.6381070709704036,7:0.6381070709704036):3.001655205593089,9:3.6397622765634923):0.12980497307153308):1.08436711640704,((2:2.4934504876925665,(6:1.1080941614199418,12:1.1080941614199418):1.3853563262726245):0.5286211707585364,(4:2.3219752288362177,10:2.3219752288362177):0.7000964296148849):1.831862707590963):5.146065633957934);
                TREE 'Tree # 7 simulated by Uniform speciation (Yule)' = ((12:5.9242576929998645,(7:1.903933411394498,11:1.903933411394498):4.020324281605367):4.075742307000133,(2:9.81121572915273,((1:4.031738766259758,((6:0.0744051889643847,10:0.0744051889643847):2.734316527529407,(3:2.6865818708597455,4:2.6865818708597455):0.1221398456340459):1.2230170497659658):5.719051985215255,((5:0.7859796988528565,9:0.7859796988528565):2.3780318497092656,8:3.164011548562122):6.58677920291289):0.06042497767771555):0.1887842708472696);
                TREE 'Tree # 8 simulated by Uniform speciation (Yule)' = (((2:4.3105571484354455,((5:2.166316798985596,1:2.166316798985596):1.8665026328845489,(((9:0.03976144687584101,8:0.03976144687584101):0.511772906188602,11:0.5515343530644431):1.3007780411548826,12:1.8523123942193258):2.180507037650819):0.2777377165653014):3.5143557897335005,10:7.824912938168946):2.1750870618310536,((4:0.9409065295517804,(7:0.4502105828006179,6:0.4502105828006179):0.49069594675116257):1.6428494097231154,3:2.583755939274896):7.416244060725104);
                TREE 'Tree # 9 simulated by Uniform speciation (Yule)' = (((11:0.3179882624839411,8:0.3179882624839411):2.6270360070710397,(1:0.19508249512113912,6:0.19508249512113912):2.7499417744338417):7.05497573044502,((4:7.371973710214715,3:7.371973710214715):2.3703833686631777,((2:6.285333053310127,((9:6.07663357307867,5:6.07663357307867):0.1765796057970613,7:6.253213178875732):0.03211987443439625):3.155381197416443,(10:9.091983017663157,12:9.091983017663157):0.34873123306341397):0.3016428281513226):0.2576429211221062);
                TREE 'Tree # 10 simulated by Uniform speciation (Yule)' = ((8:5.510639717184935,9:5.510639717184935):4.489360282815063,(((11:1.0027872608906918,(4:0.9189144699809287,(2:0.8482922713434907,10:0.8482922713434907):0.07062219863743803):0.08387279090976313):1.4070657766214212,(1:2.3011811356458574,12:2.3011811356458574):0.10867190186625546):7.04163287214505,(3:4.397130886637099,(6:4.153635258997215,(5:0.42095274625537166,7:0.42095274625537166):3.7326825127418437):0.24349562763988294):5.054355023020064):0.5485140903428354);

            END;
        """,
        """
            #NEXUS
            BEGIN TREES;
                Title Pythonidae_CytB_MLE;
                TRANSLATE
                    1 Python_molurus,
                    2 Python_sebae,
                    3 Morelia_tracyae,
                    4 Morelia_amethistina,
                    5 Morelia_nauta,
                    6 Morelia_kinghorni,
                    7 Morelia_clastolepis,
                    8 Morelia_boeleni,
                    9 Python_reticulatus,
                    10 Python_timoriensis,
                    11 Morelia_oenpelliensis,
                    12 Morelia_viridis,
                    13 Morelia_carinata,
                    14 Morelia_spilota,
                    15 Morelia_bredli,
                    16 Antaresia_maculosa,
                    17 Antaresia_childreni,
                    18 Antaresia_stimsoni,
                    19 Antaresia_perthensis,
                    20 Leiopython_albertisii,
                    21 Bothrochilus_boa,
                    22 Liasis_olivaceus,
                    23 Liasis_mackloti,
                    24 Liasis_fuscus,
                    25 Apodora_papuana,
                    26 Aspidites_ramsayi,
                    27 Aspidites_melanocephalus,
                    28 Python_brongersmai,
                    29 Python_regius;
                TREE 0 = (1:0.0779719244,((2:0.1414715009,(((((3:0.0435011998,(4:0.0305993564,((5:0.0092774432,6:0.0093145395):0.005595,7:0.005204698):0.023435):0.012223):0.025359,8:0.0863199106):0.019894,((9:0.0828549023,10:0.0963051344):0.072003,11:0.0820543043):0.002785):0.00274,((((12:0.0925974416,(13:0.0943697342,(14:0.0237557178,15:0.0357358071):0.041377):0.005225):0.004424,(16:0.1141193265,((17:0.0363195704,18:0.0188535952):0.043287,19:0.0947695442):0.019148):0.007921):0.022413,(20:0.0698883547,21:0.0811607602):0.020941):0.007439,((22:0.0449896545,(23:0.0331564496,24:0.0230286886):0.058253):0.016766,25:0.0847328612):0.008417):0.006539):0.011557,(26:0.0349772256,27:0.0577536309):0.042499):0.036177):0.016859,28:0.1147218285):0.001271,29:0.1800489093):0.0;

            END;
        """
        ]

    def check_full_dataset_taxon_references(self, dataset):
        self.check_taxon_sets(dataset)
        self.check_tree_lists(dataset)

    def check_taxon_sets(self, dataset):
        self.assertEqual(len(dataset.taxon_sets), len(self.taxon_set_names))
        for i, tax_label in enumerate(self.taxon_set_names):
            self.assertEqual(tax_label, dataset.taxon_sets[i].label)
        self.assertEqual(len(dataset.taxon_sets), len(self.taxon_set_taxon_labels))
        for i, taxon_set in enumerate(dataset.taxon_sets):
            tax_labels = self.taxon_set_taxon_labels[i]
            self.assertEqual(len(tax_labels), len(taxon_set))
            ds_tax_labels = [t.label for t in taxon_set]
            self.assertEqual(tax_labels, ds_tax_labels)

    def check_tree_lists(self, dataset):
        self.assertEqual(self.tree_list_labels, [t.label for t in dataset.tree_lists])
        for t in dataset.tree_lists:
            self.assertEqual(t.taxon_set.label, self.tree_list_taxon_set_labels[t.label])
        self.assertEqual(len(self.tree_list_strings), len(dataset.tree_lists))
        for i, tree_list_str in enumerate(self.tree_list_strings):
            ds_tlist = dataset.tree_lists[i]
            check_tlist = dendropy.TreeList.get_from_string(tree_list_str, "nexus")
            self.assertDistinctButEqual(check_tlist, ds_tlist, distinct_taxa=True)
