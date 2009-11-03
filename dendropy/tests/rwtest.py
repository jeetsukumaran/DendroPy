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
Wraps reading/writing/re-reading tests.
"""

import os
import unittest
from dendropy.utility import messaging
from dendropy import tests


def text_to_expected(text):
    """
    Takes a tab-delimited string in the form of:
        <TAXON_NAME>\t<CHARACTERS>
    and returns a list of pairs with first element the
    taxon name and the second a string of state symbols.
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

def char_array_to_expected(char_array):
    """
    Takes a `char_array` and returns a list of pairs with first element the
    taxon name and the second a string of state symbols.
    """
    data = []
    for t in char_array.taxon_set:
        data.append((t, char_array[t]))
    return data

def compare_labels(label1, label2, underscore_substitution=False):
    if underscore_substitution:
        label1 = label1.replace("_", " ")
        label2 = label2.replace("_", " ")
    assert label1 == label2, "'%s' != '%s'" % (label1, label2)

def check_char_array_parse_against_expected(char_array,
        expected_data,
        expected_type,
        underscore_substitution=False,
        logger=None):
    """
    Takes the following:
        - `char_array` is a CharacterArray,
        - `expected_type` is a the expected subtype of CharacterArray
        - `expected_data` is a list of pairs with
           the first element the taxon name and the second element the string of
           state symbols.

    Returns True if data matches, False if otherwise.
    """

    assert isinstance(char_array, expected_type), "%s != %s" % (type(char_array), expected_type)
    assert len(expected_data) == len(char_array)
    assert len(expected_data) == len(char_array.taxon_set)

    for tax_idx, (exp_label, exp_seq) in enumerate(expected_data):
        taxon = char_array.taxon_set[tax_idx]
        label = taxon.label
        compare_labels(label, exp_label, underscore_substitution)
        assert char_array[taxon] == char_array.taxon_seq_map[taxon]
        assert char_array[tax_idx] == char_array.taxon_seq_map[taxon]
        assert len(exp_seq) == len(char_array.taxon_seq_map[taxon])
        for col_idx, symbol1 in enumerate(exp_seq):
            test_state = char_array.taxon_seq_map[taxon][col_idx].value
            if char_array.taxon_seq_map[taxon][col_idx].column_type is not None:
                state_alpha = char_array.taxon_seq_map[taxon][col_idx].column_type.state_alphabet
            else:
                state_alpha = char_array.default_state_alphabet
            exp_state = state_alpha.state_for_symbol(symbol1)
            assert test_state == exp_state, "Expecting '%s' but found '%s'." % (exp_state, test_state)

def check_canonical_Pythonidae_cytb_tree_parse(reader, srcpath, logger, underscore_substitution=False):
    logger.info("Reading trees from '%s' using '%s' ..." % (srcpath, reader.__class__.__name__))
    dataset = reader.read(istream=open(srcpath, "rU"))
    expected_taxa = [
        "Aspidites_ramsayi",
        "Bothrochilus_boa",
        "Liasis_fuscus",
        "Antaresia_stimsoni",
        "Morelia_viridis",
        "Morelia_bredli",
        "Antaresia_perthensis",
        "Python_timoriensis",
        "Antaresia_maculosa",
        "Morelia_carinata",
        "Python_brongersmai",
        "Morelia_boeleni",
        "Morelia_oenpelliensis"
    ]
    logger.info("Verifying trees ...")
    assert len(dataset.taxon_sets[0]) == len(expected_taxa), \
            "Expecting %d taxa, but found %d." % (len(dataset.taxon_sets[0]), len(expected_taxa))
    if underscore_substitution:
        taxlabels = [t.label.replace(' ', '_') for t in dataset.taxon_sets[0]]
    else:
        taxlabels = [t.label for t in dataset.taxon_sets[0]]
    for et in expected_taxa:
        assert et in taxlabels, \
            "Could not find expected taxon '%s' in observed taxa: '%s'." % (et, taxlabels)
    for ot in taxlabels:
        assert ot in expected_taxa, \
            "Could not find observed taxon '%s' in expected taxa: '%s'." % (ot, taxlabels)
    assert len(dataset.tree_lists) == 1, \
            "Expecting %d tree collections, but found %d." % (1, len(dataset.tree_lists))
    assert len(dataset.tree_lists[0]) == 10, \
            "Expecting %d trees, but found %d." % (10, len(dataset.tree_lists[0]))
    for tree in dataset.tree_lists[0]:
        tree.debug_check_tree(logger)
    logger.info("Trees verified.")



###############################################################################
## Comparison between two phylogenetic data objects that have been instantiated
## separately but from the same data source

def compare_datasets(ds1, ds2, tester):
    tester.logger.info("Comparing dataset taxon sets ...")
    compare_dataset_taxon_sets(ds1, ds2, tester)

    tester.logger.info("Comparing dataset tree lists ...")
    compare_dataset_tree_lists(ds1, ds2, tester)

    tester.logger.info("Comparing dataset character arrays ...")
    compare_dataset_char_arrays(ds1, ds2, tester)

def compare_dataset_taxon_sets(ds1, ds2, tester):
    tester.assertEqual(len(ds1.taxon_sets), len(ds2.taxon_sets))
    for ts_idx, ts1 in enumerate(ds1.taxon_sets):
        ts2 = ds2.taxon_sets[ts_idx]
        tester.logger.info("Comparing taxa of taxon set %d: %d taxa vs. %d taxa" \
            % (ts_idx, len(ts1), len(ts2)))
        compare_individual_taxon_sets(ts1, ts2, tester)

def compare_individual_taxon_sets(ts1, ts2, tester):
    tester.assertEqual(len(ts1), len(ts2))
    for taxon_idx, taxon1 in enumerate(ts1):
        tester.logger.debug("Taxon %d: '%s' == '%s'" % (taxon_idx, taxon1.label, ts2[taxon_idx].label))
        tester.assertEqual(taxon1.label, ts2[taxon_idx].label)

def compare_dataset_tree_lists(ds1, ds2, tester):
    tester.assertEqual(len(ds1.tree_lists), len(ds2.tree_lists))
    for tree_list_idx, tree_list1 in enumerate(ds1.tree_lists):
        tree_list2 = ds2.tree_lists[tree_list_idx]
        compare_individual_tree_lists(tree_list1, tree_list2, tester)

def compare_individual_tree_lists(tree_list1, tree_list2, tester):
    tester.assertEqual(len(tree_list1), len(tree_list2))
    for tree_idx, tree1 in enumerate(tree_list1):
        tree2 = tree_list2[tree_idx]
        tester.logger.debug(tree1.to_newick_str())
        tree1.debug_check_tree(logger=tester.logger)
        tester.logger.debug(tree2.to_newick_str())
        tree2.debug_check_tree(logger=tester.logger)
        tree1_nodes = [nd for nd in tree1.postorder_node_iter()]
        tree2_nodes = [nd for nd in tree2.postorder_node_iter()]
        tester.assertEqual(len(tree1_nodes), len(tree2_nodes))
        for nd_idx, node1 in enumerate(tree1_nodes):
            node2 = tree2_nodes[nd_idx]
            if node1.taxon is not None:
                tester.assert_(node2.taxon is not None)
                tester.assertEqual(node1.taxon.label, node2.taxon.label)
            else:
                tester.assert_(node2.taxon is None)
            if node1.edge.length is not None:
                tester.assert_(node2.edge.length is not None)
                tester.assertAlmostEqual(node1.edge.length, node2.edge.length, 3)
            else:
                tester.assert_(node2.edge.length is None)
            tester.assertEqual(len(node1.child_nodes()), len(node2.child_nodes()))

def compare_dataset_char_arrays(ds1, ds2, tester):
    tester.assertEqual(len(ds1.char_arrays), len(ds2.char_arrays))
    for char_array_idx, char_array1 in enumerate(ds1.char_arrays):
        char_array2 = ds2.char_arrays[char_array_idx]
        compare_individual_char_arrays(char_array1, char_array2, tester)

def compare_individual_char_arrays(char_array1, char_array2, tester):
    tester.assertEqual(len(char_array1), len(char_array2))
    tester.assertEqual(len(char_array1.taxon_set), len(char_array2.taxon_set))
    for taxon_idx, taxon1 in enumerate(char_array1.taxon_set):
        tester.assertEqual(char_array1.taxon_set[taxon_idx].label,
                char_array2.taxon_set[taxon_idx].label)
        seq1 = char_array1[taxon_idx]
        seq2 = char_array2[taxon_idx]
        tester.assertEqual(len(seq1), len(seq2))
        for cell_idx, cell1 in enumerate(seq1):
            cell2 = seq2[cell_idx]
            state1 = cell1.value
            state2 = cell2.value
            tester.assertEqual(state1.symbol, state2.symbol)
            tester.assertEqual(state1.token, state2.token)
            tester.assertEqual(state1.multistate, state2.multistate)
            tester.assertEqual(state1.fundamental_symbols, state2.fundamental_symbols)

###############################################################################
## DatasetReadWriteTest -- round-tripping.

class DatasetReadWriteTest(unittest.TestCase):
    """
    Reads a Dataset, writes the data out, and re-reads the data into a new
    Dataset, checking for identical domain-level information across the first
    and second Dataset.
    """

    def setUp(self, reader_type, writer_type):
        self.reader_type = reader_type
        self.writer_type = writer_type
        self.logger = messaging.get_logger(self.__class__.__name__)

    def dataset_read_write_test(self, src):
        if hasattr(src, "name"):
            src_name = src.name
        else:
            src_name = "STDIN"
        self.logger.info("Read/Write/Re-read/Compare Test: '%s' ..." % os.path.basename(src_name))
        self.logger.info("Using Reader: '%s'" % self.reader_type.__name__)
        self.logger.info("Using Writer: '%s'" % self.writer_type.__name__)

        reader = self.reader_type()
        self.logger.info("Reading '%s' using '%s' ..." \
            % (src_name, reader.__class__.__name__))
        ds1 = reader.read(istream=src)

        writer = self.writer_type(dataset=ds1)
        outfile = open(tests.named_output_file_path(filename="%s.rwtest" % os.path.basename(src_name)), "w")
        assert outfile is not None
        self.logger.info("Writing to '%s' using '%s' ..." % (outfile.name, writer.__class__.__name__))
        writer.write(ostream=outfile)
        outfile.flush()

        reader = self.reader_type()
        infile = open(outfile.name, "r")
        self.logger.info("Re-reading '%s' using '%s' ..." \
            % (infile.name, reader.__class__.__name__))
        ds2 = reader.read(istream=infile)

        compare_datasets(ds1, ds2, self)
