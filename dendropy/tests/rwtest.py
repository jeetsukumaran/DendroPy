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
from dendropy.tests import datacompare

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

        datacompare.compare_datasets(ds1, ds2, self)
