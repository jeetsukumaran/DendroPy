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
NEXUS data read/write parse/format tests.
"""

import sys
import os
import unittest
import tempfile

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
import dendropy
from dendropy.dataio import nexus

class NexusGeneralParseCharsTest(datatest.DataObjectVerificationTestCase):

    def check_chars_against_expected(self, data_filename, expected_filename, datatype):
        reader = nexus.NexusReader()
        dataset = reader.read(stream=pathmap.char_source_stream(data_filename))
        expected_label_symbol_stream = pathmap.char_source_stream(expected_filename)
        self.assertEqual(len(dataset.char_arrays), 1)
        self.assertEqualCharArrayLabelSymbols(dataset.char_arrays[0], \
            expected_label_symbol_stream=expected_label_symbol_stream)

    def testCharParse(self):
        test_sets = [
            ["pythonidae_cytb.chars.nexus", "pythonidae_cytb.chars.txt", dendropy.DnaCharacterArray],
            ["caenophidia_mos.chars.nexus", "caenophidia_mos.chars.txt", dendropy.ProteinCharacterArray],
            ["angiosperms.chars.nexus", "angiosperms.chars.txt", dendropy.StandardCharacterArray],
        ]
        for t in test_sets:
            self.logger.info("Checking '%s' => %s" % (t[0], t[2].__name__))
            self.check_chars_against_expected(t[0], t[1], t[2])

class NexusParseStandardCharsWithMultistateTest(datatest.DataObjectVerificationTestCase):
    """
    This tests the capability of the NEXUS parser in handling "{}" and
    "()" constructs in the data. Two files are used, one in which the
    ambiguous data are marked up using "{}" and "()" constructs, and the other
    in which these are substituted by symbols representing the appropriate
    multistate. The first file is parsed, and the result character array's
    state alphabet is hacked to map the ambiguous states to the symbols used
    in the second file. The resulting label-symbol lists are then compared.
    """

    def map_multistate_to_symbols(self, char_array):
        self.assertEqual(len(char_array.state_alphabets), 1)
        sa = char_array.state_alphabets[0]
        for sae in sa:
            if sae.multistate != dendropy.StateAlphabetElement.SINGLE_STATE \
                    and sae.symbol is None:
                member_symbols = sae.fundamental_symbols
                if member_symbols == set('01'):
                    sae.symbol = 'J'
                elif member_symbols == set('023'):
                    sae.symbol = 'K'
                elif member_symbols == set('12'):
                    sae.symbol = 'L'
                elif member_symbols == set('23'):
                    sae.symbol = 'M'
                elif member_symbols == set('13'):
                    sae.symbol = 'N'
                elif member_symbols == set('02'):
                    sae.symbol = 'P'
                elif member_symbols == set('03'):
                    sae.symbol = 'Q'
                elif member_symbols == set('012'):
                    sae.symbol = 'R'
                else:
                    raise self.failureException("Unexpected multistate: %s" % member_symbols)
        return sae

    def testStandardWithMultistateInBraces(self):
        reader = nexus.NexusReader()
        dataset = reader.read(stream=pathmap.char_source_stream("apternodus.chars.nexus"))
        self.assertEqual(len(dataset.char_arrays), 1)
        self.map_multistate_to_symbols(dataset.char_arrays[0])
        expected_label_symbol_stream = pathmap.char_source_stream("apternodus.chars.hacked-for-tests.txt")
        self.assertEqualCharArrayLabelSymbols(dataset.char_arrays[0], \
            expected_label_symbol_stream = expected_label_symbol_stream)

class NexusTreeListReaderTest(datatest.DataObjectVerificationTestCase):

    def testReferenceTreeFileDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = nexus.read_tree_list(stream=pathmap.tree_source_stream("reference.trees.nexus"))
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = nexus.read_tree_list(stream=pathmap.tree_source_stream("reference.trees.nexus"), taxon_set=ref_tree_list.taxon_set)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=False,
                equal_oids=None)

class NexusTreeDocumentReaderTest(datatest.DataObjectVerificationTestCase):

    def testReferenceTreeFileDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = nexus.NexusReader()
        dataset = reader.read(stream=pathmap.tree_source_stream("reference.trees.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=True,
                equal_oids=None)

#class NexusDocumentTest(rwtest.DatasetReadWriteTest):
#
#    def setUp(self):
#        rwtest.DatasetReadWriteTest.setUp(self,
#            reader_type=nexus.NexusReader,
#            writer_type=nexus.NexusWriter)
#
#    def testReadWriteStandardAndTrees(self):
#        datafiles = [
#            "angiosperms.char_and_trees.nex",
#            "apternodus.nex",
#            "pythonidae_cytb.nex",
#            "caenophidia_mos.nex",
#            "pythonidae_cytb.nexus.tre"
#        ]
#        for datafile in datafiles:
#            self.dataset_read_write_test(open(tests.data_source_path(datafile)))

if __name__ == "__main__":
    unittest.main()
