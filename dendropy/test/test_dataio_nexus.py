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
        self.logger.info("Checking '%s' => %s" % (data_filename, datatype.__name__))
        reader = nexus.NexusReader()
        dataset = reader.read(stream=pathmap.char_source_stream(data_filename))
        expected_label_symbol_stream = pathmap.char_source_stream(expected_filename)
        self.assertEqual(len(dataset.char_arrays), 1)
        self.assertEqualCharArrayLabelSymbols(dataset.char_arrays[0], \
            expected_label_symbol_stream=expected_label_symbol_stream)

class NexusParseDnaCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("pythonidae_cytb.chars.nexus",
                "pythonidae_cytb.chars.txt",
                dendropy.DnaCharacterArray)

class NexusParseDnaCharsInterleavedTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("pythonidae_cytb.chars.interleaved.nexus",
                "pythonidae_cytb.chars.txt",
                dendropy.DnaCharacterArray)

class NexusParseProteinCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("caenophidia_mos.chars.nexus",
                "caenophidia_mos.chars.txt",
                dendropy.ProteinCharacterArray)

class NexusParseStandardCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("angiosperms.chars.nexus",
                "angiosperms.chars.txt",
                dendropy.ProteinCharacterArray)

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

    def check_parse_with_ambiguities(self, data_filename, expected_filename):
        reader = nexus.NexusReader()
        dataset = reader.read(stream=pathmap.char_source_stream(data_filename))
        self.assertEqual(len(dataset.char_arrays), 1)
        self.map_multistate_to_symbols(dataset.char_arrays[0])
        expected_label_symbol_stream = pathmap.char_source_stream(expected_filename)
        self.assertEqualCharArrayLabelSymbols(dataset.char_arrays[0], \
            expected_label_symbol_stream = expected_label_symbol_stream)

    def testStandardWithMultistateInBraces(self):
        self.check_parse_with_ambiguities("apternodus.chars.nexus", "apternodus.chars.hacked-for-tests.txt")

    def testStandardWithMultistateInBracesInterleaved(self):
        self.check_parse_with_ambiguities("apternodus.chars.interleaved.nexus", "apternodus.chars.hacked-for-tests.txt")

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

class NexusTreeListWriterTest(datatest.DataObjectVerificationTestCase):

    def testWriteTreeListDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        output_path = pathmap.named_output_path(filename="reference.trees.out.nexus", suffix_timestamp=True)
        nexus.write_tree_list(ref_tree_list, stream=open(output_path, "w"))
        t_tree_list = nexus.read_tree_list(stream=open(output_path, "rU"))
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None)

class MultiTreeSourceIterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.ref_tree_list = datagen.reference_tree_list()
        self.ref_idx = -1

    def next_ref_tree(self):
        self.ref_idx += 1
        if self.ref_idx >= len(self.ref_tree_list):
            self.ref_idx = 0
        return self.ref_tree_list[ref_idx]

class NexusOrNewickTreeSourceIterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.ref_tree_list = datagen.reference_tree_list()

    def testNexusDistinctTaxa(self):
        stream = pathmap.tree_source_stream('reference.trees.nexus')
        for idx, test_tree in enumerate(nexus.generalized_tree_source_iter(stream=stream)):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=True)

    def testNexusSameTaxa(self):
        stream = pathmap.tree_source_stream('reference.trees.nexus')
        for idx, test_tree in enumerate(nexus.generalized_tree_source_iter(stream=stream, taxon_set=self.ref_tree_list.taxon_set)):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=False)

    def testNewickDistinctTaxa(self):
        stream = pathmap.tree_source_stream('reference.trees.newick')
        for idx, test_tree in enumerate(nexus.generalized_tree_source_iter(stream=stream)):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=True, ignore_taxon_order=True)

    def testNewickSameTaxa(self):
        stream = pathmap.tree_source_stream('reference.trees.newick')
        for idx, test_tree in enumerate(nexus.generalized_tree_source_iter(stream=stream, taxon_set=self.ref_tree_list.taxon_set)):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=False)

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

class NexusDocumentReadWriteTest(datatest.DataObjectVerificationTestCase):

    def testRoundTripReference(self):
        reference_dataset = datagen.reference_single_taxonset_dataset()
        self.roundTripDataSetTest(reference_dataset, "nexus", ignore_columns=True)

    def testRoundTripProtein(self):
        s = pathmap.char_source_stream("caenophidia_mos.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="nexus")
        self.roundTripDataSetTest(d1, "nexus", ignore_columns=True)

    def testRoundTripStandard1(self):
        s = pathmap.char_source_stream("angiosperms.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="nexus")
        self.roundTripDataSetTest(d1, "nexus", ignore_columns=True)

    def testRoundTripStandard2(self):
        s = pathmap.char_source_stream("apternodus.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="nexus")
        self.roundTripDataSetTest(d1, "nexus", ignore_columns=True)

if __name__ == "__main__":
    unittest.main()
