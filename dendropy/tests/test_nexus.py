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

import dendropy
from dendropy.utility import messaging
from dendropy import tests
from dendropy.tests import rwtest
from dendropy.dataio import nexus

_LOG = messaging.get_logger(__name__)

class NexusCharTest(unittest.TestCase):

    def check_char_io(self, data_filepath, reference_filepath, datatype):
        _LOG.info("Reading ...")
        reader = nexus.NexusReader()
        data = reader.read(path=tests.data_source_path(data_filepath))
        expected = rwtest.text_to_expected(
                open(tests.data_source_path(reference_filepath), "rU").read())
        assert len(data.taxon_sets[0]) == len(expected), \
                "%d != %d" % (len(data.taxon_sets[0]), len(expected))
        assert len(data.char_arrays) == 1
        rwtest.check_char_array_parse_against_expected(data.char_arrays[0],
                expected, \
                datatype,
                underscore_substitution=True,
                logger=_LOG)

        _LOG.info("Writing ...")
        writer = nexus.NexusWriter()

    def testCharIO(self):
        test_sets = [
            ["pythonidae_cytb.nex", "pythonidae_cytb.txt", dendropy.DnaCharacterArray],
            ["caenophidia_mos.nex", "caenophidia_mos.txt", dendropy.ProteinCharacterArray],
            ["angiosperms.nex", "angiosperms.txt", dendropy.StandardCharacterArray],
#             ["apternodus.nex", "apternodus.txt", dendropy.StandardCharacterArray],
        ]
        for t in test_sets:
            _LOG.info("Checking '%s' => %s" % (t[0], t[2].__name__))
            self.check_char_io(t[0], t[1], t[2])

    def testStandardWithAmbiguities(self):
        reader = nexus.NexusReader()

        def special_alphabet_builder(char_block, symbols):
            sa = dendropy.StateAlphabet()
            sa.append(dendropy.StateAlphabetElement(symbol='0'))
            sa.append(dendropy.StateAlphabetElement(symbol='1'))
            sa.append(dendropy.StateAlphabetElement(symbol='2'))
            sa.append(dendropy.StateAlphabetElement(symbol='3'))
            sa.append(dendropy.StateAlphabetElement(symbol='4'))
            sa.append(dendropy.StateAlphabetElement(symbol='-',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='01234')))
            sa.append(dendropy.StateAlphabetElement(symbol='?',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='01234')))
            sa.append(dendropy.StateAlphabetElement(symbol='J',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='01')))
            sa.append(dendropy.StateAlphabetElement(symbol='K',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='023')))
            sa.append(dendropy.StateAlphabetElement(symbol='L',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='12')))
            sa.append(dendropy.StateAlphabetElement(symbol='M',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='23')))
            sa.append(dendropy.StateAlphabetElement(symbol='N',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='13')))
            sa.append(dendropy.StateAlphabetElement(symbol='P',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='02')))
            sa.append(dendropy.StateAlphabetElement(symbol='Q',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='03')))
            sa.append(dendropy.StateAlphabetElement(symbol='R',
                    multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                    member_states=sa.get_states(symbols='012')))
            char_block.state_alphabets = [sa]
            char_block.default_state_alphabet = char_block.state_alphabets[0]

        reader.build_state_alphabet = special_alphabet_builder
        data = reader.read(path=tests.data_source_path("apternodus.nex"))

        expected = rwtest.text_to_expected(
                open(tests.data_source_path("apternodus.hacked-for-tests.txt"), "rU").read())
        assert len(data.taxon_sets) == 1
        assert len(data.taxon_sets[0]) == len(expected), \
                "%d != %d" % (len(data.taxon_sets), len(expected))
        assert len(data.char_arrays) == 1,  \
            "Expecting 1 CharacterArray but found %d.\n" % (len(data.char_arrays))
        rwtest.check_char_array_parse_against_expected(data.char_arrays[0],
                expected, \
                dendropy.StandardCharacterArray,
                underscore_substitution=True,
                logger=_LOG)

class NexusTreeTest(unittest.TestCase):

    def testReadTreeList(self):
        rwtest.check_canonical_Pythonidae_cytb_tree_parse(reader = nexus.NexusReader(),
            srcpath=tests.data_source_path("pythonidae_cytb.nexus.tre"),
            logger=_LOG,
            underscore_substitution=True)

    def testWriteTreeList(self):
        _LOG.info("Reading in trees for NEXUS writing test")
        reader = nexus.NexusReader()
        ds1 = reader.read(path=tests.data_source_path("pythonidae_cytb.nexus.tre"))

        outfile = tempfile.NamedTemporaryFile()
        _LOG.info("Writing trees to temporary file '%s'" % outfile.name)
        writer = nexus.NexusWriter(dataset=ds1)
        writer.write(file=outfile)
        outfile.flush()

        _LOG.info("Re-reading trees")
        rwtest.check_canonical_Pythonidae_cytb_tree_parse(
            reader = nexus.NexusReader(),
            srcpath=outfile.name,
            logger=_LOG,
            underscore_substitution=True)

    def read_write_test(self, filename):
        _LOG.info("NEXUS Read/Write Tests: '%s'" % os.path.basename(filename))
        dataset_read_write_test(
            reader_type=nexus.NexusReader,
            writer_type=nexus.NexusWriter,
            file=open(tests.data_source_path(filename), "rU"),
            logger=_LOG)

class NexusDocumentTest(rwtest.DatasetReadWriteTest):

    def setUp(self):
        rwtest.DatasetReadWriteTest.setUp(self,
            reader_type=nexus.NexusReader,
            writer_type=nexus.NexusWriter)

    def testReadWriteStandardAndTrees(self):
        datafiles = [
            "angiosperms.char_and_trees.nex",
            "apternodus.nex",
            "pythonidae_cytb.nex",
            "caenophidia_mos.nex",
            "pythonidae_cytb.nexus.tre"
        ]
        for datafile in datafiles:
            self.dataset_read_write_test(open(tests.data_source_path(datafile)))

if __name__ == "__main__":
    unittest.main()
