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
NEXUS data read/write parse/format tests.
"""

import sys
import os
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.test.support import extendedtest
from dendropy import dataio
from dendropy.dataio.nexustokenizer import TooManyTaxaError
from dendropy.utility import error
from dendropy.utility.messaging import get_logger
import dendropy

_LOG = get_logger(__name__)

class NexusGeneralParseCharsTest(datatest.DataObjectVerificationTestCase):

    def check_chars_against_expected(self, data_filename, expected_filename, datatype):
        self.logger.info("Checking '%s' => %s" % (data_filename, datatype.__name__))
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.char_source_stream(data_filename))
        expected_label_symbol_stream = pathmap.char_source_stream(expected_filename)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqualCharMatrixLabelSymbols(dataset.char_matrices[0], \
            expected_label_symbol_stream=expected_label_symbol_stream)

    def check_continuous_chars_against_expected(self, data_filename, expected_filename, datatype):
        self.logger.info("Checking '%s' => %s" % (data_filename, datatype.__name__))
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.char_source_stream(data_filename))
        expected_label_symbol_stream = pathmap.char_source_stream(expected_filename)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqualCharMatrixLabelContinuousValues(dataset.char_matrices[0], \
            expected_label_symbol_stream=expected_label_symbol_stream)

class NexusParseDnaCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("pythonidae.chars.nexus",
                "pythonidae.chars.txt",
                dendropy.DnaCharacterMatrix)

class NexusParseDnaCharsInterleavedTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("pythonidae.chars.interleaved.nexus",
                "pythonidae.chars.txt",
                dendropy.DnaCharacterMatrix)

class NexusParseProteinCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("caenophidia_mos.chars.nexus",
                "caenophidia_mos.chars.txt",
                dendropy.ProteinCharacterMatrix)

class NexusParseProteinMatchCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("avian-ovomucoids.chars.nexus",
                "avian-ovomucoids.chars.txt",
                dendropy.ProteinCharacterMatrix)

class NexusParseStandardCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_chars_against_expected("angiosperms.chars.nexus",
                "angiosperms.chars.txt",
                dendropy.ProteinCharacterMatrix)

class NexusParseContinuousCharsTest(NexusGeneralParseCharsTest):

    def runTest(self):
        self.check_continuous_chars_against_expected("pythonidae_continuous.chars.nexus",
                "pythonidae_continuous.chars.txt",
                dendropy.ContinuousCharacterMatrix)

class NexusParseStandardCharsWithMultistateTest(datatest.DataObjectVerificationTestCase):
    """
    This tests the capability of the NEXUS parser in handling "{}" and
    "()" constructs in the data. Two files are used, one in which the
    ambiguous data are marked up using "{}" and "()" constructs, and the other
    in which these are substituted by symbols representing the appropriate
    multistate. The first file is parsed, and the result character matrix's
    state alphabet is hacked to map the ambiguous states to the symbols used
    in the second file. The resulting label-symbol lists are then compared.
    """

    def map_multistate_to_symbols(self, char_matrix):
        self.assertEqual(len(char_matrix.state_alphabets), 1)
        sa = char_matrix.state_alphabets[0]
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
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.char_source_stream(data_filename))
        self.assertEqual(len(dataset.char_matrices), 1)
        self.map_multistate_to_symbols(dataset.char_matrices[0])
        expected_label_symbol_stream = pathmap.char_source_stream(expected_filename)
        self.assertEqualCharMatrixLabelSymbols(dataset.char_matrices[0], \
            expected_label_symbol_stream = expected_label_symbol_stream)

    def testStandardWithMultistateInBraces(self):
        self.check_parse_with_ambiguities("apternodus.chars.nexus", "apternodus.chars.hacked-for-tests.txt")

    def testStandardWithMultistateInBracesInterleaved(self):
        self.check_parse_with_ambiguities("apternodus.chars.interleaved.nexus", "apternodus.chars.hacked-for-tests.txt")

class NexusTreeListReaderTest(datatest.DataObjectVerificationTestCase):

    def testReferenceTreeFileDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_stream(pathmap.tree_source_stream(datagen.reference_trees_filename(schema="nexus")), "nexus")
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_stream(pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-block.nexus"), "nexus", taxon_set=ref_tree_list.taxon_set)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=False,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_stream(pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-block.nexus"), "nexus")
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-block.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockNoTranslateBlockDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_stream(pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-no-translate-block.nexus"), "nexus")
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockNoTranslateBlockDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-no-translate-block.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=True,
                ignore_taxon_order=True,
                equal_oids=None)

class NexusTreeListWriterTest(datatest.DataObjectVerificationTestCase):

    def testWriteTreeListDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        output_path = pathmap.named_output_path(filename="reference.trees.out.nexus", suffix_timestamp=True)
        ref_tree_list.write_to_path(output_path, "nexus")
        t_tree_list = dendropy.TreeList.get_from_path(output_path, "nexus")
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None)

class ExtraSemiColonsNexusTest(unittest.TestCase):

    def setUp(self):
        self.str1 = """\
#NEXUS
begin trees;
    ;
    tree 1 = ((A,B),(C,D));;;
    ;
    tree 2 = ((A,(B,(C,D))));;;
    ;
    tree 3 = ((A,C),(D,B));;;
end;
"""

        self.str2 = """\
#NEXUS
begin trees;
    ;
    tree 1 = ();;
    ;
    tree 2 = ();;;
    tree 3 = ();;
end;
"""

        self.str3 = """\
#NEXUS
begin trees;
    ;
    tree 1 = ;;
    ;
    tree 2 = ;;;
    tree 3 = ;;
end;
"""

    def testStr1AsDoc(self):
        tlist = dendropy.TreeList.get_from_string(self.str1, "nexus")
        self.assertEqual(len(tlist), 3)

    def testStr1Iter(self):
        for t in dendropy.dataio.tree_source_iter(StringIO(self.str1), "nexus"):
            _LOG.info(t.as_string("newick"))

    def testStr2AsDoc(self):
        self.assertRaises(error.DataParseError, dendropy.TreeList.get_from_string, self.str2, "nexus")

#    def testStr2Iter(self):
#        for t in dendropy.dataio.tree_source_iter(StringIO(self.str2), "nexus"):
#            _LOG.info(t.as_string("newick"))

    def testStr3AsDoc(self):
#        self.assertRaises(error.DataParseError, dendropy.TreeList.get_from_string, self.str3, "nexus")
        tlist = dendropy.TreeList.get_from_string(self.str3, "nexus")
        _LOG.info(tlist.as_string("nexus"))

    def testStr3Iter(self):
        for t in dendropy.dataio.tree_source_iter(StringIO(self.str3), "nexus"):
            _LOG.info(t.as_string("newick"))

class MultiTreeSourceIterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.ref_tree_list = datagen.reference_tree_list()
        self.ref_index = 0

    def next_ref_tree(self, restart_index=0):
        t = self.ref_tree_list[self.ref_index]
        self.ref_index += 1
        if self.ref_index >= len(self.ref_tree_list):
            self.ref_index = restart_index
        return t

    def testMixedNexusAndNewickDistinctTaxa(self):
        filenames = [datagen.reference_trees_filename(schema="newick"),
                     datagen.reference_trees_filename(schema="nexus"),
                     datagen.reference_trees_filename(schema="newick"),
                     datagen.reference_trees_filename(schema="nexus")]
        filepaths = [pathmap.tree_source_path(f) for f in filenames]
        taxon_set = dendropy.TaxonSet()
        for idx, test_tree in enumerate(dataio.multi_tree_source_iter(filepaths, schema="nexus/newick", taxon_set=taxon_set)):
            self.assertDistinctButEqualTree(self.next_ref_tree(), test_tree, distinct_taxa=True, ignore_taxon_order=True)
        self.assertEqual(idx, 43)

    def testMixedNexusAndNewickSameTaxa(self):
        filenames = [datagen.reference_trees_filename(schema="newick"),
                     datagen.reference_trees_filename(schema="nexus"),
                     datagen.reference_trees_filename(schema="newick"),
                     datagen.reference_trees_filename(schema="nexus")]
        filepaths = [pathmap.tree_source_path(f) for f in filenames]
        taxon_set = self.ref_tree_list.taxon_set
        for idx, test_tree in enumerate(dataio.multi_tree_source_iter(filepaths, schema="nexus/newick", taxon_set=taxon_set)):
            self.assertDistinctButEqualTree(self.next_ref_tree(), test_tree, distinct_taxa=False, ignore_taxon_order=True)
        self.assertEqual(idx, 43)

    def testBurnIn(self):
        filenames = [datagen.reference_trees_filename(schema="newick"),
                     datagen.reference_trees_filename(schema="nexus"),
                     datagen.reference_trees_filename(schema="newick"),
                     datagen.reference_trees_filename(schema="nexus")]
        filepaths = [pathmap.tree_source_path(f) for f in filenames]
        taxon_set = self.ref_tree_list.taxon_set
        self.ref_index = 5
        for idx, test_tree in enumerate(dataio.multi_tree_source_iter(filepaths,
                schema="nexus/newick",
                taxon_set=taxon_set,
                tree_offset=5)):
            check_tree = self.next_ref_tree(restart_index=5)
            self.assertDistinctButEqualTree(check_tree, test_tree, distinct_taxa=False, ignore_taxon_order=True)
        self.assertEqual(idx, 23)

class NexusOrNewickTreeSourceIterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.ref_tree_list = datagen.reference_tree_list()

    def testNexusDistinctTaxa(self):
        stream = pathmap.tree_source_stream(datagen.reference_trees_filename(schema="nexus"))
        for idx, test_tree in enumerate(dataio.tree_source_iter(stream=stream, schema='nexus/newick')):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=True)

    def testNexusSameTaxa(self):
        stream = pathmap.tree_source_stream(datagen.reference_trees_filename(schema="nexus"))
        for idx, test_tree in enumerate(dataio.tree_source_iter(stream=stream, schema='nexus/newick', taxon_set=self.ref_tree_list.taxon_set)):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=False)

    def testNewickDistinctTaxa(self):
        stream = pathmap.tree_source_stream(datagen.reference_trees_filename(schema="newick"))
        for idx, test_tree in enumerate(dataio.tree_source_iter(stream=stream, schema='nexus/newick')):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=True, ignore_taxon_order=True)

    def testNewickSameTaxa(self):
        stream = pathmap.tree_source_stream(datagen.reference_trees_filename(schema="newick"))
        for idx, test_tree in enumerate(dataio.tree_source_iter(stream=stream, schema='nexus/newick', taxon_set=self.ref_tree_list.taxon_set)):
            self.assertDistinctButEqualTree(self.ref_tree_list[idx], test_tree, distinct_taxa=False)

class NexusTreeDocumentReaderTest(datatest.DataObjectVerificationTestCase):

    def testReferenceTreeFileDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(schema="nexus")))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus', taxon_set=ref_tree_list.taxon_set)
        dataset = reader.read(stream=pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-block.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=False,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-block.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=True,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockNoTranslateBlockSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus', taxon_set=ref_tree_list.taxon_set)
        dataset = reader.read(stream=pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-no-translate-block.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=False,
                equal_oids=None)

    def testReferenceTreeFileNoTaxaBlockNoTranslateBlockDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.tree_source_stream("pythonidae.reference-trees.no-taxa-no-translate-block.nexus"))
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                dataset.tree_lists[0],
                distinct_taxa=True,
                ignore_taxon_order=True,
                equal_oids=None)

class NexusDocumentReadWriteTest(datatest.DataObjectVerificationTestCase):

    def testRoundTripReference(self):
        reference_dataset = datagen.reference_single_taxonset_dataset()
        self.roundTripDataSetTest(reference_dataset, "nexus")

    def testRoundTripProtein(self):
        s = pathmap.char_source_stream("caenophidia_mos.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexus")

    def testRoundTripStandard1(self):
        s = pathmap.char_source_stream("angiosperms.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexus")

    def testRoundTripStandard2(self):
        s = pathmap.char_source_stream("apternodus.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexus")

    def testRoundTripContinuous(self):
        s = pathmap.char_source_stream("pythonidae_continuous.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexus")

class MesquiteNexusMultiTaxaTest(datatest.ComplexMultiTaxonSetDataVerificationTest):

    def testParseMesquiteMultiTaxa(self):
        reader = dataio.get_reader('nexus')
        dataset = reader.read(stream=pathmap.mixed_source_stream("multitaxa_mesquite.nex"))
        self.check_full_dataset_taxon_references(dataset)

    def testRoundTripMesquiteMultiTaxa(self):
        s = pathmap.mixed_source_stream("multitaxa_mesquite.nex")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexus")

class NexusInterleavedWhitespace(extendedtest.ExtendedTestCase):

    def testSpaces1(self):
        tax_labels = (" A  1   ", " B  2   ", " C 3   ", " D  4   ")
        data_rows = [ """#NEXUS  """,
                      """   begin taxa;   """,
                      """  dimensions ntax = 4; """,
                      """ taxlabels %s; """ % (" ".join([("'%s'") % x for x in tax_labels])),
                      """ end; """,
                      """begin characters;""",
                      """  dimensions nchar=8;""",
                      """  format interleave datatype=dna; """,
                      """    matrix  """,
                      """ '%s' {ACTG}   C  T [some comment] G      """ % tax_labels[0],
                      """ '%s' A [another comment]  C  T  G      """ % tax_labels[1],
                      """ '%s' [starting with a comment] A   C {A T C}  G      """ % tax_labels[2],
                      """ '%s' -   C  T  G      """ % tax_labels[3],
                      """    """,
                      """ '%s' A   C  (T  C G) G      """ % tax_labels[0],
                      """ '%s' A   C  T  G      """ % tax_labels[1],
                      """ '%s' {C  A}   ?  -  -      """ % tax_labels[2],
                      """ '%s' -   C  T  G      """ % tax_labels[3],
                      """ ; """,
                      """ end; """,
                    ]
        data_src = "\n".join(data_rows)
        data = dendropy.DnaCharacterMatrix.get_from_string(data_src, "nexus")

        expected_symbols = {
            " A  1   "  : "NCTGACBG",
            " B  2   "  : "ACTGACTG",
            " C 3   "  : "ACHGM?--",
            " D  4   "  : "-CTG-CTG",
        }

        self.assertEqual(len(data.taxon_set), len(tax_labels))
        self.assertEqual(len(data), len(tax_labels))
        for i, t in enumerate(data.taxon_set):
            self.assertEqual(t.label, tax_labels[i])
            self.assertIn(t, data)
            self.assertIs(data[i], data[t])
            s1 = data[t].symbols_as_list()
            s2 = expected_symbols[t.label]
            self.assertEqual(len(s1), len(s2))
            for j, c in enumerate(s1):
                self.assertEqual(c, s2[j])

class NexusTaxaCaseInsensitivityTest(extendedtest.ExtendedTestCase):

    def setUp(self):
        self.data_str = """\
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=5;
    TAXLABELS AAA BBB CCC DDD EEE;
END;

BEGIN CHARACTERS;
    DIMENSIONS  NCHAR=8;
    FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=. INTERLEAVE;
    MATRIX
        AAA ACGT
        BBB ACGT
        CCC ACGT
        DDD ACGT
        EEE ACGT

        aaa ACGT
        bbb ACGT
        ccc ACGT
        ddd ACGT
        eee ACGT
    ;
END;
"""
        self.tree_str = """\
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=5;
    TAXLABELS AAA BBB CCC DDD EEE;
END;

BEGIN TREES;
    TREE 1 = (AAA, (bbb, (ccc, (ddd, eee))));
END;
"""

    def testCaseInsensitiveTrees(self):
        d = dendropy.TreeList.get_from_string(self.tree_str,
                'nexus',
                case_sensitive_taxon_labels=False)
        self.assertEqual(len(d.taxon_set), 5)

    def testDefaultCaseSensitivityTrees(self):
        d = dendropy.TreeList.get_from_string(self.tree_str,
                'nexus')
        self.assertEqual(len(d.taxon_set), 5)

    def testCaseSensitiveTrees(self):
        d = dendropy.TreeList.get_from_string(self.tree_str,
                'nexus',
                case_sensitive_taxon_labels=True)
        self.assertEqual(len(d.taxon_set), 9)

    def testCaseInsensitiveChars(self):
        d = dendropy.DnaCharacterMatrix.get_from_string(self.data_str, 'nexus', case_sensitive_taxon_labels=False)
        self.assertEqual(len(d.taxon_set), 5)

    def testCaseSensitiveChars(self):
        #d = dendropy.DnaCharacterMatrix.get_from_string(self.data_str, 'nexus', case_sensitive_taxon_labels=False)
        self.assertRaises(error.DataParseError,
                dendropy.DnaCharacterMatrix.get_from_string,
                self.data_str,
                'nexus',
                case_sensitive_taxon_labels=True)

    def testDefaultCaseSensitivityChars(self):
        d = dendropy.DnaCharacterMatrix.get_from_string(self.data_str, 'nexus')
        self.assertEqual(len(d.taxon_set), 5)


class NexusTooManyTaxaTest(extendedtest.ExtendedTestCase):

    def testTooManyTaxaNonInterleaved(self):
        data_str = """\
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=5;
    TAXLABELS AAA BBB CCC DDD EEE;
END;

BEGIN CHARACTERS;
    DIMENSIONS  NCHAR=8;
    FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
    MATRIX
        AAA ACGTACGT
        BBB ACGTACGT
        CCC ACGTACGT
        DDD ACGTACGT
        EEE ACGTACGT
        FFF ACGTACGT
    ;
END;
"""
        #d = dendropy.DnaCharacterMatrix.get_from_string(data_str, 'nexus')
        self.assertRaises(TooManyTaxaError, dendropy.DnaCharacterMatrix.get_from_string, data_str, 'nexus')

    def testTooManyTaxaInterleaved(self):
        data_str = """\
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=5;
    TAXLABELS AAA BBB CCC DDD EEE;
END;

BEGIN CHARACTERS;
    DIMENSIONS  NCHAR=8;
    FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=. INTERLEAVE;
    MATRIX
        AAA ACGT
        BBB ACGT
        CCC ACGT
        DDD ACGT
        EEE ACGT

        AAA ACGT
        BBB ACGT
        CCC ACGT
        DDD ACGT
        EEE ACGT
        FFF ACGT
    ;
END;
"""
        #d = dendropy.DnaCharacterMatrix.get_from_string(data_str, 'nexus')
        self.assertRaises(TooManyTaxaError, dendropy.DnaCharacterMatrix.get_from_string, data_str, 'nexus')

class NexusInternalNodeLabelsAsNonTaxa(extendedtest.ExtendedTestCase):

    def setUp(self):
        self.trees_string = """\
#NEXUS

begin taxa;
    dimensions ntax=5;
    taxlabels aaa bbb ccc ddd eee;
end;

begin trees;
    translate
        1  aaa,
        2  bbb,
        3  ccc,
        4  ddd,
        5  eee
    ;
    tree 0 = (((1, 2)1,(3,4)1),5)1;
end;
"""

    def testParseWithKeyword(self):
        t = dendropy.TreeList.get_from_string(self.trees_string, 'nexus', suppress_internal_node_taxa=True)

class NexusCharsSubsetsTest(datatest.DataObjectVerificationTestCase):

    def verify_subsets(self, src_filename, expected_sets):
        """
        `src_filename` -- name of file containing full data and charsets
                          statement
        `expected_sets` -- dictionary with keys = label of charset, and values
                           = name of file with subset of characters correspond
                           to the charset.
        """
        _LOG.debug(src_filename)
        src_data = dendropy.DnaCharacterMatrix.get_from_path(
                pathmap.char_source_path(src_filename),
                'nexus')
        state_alphabet = src_data.default_state_alphabet
        self.assertEqual(len(src_data.character_subsets), len(expected_sets))
        for label, expected_data_file in expected_sets.items():

            _LOG.debug(label)

            self.assertTrue(label in src_data.character_subsets)
            result_subset = src_data.export_character_subset(label)
            expected_subset = dendropy.DnaCharacterMatrix.get_from_path(
                pathmap.char_source_path(expected_data_file),
                'nexus')

            # confirm subset is correct
            self.assertDistinctButEqualDiscreteCharMatrix(result_subset, expected_subset)

            # mutate new and confirm that old remains unchanged
            e1_symbols = src_data[0].symbols_as_string()
            r1 = result_subset[0]
            dummy_state = state_alphabet.state_for_symbol('A')
            for idx in range(len(r1)):
                r1[idx].value = dummy_state
            self.assertEqual(e1_symbols, src_data[0].symbols_as_string())

            # mutate old and confirm that new remains unchanged
            r2_symbols = result_subset[1].symbols_as_string()
            e2 = src_data[1]
            dummy_state = state_alphabet.state_for_symbol('A')
            for idx in range(len(e2)):
                e2[idx].value = dummy_state
            self.assertEqual(r2_symbols, result_subset[1].symbols_as_string())

    def testNonInterleaved(self):
        """
        Charsets here go through all forms of position specification.
        """
        expected_sets = {
                "coding" : "primates.chars.subsets-coding.nexus",
                "noncoding" : "primates.chars.subsets-noncoding.nexus",
                "1stpos" : "primates.chars.subsets-1stpos.nexus",
                "2ndpos" : "primates.chars.subsets-2ndpos.nexus",
                "3rdpos" : "primates.chars.subsets-3rdpos.nexus",
                }
        self.verify_subsets('primates.chars.subsets-all.nexus', expected_sets)

    def testInterleaved(self):
        """
        A bug in DendroPy resulted in the block immediately following an
        interleaved character matrix DATA or CHARACTERS block being skipped.
        This tests for it by ensuring that the ASSUMPTIONS block following an
        interleaved CHARACTERS block is parsed. A better test would approach
        the issue more directly, by checking to see if block parsing left the
        stream reader in the correct position.
        """
        expected_sets = {
                "c1" : "interleaved-charsets-c1.nex",
                "c2" : "interleaved-charsets-c2.nex",
                "c3" : "interleaved-charsets-c3.nex",
                }
        self.verify_subsets('interleaved-charsets-all.nex', expected_sets)

if __name__ == "__main__":
    unittest.main()
