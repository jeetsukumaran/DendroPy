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
Tests composition and taxon indexing of CharacterArray.
"""

import unittest
import dendropy
from dendropy.test.support.extendedtest import ExtendedTestCase
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class TestCharStruct(ExtendedTestCase):

    def setUp(self):
        self.tb1 = dendropy.TaxonSet(label="TI1")
        for i in range(1,11):
            self.tb1.new_taxon(label="T%02d" % i)
        self.cb1 = dendropy.DnaCharacterArray(taxon_set=self.tb1, label="TI1, CA1")
        for t in self.tb1:
            self.cb1.append_taxon_sequence(t, state_symbols="AAAAAAAAAA")
        self.tb2 = dendropy.TaxonSet(label="TI2")
        for i in range(1,21):
            self.tb2.new_taxon(label="T%02d" % i)
        self.cb2 = dendropy.DnaCharacterArray(taxon_set=self.tb2, label="TI2, CA2")
        for t in self.tb2:
            self.cb2.append_taxon_sequence(t, state_symbols="CCCCCCCCCC")

class TestExtendCharacters(TestCharStruct):

    def runTest(self):
        ntax_pre = len(self.cb1)
        nchars_pre = len(self.cb1.values()[0])
        self.cb1.extend_characters(self.cb2)
        self.assertEqual(len(self.cb1), ntax_pre)
        for t in self.cb1:
            self.assertEqual(len(self.cb1[t]), 20)
            self.assertEqual(self.cb1[t].values_as_string(), "AAAAAAAAAACCCCCCCCCC")

class TestExtendSequencesOverwrite(TestCharStruct):

    def runTest(self):
        self.cb1.extend(self.cb2, overwrite_existing=True)
        target_ntax = 20
        self.assertEqual(len(self.cb1), target_ntax)
        self.assertEqual(len(self.cb1.taxon_set), target_ntax)
        for t in self.tb2:
            cb_tb_labels = self.cb1.taxon_set.labels()
            self.assertIsContainedIn(t.label, cb_tb_labels)
            cb_labels = [t.label for t in self.cb1]
            self.assertIsContainedIn(t.label, cb_labels)
        for t in self.cb1:
            self.assertEqual(len(self.cb1[t]), 10)
            self.assertEqual(self.cb1[t].values_as_string(), "CCCCCCCCCC",)

class TestExtendSequencesAppend(TestCharStruct):

    def runTest(self):
        self.cb1.extend(self.cb2, append_existing=True)
        target_ntax = 20
        self.assertEqual(len(self.cb1), target_ntax)
        self.assertEqual(len(self.cb1.taxon_set), target_ntax)
        for t in self.tb2:
            cb_tb_labels = self.cb1.taxon_set.labels()
            self.assertIsContainedIn(t.label, cb_tb_labels)
            cb_labels = [t.label for t in self.cb1]
            self.assertIsContainedIn(t.label, cb_labels)
        for t in self.cb1:
            if int(t.label[-2:]) > 10:
                self.assertEqual(len(self.cb1[t]), 10)
                self.assertEqual(self.cb1[t].values_as_string(), "CCCCCCCCCC")
            else:
                self.assertEqual(len(self.cb1[t]), 20)
                self.assertEqual(self.cb1[t].values_as_string(), "AAAAAAAAAACCCCCCCCCC")

if __name__ == "__main__":
    unittest.main()

