#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests for rooting intepretation.
"""

import unittest
import warnings
from dendropy.dataio import newickreader

class RootingInterpreterLegacySupportTestCase(unittest.TestCase):
    legacy_kw = (
            ("as_unrooted", True,  "force-unrooted"),
            ("as_unrooted", False, "force-rooted"),
            ("as_rooted",   True,  "force-rooted"),
            ("as_rooted",   False, "force-unrooted"),
            ("default_as_unrooted", True,  "default-unrooted"),
            ("default_as_unrooted", False, "default-rooted"),
            ("default_as_rooted",   True,  "default-rooted"),
            ("default_as_rooted",   False, "default-unrooted"),
    )

    def test_multiple_keyword_error(self):
        with self.assertRaises(ValueError):
            nr = newickreader.NewickReader(
                    as_rooted=True,
                    as_unrooted=True)

    def test_legacy_keywords(self):
        for kwset in self.legacy_kw:
            with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
                warnings.simplefilter("always")
                # Trigger a warning.
                kwargs = {kwset[0]: kwset[1]}
                nr = newickreader.NewickReader(**kwargs)
                self.assertEqual(nr.rooting, kwset[2])

class RootingInterpreterForceUnrootedTestCase(unittest.TestCase):

    def test_with_no_token(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state(""), False)

    def test_with_no_token2(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state(), False)

    def test_with_rooted_token_lcase(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&r"), False)

    def test_with_rooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&R"), False)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&u"), False)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&U"), False)

    def test_with_meaningless_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("zz"), False)

class RootingInterpreterForceRootedTestCase(unittest.TestCase):

    def test_with_no_token(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state(""), True)

    def test_with_no_token2(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state(), True)

    def test_with_rooted_token_lcase(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&r"), True)

    def test_with_rooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&R"), True)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&u"), True)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&U"), True)

    def test_with_meaningless_token_ucase(self):
        nr = newickreader.NewickReader(rooting="force-rooted")
        self.assertIs(nr._parse_tree_rooting_state("zz"), True)

class RootingInterpreterDefaultUnrootedTestCase(unittest.TestCase):

    def test_with_no_token(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state(""), False)

    def test_with_no_token2(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state(), False)

    def test_with_rooted_token_lcase(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&r"), True)

    def test_with_rooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&R"), True)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&u"), False)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("&U"), False)

    def test_with_meaningless_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-unrooted")
        self.assertIs(nr._parse_tree_rooting_state("zz"), False)

class RootingInterpreterDefaultRootedTestCase(unittest.TestCase):

    def test_with_no_token(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state(""), True)

    def test_with_no_token2(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state(), True)

    def test_with_rooted_token_lcase(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&r"), True)

    def test_with_rooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&R"), True)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&u"), False)

    def test_with_unrooted_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state("&U"), False)

    def test_with_meaningless_token_ucase(self):
        nr = newickreader.NewickReader(rooting="default-rooted")
        self.assertIs(nr._parse_tree_rooting_state("zz"), True)

if __name__ == "__main__":
    unittest.main()
