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
Tests Tree taxon management
"""

import unittest
import collections
import dendropy
from dendropy.test.support import dendropytest

class StateAlphabetTester(object):

    def validate_state_identities(self,
            state_container,
            state_iter,
            expected_symbols,
            case_sensitive=True):
        self.assertEqual(len(state_container), len(expected_symbols))
        states = list(state_iter())
        self.assertEqual(len(states), len(expected_symbols))
        for state, symbol in zip(states, expected_symbols):
            self.assertEqual(state.symbol, symbol)
            if case_sensitive:
                self.assertNotIn(symbol.upper(), state.symbol_synonyms)
                self.assertNotIn(symbol.lower(), state.symbol_synonyms)
            else:
                if symbol.upper() != symbol:
                    self.assertIn(symbol.upper(), state.symbol_synonyms)
                if symbol.lower() != symbol:
                    self.assertIn(symbol.lower(), state.symbol_synonyms)

class DnaStateAlphabetTest(
        StateAlphabetTester,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.expected_fundamental_state_symbols = ["A", "C", "G", "T", "-"]
        self.ambiguous_symbol_mappings = collections.OrderedDict({
                "?": "ACGT-",
                "N": "ACGT" ,
                "X": "ACGT" ,
                "R": "AG"   ,
                "Y": "CT"   ,
                "M": "AC"   ,
                "W": "AT"   ,
                "S": "CG"   ,
                "K": "GT"   ,
                "V": "ACG"  ,
                "H": "ACT"  ,
                "D": "AGT"  ,
                "B": "CGT"  ,
                })
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())
        self.sa = dendropy.DNA_STATE_ALPHABET

    def test_fundamental_state_definitions(self):
        self.validate_state_identities(
                self.sa._fundamental_states,
                self.sa.fundamental_state_iter,
                self.expected_fundamental_state_symbols,
                case_sensitive=False)

if __name__ == "__main__":
    unittest.main()
