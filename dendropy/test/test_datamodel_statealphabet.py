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
            expected_denomination,
            member_state_map,
            case_sensitive=True):
        self.assertEqual(len(state_container), len(expected_symbols))
        states = list(state_iter())
        self.assertEqual(len(states), len(expected_symbols))
        for state, symbol in zip(states, expected_symbols):
            self.assertEqual(state.symbol, symbol)
            self.assertEqual(state.state_denomination, expected_denomination)
            if member_state_map:
                expected_member_state_symbols = frozenset(member_state_map[symbol])
                self.assertEqual(state.fundamental_symbols, expected_member_state_symbols)
            else:
                self.assertEqual(state.fundamental_states, set([state]))
                self.assertEqual(state.fundamental_symbols, set([state.symbol]))
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
        self.ambiguous_symbol_mappings = collections.OrderedDict()
        self.ambiguous_symbol_mappings["?"] = "ACGT-"
        self.ambiguous_symbol_mappings["N"] = "ACGT"
        self.ambiguous_symbol_mappings["R"] = "AG"
        self.ambiguous_symbol_mappings["Y"] = "CT"
        self.ambiguous_symbol_mappings["M"] = "AC"
        self.ambiguous_symbol_mappings["W"] = "AT"
        self.ambiguous_symbol_mappings["S"] = "CG"
        self.ambiguous_symbol_mappings["K"] = "GT"
        self.ambiguous_symbol_mappings["V"] = "ACG"
        self.ambiguous_symbol_mappings["H"] = "ACT"
        self.ambiguous_symbol_mappings["D"] = "AGT"
        self.ambiguous_symbol_mappings["B"] = "CGT"
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())
        self.polymorphic_symbol_mappings = collections.OrderedDict()
        self.expected_polymorphic_state_symbols = list(self.polymorphic_symbol_mappings.keys())
        self.sa = dendropy.DNA_STATE_ALPHABET

    def test_fundamental_state_definitions(self):
        self.validate_state_identities(
                state_container=self.sa._fundamental_states,
                state_iter=self.sa.fundamental_state_iter,
                expected_symbols=self.expected_fundamental_state_symbols,
                expected_denomination=self.sa.FUNDAMENTAL_STATE,
                member_state_map=None,
                case_sensitive=False)

    def test_ambiguous_state_definitions(self):
        self.validate_state_identities(
                state_container=self.sa._ambiguous_states,
                state_iter=self.sa.ambiguous_state_iter,
                expected_symbols=self.expected_ambiguous_state_symbols,
                expected_denomination=self.sa.AMBIGUOUS_STATE,
                member_state_map=self.ambiguous_symbol_mappings,
                case_sensitive=False)

    def test_polymorphic_state_definitions(self):
        self.validate_state_identities(
                state_container=self.sa._polymorphic_states,
                state_iter=self.sa.polymorphic_state_iter,
                expected_symbols=self.expected_polymorphic_state_symbols,
                expected_denomination=self.sa.POLYMORPHIC_STATE,
                member_state_map=self.polymorphic_symbol_mappings,
                case_sensitive=False)

if __name__ == "__main__":
    unittest.main()
