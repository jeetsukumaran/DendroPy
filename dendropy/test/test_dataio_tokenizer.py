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
Tests for tokenizers classes.
"""

import unittest
from dendropy.dataio import nexusprocessing
from dendropy.utility.textprocessing import StringIO

class NexusTokenizerTestCase(unittest.TestCase):
    """
    Unit tests for NexusTokenizer.
    """

    def check_tokenization(self,
            input_str,
            expected_tokens):
        src = StringIO(input_str)
        observed = []
        for token in nexusprocessing.NexusTokenizer(src=src):
            observed.append(token)
        self.assertEqual(observed, expected_tokens)

    def test_simple_string(self):
        input_str = "the    quick    brown\t\tfox \n  jumps over\t\t\n the    lazy dog"
        expected = [
                "the", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog"
                ]
        self.check_tokenization(input_str, expected)

    def test_simple_quoted_string(self):
        input_str = "the quick 'brown fox' jumps over the 'lazy dog'"
        expected = [
                "the", "quick", "brown fox", "jumps", "over", "the", "lazy dog"
                ]
        self.check_tokenization(input_str, expected)

    def test_padded_quoted_string(self):
        input_str = "the quick 'brown fox''s friend' jumps over the 'lazy dog''s colleague'"
        expected = [
                "the", "quick", "brown fox's friend", "jumps", "over", "the", "lazy dog's colleague"
                ]
        self.check_tokenization(input_str, expected)

    def test_runon_quoted_string(self):
        input_str = "'a','b','c','d','e'"
        expected = [
                "a", ",", "b", ",", "c", ",", "d", ",", "e",
                ]
        self.check_tokenization(input_str, expected)

    def test_comments(self):
        input_str = "[&R] (foo:1 [a foo object], [start of subgroup](bar:2, c:2)[end of group][][][";
        expected = [
                "(", "foo", ":","1", ",", "(", "bar", ":", "2",  ",", "c", ":", "2", ")"
                ]
        self.check_tokenization(input_str, expected)

    def test_empty(self):
        input_str = "";
        expected = []
        self.check_tokenization(input_str, expected)

    def test_captured_delimiters(self):
        input_str = "(aaa:1.00,     (b:2.18e-1,      (ccc:11, d:1e-1)   k:  3)  u:   7)    rrr:0.0;";
        expected = [
            "(",
            "aaa",
            ":",
            "1.00",
            ",",
            "(",
            "b",
            ":",
            "2.18e-1",
            ",",
            "(",
            "ccc",
            ":",
            "11",
            ",",
            "d",
            ":",
            "1e-1",
            ")",
            "k",
            ":",
            "3",
            ")",
            "u",
            ":",
            "7",
            ")",
            "rrr",
            ":",
            "0.0",
            ";"
                ]
        self.check_tokenization(input_str, expected)

    def test_comments(self):
        input_str = "([the quick]apple[brown],([fox]banjo,([jumps]cucumber[over the],[really]dogwood)[lazy]eggplant)) rhubarb[dog];";
        expected_comments = {
            "apple"    : ["the quick", "brown"],
            "banjo"    : ["fox"               ],
            "cucumber" : ["jumps", "over the" ],
            "dogwood"  : ["really"            ],
            "eggplant" : ["lazy"              ],
            "rhubarb"  : ["dog"               ],
                }
        expected_tokens = [
                "(",
                "apple",
                ",",
                "(",
                "banjo",
                ",",
                "(",
                "cucumber",
                ",",
                "dogwood",
                ")",
                "eggplant",
                ")",
                ")",
                "rhubarb",
                ";"
                ]
        src = StringIO(input_str)
        observed_tokens = []
        tk = nexusprocessing.NexusTokenizer(src=src)
        for token in tk:
            if token in expected_comments:
                expected_comment = expected_comments[token]
                observed_comment = tk.pull_captured_comments()
                self.assertEqual(expected_comment, observed_comment)
                del expected_comments[token]
            observed_tokens.append(token)
        self.assertEqual(expected_comments, {})
        self.assertEqual(observed_tokens, expected_tokens)

if __name__ == "__main__":
    unittest.main()
