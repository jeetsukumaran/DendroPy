#! /usr/bin/env python

############################################################################
##  test_dataset_rw.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
Tests dataset reading & writing
"""

import unittest
import csv

from dendropy import get_logger
import dendropy.tests

_LOG = get_logger("DatsetReadWrite")

from dendropy import datasets


class DatasetReadTest(unittest.TestCase):

    def compare_chars(self, src, format, expected):
        """Reads 'src', checks against 'expected'"""
        _LOG.info("Reading %s" % src.name)
        d = datasets.Dataset()
        d.read(src, format)
        
        taxa_block = d.taxa_blocks[0]
        char_block = d.char_blocks[0]
        
        assert len(expected) == len(char_block)
        assert len(expected) == len(taxa_block)

        for tax_idx, (exp_taxa, exp_seq) in enumerate(expected):
        
            taxon = taxa_block[tax_idx]
            label = taxon.label
            
            # ok, this is ugly, but my nexus parser does not
            # do the "_" => " " conversion in taxlabels (yet)
            # so ...
            assert ((exp_taxa == label) \
                or (exp_taxa.replace("_", " ") == label) \
                or (exp_taxa.replace(" ", "_") == label)), \
                "(Taxon #%d) %s not eq. %s" % (tax_idx, exp_taxa, label)
                
            assert len(exp_seq) == len(char_block.matrix[taxon])   
                
            for col_idx, symbol1 in enumerate(exp_seq):
                test_state = char_block.matrix[taxon][col_idx].value
                if char_block.matrix[taxon][col_idx].column_type is not None:
                    state_alpha = char_block.matrix[taxon][col_idx].column_type.state_alphabet
                else:
                    state_alpha = char_block.default_state_alphabet
                exp_state = state_alpha.state_for_symbol(symbol1)
                assert test_state == exp_state
                
    def testBasicCharactersRead(self):
        test_cases = [
            ['anolis.chars.nexus', 'nexus', 'anolis.chars.csv'],
#             ['anolis.chars.nexml', 'nexml', 'anolis.chars.csv'],
            ['primates.chars.nexus', 'nexus', 'primates.chars.csv'],
            ['primates.chars.interleaved.nexus', 'nexus', 'primates.chars.csv'],
            ['primates.chars.simple.nexus', 'nexus', 'primates.chars.csv'],
            ['primates.chars.simple.interleaved.nexus', 'nexus', 'primates.chars.csv'],
#             ['primates.chars.nexml', 'nexml', 'primates.chars.csv'],
        ]
        
        for t in test_cases:
            src = open(dendropy.tests.data_source_path(t[0]), "rU")
            expected = [ (s[0], s[1]) for s in csv.reader(open(dendropy.tests.data_source_path(t[2]), "rU"))]
            self.compare_chars(src, t[1], expected)        
        
if __name__ == "__main__":
    unittest.main()
        