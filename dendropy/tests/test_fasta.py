#! /usr/bin/env python

############################################################################
##  test_tree_io.py
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
Tests input/output of trees from files.
"""

import unittest
import datetime
import logging
import tempfile
import os
from optparse import OptionGroup
from optparse import OptionParser
from cStringIO import StringIO

from dendropy import get_logger
from dendropy.datasets import Dataset
from dendropy import get_logging_level

import dendropy.tests
_LOG = get_logger("TreeParsingAndWriting")

from dendropy.tests import data_source_path, data_target_path
from dendropy.datasets import Dataset

### MODULES THAT WE ARE TESTING ###
from dendropy import nexus
# from dendropy import nexml
### MODULES THAT WE ARE TESTING ###
class TestFasta(unittest.TestCase):

    def testAsStrReading(self):
        fp = data_source_path("bad_names.fasta")
        fileobj = open(fp, 'rU')
        dataset = Dataset()
        dataset.read(fileobj, format='DNAFasta', row_type='str')
        fileobj.close()
        
        taxa = dataset.taxa_blocks[0]
        label = [i.label for i in taxa]
        expected = ['a Bad name', 'another', 'a Badn,ame', 'a  nothe++-_=+r', 'an!@#$o^&*()}{_ther']
        self.assertEquals(label, expected)
        tree_str = """('a Bad name','a  nothe++-_=+r',('an!@#"$o^&*()}{_ther',(another,'a Badn,ame')))"""
        tree_stream = StringIO(tree_str)
        dataset.read_trees(tree_stream, format="NEWICK")
        tree = dataset.trees_blocks[0][0]
        self.assertEquals(tree.compose_newick(), tree_str)
        
    def testAsStrReadingAndWriting(self):
        fp = data_source_path("bad_names.fasta")
        fileobj = open(fp, 'rU')
        dataset = Dataset()
        dataset.read(fileobj, format='DNAFasta', row_type='str')
        fileobj.close()
        op = tempfile.TemporaryFile()
        dataset.write(op, format="FASTA")        

if __name__ == "__main__":
    unittest.main()


    #compare_heavy(nexus.iterate_over_trees, "*.newick.tre")
    #compare_heavy(nexus.tree_iter, "*.newick.tre")
    
