#! /usr/bin/env python

############################################################################
##  test_misc.py
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
This module collects "quick-and-dirty" (general single-function-scope)
tests. As bugs are discovered during general usage of DendroPy, a simple
function that replicates the error should be immediately be added with a
unique name prefixed by the string "test" (e.g. "testFileWriting") to
the class `TestMiscellaneous`. Docstrings providing bug reporter, test author,
description etc. would be greatly appreciated. 

Full local path to test data for reading is available via:

    dendropy.tests.data_source_path(<filename>)
    
While paths to files for testing of output/writing is available via:

    dendropy.tests.data_target_path(<filename>)
    
Or, alternatively using tempfile.    

"""

import unittest
import tempfile

from dendropy import get_logger
_LOG = get_logger("GeneralTests")

import dendropy.tests
from dendropy.datasets import Dataset

### MODULES THAT WE ARE TESTING ###
from dendropy import nexus
### MODULES THAT WE ARE TESTING ###

class TestMiscellaneous(unittest.TestCase):
    
    def test_miscellaneous(self):
        """
        Sample miscellaneous test.
        Reading/writing of FASTA files read with `row_type`=`str`. Reproduces 
        bug condition discovered by Jiaye 2009-10-22. 
        """
        _LOG.info("Reading/writing of FASTA files read with `row_type`=`str`")
        fp = dendropy.tests.data_source_path("bad_names.fasta")
        fileobj = open(fp, 'rU')
        dataset = Dataset()
        dataset.read(fileobj, format='DNAFasta', row_type='str')
        fileobj.close()
        op = tempfile.TemporaryFile()
        dataset.write(op, format="FASTA")          

if __name__ == "__main__":
    unittest.main()