#! /usr/bin/env python

############################################################################
##  test_char_io.py
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
Tests input/output of characters from files.
"""

import unittest
import datetime
import logging
import tempfile
import os
from optparse import OptionGroup
from optparse import OptionParser
import StringIO

from dendropy import get_logger
from dendropy import get_logging_level

import dendropy.tests
_LOG = get_logger("Character Parsing and Writing")


### MODULES THAT WE ARE TESTING ###
from dendropy import nexus
from dendropy import nexml
### MODULES THAT WE ARE TESTING ###

class CharIOTest(unittest.TestCase):
    
    def setUp(self):
        self.formats = ["nexus",
                        "nexml",
                       ]
        self.readers = {"nexus": nexus.NexusReader,
                        "nexml": nexml.NexmlReader,
                       }        
        self.writers = {"nexus": nexus.NexusWriter,
                        "nexml": nexml.NexmlWriter,
                       }
        self.char_data = ["anolis",
                          "primates",
                          #"terrarana",
                         ]
        self.char_files = {}
        self.char_check_files = {}
        for format in self.readers:
            self.char_files[format] = []
            for cd in self.char_data:
                fpath = dendropy.tests.data_source_path(cd + ".chars." + format)
                if os.path.exists(fpath):
                    self.char_files[format].append(fpath)
                    self.char_check_files[fpath] = dendropy.tests.data_source_path(cd + ".chars." + "csv")
                
    def readCheckFile(self, fpath):
        data = {}
        flines = open(fpath, 'r')
        for line in flines:
            line = line.replace('\n','').strip()
            if line:
                taxon, seq = line.split(",")
                data[taxon] = seq
        return data            
                
    def testCharParse(self):
        for format in self.formats:
            _LOG.info('\n[Testing %s format parsing: <%s>, <%s>]'  % (format.upper(),
                                                                          self.readers[format].__name__,
                                                                          self.writers[format].__name__))
            for fpath in self.char_files[format]:
                _LOG.info("\nDATA FILE: \"%s\"" % os.path.basename(fpath))
                reader = self.readers[format]()
                writer = self.writers[format]()
                dataset = reader.read_dataset(open(fpath, 'r'))
                _LOG.info("%d taxa block(s), %d character block(s), %d trees block(s)" % (len(dataset.char_blocks),
                                                                                          len(dataset.taxa_blocks),
                                                                                          len(dataset.trees_blocks)))
                                                                                          
                _LOG.info("\nTaxa Blocks:")                                                                                          
                for tb_idx, tb in enumerate(dataset.taxa_blocks):
                    _LOG.info("-- Taxa Block %d (%s \"%s\"): %d taxa" % (tb_idx+1, tb.oid, str(tb.label), len(tb)))
                    _LOG.debug("\n".join([("      " + str(t)) for t in tb]))
                    _LOG.debug("")
                    
                _LOG.info("\nCharacter Blocks:")                     
                for cb_idx, cb in enumerate(dataset.char_blocks):
                    ntax = len(cb.matrix)
                    if cb.matrix:
                        nchar = max([len(v) for v in cb.matrix.values()])
                    _LOG.info("-- Character Block %d (%s \"%s\"): %s with %d taxa and %d characters" % (cb_idx+1, 
                                                                                            cb.oid, 
                                                                                            str(cb.label), 
                                                                                            cb.__class__.__name__,
                                                                                            ntax, 
                                                                                            nchar))
                    check_data = self.readCheckFile(self.char_check_files[fpath])
                    
                    # check that all taxa in original is in character matrix
                    cb_taxa = [str(t) for t in cb.matrix]
                    for t in check_data:
                        assert t in cb_taxa, \
                            'Taxon "%s" expected but not found in Character Block %d (%s)' % (t, cb_idx+1, cb.oid)
                    
                    # and vice versa
                    for t in cb_taxa:
                        assert t in check_data, \
                            'Taxon "%s" found in Character Block %d (%s) is unexpected' % (t, cb_idx+1, cb.oid)

                    # check for sequence correspondence
                    _LOG.debug("")                    
                    for taxon in cb.matrix:
                        tax_label = str(taxon)
                        seq = ''.join(cb.matrix[taxon].values_as_string_list())                        
                        _LOG.debug("%s          %s" % (tax_label, seq))
                        _LOG.debug("")
                        check_seq = check_data[tax_label]
                        assert seq == check_seq, \
                                    '\n\nSequences do not match for taxon "%s":\n\nEXPECTED >>>\n%s\n<<< EXPECTED\nFOUND >>>\n%s<<< FOUND' % (tax_label, check_seq, seq)


if __name__ == "__main__":
    unittest.main()
