#! /usr/bin/env python

############################################################################
##  phylip.py
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
This module wraps routines needed for reading and writing data in
PHYLIP format.

*** VERY CRUDE, AND CURRENTLY ONLY SUPPORTS WRITING IN RELAXED UNINTERLEAVED FORMAT ***
"""

from dendropy import datasets

class PhylipWriter(datasets.Writer):
    "Implements the DataWriter interface for handling PHYLIP files."
    
    def __init__(self):
        "Calls the base class constructor."
        datasets.Writer.__init__(self)
        
    def write_dataset(self, dataset, dest):
        "Writes dataset to a full PHYLIP document."
        char_block = dataset.char_blocks[0]
        n_seqs = len(char_block.matrix)
        n_sites = len(char_block.matrix.values()[0])
        dest.write("%d %d\n" % (n_seqs, n_sites))
        maxlen = max([len(str(taxon)) for taxon in char_block.matrix])
        for taxon in char_block.matrix:
            dest.write("%s        %s\n" % ( str(taxon).ljust(maxlen), str(char_block.matrix[taxon]).replace(' ', '')))
                    
        
