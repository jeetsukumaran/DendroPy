#! /usr/bin/env python

############################################################################
##  fasta.py
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
FASTA format.

*** VERY CRUDE, AND CURRENTLY ONLY SUPPORTS WRITING  ***
"""
import textwrap
from dendropy import datasets

class FastaWriter(datasets.Writer):
    "Implements the DataWriter interface for handling PHYLIP files."
    
    def __init__(self):
        "Calls the base class constructor."
        datasets.Writer.__init__(self)
        
    def write_dataset(self, dataset, dest):
        "Writes dataset to a full PHYLIP document."
        for char_block in dataset.char_blocks:
            for taxon in char_block.matrix:
                dest.write(">%s\n" % str(taxon))
                seqs = char_block.matrix[taxon].values_as_string(sep='')
                seqsf = textwrap.fill(seqs, width=79, break_long_words=True)
                dest.write("%s\n\n" % seqsf)
