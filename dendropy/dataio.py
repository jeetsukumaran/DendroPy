#! /usr/bin/env python

############################################################################
##  dataio.py
##
##  Part of the DendroPy phylogenetic computation library.
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
Convenience packaging around readers/writers.
"""

import sys
import StringIO

from dendropy.datasets import Dataset
from dendropy.trees import TreesBlock
from dendropy import nexus
from dendropy import nexml
from dendropy import fasta
from dendropy import phylip

############################################################################
## File Formats

NEXUS='NEXUS'
NEWICK='NEWICK'
NEXML='NEXML'
FASTA='FASTA'
PHYLIP='PHYLIP'
FORMATS = [NEXUS, NEXML, NEWICK, FASTA, PHYLIP]

READERS = {
    NEXUS: nexus.NexusReader,
    NEWICK: nexus.NewickReader,
    NEXML: nexml.NexmlReader,
}

WRITERS = {
    NEXUS: nexus.NexusWriter,
    NEWICK: nexus.NewickWriter,
    NEXML: nexml.NexmlWriter,
    FASTA: fasta.FastaWriter,
    PHYLIP: phylip.PhylipWriter,
}

############################################################################
## Wrappers (Reading/Parsing)

def get_dataset(format, file=None, string=None):
    """
    Returns a Dataset object parsed from the source, where:
        `format` - file format specification
        `file`   - can either be a file descriptor object/handle opened 
                   for reading or a string indicating a filepath that 
                   can be opened for reading using open().
        `string` - a string containing the data to be parsed.
    Either `file` or `string` must be given. If both are given, `file` is used.                
    """
    reader = get_reader(format)
    return reader.read_dataset(source_file_handle(file=file, string=string))
    
def get_trees(format, file=None, string=None):
    """
    Returns a *list* of TreesBlock objects parsed from the source, where:
        `format` - file format specification
        `file`   - can either be a file descriptor object/handle opened 
                   for reading or a string indicating a filepath that 
                   can be opened for reading using open().
        `string` - a string containing the data to be parsed.
    Either `file` or `string` must be given. If both are given, `file` is used.                
    """
    reader = get_reader(format)
    return reader.read_trees(source_file_handle(file=file, string=string))
    
def get_nexus(file=None, string=None):
    """
    Returns a Dataset object parsed from a NEXUS or NEWICK source.
        `file`   - can either be a file descriptor object/handle opened 
                   for reading or a string indicating a filepath that 
                   can be opened for reading using open().
        `string` - a string containing the data to be parsed.
    Either `file` or `string` must be given. If both are given, `file` is used.                
    """
    return nexus.read_dataset(source_file_handle(file=file, string=string))    
    
############################################################################
## Wrappers (Writing)    

def store_dataset(format, dataset, dest=None):
    """
    Writes the Dataset object `dataset` using `writer` (a DatasetWriter or 
    derived object) to `dest`. If `dest` is a string, then it is assumed to be
    a path name, and open() is used to construct an output stream handle from it.
    If `dest` is not given, then the dataset is written to a string and a string 
    is returned.
    """
    if dest is None:
        dest = StringIO.StringIO()
    if isinstance(dest, str):
        dest = open(dest, "w")
    writer.write_dataset(dataset, dest)
    if hasattr(dest, "getvalue"):
        return dest.getvalue()

def store_trees(format, trees, dest=None):
    """
    Writes the list of trees `trees` to `dest` using writer.
    """
    if isinstance(trees, trees_block):
        trees_block = trees
    else:
        trees_block = TreesBlock()
        for tree in trees:
            trees_block.append(trees)
        trees_block.normalize_taxa()
    dataset = Dataset()
    dataset.add_trees_block(trees_block=trees_block)
    write_dataset(dataset=dataset,
                  writer=writer,
                  dest=dest)

############################################################################
## Helpers

def source_file_handle(file=None, string=None):
    """
    Construct an appropriate file handle (i.e. something that supports read()
    operations) based on the given arguments.
    """
    if file is None and string is None:
        raise Exception("File or string source must be specified.")            
    if isinstance(file, str):
        file = open(file, "r")
    return file    
    
def get_reader(format):
    """
    Return reader of the appropriate format.
    """
    format = format.upper()
    if format not in READERS:
        raise Exception('Unrecognized format specificiation "%s", must be one of: %s' % (format,
                                                                                         ", ".join([('"'+f+'"') for f in FORMATS]),
                                                                                         ))
    return READERS[format]()            
    