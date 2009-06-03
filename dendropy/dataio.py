#! /usr/bin/env python

############################################################################
##  dataio.py
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
Convenience packaging around readers/writers.
"""

import sys
from cStringIO import StringIO

from dendropy.datasets import Dataset
from dendropy.trees import TreesBlock
from dendropy import nexus
from dendropy.nexus import RootingInterpretation
from dendropy import nexml
from dendropy import fasta
from dendropy import phylip
from dendropy import get_logger, deprecation
_LOG = get_logger('dataio')

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
   
def dataset_from_file(file_obj, format, **kwargs):
    """
    Returns a Dataset object parsed from the source, where:
        `file_obj`   - can either be a file descriptor object/handle opened 
                   for reading or a string indicating a filepath that 
                   can be opened for reading using open().     
        `format` - file format specification               
    """
    deprecation("'dataio.dataset_from_file()' is deprecated; use 'read()' method of a Dataset object instead", logger_obj=_LOG)
    reader = get_reader(format)
    return reader.read_dataset(source_file_handle(file_obj=file_obj), **kwargs)
    
def dataset_from_string(string, format, **kwargs):
    """
    Returns a Dataset object parsed from the source, where:
        `string` - a string containing the data to be parsed.   
        `format` - file format specification               
    """
    deprecation("'dataio.dataset_from_string()' is deprecated; use 'from_string()' method of a Dataset object instead", logger_obj=_LOG)
    reader = get_reader(format)
    return reader.read_dataset(source_file_handle(string=string), **kwargs)    
    
def trees_from_file(file_obj, format, encode_splits=False, rooted=RootingInterpretation.UNKNOWN_DEF_UNROOTED, **kwargs):
    """
    Returns a *list* of TreesBlock objects parsed from the source, where:
        `file_obj`   - can either be a file descriptor object/handle opened 
                   for reading or a string indicating a filepath that 
                   can be opened for reading using open().    
        `format` - file format specification               
    """
    deprecation("'dataio.trees_from_file()' is deprecated: use 'read_trees()' method of a Dataset object instead", logger_obj=_LOG)
    reader = get_reader(format)
    return reader.read_trees(source_file_handle(file_obj=file_obj), encode_splits=encode_splits, rooted=rooted, **kwargs)
    
def trees_from_string(string, format, encode_splits=False, rooted=RootingInterpretation.UNKNOWN_DEF_UNROOTED, **kwargs):
    """
    Returns a *list* of TreesBlock objects parsed from the source, where:
        `string` - a string containing the data to be parsed.   
        `format` - file format specification               
    """
    deprecation("'dataio.trees_from_string()' is deprecated: use 'trees_from_string()' method of a Dataset object instead", logger_obj=_LOG)
    reader = get_reader(format)
    return reader.read_trees(source_file_handle(string=string), encode_splits=encode_splits, rooted=rooted, **kwargs)
    
def from_nexus(file_obj=None, string=None, **kwargs):
    """
    Returns a Dataset object parsed from a NEXUS or NEWICK source.
        `file_obj`   - can either be a file descriptor object/handle opened 
                   for reading or a string indicating a filepath that 
                   can be opened for reading using open().
        `string` - a string containing the data to be parsed.
    Either `file_obj` or `string` must be given. If both are given, `file_obj` is used.                
    """
    deprecation("'dataio.from_nexus() is deprecated: use 'read()' method of a Dataset object instead", logger_obj=_LOG)
    return nexus.read_dataset(source_file_handle(file_obj=file_obj, string=string), **kwargs)    
    
############################################################################
## Wrappers (Writing)    

def store_dataset(dataset, format, dest=None):
    """
    Writes the Dataset object `dataset` using `writer` (a DatasetWriter or 
    derived object) to `dest`. If `dest` is a string, then it is assumed to be
    a path name, and open() is used to construct an output stream handle from it.
    If `dest` is not given, then the dataset is written to a string and a string 
    is returned.
    """
    deprecation("'dataio.store_dataset()' is deprecated: use 'write()' method of a Dataset object instead", logger_obj=_LOG)
    writer = get_writer(format)
    if dest is None:
        dest = StringIO()
    if isinstance(dest, str):
        dest = open(dest, "w")
    writer.write_dataset(dataset, dest)
    if hasattr(dest, "getvalue"):
        return dest.getvalue()

def store_trees(trees_collection, format, dest=None):
    "Writes the list of trees `trees` to `dest` using writer."
    deprecation("'dataio.store_trees()' is deprecated: use 'write()' method of a Dataset object instead", logger_obj=_LOG)
    if isinstance(trees_collection, TreesBlock):
        trees_block = trees_collection
    else:
        trees_block = TreesBlock()
        for tree in trees_collection:
            trees_block.append(tree)
        trees_block.normalize_taxa()
    dataset = Dataset()
    dataset.add_trees_block(trees_block=trees_block)
    return store_dataset(dataset=dataset,
        format=format,                  
        dest=dest)
                  
def store_chars(char_block, format, dest=None):
    "Writes the CharacterBlock `char_block` to `dest` using writer."
    deprecation("'dataio.store_chars()' is deprecated: use 'write()' method of a Dataset object instead", logger_obj=_LOG)
    dataset = Dataset()
    dataset.add_char_block(char_block=char_block)
    return store_dataset(dataset=dataset,
        format=format,                  
        dest=dest)                  

############################################################################
## Helpers

def source_file_handle(file_obj=None, string=None):
    """
    Construct an appropriate file_obj handle (i.e. something that supports read()
    operations) based on the given arguments.
    """
    if file_obj is None and string is None:
        raise Exception("File or string source must be specified.")            
    if file_obj is not None:        
        if isinstance(file_obj, str):
            file_obj = open(file_obj, "r")        
        return file_obj
    else:
        return StringIO(string)
        
def get_writer(format):
    "Return reader of the appropriate format."
    format = format.upper()
    if format not in WRITERS:
        raise Exception('Unrecognized format specificiation "%s", ' \
            'must be one of: %s' % (format,
             ", ".join([('"'+f+'"') for f in WRITERS]),
             ))
    return WRITERS[format]()      
    
def get_reader(format):
    "Return reader of the appropriate format."
    format = format.upper()
    if format not in READERS:
        raise Exception('Unrecognized format specificiation "%s", ' \
            'must be one of: %s' % (format,
             ", ".join([('"'+f+'"') for f in READERS]),
             ))
    return READERS[format]()            
    

def trees_from_newick(nl, taxa_block=None, **kwargs):
    """Takes an iterable list of newick strings (or files with just newick strings
    in them.
    """
    reader = get_reader(NEWICK)
    if taxa_block is not None:
        dataset = Dataset(taxa_blocks=[taxa_block])
    else:
        dataset = Dataset()
    for t in nl:
        f = t
        if isinstance(t, str):
            f = StringIO(t)
        reader.read_dataset(file_obj=f, dataset=dataset, **kwargs)
    return dataset



class MultiFileTreeIterator(object):
    def __init__(self, sources=[], 
                       core_iterator=None, 
                       taxa_block=None, 
                       dataset=None, 
                       format=None, 
                       from_index=0, 
                       progress_func=None,
                       **kwargs):
        """An iterable collection of trees from multiple sources
            `sources` is as list of tree sources each can be either a file path (a str) 
                or a file-like object
            Either
                `core_iterator` or `dataset` must be specified as the source of the iterator
                If the `dataset` is used, the the `format` must be specified
             
            `from_index` can be used to skip the first `from_index` trees from _EACH_ file
                this is useful if you want to discard a certain number of trees from the beginning of each run
                as burnin.
               
            If `progress_func` is specified verbose messages will be sent to it for every tree processed.
        """
        if dataset is None:
            self.dataset = Dataset()
        else:
            self.dataset = dataset
        self.taxa_block = taxa_block
        self.format = format        
        if core_iterator is None:
            if self.format is None:
                raise ValueError("Either 'core_iterator' or 'format' flags must be used")
            self.using_data_it = True
        else:
            self.using_data_it = False
            self._core_iterator = core_iterator
        self.progress_func = progress_func
        self.sources = sources
        self.total_trees_read = 0
        self.total_trees_ignored = 0
        self.total_num_sources_read = 0
        self.from_index = from_index
        self.iterator_kwargs = kwargs

    def __iter__(self):
        si = self.from_index
        tb = self.taxa_block
        progress_func = self.progress_func
        self.curr_trees_read = 0
        self.curr_trees_ignored = 0
        self.curr_num_sources_read = 0
        for source_ind, tree_source in enumerate(self.sources):
            if isinstance(tree_source, str):
                fo = open(tree_source, "rU")
            else:
                fo = tree_source
            if progress_func:
                current_file_note = "Tree file %d of %d: " % (source_ind + 1, len(self.sources))
            self.curr_num_sources_read += 1
            self.total_num_sources_read += 1
            for n, tree in enumerate(self._raw_iter(fo, tb)):
                if (not si) or (n >= si):
                    if tb is None:
                        tb = tree.taxa_block
                    self.total_trees_read += 1
                    self.curr_trees_read += 1
                    if progress_func:
                        progress_func("%sProcessing tree %d" % (current_file_note, (n+1)))
                    yield tree
                else:
                    self.total_trees_ignored += 1
                    self.curr_trees_ignored += 1
                    if progress_func:
                        progress_func("%sSkipping tree %d (# to skip=%d)" % (current_file_note, (n+1), si))

    def _raw_iter(self, fo, tb):
        if self.using_data_it:
            for tree in self.dataset.iterate_over_trees(fo, taxa_block=tb, format=self.format, **(self.iterator_kwargs)):
                yield tree
        else:
            for tree in self._core_iterator(fo, taxa_block=tb, file_format=self.format, **(self.iterator_kwargs)):
                yield tree
        
        



