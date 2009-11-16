#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Infrastructure for phylogenetic data object serialization and deserialization.
Provides support for reading/parsing and formatting/writing phylogenetic data
in various formats.
"""

import sys
import os

from dendropy.utility import error
from dendropy.utility import iosys
from dendropy.utility.containers import OrderedCaselessDict

###############################################################################
## _DataFormat

class _DataFormat(object):

    def __init__(self,
                 name,
                 reader_type=None,
                 writer_type=None,
                 tree_source_iter=None):
        self.name = name
        self.reader_type = reader_type
        self.writer_type = writer_type
        self.tree_source_iter = tree_source_iter

    def has_reader(self):
        return self.reader_type is not None

    def has_writer(self):
        return self.writer_type is not None

    def has_tree_source_iter(self):
        return self.tree_source_iter is not None

    def get_reader(self, **kwargs):
        if self.reader_type is None:
            raise error.UnsupportedFormatError("Reading is not currently supported for data format '%s'" % self.name)
        return self.reader_type(**kwargs)

    def get_writer(self, **kwargs):
        if self.writer_type is None:
            raise error.UnsupportedFormatError("Writing is not currently supported for data format '%s'" % self.name)
        return self.writer_type(**kwargs)

    def get_tree_source_iter(self, stream, **kwargs):
        if self.tree_source_iter is None:
            raise error.UnsupportedFormatError("Iteration over source trees not currently supported for data format '%s'" % self.name)
        return self.tree_source_iter(stream, **kwargs)

###############################################################################
## _DataFormatRegistry

class _DataFormatRegistry(object):

    def __init__(self):
        self.formats = OrderedCaselessDict()

    def add_format(self, data_format):
        self.formats[data_format.name] = data_format

    def remove_format(self, data_format):
        del(self.formats[data_format.name])

    def add(self, name, reader_type=None, writer_type=None, tree_source_iter=None):
        self.formats[name] = _DataFormat(name,
                reader_type=reader_type,
                writer_type=writer_type,
                tree_source_iter=tree_source_iter)

    def remove(self, name):
        del(self.formats[name])

    def get_reader(self, name, **kwargs):
        ## hack to avoid confusing users ##
        if name.lower() == "fasta":
            raise error.UnsupportedFormatError("FASTA format needs to be specified as 'dnafasta', 'rnafasta', or 'proteinfasta'")
        if name not in self.formats:
            raise error.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].get_reader(**kwargs)

    def get_writer(self, name, **kwargs):
        if name not in self.formats:
            raise error.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].get_writer(**kwargs)

    def tree_source_iter(self, stream, name, **kwargs):
        if name not in self.formats:
            raise error.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].get_tree_source_iter(stream, **kwargs)

###############################################################################
## Client Code Interface

_GLOBAL_DATA_FORMAT_REGISTRY = _DataFormatRegistry()

def register(format, reader, writer, tree_source_iter):
    _GLOBAL_DATA_FORMAT_REGISTRY.add(format, reader, writer, tree_source_iter)

def get_reader(format, **kwargs):
    """
    Returns a reader object of the appropriate format-handler as specified by
    `format`.

    `format` is a string that is name of one of the registered data
    formats, such as `nexus`, `newick`, etc, for which a specialized
    reader is available. If this is not implemented for the format
    specified, then a `UnsupportedFormatError` is raised.

    The following keyword arguments are recognized:

        - `dataset`: All data read from the source will be instantiated
                as objects within this `DataSet` object.
        - `taxon_set`: A`TaxonSet` object. If given, results in all the
                taxa being accessioned into a single `TaxonSet` (and all
                TaxonSetLinked objects instantiated being associated with that
                `TaxonSet` object), even if multiple taxon collection
                definitions are encountered in the source.
        - `exclude_trees`: Trees in the source will be skipped.
        - `exclude_chars`: Characters in the source will be skipped.
        - `encode_splits`: Specifies whether or not splits will be
                automatically-encoded upon a tree being read.

    Other keywords may be implemented by specific readers (e.g. NexusReader,
    NewickReader). Refer to their documentation for details.
    """
    return _GLOBAL_DATA_FORMAT_REGISTRY.get_reader(format, **kwargs)

def get_writer(format, **kwargs):
    """
    Returns a writer object of the appropriate format as specified by `format`.

    `format` is a string that is name of one of the registered data
    formats, such as `nexus`, `newick`, etc, for which a specialized
    writer is available. If this is not implemented for the format
    specified, then a `UnsupportedFormatError` is raised.

    The following keyword arguments are recognized:

        - `dataset`: A `DataSet` object that will be the default source of the
                data to be written.
        - `exclude_trees`: Trees in the `DataSet` or `TaxonDomain` will not be
                written.
        - `exclude_chars`: Characters in the `DataSet` or `TaxonDomain` will
                not written.
    """
    return _GLOBAL_DATA_FORMAT_REGISTRY.get_writer(format, **kwargs)

def tree_source_iter(stream, format, **kwargs):
    """
    Returns an iterator over trees in `format`-formatted data
    in from file-like object source `stream`. Keyword arguments
    are passed to format-specialized implementation of the iterator
    invoked.

    Keyword arguments accepted (handled here):

        - `from_index` 0-based index specifying first tree to actually return
           (raises KeyError if >= #trees)

    Keyword arguments that should be handled by implementing Readers:

        - `taxon_set` specifies the `TaxonSet` object to be attached to the
           trees parsed and manage their taxa. If not specified, then a
           (single) new `TaxonSet` object will be created and for all the
           `Tree` objects.
        - `encode_splits` specifies whether or not split bitmasks will be
           calculated and attached to the edges.
        - `finish_node_func` is a function that will be applied to each node
           after it has been constructed.
        - `edge_len_type` specifies the type of the edge lengths (int or float)

    """
    if "from_index" in kwargs:
        from_index = kwargs["from_index"]
        del(kwargs["from_index"])
    else:
        from_index = 0
    if "write_progress" in kwargs:
        write_progress = kwargs["write_progress"]
        del(kwargs["write_progress"])
    else:
        progress_writer = None
    tree_iter = _GLOBAL_DATA_FORMAT_REGISTRY.tree_source_iter(stream, format, **kwargs)
    for count, t in enumerate(tree_iter):
        if count >= from_index and t is not None:
            if write_progress is not None:
                write_progress("Processing tree at index %d" % count)
            count += 1
            yield t
        else:
            if write_progress is not None:
                write_progress("Skipping tree at index %d" % count)
    if count < from_index:
        raise KeyError("0-based index out of bounds: %d (trees=%d, from_index=[0, %d])" % (from_index, count, count-1))

def multi_tree_source_iter(sources, format, **kwargs):
    """
    Iterates over trees from multiple sources, which may be given as file-like
    objects or filepaths (strings). Note that unless a TaxonSet object is
    explicitly passed using the 'taxon_set' keyword argument, the trees in each
    file will be associated with their own distinct, independent taxon set.
    """
#    if "taxon_set" not in kwargs:
#        kwargs["taxon_set"] = TaxonSet()
    if "write_progress" in kwargs:
        write_progress = kwargs["write_progress"]
        del(kwargs["write_progress"])
    else:
        write_progress = None
    num_sources = len(sources)
    for i, s in enumerate(sources):
        if isinstance(s, str):
            src = open(s, "rU")
        else:
            src = s
        if write_progress is not None:
            write_subprogress = lambda x: write_progress("Tree source %d of %d: %s\n"
                    % (i+1, num_sources, str(x)))
        else:
            write_subprogress = None
        for t in tree_source_iter(src, format, write_progress=write_subprogress, **kwargs):
            yield t
