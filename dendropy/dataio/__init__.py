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

from dendropy.utility import errors
from dendropy.utility import iosys
from dendropy.utility.containers import OrderedCaselessDict
from dendropy.dataio import newick
from dendropy.dataio import nexus
from dendropy.dataio import fasta
from dendropy.dataio import phylip
from dendropy.dataio import nexml

###############################################################################
## Client code interface

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
                as objects within this `Dataset` object.
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
    return ioregister.get_reader(format, **kwargs)

def get_writer(format, **kwargs):
    """
    Returns a writer object of the appropriate format as specified by `format`.

    `format` is a string that is name of one of the registered data
    formats, such as `nexus`, `newick`, etc, for which a specialized
    writer is available. If this is not implemented for the format
    specified, then a `UnsupportedFormatError` is raised.

    The following keyword arguments are recognized:

        - `dataset`: A `Dataset` object that will be the default source of the
                data to be written.
        - `exclude_trees`: Trees in the `Dataset` or `TaxonDomain` will not be
                written.
        - `exclude_chars`: Characters in the `Dataset` or `TaxonDomain` will
                not written.
    """
    return ioregister.get_writer(format, **kwargs)

def tree_source_iter(istream, format, **kwargs):
    """
    Returns an iterator over trees in `format`-formatted data
    in from file-like object source `stream`. Keyword arguments
    are passed to format-specialized implementation of the iterator
    invoked.
    """
    return ioregister.tree_source_iter(istream, format, **kwargs)

def write_tree_list(tree_list, format, ostream, **kwargs):
    """
    Writes `tree_list`, a `TreeList` object in `format` format to
    a destination given by file-like object `ostream`.

    `format` is a string that is name of one of the registered data
    formats, such as `nexus`, `newick`, etc, for which a specialized
    tree list writer is available. If this is not implemented for the format
    specified, then a `UnsupportedFormatError` is raised.

    Additionally, for some formats, the following keywords are recognized:

        - `edge_lengths` : if False, edges will not write edge lengths [True]
        - `internal_labels` : if False, internal labels will not be written [True]
    """
    return ioregister.write_tree_list(tree_list, ostream, format, **kwargs)


###############################################################################
## Under the hood ...

class _DataFormat(object):

    def __init__(self,
                 name,
                 reader_type=None,
                 writer_type=None,
                 tree_source_iter=None,
                 tree_list_writer=None):
        self.name = name
        self.reader_type = reader_type
        self.writer_type = writer_type
        self.tree_source_iter = tree_source_iter
        self.tree_list_writer = tree_list_writer

    def has_reader(self):
        return self.reader_type is not None

    def has_writer(self):
        return self.writer_type is not None

    def has_tree_source_iter(self):
        return self.tree_source_iter is not None

    def has_tree_list_writer(self):
        return self.tree_list_writer is not None

    def get_reader(self, **kwargs):
        if self.reader_type is None:
            raise errors.UnsupportedFormatError("Reading is not currently supported for data format '%s'" % self.name)
        return self.reader_type(**kwargs)

    def get_writer(self, **kwargs):
        if self.writer_type is None:
            raise errors.UnsupportedFormatError("Writing is not currently supported for data format '%s'" % self.name)
        return self.writer_type(**kwargs)

    def get_tree_source_iter(self, istream, **kwargs):
        if self.tree_source_iter is None:
            raise errors.UnsupportedFormatError("Iteration over source trees not currently supported for data format '%s'" % self.name)
        return self.tree_source_iter(istream, **kwargs)

    def write_tree_list(self, tree_list, ostream, **kwargs):
        if self.tree_list_writer is None:
            raise errors.UnsupportedFormatError("Writing of stand-alone tree lists is not currently supported for data format '%s'" % self.name)
        self.tree_list_writer(tree_list, ostream, **kwargs)

class _DataFormatRegister(object):

    def __init__(self):
        self.formats = OrderedCaselessDict()

    def add_format(self, data_format):
        self.formats[data_format.name] = data_format

    def remove_format(self, data_format):
        del(self.formats[data_format.name])

    def add(self, name, reader_type=None, writer_type=None, tree_source_iter=None, tree_list_writer=None):
        self.formats[name] = _DataFormat(name,
                reader_type=reader_type,
                writer_type=writer_type,
                tree_source_iter=tree_source_iter,
                tree_list_writer=tree_list_writer)

    def remove(self, name):
        del(self.formats[name])

    def get_reader(self, name, **kwargs):
        if name not in self.formats:
            raise errors.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].get_reader(**kwargs)

    def get_writer(self, name, **kwargs):
        if name not in self.formats:
            raise errors.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].get_writer(**kwargs)

    def tree_source_iter(self, istream, name, **kwargs):
        if name not in self.formats:
            raise errors.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].get_tree_source_iter(istream, **kwargs)

    def write_tree_list(self, tree_list, name, ostream, **kwargs):
        if name not in self.formats:
            raise errors.UnsupportedFormatError("Format '%s' is not a recognized data format name" % name)
        return self.formats[name].write_tree_list(tree_list, ostream, **kwargs)

ioregister = _DataFormatRegister()
ioregister.add("newick", newick.NewickReader, newick.NewickWriter, newick.tree_source_iter, newick.write_tree_list)
ioregister.add("nexus", nexus.NexusReader, nexus.NexusWriter, nexus.tree_source_iter, None)
ioregister.add("fasta", None, fasta.FastaWriter, None, None)
ioregister.add("dnafasta", fasta.DNAFastaReader, fasta.FastaWriter, None, None)
ioregister.add("rnafasta", fasta.RNAFastaReader, fasta.FastaWriter, None, None)
ioregister.add("proteinfasta", fasta.ProteinFastaReader, fasta.FastaWriter, None, None)
ioregister.add("phylip", None, phylip.PhylipWriter, None, None)
ioregister.add("nexml", nexml.NexmlReader, nexml.NexmlWriter, None, None)
