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
Maps data formats to their specific parsers/writers.
"""

from dendropy.utility import error
from dendropy.utility import iosys
from dendropy.utility.containers import OrderedCaselessDict

###############################################################################
## DataFormat

class DataFormat(object):

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
## DataFormatRegistry

class DataFormatRegistry(object):

    def __init__(self):
        self.formats = OrderedCaselessDict()

    def add_format(self, data_format):
        self.formats[data_format.name] = data_format

    def remove_format(self, data_format):
        del(self.formats[data_format.name])

    def add(self, name, reader_type=None, writer_type=None, tree_source_iter=None):
        self.formats[name] = DataFormat(name,
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


