#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Maps data formats to their specific parsers/writers.
"""

from dendropy.utility import error
from dendropy.utility import iosys
from dendropy.utility.containers import OrderedCaselessDict

###############################################################################
## DataSchema

class DataSchema(object):

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
            raise error.UnsupportedSchemaError("Reading is not currently supported for data schema '%s'" % self.name)
        return self.reader_type(**kwargs)

    def get_writer(self, **kwargs):
        if self.writer_type is None:
            raise error.UnsupportedSchemaError("Writing is not currently supported for data schema '%s'" % self.name)
        return self.writer_type(**kwargs)

    def get_tree_source_iter(self, stream, **kwargs):
        if self.tree_source_iter is None:
            raise error.UnsupportedSchemaError("Iteration over source trees not currently supported for data schema '%s'" % self.name)
        return self.tree_source_iter(stream, **kwargs)

###############################################################################
## DataSchemaRegistry

class DataSchemaRegistry(object):

    def __init__(self):
        self.formats = OrderedCaselessDict()

    def add_format(self, data_format):
        self.formats[data_format.name] = data_format

    def remove_format(self, data_format):
        del(self.formats[data_format.name])

    def add(self, name, reader_type=None, writer_type=None, tree_source_iter=None):
        self.formats[name] = DataSchema(name,
                reader_type=reader_type,
                writer_type=writer_type,
                tree_source_iter=tree_source_iter)

    def remove(self, name):
        del(self.formats[name])

    def get_reader(self, name, **kwargs):
        ## hack to avoid confusing users ##
#        if name.lower() == "fasta":
#            raise error.UnsupportedSchemaError("FASTA data needs to be specified as 'dnafasta', 'rnafasta', or 'proteinfasta'")
        if name not in self.formats:
            raise error.UnsupportedSchemaError("'%s' is not a recognized data schema name" % name)
        return self.formats[name].get_reader(**kwargs)

    def get_writer(self, name, **kwargs):
        if name not in self.formats:
            raise error.UnsupportedSchemaError("'%s' is not a recognized data schema name" % name)
        return self.formats[name].get_writer(**kwargs)

    def tree_source_iter(self, stream, name, **kwargs):
        if name not in self.formats:
            raise error.UnsupportedSchemaError("'%s' is not a recognized data schema name" % name)
        return self.formats[name].get_tree_source_iter(stream, **kwargs)


