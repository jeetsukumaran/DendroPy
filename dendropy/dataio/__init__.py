#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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


import collections
from dendropy.dataio import newickreader
from dendropy.utility import container

_IOServices = collections.namedtuple(
        "_IOServices",
        ["reader", "writer", "tree_iterator"]
        )

_IO_SERVICE_REGISTRY = container.CaseInsensitiveDict()
_IO_SERVICE_REGISTRY["newick"] = _IOServices(newickreader.NewickReader, None, None)

def get_reader(schema, **kwargs):
    try:
        reader_type =_IO_SERVICE_REGISTRY[schema].reader
        if reader_type is None:
            raise KeyError
        reader = reader_type(**kwargs)
        return reader
    except KeyError:
        raise NotImplementedError("'{}' is not a supported data reading schema".format(schema))

def get_writer(
        schema,
        **kwargs):
    pass

