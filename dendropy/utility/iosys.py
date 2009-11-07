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
Base classes for all readers/parsers and formatters/writers.
"""

import os
from cStringIO import StringIO
from dendropy.utility import errors

###############################################################################
## KEYWORD ARGUMENT PROCESSING

def get_format_from_kwargs(kwdict):
    """
    Extract and return format term from keyword dictionary, deleting term if
    encountered.
    """
    if "format" in kwdict:
        format = kwdict["format"]
        del(kwdict["format"])
        return format
    else:
        return None

def require_format_from_kwargs(kwdict):
    """
    As with `get_format_arg`, but raises Exception if not specified.
    """
    format = get_format_from_kwargs(kwdict)
    if format is None:
        raise errors.UnspecifiedFormatError("Must specify `format`.")
    return format


###############################################################################
## IOService

class IOService(object):
    """
    Base class for all readers/writers.
    """

    def __init__(self, **kwargs):
        """
        The following keyword arguments are recognized:

        - `dataset`:  A `DataSet` object. For input clients, all data read
                from the source will be instantiated as objects within
                this `DataSet` object. For output clients, the `DataSet`
                object will be the source of the data to be written.
        - `taxon_set`: A`TaxonSet` object. For input clients, results in all the
                taxa being accessioned into a single `TaxonSet` (and all
                TaxonSetLinked objects instantiated being associated with that
                `TaxonSet` object), even if multiple taxon collection
                definitions are encountered in the source. For output clients,
                only data associated with the specified `TaxonSet` object
                will be written.
        - `exclude_trees`: Trees in the source will be skipped.
        - `exclude_chars`: Characters in the source will be skipped.
        """
        self.dataset = kwargs.get("dataset", None)
        self.bound_taxon_set = kwargs.get("taxon_set", None)
        self.exclude_trees = kwargs.get("exclude_trees", False)
        self.exclude_chars = kwargs.get("exclude_chars", False)

###############################################################################
## DataReader

class DataReader(IOService):
    """
    Instantiates DendroPy phylogenetic data objects based data given in
    input sources. Abstract class, to be implemented by derived classes
    specializing in particular data formats.
    """

    def __init__(self, **kwargs):
        """
        Instantiates a `DataReader` object. See `IOClient` for details on
        keyword arguments recognized. In addition, the keyword `encode_splits`
        specifies whether or not splits will be automatically-encoded upon
        a tree being read.
        """
        IOService.__init__(self, **kwargs)
        self.encode_splits = kwargs.get("encode_splits", False)

    def read(self, stream, **kwargs):
        """
        Reads data from the file-like object `stream`, and populates
        and returns the bound `DataSet` object or a new `DataSet` object
        if none is bound.
        """
        raise NotImplementedError

###############################################################################
## DataReader

class DataWriter(IOService):
    """
    Writes DendroPy phylogenetic data object. Abstract class, to be
    implemented by derived classes specializing in particular data
    formats.
    """

    def __init__(self, **kwargs):
        """
        Instantiates a `DataWriter` object. See `IOClient` for details on
        keyword arguments recognized.
        """
        IOService.__init__(self, **kwargs)

    def write(self, stream, **kwargs):
        """
        Writes data in the bound `DataSet` object to a destination given
        by the file-like object `stream`.
        """
        raise NotImplementedError

###############################################################################
## Readable

class Readable(object):
    """
    Data object that can be instantiated using a `DataReader` service.
    """

    def __init__(self, *args, **kwargs):
        pass

    def process_source_kwargs(self, **kwargs):
        if "stream" in kwargs:
            stream = kwargs["stream"]
            del(kwargs["stream"])
            format = require_format_from_kwargs(kwargs)
            self.read(stream=stream, format=format, **kwargs)

    def read(self, stream, format, **kwargs):
        """
        Populates/constructs objects of this type from `format`-formatted
        data in the file-like object source `stream`.
        """
        raise NotImplementedError

    def read_from_file(self, fileobj, format, **kwargs):
        """
        Reads from file (exactly equivalent to just `read()`, provided
        here as a separate method for completeness.
        """
        return self.read(stream=fileobj, format=format, **kwargs)

    def read_from_path(self, filepath, format, **kwargs):
        """
        Reads from file specified by `filepath`.
        """
        f = os.expandvars(os.expanduser(filepath))
        return self.read(stream=f, format=format, **kwargs)

    def read_from_string(self, src_str, format, **kwargs):
        """
        Reads a string object.
        """
        s = StringIO(src_str)
        return self.read(stream=s, format=format, **kwargs)

###############################################################################
## Writeable

class Writeable(object):
    """
    Data object that can be instantiated using a `DataReader` service.
    """

    def __init__(self, *args, **kwargs):
        pass

    def write(self, stream, format, **kwargs):
        """
        Writes the object to the file-like object `stream` in `format`
        format.
        """
        raise NotImplementedError

    def write_to_file(self, fileobj, format, **kwargs):
        """
        Writes to file-like object `fileobj`.
        """
        return self.write(stream=fileobj, format=format, **kwargs)

    def write_to_path(self, filepath, format, **kwargs):
        """
        Writes to file specified by `filepath`.
        """
        f = os.expandvars(os.expanduser(filepath))
        return self.write(stream=f, format=format, **kwargs)

    def as_string(self, format, **kwargs):
        """
        Composes and returns string representation of the data.
        """
        s = StringIO(src_str)
        self.write(stream=s, format=format, **kwargs)
        return s.getvalue()

