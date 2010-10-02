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
Base classes for all readers/parsers and formatters/writers.
"""

import os
from cStringIO import StringIO
from dendropy.utility import error

###############################################################################
## KEYWORD ARGUMENT PROCESSING

def get_format_from_kwargs(kwdict):
    """
    Extract and return schema term from keyword dictionary, deleting term if
    encountered.
    """
    if "schema" in kwdict:
        schema = kwdict["schema"]
        del(kwdict["schema"])
        return schema
    else:
        return None

def require_format_from_kwargs(kwdict):
    """
    As with `get_format_arg`, but raises Exception if not specified.
    """
    schema = get_format_from_kwargs(kwdict)
    if schema is None:
        raise error.UnspecifiedSchemaError("Must specify `schema`.")
    return schema

def extract_kwarg(kwdict, kw, default=None):
    if kw in kwdict:
        kwarg = kwdict[kw]
        del(kwdict)
        return kwarg
    else:
        return default

###############################################################################
## IOService

class IOService(object):
    """
    Base class for all readers/writers.
    """

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keyword arguments:

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
        self.dataset = extract_kwarg(kwargs, "dataset", None)
        self.attached_taxon_set = extract_kwarg(kwargs, "taxon_set", None)
        self.exclude_trees = extract_kwarg(kwargs, "exclude_trees", False)
        self.exclude_chars = extract_kwarg(kwargs, "exclude_chars", False)

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
        __init__ creates a `DataReader` object. See `IOClient` for details on
        keyword arguments recognized. In addition, the keyword `encode_splits`
        specifies whether or not splits will be automatically-encoded upon
        a tree being read.
        """
        IOService.__init__(self, **kwargs)
        self.encode_splits = kwargs.get("encode_splits", False)

    def read(self, stream):
        """
        Reads data from the file-like object `stream`, and populates
        and returns the attached `DataSet` object or a new `DataSet` object
        if none is attached.
        """
        raise NotImplementedError

    def get_default_taxon_set(self, **kwargs):
        """
        Returns an appropriate TaxonSet object, based on current settings.
        """
        if self.dataset is None:
            raise TypeError("'dataset' is not defined")
        if "taxon_set" in kwargs:
            self.attached_taxon_set = kwargs["taxon_set"]
        if self.dataset.attached_taxon_set is not None:
            if self.attached_taxon_set is not None \
                    and self.dataset.attached_taxon_set is not self.attached_taxon_set:
                raise TypeError("DataSet is attached to different TaxonSet than that specified by 'taxon_set'")
            self.attached_taxon_set = self.dataset.attached_taxon_set
        if self.attached_taxon_set is not None:
            if self.attached_taxon_set not in self.dataset.taxon_sets:
                self.dataset.add(self.attached_taxon_set)
            return self.attached_taxon_set
        else:
            return self.dataset.new_taxon_set(**kwargs)

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
        __init__ creates a `DataWriter` object. See `IOClient` for details on
        keyword arguments recognized.
        """
        IOService.__init__(self, **kwargs)

    def write(self, stream):
        """
        Writes data in the attached `DataSet` object to a destination given
        by the file-like object `stream`.
        """
        raise NotImplementedError

###############################################################################
## Readable

class Readable(object):
    """
    Data object that can be instantiated using a `DataReader` service.
    """

    def get_from_stream(cls, src, schema, **kwargs):
        """
        Factory method to return new object of this class from file-like
        object `src`.
        """
        readable = cls(**kwargs)
        readable.read_from_stream(src, schema, **kwargs)
        return readable
    get_from_stream = classmethod(get_from_stream)

    def get_from_path(cls, src, schema, **kwargs):
        """
        Factory method to return new object of this class from file
        specified by string `src`.
        """
        readable = cls(**kwargs)
        readable.read_from_path(src, schema, **kwargs)
        return readable
    get_from_path = classmethod(get_from_path)

    def get_from_string(cls, src, schema, **kwargs):
        """
        Factory method to return new object of this class from string `src`.
        """
        readable = cls(**kwargs)
        readable.read_from_string(src, schema, **kwargs)
        return readable
    get_from_string = classmethod(get_from_string)

    def __init__(self, *args, **kwargs):
        pass

    def process_source_kwargs(self, **kwargs):
        """
        If `stream` is specified, then process_source_kwargs:

            # checks that `schema` keyword argument is specified, and then
            # calls self.read()
        
        If `stream` is not specified then `source_string` and `source_filepath`
            kwarg arguments are checked. The effect of thes is similar to calling
            read_from_string and read_from_path
        """
        if "stream" in kwargs:
            stream = kwargs["stream"]
            del(kwargs["stream"])
            schema = require_format_from_kwargs(kwargs)
            self.read(stream=stream, schema=schema, **kwargs)
        elif "source_string" in kwargs:
            as_str = kwargs["source_string"]
            del(kwargs["source_string"])
            kwargs["stream"] = StringIO(as_str)
            return self.process_source_kwargs(**kwargs)
        elif "source_file" in kwargs:
            fp = kwargs["source_filepath"]
            fo = open(os.path.expandvars(os.path.expanduser(filepath)), "rU")
            del(kwargs["source_filepath"])
            kwargs["stream"] = fo
            return self.process_source_kwargs(**kwargs)

    def read(self, stream, schema, **kwargs):
        """
        Populates/constructs objects of this type from `schema`-formatted
        data in the file-like object source `stream`.
        """
        raise NotImplementedError

    def read_from_stream(self, fileobj, schema, **kwargs):
        """
        Reads from file (exactly equivalent to just `read()`, provided
        here as a separate method for completeness.
        """
        return self.read(stream=fileobj, schema=schema, **kwargs)

    def read_from_path(self, filepath, schema, **kwargs):
        """
        Reads from file specified by `filepath`.
        """
        f = open(os.path.expandvars(os.path.expanduser(filepath)), "rU")
        return self.read(stream=f, schema=schema, **kwargs)

    def read_from_string(self, src_str, schema, **kwargs):
        """
        Reads a string object.
        """
        s = StringIO(src_str)
        return self.read(stream=s, schema=schema, **kwargs)

###############################################################################
## Writeable

class Writeable(object):
    """
    Data object that can be instantiated using a `DataReader` service.
    """

    def __init__(self, *args, **kwargs):
        pass

    def write(self, stream, schema, **kwargs):
        """
        Writes the object to the file-like object `stream` in `schema`
        schema.
        """
        raise NotImplementedError

    def write_to_stream(self, dest, schema, **kwargs):
        """
        Writes to file-like object `dest`.
        """
        return self.write(stream=dest, schema=schema, **kwargs)

    def write_to_path(self, dest, schema, **kwargs):
        """
        Writes to file specified by `dest`.
        """
        f = open(os.path.expandvars(os.path.expanduser(dest)), "w")
        return self.write(stream=f, schema=schema, **kwargs)

    def as_string(self, schema, **kwargs):
        """
        Composes and returns string representation of the data.
        """
        s = StringIO()
        self.write(stream=s, schema=schema, **kwargs)
        return s.getvalue()

