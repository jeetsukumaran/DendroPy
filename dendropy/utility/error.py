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
Exceptions and error.
"""
import sys

class DataError(Exception):

    def __init__(self, message=None, row=None, column=None, filename=None, stream=None):
        Exception.__init__(self)
        self.row = row
        self.column = column
        self.msg = message
        self.filename = None
        self.decorate_with_name(filename=filename, stream=stream)

    def decorate_with_name(self, filename=None, stream=None):
        if filename is not None:
            self.filename = filename
        if stream is not None:
            try:
                self.filename = stream.name
            except AttributeError:
                pass

    def __str__(self):
        f, l, c = "", "", ""
        if self.filename:
            f =  ' "%s"' % self.filename
        if self.row is not None:
            l =  " on line %d" % self.row
        if self.column is not None:
            c =  " at column %d" % self.column
        return 'Error parsing data source%s%s%s: %s' % (f, l, c, self.msg)

class DataParseError(DataError):

    def __init__(self, message=None, row=None, column=None, filename=None, stream=None):
        DataError.__init__(self, message=message, row=row, column=column, filename=filename, stream=stream)

class UnsupportedSchemaError(NotImplementedError):

    def __init__(self, *args, **kwargs):
        NotImplementedError.__init__(self, *args, **kwargs)

class UnspecifiedSchemaError(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class UnspecifiedSourceError(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class TooManyArgumentsError(TypeError):

    def __init__(self, message=None, func_name=None, max_args=None, args=None):
        if message is None and (func_name is not None and max_args):
            message = "%s() takes a maximum of %d arguments (%d given)" % (func_name, max_args, len(args))
        TypeError.__init__(self, message)

class InvalidArgumentValueError(ValueError):

    def __init__(self, message=None, func_name=None, arg=None):
        if message is None and (func_name is not None and arg is not None):
            message = "%s() does not accept objects of type '%s' as an argument" % (func_name, arg.__class__.__name__)
        ValueError.__init__(self, message)

class MultipleInitializationSourceError(TypeError):
    def __init__(self, message=None, class_name=None, arg=None):
        if message is None and (class_name is not None and arg is not None):
            message = "%s() does not accept data 'stream' or 'schema' arguments when initializing with another object" % (class_name)
        TypeError.__init__(self, message)


