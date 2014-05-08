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

"""
Exceptions and errors.
"""

import sys
import warnings
import inspect

def get_calling_code_info(stack_level):
    frame = inspect.stack()[stacklevel]
    filename = inspect.getfile(frame[0])
    lineno = inspect.getlineno(frame[0])
    return filename, lineno

def dump_stack(out=None):
    if out is None:
        out = sys.stderr
    for frame, filename, line_num, func, source_code, source_index in inspect.stack()[2:]:
        if source_code is None:
            out.write("{}: {}\n".format(filename, line_num))
        else:
            out.write("{}: {}: {}\n".format(filename, line_num, source_code[source_index].strip()))

# def deprecation_alert(message,
#         logger_obj=None,
#         stacklevel=2,
#         force_warning=False):
#     with warnings.catch_warnings() as w:
#         if force_warning:
#             frame = inspect.stack()[stacklevel]
#             warnings.warn_explicit(
#                     message=message,
#                     category=UserWarning,
#                     filename=inspect.getfile(frame[0]),
#                     lineno=inspect.getlineno(frame[0]))
#         else:
#             warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)
#     if logger_obj:
#         logger_obj.warning(message)

class DataError(Exception):

    def __init__(self,
            message=None,
            line_num=None,
            col_num=None,
            filename=None,
            stream=None):
        Exception.__init__(self)
        self.line_num = line_num
        self.col_num = col_num
        self.msg = message
        self.filename = None
        self.decorate_with_name(filename=filename, stream=stream)

    def decorate_with_name(self,
            filename=None,
            stream=None):
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
            f =  " '{}'".format(self.filename)
        if self.line_num is not None:
            l =  " on line {}".format(self.line_num)
        if self.col_num is not None:
            c =  " at column {}".format(self.col_num)
        return "Error parsing data source{}{}{}: {}".format(f, l, c, self.msg)

class DataParseError(DataError):

    def __init__(self,
            message=None,
            line_num=None,
            col_num=None,
            filename=None,
            stream=None):
        DataError.__init__(self,
                message=message,
                line_num=line_num,
                col_num=col_num,
                filename=filename,
                stream=stream)

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

    def __init__(self,
            message=None,
            func_name=None,
            max_args=None,
            args=None):
        if message is None and (func_name is not None and max_args):
            message = "{}() takes a maximum of {} arguments ({} given)".format(func_name, max_args, len(args))
        TypeError.__init__(self, message)

class InvalidArgumentValueError(ValueError):

    def __init__(self, message=None, func_name=None, arg=None):
        if message is None and (func_name is not None and arg is not None):
            message = "{}() does not accept objects of type '{}' as an argument".format(func_name, arg.__class__.__name__)
        ValueError.__init__(self, message)

class MultipleInitializationSourceError(TypeError):
    def __init__(self, message=None, class_name=None, arg=None):
        if message is None and (class_name is not None and arg is not None):
            message = "{}() does not accept data 'stream' or 'schema' arguments when initializing with another object".format(class_name)
        TypeError.__init__(self, message)

class TaxonNamespaceError(ValueError):

    def __init__(self, message=None):
        # common = "TaxonNamespace references not identical"
        common = "TaxonNamespace references not identical"
        if message is not None:
            message = common + ": " + message
        else:
            message = common
        ValueError.__init__(self, message)

