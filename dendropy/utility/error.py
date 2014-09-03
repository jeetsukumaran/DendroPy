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
Errors, exceptions, warnings, etc.
"""

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import re
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

class CriticalDeprecationWarning(DeprecationWarning):
    pass

def critical_deprecation_alert(message, stacklevel=4):
    stack = inspect.stack()
    frame = None
    while frame is None and stacklevel > 0:
        try:
            frame, filename, line_num, func, source_code, source_index = stack[stacklevel]
        except IndexError:
            stacklevel -= 1
    warnings.simplefilter('once')
    warnings.warn_explicit(
            message=message,
            category=CriticalDeprecationWarning,
            filename=inspect.getfile(frame),
            lineno=inspect.getlineno(frame),
            )

def dendropy_construct_migration_warning(
        old_construct,
        new_construct,
        stacklevel=5):
    message = (
        """Instead of:\n"""
        """    {}\n"""
        """Use:\n"""
        """    {}\n""").format(
                old_construct, new_construct)
    critical_deprecation_alert(message, stacklevel)

def dendropy_module_migration_warning(
        old_module_path,
        new_module_path,
        stacklevel=5):
    x1, x2 = old_module_path.rsplit(".", 1)
    y1, y2 = new_module_path.rsplit(".", 1)
    message = (
        """\nThis module has been moved to:\n"""
        """    {}\n"""
        """Instead of e.g.:\n"""
        """    from {} import {}\n"""
        """Use:\n"""
        """    from {} import {}\n""").format(
                new_module_path, x1, x2, y1, y2)
    critical_deprecation_alert(message, stacklevel)

def dendropy_migration_warning(
        old_api_construct,
        new_api_construct,
        old_api_construct_detection_pattern,
        start_offset=3,
        stack_dump_out=None):
    pattern_to_match = re.compile(old_api_construct_detection_pattern)
    found = False
    stacklevel = 0
    stack = inspect.stack()
    for level, (frame, filename, line_num, func, source_code, source_index) in enumerate(stack[start_offset:]):
        if source_code and len(source_code) >= 1:
            source_line = source_code[0].replace("\n","").strip()
        else:
            source_line = ""
        if found:
            leader = "   "
        elif pattern_to_match.match(source_line):
            found = True
            stacklevel = level + start_offset
            leader = ">> "
        else:
            leader = "   "
        if stack_dump_out:
            stack_dump_out.write('{}File "{}", line {}:\n{}    {}\n'.format(
                leader,
                filename,
                line_num,
                "    ",
                source_line))
    if not found:
        stacklevel = start_offset + 1
    frame, filename, line_num, func, source_code, source_index = stack[stacklevel]
    message = "'{}' has been deprecated since DendroPy 4.0.0 and will no longer be supported in future releases: instead of '{}', use '{}'".format(old_api_construct, old_api_construct, new_api_construct)
    warnings.simplefilter('once')
    warnings.warn_explicit(
            message=message,
            category=CriticalDeprecationWarning,
            filename=inspect.getfile(frame),
            lineno=inspect.getlineno(frame),
            )

class ImmutableTaxonNamespaceError(TypeError):
    def __init__(self, message):
        TypeError.__init__(self, message)

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

class TaxonNamespaceReconstructionError(ValueError):

    def __init__(self, message=None):
        ValueError.__init__(self, message)

def process_attached_taxon_namespace_directives(kwargs_dict):
    """
    The following idioms are supported:

        `taxon_namespace=tns`
            Attach `tns` as the bound (single, unified) taxonomic namespace
            reference for all objects.
        `attached_taxon_namespace=tns`
            Attach `tns` as the bound (single, unified) taxonomic namespace
            reference for all objects.
        `attach_taxon_namespace=True, attached_taxon_namespace=tns`
            Attach `tns` as the bound (single, unified) taxonomic namespace
            reference for all objects.
        `attach_taxon_namespace=True`
            Create a *new* :class:`TaxonNamespace` and set it as the bound
            (single, unified) taxonomic namespace reference for all
            objects.
    """
    deprecated_kw = [
            "taxon_namespace",
            "attach_taxon_namespace",
            "attached_taxon_namespace",
            "taxon_set",
            "attach_taxon_set",
            "attached_taxon_set",
            ]
    for kw in deprecated_kw:
        if kw in kwargs_dict:
            raise TypeError("'{}' is no longer supported as a keyword argument. Use the instance method 'attach_taxon_namespace()' of the data object instead".format(kw))
    taxon_namespace = None
    attach_taxon_namespace = False
    if ( ("taxon_set" in kwargs_dict or "taxon_namespace" in kwargs_dict)
            and ("attached_taxon_set" in kwargs_dict or "attached_taxon_namespace" in kwargs_dict)
            ):
        raise TypeError("Cannot specify both 'taxon_namespace'/'taxon_set' and 'attached_taxon_namespace'/'attached_taxon_set' together")
    if "taxon_set" in kwargs_dict:
        if "taxon_namespace" in kwargs_dict:
            raise TypeError("Both 'taxon_namespace' and 'taxon_set' cannot be specified simultaneously: use 'taxon_namespace' ('taxon_set' is only supported for legacy reasons)")
        kwargs_dict["taxon_namespace"] = kwargs_dict["taxon_set"]
        del kwargs_dict["taxon_set"]
    if "attached_taxon_set" in kwargs_dict:
        if "attached_taxon_namespace" in kwargs_dict:
            raise TypeError("Both 'attached_taxon_namespace' and 'attached_taxon_set' cannot be specified simultaneously: use 'attached_taxon_namespace' ('attached_taxon_set' is only supported for legacy reasons)")
        kwargs_dict["attached_taxon_namespace"] = kwargs_dict["attached_taxon_set"]
        del kwargs_dict["attached_taxon_set"]
    if "taxon_namespace" in kwargs_dict:
        taxon_namespace = kwargs_dict.pop("taxon_namespace", None)
        attach_taxon_namespace = True
    elif "attached_taxon_namespace" in kwargs_dict:
        taxon_namespace = kwargs_dict["attached_taxon_namespace"]
        if not isinstance(taxon_namespace, TaxonNamespace):
            raise TypeError("'attached_taxon_namespace' argument must be an instance of TaxonNamespace")
        attach_taxon_namespace = True
    else:
        taxon_namespace = None
        attach_taxon_namespace = kwargs_dict.get("attach_taxon_namespace", False)
    kwargs_dict.pop("taxon_namespace", None)
    kwargs_dict.pop("attach_taxon_namespace", None)
    kwargs_dict.pop("attached_taxon_namespace", None)
    return (attach_taxon_namespace, taxon_namespace)
