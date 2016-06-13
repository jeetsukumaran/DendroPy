#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
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

import sys
import re
import warnings
import inspect
import subprocess

class ImmutableTaxonNamespaceError(TypeError):
    def __init__(self, message):
        TypeError.__init__(self, message)

class DataParseError(Exception):

    def __init__(self,
            message=None,
            line_num=None,
            col_num=None,
            filename=None,
            stream=None):
        Exception.__init__(self)
        self.line_num = line_num
        self.col_num = col_num
        self.message = message
        self.stream = stream
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
        return "Error parsing data source{}{}{}: {}".format(f, l, c, self.message)

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

class TaxonNamespaceIdentityError(ValueError):
    def __init__(self, o1, o2):
        message = "Non-identical taxon namespace references: {} is not {}".format(
                "<TaxonNamespace object at {}>".format(str(hex(id((o1.taxon_namespace))))),
                "<TaxonNamespace object at {}>".format(str(hex(id((o2.taxon_namespace))))),
                )
        ValueError.__init__(self, message)

class SingleTaxonAssemblageException(ValueError):
    def __init__(self, message=None):
        ValueError.__init__(self, message)

class NullAssemblageException(ValueError):
    def __init__(self, message=None):
        ValueError.__init__(self, message)

class MixedRootingError(ValueError):
    def __init__(self, message=None):
        ValueError.__init__(self, message)

class TaxonNamespaceReconstructionError(ValueError):
    def __init__(self, message=None):
        ValueError.__init__(self, message)

class UltrametricityError(ValueError):
    def __init__(self, message=None):
        ValueError.__init__(self, message)

class ProcessFailedException(Exception):
    """Exception to be raised when branching process results in all lineages going extinct."""
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class TreeSimTotalExtinctionException(ProcessFailedException):
    """Exception to be raised when branching process results in all lineages going extinct."""
    def __init__(self, *args, **kwargs):
        ProcessFailedException.__init__(self, *args, **kwargs)

class NullLeafSetException(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class SeedNodeDeletionException(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class InvalidMultispeciesCoalescentStructureError(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class ExternalServiceError(Exception):

    def __init__(self,
            service_name,
            invocation_command,
            service_input,
            returncode,
            stdout,
            stderr):
        self.service_name = service_name
        self.invocation_command = invocation_command
        self.service_input = service_input
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.message = self.compose_message()

    def __str__(self):
        return self.compose_message()

    def compose_message(self):
        parts = []
        parts.append("- Service process {} exited with return code: {}".format(self.service_name, self.returncode))
        parts.append("- Service invoked with command:")
        parts.append("   {}".format(self.invocation_command))
        parts.append("<<<<<<< (SERVICE STANDARD INPUT)")
        parts.append(self.service_input)
        parts.append(">>>>>>>")
        parts.append("<<<<<<< (SERVICE STANDARD OUTPUT)")
        parts.append(self.stdout)
        parts.append(">>>>>>>")
        parts.append("<<<<<<< (SERVICE STANDARD ERROR)")
        parts.append(self.stderr)
        parts.append(">>>>>>>")
        return "\n".join(parts)
