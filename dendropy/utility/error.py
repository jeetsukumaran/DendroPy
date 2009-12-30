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
Exceptions and error.
"""
import sys

class DataSyntaxError(Exception):

    def __init__(self, row=None, column=None, message=None):
        Exception.__init__(self)
        self.row = row
        self.column = column
        self.msg = message

    def __str__(self):
        if self.row is None:
            t = ""
        else:
            t =  "Line %d in data source: " % self.row
        return '%s%s' % (t, self.msg)

class UnsupportedFormatError(NotImplementedError):

    def __init__(self, *args, **kwargs):
        NotImplementedError.__init__(self, *args, **kwargs)

class UnspecifiedFormatError(Exception):

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


