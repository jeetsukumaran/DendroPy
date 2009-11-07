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

class DataFormatError(Exception):

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

    def __init__(self, func_name, max_args, args):
        TypeError.__init__(self,
            "%s() takes a maximum of %d arguments (%d given)" \
                % (func_name, max_args, len(args)))

class InvalidArgumentTypeError(TypeError):

    def __init__(self, func_name, arg):
        TypeError.__init__(self,
            "%s() does not accept objects of type '%s' as an argument" \
                % (func_name, arg.__class__.__name__))

class MultipleInitializationSourceError(TypeError):
    def __init__(self, class_name, arg):
        TypeError.__init__(self,
            "%s() does not accept data 'stream' or 'format' arguments when initializing with a '%s' object" \
                % (class_name, arg.__class__.__name__))
