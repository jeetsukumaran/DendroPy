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
Path mapping for various test resources.
"""

import os
from dendropy.utility import texttools
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

try:
    import pkg_resources
    _LOG.info("Using local pkg_resources path mapping")
    TESTS_DIR = pkg_resources.resource_filename("dendropy", "test")
    SCRIPTS_DIR = pkg_resources.resource_filename("dendropy", "scripts")
except:
    _LOG.info("Using local filesystem path mapping")
    TESTS_DIR = os.path.dirname(__file__)
    SCRIPTS_DIR = os.path.join(os.path.pardir, "scripts")

TESTS_DATA_DIR = os.path.join(TESTS_DIR, "data")
TESTS_OUTPUT_DIR = os.path.join(TESTS_DIR, "output")

def tree_source_stream(filename):
    return open(tree_source_path(filename), "rU")

def tree_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(TESTS_DATA_DIR, "trees", filename)

def char_source_stream(filename):
    return open(char_source_path(filename), "rU")

def char_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(TESTS_DATA_DIR, "chars", filename)

def mixed_source_stream(filename):
    return open(mixed_source_path(filename), "rU")

def mixed_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(TESTS_DATA_DIR, "mixed", filename)

def data_source_stream(filename):
    return open(data_source_path(filename), "rU")

def data_source_path(filename=None):
    if filename is None:
        filename = ""
    elif isinstance(filename, list):
        filename = os.path.sep.join(filename)
    return os.path.join(TESTS_DATA_DIR, filename)

def named_output_stream(filename=None, suffix_timestamp=True):
    return open(named_output_path(filename=filename, suffix_timestamp=suffix_timestamp), "w")

def named_output_path(filename=None, suffix_timestamp=True):
    if filename is None:
        filename = ""
    else:
        if isinstance(filename, list):
            filename = os.path.sep.join(filename)
        if suffix_timestamp:
            filename = "%s.%s" % (filename, texttools.pretty_timestamp(style=1))
    if not os.path.exists(TESTS_OUTPUT_DIR):
        os.makedirs(TESTS_OUTPUT_DIR)
    return os.path.join(TESTS_OUTPUT_DIR, filename)

def scripts_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(SCRIPTS_DIR, filename)
