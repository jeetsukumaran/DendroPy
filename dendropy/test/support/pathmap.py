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
Path mapping for various test resources.
"""

import os
from dendropy.utility import textutils
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

try:
    import pkg_resources
    TESTS_DIR = pkg_resources.resource_filename("dendropy", "test")
    SCRIPTS_DIR = pkg_resources.resource_filename("dendropy", os.path.join(os.pardir, "scripts"))
    _LOG.info("using pkg_resources path mapping")
except:
    LOCAL_DIR = os.path.dirname(__file__)
    TESTS_DIR = os.path.join(LOCAL_DIR, os.path.pardir)
    PACKAGE_DIR = os.path.join(TESTS_DIR, os.path.pardir)
    SCRIPTS_DIR = os.path.join(PACKAGE_DIR, os.path.pardir, "scripts")
    _LOG.info("using local filesystem path mapping")


TESTS_DATA_DIR = os.path.join(TESTS_DIR, "data")
TESTS_OUTPUT_DIR = os.path.join(TESTS_DIR, "output")
TESTS_COVERAGE_DIR = os.path.join(TESTS_DIR, "coverage")
TESTS_COVERAGE_REPORT_DIR = os.path.join(TESTS_COVERAGE_DIR, "report")
TESTS_COVERAGE_SOURCE_DIR = os.path.join(TESTS_COVERAGE_DIR, "source")

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
            filename = "%s.%s" % (filename, textutils.pretty_timestamp(style=1))
    if not os.path.exists(TESTS_OUTPUT_DIR):
        os.makedirs(TESTS_OUTPUT_DIR)
    return os.path.join(TESTS_OUTPUT_DIR, filename)

def script_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(SCRIPTS_DIR, filename)
