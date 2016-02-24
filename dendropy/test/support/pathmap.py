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
Path mapping for various test resources.
"""

import os
import sys
import tempfile
from dendropy.utility import messaging
from dendropy.utility.timeprocessing import pretty_timestamp
_LOG = messaging.get_logger(__name__)

try:
    import pkg_resources
    # TESTS_DIR = pkg_resources.resource_filename("dendropy", os.path.join(os.pardir, "tests"))
    TESTS_DIR = pkg_resources.resource_filename("dendropy", "test")
    APPLICATIONS_DIR = pkg_resources.resource_filename("dendropy", os.path.join(os.pardir, "applications"))
    _LOG.info("using pkg_resources path mapping")
except:
    LOCAL_DIR = os.path.dirname(__file__)
    TESTS_DIR = os.path.join(LOCAL_DIR, os.path.pardir)
    PACKAGE_DIR = os.path.join(TESTS_DIR, os.path.pardir)
    APPLICATIONS_DIR = os.path.join(PACKAGE_DIR, os.path.pardir, "applications")
    _LOG.info("using local filesystem path mapping")
TESTS_DATA_DIR = os.path.join(TESTS_DIR, "data")
TESTS_OUTPUT_DIR = os.path.join(TESTS_DIR, "output")
TESTS_COVERAGE_DIR = os.path.join(TESTS_DIR, "coverage")
TESTS_COVERAGE_REPORT_DIR = os.path.join(TESTS_COVERAGE_DIR, "report")
TESTS_COVERAGE_SOURCE_DIR = os.path.join(TESTS_COVERAGE_DIR, "source")

def tree_source_stream(filename):
    if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
        return open(tree_source_path(filename), "rU")
    else:
        return open(tree_source_path(filename), "r")

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

def splits_source_stream(filename):
    if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
        return open(splits_source_path(filename), "rU")
    else:
        return open(splits_source_path(filename), "r")

def splits_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(TESTS_DATA_DIR, "splits", filename)

def other_source_stream(filename):
    if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
        return open(other_source_path(filename), "rU")
    else:
        return open(other_source_path(filename), "r")

def other_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(TESTS_DATA_DIR, "other", filename)

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
            filename = "%s.%s" % (filename, pretty_timestamp(style=1))
    if not os.path.exists(TESTS_OUTPUT_DIR):
        os.makedirs(TESTS_OUTPUT_DIR)
    return os.path.join(TESTS_OUTPUT_DIR, filename)

def application_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(APPLICATIONS_DIR, filename)

class SandboxedFile(object):

    def __init__(self, mode="w"):
        self.mode = mode
        self.fileobj = None
        self.filepath = None

    def __enter__(self):
        self.fileobj = tempfile.NamedTemporaryFile(
                mode=self.mode,
                delete=False)
        self.filepath = self.fileobj.name
        return self.fileobj

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            self.fileobj.flush()
            self.fileobj.close()
        except ValueError:
            # If client code closes:
            #   ValueError: I/O operation on closed file.
            pass
        try:
            os.remove(self.filepath)
        except OSError:
            pass
