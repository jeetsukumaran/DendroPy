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

import sys
import os
import argparse
import collections
import unittest
from dendropy.utility import metavar
from dendropy.utility import messaging
sys.path.insert(0, os.path.dirname(__file__))
from dendropy.test.support import dendropytest

def main():
    group_names = (
        ("@all"       , ".*"),
        ("@datamodel" , ".*_datamodel_.*"),
        ("@dataio"    , ".*_dataio_.*"),
        ("@newick"    , ".*_newick_.*"),
        ("@tree"      , ".*_tree_.*"),
        )
    test_group_patterns = collections.OrderedDict(group_names)
    test_group_names = list(test_group_patterns)
    parser = argparse.ArgumentParser()
    parser.add_argument("test_names",
            metavar="TEST",
            nargs="*",
            help= "Name of test(s) to run. These can be (dot-)qualified module, test"
            "case, or test name (e.g., 'test_module', 'test_module.TestCase1',"
            "'test_module.TestCase1.test1') or special pre-defined groups of"
            "tests (e.g., '@datamodel', '@dataio'). Type '--help-testgroups' for"
            "a list of available groups.")
    parser.add_argument("--help-testgroups",
            action="store_true",
            default=False,
            help="Show list of available test groups and exit.")
    parser.add_argument("--list-only",
            action="store_true",
            default=False,
            help="Do not actually run tests: just print list of test module names and exit.")
    parser.add_argument("-v", "--verbosity",
            default=3,
            type=int,
            help="Messaging noisiness (default: %(default)s)")
    parser.add_argument("--logging-level",
            default=os.environ.get(metavar.LOGGING_LEVEL_ENVAR, "NOTSET"),
            choices=["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            help="Test logging level (default: '%(default)s')")
    parser.add_argument("-f", "--fail-fast",
            action="store_true",
            default=False,
            help="Stop the test run on the first error or failure.")
    parser.add_argument("-I", "--fail-incomplete",
            action="store_true",
            default=False,
            help="Fail incomplete or partially-complete test stubs.")
    args = parser.parse_args()

    if args.help_testgroups:
        out = sys.stdout
        out.write("Available special test groups:\n")
        for name in test_group_names:
            out.write("  - {}\n".format(name))
        sys.exit(0)

    # Set logging level:
    os.environ[metavar.LOGGING_LEVEL_ENVAR] = args.logging_level
    _LOG = messaging.get_logger("dendropy")

    # Set test specifications
    if args.fail_incomplete:
        os.environ[metavar.FAIL_INCOMPLETE_TESTS_ENVAR] = "1"

    # get test modules
    test_names = []
    filter_patterns = []
    for name in args.test_names:
        if name is None:
            continue
        if name.startswith("@"):
            try:
                filter_patterns.append(test_group_patterns[name])
            except KeyError:
                sys.exit("Unrecognized test group name '{}'. Accepted names: {}".format(name, test_group_names))
        else:
            name = name.replace(os.sep, ".")
            if name.endswith(".py"):
                name = name[:-3]
            if not name.startswith("dendropy.test."):
                if name.startswith("test."):
                    name = "dendropy." + name
                else:
                    name = "dendropy.test." + name
            test_names.append(name)

    if not test_names and not filter_patterns:
        test_names = dendropytest.discover_test_module_paths() # get all
    if filter_patterns:
        test_names.extend(dendropytest.discover_test_module_paths(filter_patterns))
    test_names = sorted(set(test_names))

    # 0: nothing
    # 1: errors and mishaps only + 0
    # 2: warnings + 1
    # 3: general messages + 2
    if args.verbosity >= 3 or args.list_only:
        if args.list_only:
            out = sys.stdout
        else:
            out = sys.stderr
        out.write("DendroPy tests to be run:\n")
        for mp in test_names:
            out.write(" + {}\n".format(mp))

    if args.list_only:
        sys.exit(0)

    tests = unittest.defaultTestLoader.loadTestsFromNames(test_names)
    test_suite = unittest.TestSuite(tests)
    test_runner = unittest.TextTestRunner(verbosity=args.verbosity, failfast=args.fail_fast)
    test_runner.run(test_suite)

if __name__ == '__main__':
    main()


