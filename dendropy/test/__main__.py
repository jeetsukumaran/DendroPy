#! /usr/bin/env python

###############################################################################
##
##  Copyright 2012 Jeet Sukumaran.
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

import sys
import os
import argparse
import collections
from dendropy.utility import messaging
from dendropy import test

def main():
    names = (
        (None        , ".*"),
        ("all"       , ".*"),
        ("datamodel" , ".*_datamodel_.*"),
        ("dataio"    , ".*_dataio_.*"),
        ("newick"    , ".*_newick_.*"),
        ("tree"      , ".*_tree_.*"),
        )
    test_group_patterns = collections.OrderedDict(names)
    parser = argparse.ArgumentParser()
    parser.add_argument("test_groups",
            metavar="GROUP",
            action="append",
            nargs="?",
            choices=list(test_group_patterns.keys()),
            help="Groups of tests to run. May be specified multiple times. Defaults to 'all'.")
    parser.add_argument("--list-only",
            action="store_true",
            default=False,
            help="Do not actually run tests: just print list of test module names and exit.")
    parser.add_argument("-v", "--verbosity",
            default=3,
            type=int,
            help="Messaging noisiness (default: %(default)s)")
    parser.add_argument("--logging-level",
            default=os.environ.get(messaging._LOGGING_LEVEL_ENVAR, "NOTSET"),
            choices=["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            help="Test logging level (default: '%(default)s')")
    args = parser.parse_args()

    # Set logging level:
    os.environ[messaging._LOGGING_LEVEL_ENVAR] = args.logging_level
    _LOG = messaging.get_logger("dendropy")

    # get test modules
    filter_patterns = [test_group_patterns[n] for n in args.test_groups]
    test_modules = test.get_test_module_names(filter_patterns)
    test_modules = sorted(set(test_modules))

    # 0: nothing
    # 1: errors and mishaps only + 0
    # 2: warnings + 1
    # 3: general messages + 2
    if args.verbosity >= 3 or args.list_only:
        sys.stderr.write("DendroPy test modules to be run:\n")
        for mp in test_modules:
            sys.stderr.write(" + {}\n".format(mp))

    if args.list_only:
        sys.exit(0)

    test_suite = test.get_test_suite(test_modules)
    test.run(test_suite, args.verbosity)

if __name__ == '__main__':
    main()


