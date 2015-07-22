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
Support for coverage analysis.
"""

import unittest
import shutil
import sys
from optparse import OptionParser
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

DENDROPY_COVERAGE_ANALYSIS_AVAILABLE = False
try:
    from setuptools import Command
except ImportError:
    _LOG.warn("setuptools.Command could not be imported: setuptools extensions not available")
else:
    try:
        import coverage
    except ImportError:
        _LOG.warn("coverage could not be imported: test coverage analysis not available")
    else:
        _LOG.info("coverage imported successfully: test coverage analysis available")
        DENDROPY_COVERAGE_ANALYSIS_AVAILABLE = True

        from dendropy.test.support.dendropytest import get_test_suite
        from dendropy.test.support import pathmap

        class CoverageAnalysis(Command):
            """
            Code coverage analysis command.
            """

            description = "run test coverage analysis"
            user_options = [
                ('erase', None, "remove all existing coverage results"),
                ('branch', 'b', 'measure branch coverage in addition to statement coverage'),
                ('test-module=', 't', "explicitly specify a module to test (e.g. 'dendropy.test.test_containers')"),
                ('no-annotate', None, "do not create annotated source code files"),
                ('no-html', None, "do not create HTML report files"),
            ]

            def initialize_options(self):
                """
                Initialize options to default values.
                """
                self.test_module = None
                self.branch = False
                self.erase = False
                self.no_annotate = False
                self.no_html = False
                self.omit_prefixes = ['dendropy/test']

            def finalize_options(self):
                pass

            def run(self):
                """
                Main command implementation.
                """

                if self.erase:
                    _LOG.warn("removing coverage results directory: '%s'" % pathmap.TESTS_COVERAGE_DIR)
                    try:
                        shutil.rmtree(pathmap.TESTS_COVERAGE_DIR)
                    except:
                        pass
                else:
                    _LOG.info("running coverage analysis ...")
                    if self.test_module is None:
                        test_suite = get_test_suite()
                    else:
                        test_suite = get_test_suite([self.test_module])
                    runner = unittest.TextTestRunner()
                    cov = coverage.coverage(branch=self.branch)
                    cov.start()
                    runner.run(test_suite)
                    cov.stop()
                    if not self.no_annotate:
                        cov.annotate(omit_prefixes=self.omit_prefixes,
                                directory=pathmap.TESTS_COVERAGE_SOURCE_DIR)
                    if not self.no_html:
                        cov.html_report(omit_prefixes=self.omit_prefixes,
                                directory=pathmap.TESTS_COVERAGE_REPORT_DIR)
                    cov.report(omit_prefixes=self.omit_prefixes)

#if __name__ == "__main__":
#    if DENDROPY_COVERAGE_ANALYSIS_AVAILABLE:
#        parser = OptionParser(add_help_option=True)
#        parser.add_option('--erase', dest='erase', action="store_true", default=False, help="remove all existing coverage results")
#        parser.add_option('--branch', '-b', dest='branch', action="store_true", default=False, help='measure branch coverage in addition to statement coverage')
#        parser.add_option('--test-file', '-t', dest='test_module', default=None, help="explicitly specify a module to test (e.g. 'dendropy.test.test_containers')")
#        parser.add_option('--no-annotate', dest='no_annotate', action="store_true", default=False, help="do not create annotated source code files"),
#        parser.add_option('--no-html', dest='no_html', action="store_true", default=False, help="do not create HTML report files"),
#        (opts, args) = parser.parse_args()
#        cov = CoverageAnalysis()
#        cov.erase = opts.erase
#        cov.branch = opts.branch
#        cov.test_module = opts.test_module
#        cov.no_annotate = opts.no_annotate
#        cov.no_html = opt.no_html
#        cov.run()
#    else:
#        sys.stderr.write("Coverage command extension not available: either setuptools or coverage or both could not be imported.\n")
#        sys.exit(1)
