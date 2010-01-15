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
Support for coverage analysis.
"""

import unittest
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

try:
    from setuptools import Command
    DENDROPY_COVERAGE_ANALYSIS_AVAILABLE = False
    try:
        import coverage
        _LOG.info("coverage imported successfully: test coverage analysis available")
        DENDROPY_COVERAGE_ANALYSIS_AVAILABLE = True

        from dendropy.test import get_test_file_names, get_test_suite
        from dendropy.test.support import pathmap

        class CoverageAnalysis(Command):
            """
            Code coverage analysis command.
            """

            description = "run test coverage analysis"
            user_options = [
                ('test-file=', 't', "explicitly specify a module to test (e.g. 'dendropy.test.test_containers')"),
#                ('annotate', 'a', 'annotate source files with execution information'),
#                ('branch', 'b', 'meassure branch coverage in addition to statement coverage'),
#                ('directory=', 'd', 'write output files to DIR')
            ]

            def initialize_options(self):
                """init options"""
                self.test_file = None
                self.annotate = False
                self.branch = False
                self.directory = None
                self.omit_prefixes = ['dendropy/test']

            def finalize_options(self):
                """finalize options"""
                pass

            def run(self):
                """runner"""
                self.test_file = 'dendropy.test.test_containers'
                if self.test_file is None:
                    test_suite = get_test_suite()
                else:
                    test_suite = get_test_suite([self.test_file])
                runner = unittest.TextTestRunner()
                cov = coverage.coverage()
                cov.start()
                runner.run(test_suite)
                cov.stop()
                cov.annotate(omit_prefixes=self.omit_prefixes,
                        directory=pathmap.TESTS_COVERAGE_SOURCE_DIR)
                cov.report(omit_prefixes=self.omit_prefixes)

    except ImportError:
        _LOG.warn("coverage could not be imported: test coverage analysis not available")
        DENDROPY_COVERAGE_ANALYSIS_AVAILABLE = False

except ImportError:
    _LOG.warn("setuptools.Command could not be imported: setuptools extensions not available")
