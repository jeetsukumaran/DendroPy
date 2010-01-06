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
Extension to the basic unittest TestCase.
"""

import unittest
from dendropy.utility import messaging

# Defining this here means that unittest will exclude all lines from this
# module in the traceback report when an assertion fails, allowing
# for the starting point of the traceback to be the point where the assertion
# was made, rather than the point where an exception was raised because
# the assertion was false.
__unittest = True

class ExtendedTestCase(unittest.TestCase):
    """
    Extends unittest.TestCase with various new assertion tests.
    """

    def _get_logger(self):
        if not hasattr(self, "_logger") or self._logger is None:
            self._logger = messaging.get_logger(self.__class__.__name__)
        return self._logger

    def _set_logger(self, logger):
        self._logger = logger

    logger = property(_get_logger, _set_logger)

    def assertSame(self, obj1, obj2, message=None):
        if message is None:
            message = "Object %s is not same as object %s: %s vs. %s" % (id(obj1), id(obj2), obj1, obj2)
        self.assertTrue(obj1 is obj2, message)

    def assertNotSame(self, obj1, obj2, message=None):
        if message is None:
            message = "Object %s is same as object %s: %s vs. %s" % (id(obj1), id(obj2), obj1, obj2)
        self.assertTrue(obj1 is not obj2, message)

    def assertContained(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is not in: %s" % (obj1, obj2)
        self.assertTrue(obj1 in obj2, message)

    def assertNotContained(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is in: %s" % (obj1, obj2)
        self.assertTrue(obj1 not in obj2, message)
