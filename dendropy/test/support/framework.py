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
Extensions of the unittest framework.
"""

import unittest

# defining this here means that unittest will exclude all lines from this
# module in the traceback report when an assertion is raised, allowing
# for the point of failure to be the point where the assertion statement
# was made, rather than the point where an exception was raised because
# the assertion as false
__unittest = True

class DendropyTestCase(unittest.TestCase):
    """
    Extends unittest.TestCase with various new assertion tests.
    """

    def assertIsSame(self, obj1, obj2, message=None):
        if message is None:
            message = "%s (%d) is not %s (%d)" % (obj1, id(obj1), obj2, id(obj2))
        self.assertTrue(obj1 is obj2, message)

    def assertIsNotSame(self, obj1, obj2, message=None):
        if message is None:
            message = "%s (%d) is %s (%d)" % (obj1, id(obj1), obj2, id(obj2))
        self.assertTrue(obj1 is not obj2, message)

    def assertIsContainedIn(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is not in: %s" % (obj1, obj2)
        self.assertTrue(obj1 in obj2, message)

    def assertIsNotContainedIn(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is in: %s" % (obj1, obj2)
        self.assertTrue(obj1 in obj2, message)

class DataVerificationTestCase(unittest.TestCase):
    """
    Extends unittest.TestCase with various new assertion tests.
    """
    pass

