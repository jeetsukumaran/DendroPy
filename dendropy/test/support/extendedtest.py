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

    def assertIs(self, obj1, obj2, message=None):
        if message is None:
            message = "Object %s is not same as object %s: %s vs. %s" % (id(obj1), id(obj2), obj1, obj2)
        self.assertTrue(obj1 is obj2, message)

    def assertIsNot(self, obj1, obj2, message=None):
        if message is None:
            message = "Object %s is same as object %s: %s vs. %s" % (id(obj1), id(obj2), obj1, obj2)
        self.assertTrue(obj1 is not obj2, message)

    def assertIn(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is not in: %s" % (obj1, obj2)
        self.assertTrue(obj1 in obj2, message)

    def assertNotIn(self, obj1, obj2, message=None):
        if message is None:
            message = "%s is in: %s" % (obj1, obj2)
        self.assertTrue(obj1 not in obj2, message)
