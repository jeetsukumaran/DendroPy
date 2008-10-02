#! /usr/bin/env python

############################################################################
##  test_char_io.py
##
##  Part of the DendroPy phylogenetic computation library.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Tests input/output of characters from files.
"""

import unittest
import datetime
import logging
import tempfile
import os
from optparse import OptionGroup
from optparse import OptionParser
import StringIO

from dendropy import get_logger
from dendropy import get_logging_level

import dendropy.tests
_LOG = get_logger("Character Parsing and Writing")


### MODULES THAT WE ARE TESTING ###
from dendropy import nexus
from dendropy import nexml
### MODULES THAT WE ARE TESTING ###

class CharIOTest(unittest.TestCase):
    pass


import sys
if __name__ == "__main__":
    unittest.main()
