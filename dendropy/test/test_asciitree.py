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
Tests of ascii tree plots.
"""

import unittest
from cStringIO import StringIO
from dendropy.test.support.datagen import reference_tree_list
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

class TreePHGammTest(unittest.TestCase):

    def setUp(self):
        self.tree = reference_tree_list()[0]

    def test_plot_by_depth(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='depth'))

    def test_plot_by_level(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='level'))

    def test_plot_by_age(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='age'))

    def test_plot_by_length(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='length'))

if __name__ == "__main__":
    unittest.main()
