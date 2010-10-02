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
