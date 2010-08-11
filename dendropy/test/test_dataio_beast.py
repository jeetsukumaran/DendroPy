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
BEAST data read/write parse/format tests.
"""

import sys
import os
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.test.support import extendedtest
from dendropy import dataio
from dendropy.utility import error
from dendropy.utility.messaging import get_logger
import dendropy

_LOG = get_logger(__name__)

class BeastSummaryTreeTests(datatest.DataObjectVerificationTestCase):

    def testReadSummaryTree(self):
        tree = dendropy.Tree.get_from_stream(pathmap.tree_source_stream('pythonidae.beast-summary.tre'), 'beast-summary-tree')
        #for nd in tree:
        #    print('---')
        #    print('Node: %s' % nd.oid)
        #    #print nd.comments[0]
        #    print('height          = %s' % nd.height)
        #    print('height_median   = %s' % nd.height_median)
        #    print('height_95hpd    = %s' % nd.height_95hpd)
        #    print('height_range    = %s' % nd.height_range)
        #    print('length          = %s' % nd.length)
        #    print('length_median   = %s' % nd.length_median)
        #    print('length_95hpd    = %s' % nd.length_95hpd)
        #    print('length_range    = %s' % nd.length_range)
        #    print('posterior       = %s' % nd.posterior)

    def testSummaryTreeIter(self):
        stream = pathmap.tree_source_stream('pythonidae.beast-summary.tre')
        for tree in dendropy.tree_source_iter(stream, 'beast-summary-tree'):
            pass
            #for nd in tree:
                #print('---')
                #print('Node: %s' % nd.oid)
                ##print nd.comments[0]
                #print('height          = %s' % nd.height)
                #print('height_median   = %s' % nd.height_median)
                #print('height_95hpd    = %s' % nd.height_95hpd)
                #print('height_range    = %s' % nd.height_range)
                #print('length          = %s' % nd.length)
                #print('length_median   = %s' % nd.length_median)
                #print('length_95hpd    = %s' % nd.length_95hpd)
                #print('length_range    = %s' % nd.length_range)
                #print('posterior       = %s' % nd.posterior)

if __name__ == "__main__":
    unittest.main()

