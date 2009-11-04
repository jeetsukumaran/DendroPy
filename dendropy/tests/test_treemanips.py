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
Tests of tree structural manipulations.
"""

from cStringIO import StringIO
import unittest

import dendropy
from dendropy import treecalc
from dendropy import treemanip
from dendropy import splitcalc
from dendropy.tests import services
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class RandomlyRotateTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        trees = dendropy.TreeList(StringIO(n+n), format="newick")
        ref = trees[0]
        changing = trees[1]
        rng = services.KnownRandom()
        splitcalc.encode_splits(ref)
        splitcalc.encode_splits(changing)
        orig_root = changing.seed_node
        for i in xrange(50):
            treemanip.randomly_rotate(changing, rng=rng)
            self.assertNotEqual(str(changing), n)
            self.assertEqual(orig_root, changing.seed_node)
            changing.debug_check_tree(logger_obj=_LOG, splits=True)
            if treecalc.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

class RandomlyReorientTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        k = dendropy.TreeList(StringIO(n), format="newick")[0]
        trees = dendropy.TreeList(StringIO(n+n), format="newick", encode_splits=True, taxon_set=k.taxon_set)
        ref = trees[0]
        changing = trees[1]
        rng = services.KnownRandom()
        for i in xrange(50):
            treemanip.randomly_reorient_tree(changing, rng=rng, splits=True)
            self.assertNotEqual(str(changing), n)
            changing.debug_check_tree(logger_obj=_LOG, splits=True)
            if treecalc.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

if __name__ == "__main__":
    unittest.main()
