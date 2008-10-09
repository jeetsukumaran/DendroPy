#! /usr/bin/env python

############################################################################
##  test_char_gen.py
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
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Test character generation
"""

import unittest
import subprocess
import re
import os
import sys

from dendropy import dataio
from dendropy import get_logger
import dendropy.tests
from dendropy.tests import paup

_LOG = get_logger("Splits")

from dendropy import chargen

model_tree_string = """
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=5;
    TAXLABELS
        A
        B
        C
        D
        E
  ;
END;
begin trees;
    tree true=(A:0.25,(B:0.25,(C:0.25,(D:0.25,E:0.25):0.25):0.25):0.25):0.25;
end;
"""

class CharGenTest(unittest.TestCase):
    
    def testCharGen(self):
        source_ds = dataio.get_nexus(string=model_tree_string)
        tree_model = source_ds.trees_blocks[0][0]
        output_ds = chargen.generate_hky_dataset(10000, tree_model=tree_model)
        tb = output_ds.taxa_blocks[0]
        cb = output_ds.char_blocks[0]
        cb_tb = output_ds.char_blocks[0].taxa_block
#         _LOG.info("--      Taxa in Taxa Block: %s" % (" ".join([str(t) for t in tb]))) 
#         _LOG.info("--Taxa in Characters Block: %s" % (" ".join([str(t) for t in cb]))) 
#         _LOG.info("--Taxa in Charac. Block TB: %s" % (" ".join([str(t) for t in cb_tb])))        
#         _LOG.info("\n--Sequences:")
#         for t in cb:
#             _LOG.info("\n%s:      %s" % (str(t), cb[t].values_as_string()))
            
        _LOG.debug(dataio.store_dataset(dataset=output_ds, format='nexus'))
                    
if __name__ == "__main__":
    unittest.main()

# lset nst=6 rmatrix=estimate basefreq=equal rates=equal pinvar=0
# lscore / userbrlen
# Likelihood scores of tree(s) in memory:
#   Likelihood settings:
#     Number of substitution types  = 6
#     Substitution rate-matrix parameters estimated via ML
#     Assumed nucleotide frequencies (set by user):
#       A=0.25000  C=0.25000  G=0.25000  T=0.25000
#     Among-site rate variation:
#       Assumed proportion of invariable sites  = none
#       Distribution of rates at variable sites = equal
#     These settings correspond to a submodel of the GTR model
#     Number of distinct data patterns under this model = 954
#     Molecular clock not enforced
#     Branch lengths constrained to user-input values
#     -ln L (unconstrained) = 58834.70584
# 
# Tree              1
# -------------------
# -ln L   59414.02245
# Rate matrix R:
#   AC        0.95326
#   AG        0.96850
#   AT        0.97641
#   CG        1.00479
#   CT        0.98176
#   GT        1.00000
