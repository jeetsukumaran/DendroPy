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
Annotations testing.
"""

import sys
import os
import unittest
import tempfile

from dendropy.test.support import pathmap
from dendropy.test.support import nexmlvalidator

class NexmlValidatorTest(unittest.TestCase):

    def testValidatorGoodXml1(self):
        s = pathmap.tree_source_path("pythonidae.annotated.nexml")
        nexmlvalidator.validate_nexml(s)

if __name__ == "__main__":
    unittest.main()


