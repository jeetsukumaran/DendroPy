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
Low-level data structures, classes, constants, etc. that are used
by the DendroPy library but generally not directly by client code.
"""

import os
import random

###############################################################################
## USER-SPECIFIC

# global debugging flag
if "DENDROPY_DEBUG" in os.environ:
    if os.environ["DENDROPY_DEBUG"] \
        and os.environ["DENDROPY_DEBUG"].lower()[0] in ["1", "t", "y", "d"]:
        GLOBAL_DEBUG = True
    else:
        GLOBAL_DEBUG = False
else:
    GLOBAL_DEBUG = False

_user_ini_checked = False
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.dendropy/startup.py")
    if os.path.exists(p):
        execfile(p)
    del p

###############################################################################
## GLOBAL RANDOM NUMBER GENERATOR

GLOBAL_RNG = random.Random()
