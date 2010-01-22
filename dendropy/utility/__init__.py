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
