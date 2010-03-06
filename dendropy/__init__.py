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
Imports into the `dendropy` namespace all fundamental
classes and methods for instantiating objects in the
`dendropy.dataobject` subpackage to for usage by client code.
"""

import sys
import os
from dendropy.utility import vcsinfo
from dendropy.dataobject.base import *
from dendropy.dataobject.taxon import *
from dendropy.dataobject.tree import *
from dendropy.dataobject.char import *
from dendropy.dataobject.dataset import *

from dendropy.dataio import get_reader, get_writer, tree_source_iter, multi_tree_source_iter
#from dendropy.interop import paup
#from dendropy.interop import ape

###############################################################################
## PACKAGE METADATA

__project__ = "DendroPy"
__version__ = "3.1.3"
try:
    __source_path__ = os.path.dirname(os.path.abspath(__file__))
    __revision__ = vcsinfo.Revision(repo_path=__source_path__)
except OSError:
    __source_path__ = None
    __revision__ = vcsinfo.Revision(repo_path=None)
__author__ = "Jeet Sukumaran and Mark T. Holder"
__copyright__ = "Copyright 2009 Jeet Sukumaran and Mark T. Holder."
__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""
PACKAGE_VERSION = __version__ # for backwards compatibility (with sate)

def description():
    if __revision__.is_available:
        revision_text = " (%s)" % str(__revision__)
    else:
        revision_text = ""
    return "%s %s%s" % (__project__, __version__, revision_text)

if __name__ == "__main__":
    sys.stdout.write("%s\n" % description())


