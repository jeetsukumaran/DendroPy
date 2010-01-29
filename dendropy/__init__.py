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

__version__ = "3.1.0"
PROJECT_NAME = "DendroPy"
PROJECT_VERSION = __version__
PACKAGE_VERSION = PROJECT_VERSION # for backwards compatibility (with sate)
PROJECT_AUTHOR = "Jeet Sukumaran and Mark T. Holder"
PROJECT_COPYRIGHT = "Copyright 2009 Jeet Sukumaran and Mark T. Holder."
PROJECT_LICENSE = """
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

try:
    source_path = os.path.dirname(os.path.abspath(__file__))
    revision = vcsinfo.Revision(repo_path=source_path)
except OSError:
    source_path = None
    revision = vcsinfo.Revision(repo_path=None)

release = __version__
def version_info():
    if revision.is_available:
        revision_text = " (%s)" % (revision)
    else:
        revision_text = ""
    return "%s %s%s" % (PROJECT_NAME, release, revision_text)

if __name__ == "__main__":
    sys.stdout.write("%s\n" % version_info())


