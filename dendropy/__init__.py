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
Imports into the `dendropy` namespace all fundamental
classes and methods for instantiating objects in the
`dendropy.dataobject` subpackage to for usage by client code.
"""

import sys
import os

###############################################################################
## Populate the 'dendropy' namespace

from dendropy.utility import messaging
from dendropy.utility import vcsinfo
from dendropy.dataobject.base import *
from dendropy.dataobject.taxon import *
from dendropy.dataobject.tree import *
from dendropy.dataobject.char import *
from dendropy.dataobject.dataset import *

from dendropy.dataio import get_reader, get_writer, tree_source_iter, multi_tree_source_iter

from dendropy.utility import probability
from dendropy import coalescent
from dendropy import continuous
from dendropy import treecalc
from dendropy import popgensim
from dendropy import popgenstat
from dendropy import reconcile
from dendropy import seqmodel
from dendropy import seqsim
from dendropy import treecalc
from dendropy import treemanip
from dendropy import treesim
from dendropy import treesplit
from dendropy import treesum

###############################################################################
## PACKAGE METADATA

__project__ = "DendroPy"
__version__ = "3.10.0"

try:
    try:
        __homedir__ = __path__[0]
    except AttributeError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except IndexError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
except OSError:
    __homedir__ = None
except:
    __homedir__ = None

__revision__ = vcsinfo.Revision(repo_path=__homedir__)
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


