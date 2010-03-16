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
Wrappers for interacting with the ETE library.
"""

import tempfile
import re
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)
import dendropy

DENDROPY_ETE_INTEROPERABILITY = False
try:
    import ete2
    DENDROPY_ETE_INTEROPERABILITY = True
except ImportError:
    _LOG.warn("ete2 not installed: ETE interoperability not available")
else:

    def as_ete_object(o):
        if isinstance(o, ete2.Tree):
            return o
        elif isinstance(o, dendropy.Tree) or isinstance(o, dendropy.Node):
            s = o.as_newick_string() + ";"
            _LOG.debug(s)
            return ete2.Tree(s)
        elif isinstance(o, list) or isinstance(o, dendropy.TreeList):
            return [as_ete_object(t) for t in o]
        else:
            raise ValueError("Object of type '%s' has not native ete2 representation" % type(o))

    def as_dendropy_object(o, taxon_set=None):
        if isinstance(o, dendropy.Tree) or isinstance(o, dendropy.Node) or isinstance(o, dendropy.TreeList):
            return o
        elif isinstance(o, ete2.Tree):
            s = o.write()
            _LOG.debug(s)
            return dendropy.Tree.get_from_string(s, 'newick', taxon_set=taxon_set)
        elif isinstance(o, list) or isinstance(o, dendropy.TreeList):
            return dendropy.TreeList([as_dendropy_object(t, taxon_set=taxon_set) for t in o], taxon_set=taxon_set)
        else:
            raise ValueError("Object of type '%s' has not native DendroPy representation" % type(o))

    def show(o):
        ete_o = as_ete_object(o)
        ete_o.show()

