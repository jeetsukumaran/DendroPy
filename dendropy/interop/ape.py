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
Wrappers for interacting with the APE library for R.
"""

import tempfile
from rpy2 import robjects
import dendropy

_R = robjects.r
_R('library(ape)')

def as_ape_tree(t):
    """
    Returns `t` as an ape object. If `t` is a TreeList or list of Trees,
    then the ape object is `multiPhylo` list of `phylo` objects.
    If `t` is a Tree, then the ape object is a `phylo` object.
    """
    kwargs = {}
    if isinstance(t, list):
        kwargs['keep.multi'] = True
        if isinstance(t, TreeList):
            text = t.as_string("newick", spaces_to_underscore=True)
        else:
            text = ";\n".join([i.as_string("newick", spaces_to_underscore=True) for i in t])
    else:
        kwargs['keep.multi'] = False
        text = t.as_string(format="newick", spaces_to_underscore=True)
    t = _R['read.tree'](text=text, **kwargs)
    return t

def as_dendropy_tree(t, taxon_set=None):
    """
    Returns a DendroPy object corresponding to the ape object `t`. If `t` is
    a single tree (i.e., `phylo`), then a DendroPy Tree is returned. If `t` is
    a list of trees (i.e., a `multiPhylo` object, or list of `phylo` objects),
    then a DendroPy TreeList is returned.
    """
    f = tempfile.NamedTemporaryFile()
    _R['write.nexus'](t, file=f.name)
    if t.rclass[0] == "multiPhylo":
        return dendropy.TreeList.get_from_path(f.name, "nexus", taxon_set=taxon_set)
    else:
        return dendropy.Tree.get_from_path(f.name, "nexus", taxon_set=taxon_set)
