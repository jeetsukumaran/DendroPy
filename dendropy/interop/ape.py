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

def as_ape_object(o):
    """
    Returns `o` as an ape object.
    """
    kwargs = {}
    if isinstance(o, dendropy.TreeList):
        kwargs['keep.multi'] = True
        text = o.as_string("newick", spaces_to_underscore=True)
        return _R['read.tree'](text=text, **kwargs)
    elif isinstance(o, dendropy.Tree):
        kwargs['keep.multi'] = False
        text = o.as_string("newick", spaces_to_underscore=True)
        return _R['read.tree'](text=text, **kwargs)
    elif isinstance(o, dendropy.CharacterArray):
        f = tempfile.NamedTemporaryFile()
        o.write_to_stream(f, "nexus", simple=True, spaces_to_underscore=True)
        f.flush()
        return _R['read.nexus.data'](f.name)
    else:
        return robjects.default_py2ri(o)

def as_dendropy_object(o, taxon_set=None):
    """
    Returns a DendroPy object corresponding to the ape object `o`. If `o` is
    a single tree (i.e., `phylo`), then a DendroPy Tree is returned. If `o` is
    a list of trees (i.e., a `multiPhylo` object, or list of `phylo` objects),
    then a DendroPy TreeList is returned.
    """
    if o.rclass[0] == "multiPhylo":
        f = tempfile.NamedTemporaryFile()
        _R['write.nexus'](o, file=f.name)
        return dendropy.TreeList.get_from_path(f.name, "nexus", taxon_set=taxon_set)
    elif o.rclass[0] == "phylo":
        f = tempfile.NamedTemporaryFile()
        _R['write.nexus'](o, file=f.name)
        return dendropy.Tree.get_from_path(f.name, "nexus", taxon_set=taxon_set)
    elif o.rclass[0] == "list":
        f = tempfile.NamedTemporaryFile()
        _R['write.nexus.data'](o, file=f.name)
#        print open(f.name, "r").read()
        d = dendropy.DataSet.get_from_path(f.name, "nexus", taxon_set=taxon_set)
        if len(d.char_arrays) == 0:
            raise ValueError("No character data found")
        elif len(d.char_arrays) == 1:
            return d.char_arrays[0]
        else:
            raise ValueError("Multiple character matrices returned")
    else:
        return robjects.default_ri2py(o)
