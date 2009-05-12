#! /usr/bin/env python

############################################################################
##  high_level.py
##
##  Part of the DendroPy library for phylogenetic computing.
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
Contains some functions that rely heavily on the rest of dendropy, and perform
    helpful, but specialized tasks.  These functions should not be called by
    the rest of the library, but are may be called in bundled scripts.
"""

def collapse_short_edges_on_trees(tree_iterator, **kwargs):
    print "inside collapse_short_edges_on_trees"
