#! /usr/bin/env python

############################################################################
##  treemanip.py
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
Tree manipulators
"""

import sys
import copy
import math

from dendropy import trees
from dendropy import get_logger
_LOG = get_logger('dendropy.treemanip')


def collapse_edge(edge):
    '''Inserts all children of the head_node of edge as children of the
    tail_node of edge in the same place in the child_node list that head_node
    had occupied.

    The edge length and head_node will no longer be part of the tree.'''
    to_del = edge.head_node
    p = edge.tail_node
    if not p:
        return
    children = to_del.child_nodes()
    if not children:
        raise ValueError('collapse_edge called with a terminal.')
    pos = p.child_nodes().index(to_del)
    p.remove_child(to_del)
    for child in children:
        p.add_child(child, pos=pos)
        pos += 1

def collapse_clade(node):
    """Collapses all internal edges that are descendants of node."""
    if node.is_leaf():
        return
    leaves = [i for i in trees.Node.leaf_iter(node)]
    node.set_children(leaves)
