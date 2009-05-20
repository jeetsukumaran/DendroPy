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
    
def prune_taxa(tree, taxa):
    """Removes terminal edges associated with taxa in `taxa` from `tree`."""
    for taxon in taxa:
        nd = tree.find_node(lambda x: x.taxon == taxon)
        if nd is not None:
            nd.edge.tail_node.remove_child(nd)
    # clean up dead leaves            
    for nd in tree.postorder_node_iter():
        if nd.taxon is None and len(nd.child_nodes()) == 0:
            dnd = nd
            while dnd.taxon is None and len(dnd.child_nodes()) == 0:
                new_dnd = dnd.parent_node
                new_dnd.remove_child(dnd)
                dnd = new_dnd
    # remove outdegree 1 nodes                
    for nd in tree.postorder_node_iter():
        children = nd.child_nodes()
        if nd.parent_node is not None and len(children) == 1:
            nd.parent_node.add_child(children[0])
            if nd.edge.length is not None:
                if children[0].edge.length is None:
                    children[0].edge.length = nd.edge.length
                else:                    
                    children[0].edge.length += nd.edge.length
            nd.parent_node.remove_child(nd)
    return tree
    
def retain_taxa(tree, taxa):
    """Removes all taxa *not* in `taxa` from the tree."""
    to_prune = [t for t in tree.taxa_block if t not in taxa]
    return prune_taxa(tree, to_prune)
