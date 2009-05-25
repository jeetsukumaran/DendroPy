#! /usr/bin/env python

############################################################################
##  treestruct.py
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
_LOG = get_logger('dendropy.treestruct')

def mrca(start_node, split, taxa_mask):
    """Returns the shallowest node in the tree (the node furthest from 
    `start_node`) that has all of the taxa that are specified in `split` or
    None if no appropriate node is found.

    Assumes that edges on tree have been decorated with splits.
    
    It is possible that split is not compatible with the subtree that is 
        returned! (compatibility tests are not fully performed).
        
    This function is used to find the "insertion point" for a new split via a
        root to tip search.
    """
    if (start_node.edge.clade_mask & split) != split:
        return None
    curr_node = start_node
    last_match = start_node
    nd_source = iter(start_node.child_nodes())
    try:
        while True:
            cm = curr_node.edge.clade_mask
            cms = (cm & split)
            if cms:
                # for at least one taxon cm has 1 and split has 1
                if cms == split:
                    # curr_node has all of the 1's that split has
                    if cm == split:
                        return curr_node
                    last_match = curr_node
                    nd_source = iter(curr_node.child_nodes())
                else:
                    # we have reached a child that has some, but not all of the
                    #   required taxa as descendants, so we return the last_match
                    return last_match
            curr_node = nd_source.next()
    except StopIteration:
        # we shouldn't reach this if all of the descendants are properly
        #   decorated with clade_mask attributes, but there may be some hacky
        #   context in which we want to allow the function to be called with
        #   leaves that have not been encoded with clade_masks.
        return last_match

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
            while dnd.taxon is None and dnd.parent_node is not None and len(dnd.child_nodes()) == 0:
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
