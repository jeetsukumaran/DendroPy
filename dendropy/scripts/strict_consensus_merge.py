#!/usr/bin/env python
import sys
import copy
import logging

from dendropy import dataio
from dendropy.splits import encode_splits, split_to_list, count_bits, lowest_bit_only
from dendropy import get_logger
from dendropy.treemanip import collapse_clade, collapse_edge
from dendropy.dataio import trees_from_newick
from dendropy.trees import debug_check_tree, format_split

_LOG = get_logger('scripts.strict_consensus_merge')
verbose = False
IS_DEBUG_LOGGING = _LOG.isEnabledFor(logging.DEBUG)

from dendropy.utils import NormalizedBitmaskDict
g_taxa_block = None # temp hack
def reroot_on_lowest_common_index_path(t, common_mask):
    """This operation is only for unrooted trees that are being merged using
    SCM. The path the separates the lowest index taxon in the leaf set 
    intersection is placed as the first descendant path of the "seed_node" for
    the tree.
    This assures that all split representations are oriented in the same way 
    for subsequent operations.
    The mask most contain more that 2 bits (there must be an internal node in
    the tree that is has degree > 2  wrt the common leafset).
    """
    l = lowest_bit_only(common_mask)
    assert(l > 0)
    assert(count_bits(common_mask) > 2)
    # start at the lowest leaf in common.
    curr_n = t.split_edges[l].head_node
    # walk back toward the root until we find a node that has another bit
    p = curr_n.parent_node
    while p:
        if (p.edge.clade_mask & common_mask) != l:
            break
        curr_n = p
        p = curr_n.parent_node

    without_lowest = common_mask^l

    taxa_mask = t.seed_node.edge.clade_mask
    if (curr_n.edge.clade_mask & common_mask) == l:
        # we did not make it to the root.  Make curr_n, the first_child of the root
        t.to_outgroup_position(curr_n, flip_splits=True, suppress_deg_two=True)
        avoid = curr_n
        nd_source = iter(t.seed_node.child_nodes())
        
        try:
            while True:
                curr_n = nd_source.next()
                if curr_n is not avoid:
                    cm = (curr_n.edge.clade_mask & without_lowest)
                    if cm:
                        if cm == without_lowest:
                            r = t.seed_node
                            assert curr_n.parent_node is r
                            t.reroot_at(curr_n, flip_splits=True, suppress_deg_two=True)
                            t.to_outgroup_position(r, flip_splits=True, suppress_deg_two=True)
                            nd_source = iter(curr_n.child_nodes())
                            avoid = r
                        else:
                            return
        except StopIteration:
            assert False
            return
    # we hit the root, now we walk up the tree, to find the a relevant internal
    lowest_on_path_to_l = curr_n
    comp_mask = (~common_mask) & taxa_mask
    children = curr_n.child_nodes()
    assert(len(children) > 1)
    nd_source = iter(children)
    try:
        while True:
            c = nd_source.next()
            cm = c.edge.clade_mask
            masked_cm = cm & common_mask
            if masked_cm:
                if masked_cm == without_lowest:
                    curr_n = c
                    children = curr_n.child_nodes()
                    assert(len(children) > 1)
                    nd_source = iter(children)
                else:
                    break
    except StopIteration:
        raise AssertionError("Reaching supposedly unreachable code")
    
    if curr_n is not t.seed_node:
        # We have found the first relevant internal node, we want to make it 
        #   the root.  We can do this by putting one of its children into the
        #   "outgroup position" and then putting the path to lowest commond
        #   leaf in the outgroup position (this last operation is just a 
        #   rearrangement of the order of children in the root.
        children = curr_n.child_nodes()
        assert(len(children) > 1)
        p = curr_n.parent
        t.to_outgroup_position(children[0], flip_splits=True, suppress_deg_two=True)
        t.to_outgroup_position(p, flip_splits=True, suppress_deg_two=True)
    else:
        # if the root first relevant, node then we just make the path leading
        #   to the lowest index node the first child of the root
        t.to_outgroup_position(lowest_on_path_to_l, flip_splits=True, suppress_deg_two=True)

def _collapse_paths_not_found(f, s, other_dict=None):
    to_del = []
    for masked_split, path in f.iteritems():
        if masked_split not in s:
            for edge in path:
                if other_dict:
                    del other_dict[edge.clade_mask]
                collapse_edge(edge)
            to_del.append(masked_split)
    for k in to_del:
        del f[k]


def add_to_scm(to_modify, to_consume, rooted=False, taxa_block=None):
    """Adds the tree `to_consume` to the tree `to_modify` in a strict consensus
    merge operation.  Both trees must have had encode_splits called on them."""
    if rooted:
        raise NotImplementedError("rooted form of add_to_scm not implemented")
    to_mod_root = to_modify.seed_node
    to_mod_split = to_mod_root.edge.clade_mask

    to_consume_root = to_consume.seed_node
    to_consume_split = to_consume_root.edge.clade_mask

    leaf_intersection = to_mod_split & to_consume_split
    if IS_DEBUG_LOGGING:
        _LOG.debug("add_to_scm:\n  %s\n  + %s\n%s" % (str(to_modify), str(to_consume), format_split(leaf_intersection, taxa=taxa_block)))

    n_common_leaves = count_bits(leaf_intersection)
    if n_common_leaves < 2:
        _LOG.error('trees must have at least 2 common leaves')
        raise ValueError('trees must have at least 2 common leaves')
    if n_common_leaves == 2: 
        # SCM with 2 leaves in common results in a polytomy
        collapse_clade(to_mod_root)
        collapse_clade(to_consume_root)
        leaves_to_steal = [c for c in to_consume_root.child_nodes() if not (leaf_intersection & c.edge.clade_mask)]
        for leaf in leaves_to_steal:
            to_mod_root.add_child(leaf)
            to_mod_root.edge.clade_mask |= leaf.edge.clade_mask
        to_modify.split_edges = {to_mod_root.edge.clade_mask : to_mod_root.edge}
        for child in to_mod_root.child_nodes():
            to_modify.split_edges[child.edge.clade_mask] = child.edge
        return
    
    # at least 3 leaves in common
    tmse = to_modify.split_edges
    
    to_mod_relevant_splits = {}
    to_consume_relevant_splits = {}
    if not rooted:
        if IS_DEBUG_LOGGING:
            debug_check_tree(to_modify, splits=True, logger_obj=_LOG)
            debug_check_tree(to_consume, splits=True, logger_obj=_LOG)

        reroot_on_lowest_common_index_path(to_modify, leaf_intersection)
        reroot_on_lowest_common_index_path(to_consume, leaf_intersection)

        if IS_DEBUG_LOGGING:
            debug_check_tree(to_modify, splits=True, logger_obj=_LOG)
            debug_check_tree(to_consume, splits=True, logger_obj=_LOG)

        to_mod_root = to_modify.seed_node
        assert(to_mod_root.edge.clade_mask == to_mod_split)
        to_consume_root = to_consume.seed_node
        assert(to_consume_root.edge.clade_mask == to_consume_split)

    for s, e in tmse.iteritems():
        s = e.clade_mask
        masked = s & leaf_intersection
        if masked and masked != leaf_intersection:
            e_list = to_mod_relevant_splits.setdefault(masked, [])
            e_list.append((s, e))

    for s, e in to_consume.split_edges.iteritems():
        s = e.clade_mask
        masked = s & leaf_intersection
        if masked and masked != leaf_intersection:
            e_list = to_consume_relevant_splits.setdefault(masked, [])
            e_list.append((s, e))

    # Because each of these paths radiates away from the root (none of the paths
    #   cross the root), the clade_masks for deeper edges will be supersets
    #   of the clade_masks for shallower nodes.  Thus if we reverse sort we
    #   get the edges in the order root->tip
    for split, path in to_mod_relevant_splits.iteritems():
        path.sort(reverse=True)
        t = [i[1] for i in path]
        del path[:]
        path.extend(t)
    for split, path in to_consume_relevant_splits.iteritems():
        path.sort(reverse=True)
        t = [i[1] for i in path]
        del path[:]
        path.extend(t)
    if IS_DEBUG_LOGGING:
        debug_check_tree(to_modify, splits=True, logger_obj=_LOG)
        debug_check_tree(to_consume, splits=True, logger_obj=_LOG)
    # first we'll collapse all paths in the common leafset in to_modify that 
    #   are not in to_consume
    _collapse_paths_not_found(to_mod_relevant_splits, to_consume_relevant_splits, tmse)
    # Now we'll collapse all paths in the common leafset in to_consume that 
    #   are not in to_modify
    _collapse_paths_not_found(to_consume_relevant_splits, to_mod_relevant_splits)


    # first we'll deal with subtrees that are:
    #       - not in the leaf intersection set, and
    #       - attached to "relevant" nodes
    # We simply move these subtrees from the to_consume tree to the appropriate
    #   node in to_modify
    to_steal = [i for i in to_consume_root.child_nodes() if (i.edge.clade_mask & leaf_intersection) == 0]
    for child in to_steal:
        to_mod_root.add_child(child)
        to_mod_root.edge.clade_mask |= child.edge.clade_mask
        
    for masked_split, to_consume_path in to_consume_relevant_splits.iteritems():
        to_mod_path = to_mod_relevant_splits.get(masked_split)
        if True: #to_mod_path is None:
            if IS_DEBUG_LOGGING:
                _LOG.debug("%s = mask" % format_split(leaf_intersection, taxa=taxa_block))
                _LOG.debug("%s = masked" % format_split(masked_split, taxa=taxa_block))
                _LOG.debug("%s = raw" % format_split(to_consume_path[-1].clade_mask, taxa=taxa_block))
                for k, v in to_mod_relevant_splits.iteritems():
                    _LOG.debug("%s in to_mod_relevant_splits" % format_split(k, taxa=taxa_block))
                
            assert to_mod_path is not None
        to_mod_head = to_mod_path[-1].head_node
        to_mod_head_edge = to_mod_head.edge
        to_consume_head = to_consume_path[-1].head_node
        for child in to_consume_head.child_nodes():
            if (child.edge.clade_mask & leaf_intersection) == 0:
                # child is the root of a subtree that has no children in the leaf_intersection
                to_mod_head.add_child(child)
                to_mod_head_edge.clade_mask |= child.edge.clade_mask
        if len(to_consume_path) > 1:
            if len(to_mod_path) > 1:
                # collision
                for edge in to_mod_path[1:-1]:
                    collapse_edge(edge)
                mid_node = to_mod_path[0].head_node
                for edge in to_consume_path[1:]:
                    p = edge.tail_node
                    avoid = edge.head_node
                    for child in p.child_nodes():
                        if child is not avoid:
                            mid_node.add_child(child)
                            mid_node.clade_mask |= child.edge.clade_mask
            else:
                # we have to move the subtrees from to_consume to to_modify
                to_mod_edge = to_mod_path[0]
                to_mod_tail, to_mod_head = to_mod_edge.tail_node, to_mod_edge.head_node
                deepest_edge_to_move = to_consume_path[0]
                deepest_node_to_move = deepest_edge_to_move.head_node
                tipmost_edge_to_move = to_consume_path[-1]
                tipmost_node_to_move = tipmost_edge_to_move.tail_node
                prev_head = tipmost_edge_to_move.head_node
                
                to_mod_tail.add_child(deepest_node_to_move)
                to_mod_tail.remove_child(to_mod_head)
                tipmost_node_to_move.add_child(to_mod_head)
                tipmost_node_to_move.remove_child(prev_head)
    encode_splits(to_modify, taxa_block=taxa_block)
                
    
def strict_consensus_merge(trees_to_merge, taxa_block, copy_trees=False, rooted=False):
    """Returns a tree that is the strict consensus merger of the input trees.
    
    If copy_trees is True then the trees will be copied before the merger 
    operation, if the `copy_trees` is False then the input trees will be 
    destroyed by the operation (and the modified first tree will be returned).
    """
    global g_taxa_block
    g_taxa_block = taxa_block
    if copy_trees:
        tree_list = [copy.copy(i) for i in trees_to_merge]
    else:
        tree_list = list(trees_to_merge)
        del trees_to_merge[1:]
    nTrees = len(tree_list)
    _LOG.debug('%d Trees to merge:\n%s\n' % (nTrees, '\n'.join([str(i) for i in tree_list])))
    if nTrees < 2:
        return tree_list[0]
    tree_iter = iter(tree_list)
    to_modify = tree_iter.next()
    to_modify.taxa_block = taxa_block
    
    if rooted:
        raise NotImplementedError("Rooted SCM is not implemented")
    else:
        to_modify.deroot()
    encode_splits(to_modify, taxa_block=taxa_block)
    if IS_DEBUG_LOGGING:
        assert to_modify._debug_tree_is_valid(splits=False)
    for to_consume in tree_iter:
        to_consume.taxa_block = taxa_block
        if not rooted:
            to_consume.deroot()
        encode_splits(to_consume, taxa_block=taxa_block)
        if IS_DEBUG_LOGGING:
            assert to_consume._debug_tree_is_valid(splits=True)
        add_to_scm(to_modify, to_consume, rooted, taxa_block=taxa_block)
        if IS_DEBUG_LOGGING:
            assert to_modify._debug_tree_is_valid(splits=False)

    return to_modify


if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) > 0:
        newick = args
    else:
        newick = sys.stdin.readlines()
    dataset = trees_from_newick(newick)
    trees = [i[0] for i in dataset.trees_blocks]
    o = strict_consensus_merge(trees, taxa_block=dataset.taxa_blocks[0])
    sys.stdout.write("%s;\n" % str(o))


