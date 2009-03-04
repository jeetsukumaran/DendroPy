#! /usr/bin/env python

############################################################################
##  treesum.py
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
Tree summarization.
"""

import sys
from dendropy import splits
from dendropy import taxa
from dendropy import trees
from dendropy import treegen
from dendropy import get_logger
_LOG = get_logger("dendropy.treesum")


def shallowest_containing_node(start_node, split, taxa_mask):
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

class TreeSummarizer(object):
    "Summarizes a distribution of trees."

    def __init__(self):
        "Initializes settings."
        self.support_as_labels = True
        self.support_as_percentages = False
        self.support_label_decimals = 2
        self.verbose=False
        self.write_message = None# sys.stderr.write
        self.progress_message_prefix = None
        self.progress_message_suffix = None
        self.ignore_node_ages = True
        self.total_trees_counted = 0

    def send_progress_message(self, msg):
        "Writes progress message."
        if self.verbose and self.write_message:
            if self.progress_message_prefix:
                prefix = self.progress_message_prefix
            else:
                prefix = ""
            if self.progress_message_suffix:
                suffix = self.progress_message_suffix
            else:
                suffix = ""
            self.write_message("%s%s%s" % (prefix, msg, suffix))

    def compose_support_label(self, split_support_freq):
        "Returns an appropriately composed and formatted support label."
        if self.support_as_percentages:
            if self.support_label_decimals <= 0:
                support_label = "%d" % round(split_support_freq * 100, 0)
            else:
                support_label_template = "%%0.%df" % self.support_label_decimals
                support_label = support_label_template % round(split_support_freq * 100,
                    self.support_label_decimals)
        else:
            if self.support_label_decimals <= 0:
                support_label_decimals = 2
            else:
                support_label_decimals = self.support_label_decimals
            support_label_template = "%%0.%df" % support_label_decimals
            support_label = support_label_template % round(split_support_freq,
                support_label_decimals)
        return support_label

    def map_split_support_to_node(self, node, split_support):
        "Appropriately sets up a node."
        if self.support_as_labels:
            node.label = self.compose_support_label(split_support)
        else:
            if self.support_as_percentages:
                node.edge.length = split_support * 100
            else:
                node.edge.length = split_support
        return node

    def map_split_support_to_tree(self, tree, split_distribution):
        "Maps splits support to the given tree."
        split_frequencies = split_distribution.split_frequencies
        tree.normalize_taxa(taxa_block=split_distribution.taxa_block)
        assert tree.taxa_block is split_distribution.taxa_block
        splits.encode_splits(tree)
        for split in tree.split_edges:
            if split in split_frequencies:
                split_support = split_frequencies[split]
            else:
                split_support = 0.0
            self.map_split_support_to_node(tree.split_edges[split].head_node, split_support)
        return tree

    def tree_from_splits(self,
                         split_distribution,
                         min_freq=0.5,
                         include_edge_lengths=True):
        "Returns a consensus tree from splits in `split_distribution`."
        leaf_to_root_search = True
        
        taxa_block = split_distribution.taxa_block
        con_tree = treegen.star_tree(taxa_block)
        split_freqs = split_distribution.split_frequencies
        taxa_mask = taxa_block.all_taxa_bitmask()
        splits.encode_splits(con_tree)
        leaves = con_tree.leaf_nodes()
        
        if leaf_to_root_search:
            to_leaf_dict = {}
            for leaf in leaves:
                to_leaf_dict[leaf.edge.clade_mask] = leaf
        include_edge_lengths = self.support_as_labels and include_edge_lengths
        unrooted = split_distribution.unrooted
        
        to_try_to_add = []
        for s, f in split_freqs.iteritems():
            if (f > min_freq):
                m = s & taxa_mask
                if (m != taxa_mask) and ((m-1) & m): # if not root (i.e., all "1's") and not singleton (i.e., one "1")
                    if unrooted:
                        c = (~m) & taxa_mask
                        if (c-1) & c: # not singleton (i.e., one "0")
                            if 1 & m:
                                k = c
                            else:
                                k = m
                            to_try_to_add.append((f, k, m))
                    else:
                        to_try_to_add.append((f, m, m))
        to_try_to_add.sort(reverse=True)

        root = con_tree.seed_node
        root_edge = root.edge
        # Now when we add splits in order, we will do a greedy, extended majority-rule consensus tree
        for freq, split_to_add, split_in_dict in to_try_to_add:
            if (split_to_add & root_edge.clade_mask) != split_to_add:
                continue
            elif leaf_to_root_search:
                lb = splits.lowest_bit_only(split_to_add)
                one_leaf = to_leaf_dict[lb]
                parent_node = one_leaf
                while (split_to_add & parent_node.edge.clade_mask) != split_to_add:
                    parent_node = parent_node.parent_node
            else:
                parent_node = shallowest_containing_node(start_node=con_tree.seed_node,
                                                      split=split_to_add,
                                                      taxa_mask=taxa_mask)
            if parent_node is None or parent_node.edge.clade_mask == split_to_add:
                continue # split is not in tree, or already in tree.
            new_node = trees.Node()
            self.map_split_support_to_node(node=new_node, split_support=freq)
            new_node_children = []
            new_edge = new_node.edge
            new_edge.clade_mask = 0
            for child in parent_node.child_nodes():
                # might need to modify the following if rooted splits
                # are used
                cecm = child.edge.clade_mask
                if (cecm & split_to_add ):
                    assert cecm != split_to_add
                    new_edge.clade_mask |= cecm
                    new_node_children.append(child)
            # Check to see if we have accumulated all of the bits that we
            #   needed, but none that we don't need.
            if new_edge.clade_mask == split_to_add:
                if include_edge_lengths:
                    elen = split_distribution.split_edge_lengths[split_in_dict]
                    if len(elen) > 0:
                        new_edge.length = float(sum(elen)) / len(elen)
                    else:
                        new_edge.length = None
                for child in new_node_children:
                    parent_node.remove_child(child)
                    new_node.add_child(child)
                parent_node.add_child(new_node)
                con_tree.split_edges[split_to_add] = new_edge

        ## here we add the support values and/or edge lengths for the terminal taxa ##
        for node in leaves:
            if unrooted:
                split = con_tree.split_edges.normalize_key(node.edge.clade_mask)
            else:
                split = node.edge.clade_mask
            self.map_split_support_to_node(node, 1.0)
            if include_edge_lengths:
                elen = split_distribution.split_edge_lengths.get(split, [0.0])
                if len(elen) > 0:
                    node.edge.length = float(sum(elen)) / len(elen)
                else:
                    node.edge.length = None
        return con_tree

    def count_splits_on_trees(self, tree_iterator, split_distribution=None, trees_splits_encoded=False):
        """
        Given a list of trees file, a SplitsDistribution object (a new one, or,
        if passed as an argument) is returned collating the split data in the files.
        """
        if split_distribution is None:
            split_distribution = splits.SplitDistribution()
        taxa_block = split_distribution.taxa_block
        for tree_idx, tree in enumerate(tree_iterator):
            self.total_trees_counted += 1
            if taxa_block is None:
                assert(split_distribution.taxa_block is None)
                split_distribution.taxa_block = tree.taxa_block
                taxa_block = tree.taxa_block
            else:
                assert(taxa_block is tree.taxa_block)
            if not trees_splits_encoded:
                splits.encode_splits(tree)
            split_distribution.count_splits_on_tree(tree)
        return split_distribution

