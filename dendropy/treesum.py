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


def deepest_containing_node(start_node, split, taxa_mask):
    """
    Assumes:
        - start_node is a node on a tree
        - edges on tree have been decorated with splits
        - start_node is compatible with split (root if all else fails)
    Returns:
        - deepest node in tree (i.e., furthest from root) that is
          compatible with split (i.e, includes all taxa in split) yet
          not equal to the split.
    """
    split = split & taxa_mask
    curr_node = start_node
    if (start_node.edge.clade_mask & split) != split:
        return None
    last_match = start_node
    nd_source = iter(start_node.child_nodes())
    try:
        while True:
            cm = curr_node.edge.clade_mask & taxa_mask
            cms = (cm & split)
            if cms:
                # for at least one taxon cm has 1 and split has 1
                if cms == split:
                    # curr_node has all of the 1's that split has
                    if cm == split:
                        return last_match
                    last_match = curr_node
                    nd_source = iter(curr_node.child_nodes())
                else:
                    return last_match
            curr_node = nd_source.next()
    except StopIteration:
        return last_match

class TreeSummarizer(object):
    "Summarizes a distribution of trees."

    def __init__(self):
        "Initializes settings."
        self.burnin = None
        self.support_as_labels = True
        self.support_as_percentages = False
        self.support_label_decimals = 2
        self.verbose=False
        self.write_message = None# sys.stderr.write
        self.progress_message_prefix = None
        self.progress_message_suffix = None
        self.ignore_node_ages = True
        self.total_trees_read = 0
        self.total_trees_counted = 0
        self.total_trees_ignored = 0

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
        splits.encode_splits(tree, taxa_block=split_distribution.taxa_block)
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
        taxa_block = split_distribution.taxa_block
        con_tree = treegen.star_tree(taxa_block)
        split_freqs = split_distribution.split_frequencies
        taxa_mask = taxa_block.all_taxa_bitmask()
        splits.encode_splits(con_tree, taxa_block)
        #start_leaf = con_tree.split_edges[0x01].head_node
        include_edge_lengths = self.support_as_labels and include_edge_lengths
        unrooted = split_distribution.unrooted
        to_try_to_add = []
        for s, f in split_freqs.iteritems():
            if (f > min_freq):
                m = s & taxa_mask
                if (m != taxa_mask) and ((m-1) & m): # if not root (i.e., all "1's") and not singleton (i.e., one "1")
                    if unrooted:
                        c = m ^ taxa_mask
                        if (c-1) & c: # not singleton (i.e., one "0")
                            if 1 & m:
                                k = m ^ taxa_mask
                            to_try_to_add.append((f, k, m))
                    else:
                        to_try_to_add.append((f, m, m))
        to_try_to_add.sort(reverse=True)

        # Now when we add splits in order, we will do a greedy, extended majority-rule consensus tree
        for freq, split_to_add, split_in_dict in to_try_to_add:
            parent_node = deepest_containing_node(start_node=con_tree.seed_node,
                                                  split=split_to_add,
                                                  taxa_mask=taxa_mask)
            if parent_node is None or parent_node.edge.clade_mask == split_to_add:
                #or (split_distribution.unrooted and parent_node.edge.clade_mask == (split ^ taxa_mask)):
                pass
            else:
                new_node = trees.Node()
                self.map_split_support_to_node(node=new_node, split_support=freq)
                if include_edge_lengths:
                    elen = split_distribution.split_edge_lengths[split_in_dict]
                    new_node.edge.length = float(sum(elen)) / len(elen)
                new_node_children = []
                new_edge = new_node.edge
                for child in parent_node.child_nodes():
                    # might need to modify the following if rooted splits
                    # are used
                    if (child.edge.clade_mask & split_to_add ):
                        assert child.edge.clade_mask != split_to_add
                        new_node_children.append(child)
                new_edge.clade_mask = 0
                for child in new_node_children:
                    parent_node.remove_child(child)
                    new_node.add_child(child)
                    new_edge.clade_mask |= child.edge.clade_mask
                con_tree.split_edges[new_edge.clade_mask] = new_edge.clade_mask
                parent_node.add_child(new_node)

        ## here we add the support values and/or edge lengths for the terminal taxa ##
        for node in con_tree.leaf_nodes():
            if unrooted:
                split = con_tree.split_edges.normalize_key(node.edge.clade_mask)
            else:
                split = node.edge.clade_mask
            self.map_split_support_to_node(node, 1.0)
            if include_edge_lengths:
                elen = split_distribution.split_edge_lengths.get(split, [0.0])
                node.edge.length = float(sum(elen)) / len(elen)
        return con_tree

    def count_splits_on_trees(self, tree_files, tree_iterator, split_distribution=None, taxa_block=None):
        """
        Given a list of trees file, a SplitsDistribution object (a new one, or,
        if passed as an argument) is returned collating the split data in the files.
        """
        if split_distribution is None:
            split_distribution = splits.SplitDistribution()
        split_distribution.ignore_node_ages = self.ignore_node_ages
        if taxa_block is not None:
            split_distribution.taxa_block = taxa_block
#         if not isinstance(tree_files, list):
#             tree_files = [tree_files]
        total_tree_files = len(tree_files)
        for tree_file_idx, tree_file in enumerate(tree_files):
            file_trees_read = 0
            if isinstance(tree_file, str):
                tree_file_obj = open(tree_file, "r")
            else:
                tree_file_obj = tree_file
            current_file_note = "Tree file %d of %d: " % (tree_file_idx+1, total_tree_files)
            for tree_idx, tree in enumerate(tree_iterator(file_obj=open(tree_file, "r"))):
                self.total_trees_read += 1
                file_trees_read += 1
                if not self.burnin or file_trees_read > self.burnin:
                    self.total_trees_counted += 1
                    self.send_progress_message("%sCounting splits in tree %d" % (current_file_note, (tree_idx+1)))
                    split_distribution.count_splits_on_tree(tree)
                else:
                    self.total_trees_ignored += 1
                    self.send_progress_message("%sSkipping tree %d (burn-in=%d)" % (current_file_note, (tree_idx+1), self.burnin))
        return split_distribution
