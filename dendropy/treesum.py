#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tree summarization and consensus tree building.
"""
import logging

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

import sys

import dendropy
from dendropy import treesplit
from dendropy import dataobject
from dendropy import treecalc
from dendropy import treesim
from dendropy.utility.containers import OrderedDict

# the following imports are for the strict consensus merger
_IS_DEBUG_LOGGING = _LOG.isEnabledFor(logging.DEBUG)
from dendropy.treesplit import encode_splits, count_bits, lowest_bit_only
from dendropy.treemanip import collapse_clade, collapse_edge
from dendropy.dataobject.tree import format_split
from dendropy.utility.containers import NormalizedBitmaskDict
from dendropy.utility.statistics import mean_and_sample_variance

##############################################################################
## TreeSummarizer

class TreeSummarizer(object):
    "Summarizes a distribution of trees."

    def __init__(self, **kwargs):
        """
        __init__ kwargs:

            - `support_as_labels` (boolean)
            - `support_as_percentages` (boolean)
            - `support_label_decimals` (integer)
        """
        self.support_as_labels = kwargs.get("support_as_labels", True)
        self.support_as_percentages = kwargs.get("support_as_percentages", False)
        self.default_support_label_decimals = 4
        self.support_label_decimals = kwargs.get("support_label_decimals", self.default_support_label_decimals)
        self.total_trees_counted = 0
        self.weighted_splits = False

    def tree_from_splits(self,
            split_distribution,
            min_freq=0.5,
            include_edge_lengths=True):
        """Returns a consensus tree from splits in `split_distribution`.

        If include_edge_length_var is True, then the sample variance of the
            edge length will also be calculated and will be stored as
            a length_var attribute.
        """
        taxon_set = split_distribution.taxon_set
        taxa_mask = taxon_set.all_taxa_bitmask()
        if self.weighted_splits:
            split_freqs = split_distribution.weighted_split_frequencies
        else:
            split_freqs = split_distribution.split_frequencies
        is_rooted = split_distribution.is_rooted
        include_edge_lengths = self.support_as_labels and include_edge_lengths

        to_try_to_add = []
        _almost_one = lambda x: abs(x - 1.0) <= 0.0000001
        for s, freq in split_freqs.iteritems():
            if (min_freq is None) or (freq > min_freq) or (_almost_one(min_freq) and _almost_one(freq)):
                to_try_to_add.append((freq, s))
        to_try_to_add.sort(reverse=True)
        splits_for_tree = [i[1] for i in to_try_to_add]

        #to_try_to_add = []
        #_almost_one = lambda x: abs(x - 1.0) <= 0.0000001
        #for s, freq in split_freqs.iteritems():
        #    if (min_freq is None) or (freq > min_freq) or (_almost_one(min_freq) and _almost_one(freq)):
        #        m = s & taxa_mask
        #        if (m != taxa_mask) and ((m-1) & m): # if not root (i.e., all "1's") and not singleton (i.e., one "1")
        #            if not is_rooted:
        #                c = (~m) & taxa_mask
        #                if (c-1) & c: # not singleton (i.e., one "0")
        #                    if 1 & m:
        #                        k = c
        #                    else:
        #                        k = m
        #                    to_try_to_add.append((freq, k, m))
        #            else:
        #                to_try_to_add.append((freq, m, m))
        #to_try_to_add.sort(reverse=True)

        #splits_for_tree = []
        #split_edge_lengths = {}
        #new_split_freq_map = {}
        #for item in to_try_to_add:
        #    splits_for_tree.append(item[1])
        #    split_edge_lengths[item[1]] = split_distribution.split_edge_lengths[item[2]]
        #    new_split_freq_map[item[1]] = item[0]
        #if not include_edge_lengths:
        #    split_edge_lengths = None

        con_tree = treesplit.tree_from_splits(splits=splits_for_tree,
                taxon_set=taxon_set,
                is_rooted=is_rooted)
        treesplit.encode_splits(con_tree)

        if include_edge_lengths:
            split_edge_lengths = {}
            for split, edges in split_distribution.split_edge_lengths.items():
                if len(edges) > 0:
                    mean, var = mean_and_sample_variance(edges)
                    elen = mean
                else:
                    elen = None
                split_edge_lengths[split] = elen
        else:
            split_edge_lengths = None

        for node in con_tree.postorder_node_iter():
            split = node.edge.split_bitmask
            if split in split_freqs:
                self.map_split_support_to_node(node=node, split_support=split_freqs[split])
            if include_edge_lengths and split in split_distribution.split_edge_lengths:
                edges = split_distribution.split_edge_lengths[split]
                if len(edges) > 0:
                    mean, var = mean_and_sample_variance(edges)
                    elen = mean
                else:
                    elen = None
                node.edge.length = elen

        return con_tree

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
                support_label_decimals = self.default_support_label_decimals
            else:
                support_label_decimals = self.support_label_decimals
            support_label_template = "%%0.%df" % support_label_decimals
            support_label = support_label_template % round(split_support_freq,
                support_label_decimals)
        return support_label

    def map_split_support_to_node(self,
            node,
            split_support,
            attr_name="support"):
        "Appropriately sets up a node."
        if self.support_as_percentages:
            support_value = split_support * 100
        else:
            support_value = split_support
        if self.support_as_labels:
            node.label = self.compose_support_label(split_support)
        else:
            node.edge.length = support_value
        if attr_name:
            setattr(node, attr_name, support_value)
            node.annotate(attr_name)
        return node

    def map_split_support_to_tree(self, tree, split_distribution):
        "Maps splits support to the given tree."
        if self.weighted_splits:
            split_freqs = split_distribution.weighted_split_frequencies
        else:
            split_freqs = split_distribution.split_frequencies
        tree.reindex_taxa(taxon_set=split_distribution.taxon_set)
        assert tree.taxon_set is split_distribution.taxon_set
        treesplit.encode_splits(tree)
        for split in tree.split_edges:
            if split in split_freqs:
                split_support = split_freqs[split]
            else:
                split_support = 0.0
            self.map_split_support_to_node(tree.split_edges[split].head_node, split_support)
        return tree

    def annotate_nodes_and_edges(self,
            tree,
            split_distribution):
        """
        Summarizes edge length and age information in `split_distribution` for
        each node on target tree `tree`.
        This will result in each node in `tree` being decorated with the following attributes:
            `age_mean`,
            `age_median`,
            `age_sd`,
            `age_hpd95`,
            `age_range`,
        And each edge in `tree` being decorated with the following attributes:
            `length_mean`,
            `length_median`,
            `length_sd`,
            `length_hpd95`,
            `length_range`,
        These attributes will be added to the annotations dictionary to be persisted.
        """
        assert tree.taxon_set is split_distribution.taxon_set
        if not hasattr(tree, "split_edges"):
            tree.update_splits()
        split_edge_length_summaries = split_distribution.split_edge_length_summaries
        split_node_age_summaries = split_distribution.split_node_age_summaries
        fields = ['mean', 'median', 'sd', 'hpd95', 'quant_5_95', 'range']
        for edge in tree.preorder_edge_iter():
            split = edge.split_bitmask
            nd = edge.head_node
            for summary_name, summary_target, summary_src in [ ('length', edge, split_edge_length_summaries),
                                               ('age', nd, split_node_age_summaries) ]:
                if split in summary_src:
                    summary = summary_src[split]
                    for field in fields:
                        attr_name = summary_name + "_" + field
                        setattr(summary_target, attr_name, summary[field])
                        summary_target.annotate(attr_name)
                else:
                    for field in fields:
                        attr_name = summary_name + "_" + field
                        setattr(summary_target, attr_name, None)
                        #nd.annotate(attr_name)

    def summarize_node_ages_on_tree(self,
            tree,
            split_distribution,
            set_edge_lengths=True,
            collapse_negative_edges=False,
            allow_negative_edges=False,
            summarization_func=None):
        """
        Sets the `age` attribute of nodes on `tree` (a `Tree` object) to the
        result of `summarization_func` applied to the vector of ages of the
        same node on the input trees (in `split_distribution`, a
        `SplitDistribution` object) being summarized.
        `summarization_func` should take an iterable of floats, and return a float. If `None`, it
        defaults to calculating the mean (`lambda x: float(sum(x))/len(x)`).
        If `set_edge_lengths` is `True`, then edge lengths will be set to so that the actual node ages
        correspond to the `age` attribute value.
        If `collapse_negative_edges` is True, then edge lengths with negative values will be set to 0.
        If `allow_negative_edges` is True, then no error will be raised if edges have negative lengths.
        """
        if summarization_func is None:
            summarization_func = lambda x: float(sum(x))/len(x)
        if not hasattr(tree, "split_edges"):
            tree.update_splits()
        #'height',
        #'height_median',
        #'height_95hpd',
        #'height_range',
        #'length',
        #'length_median',
        #'length_95hpd',
        #'length_range',
        #for split, edge in tree.split_edges.items():
        for edge in tree.preorder_edge_iter():
            split = edge.split_bitmask
            nd = edge.head_node
            if split in split_distribution.split_node_ages:
                ages = split_distribution.split_node_ages[split]
                nd.age = summarization_func(ages)
            else:
                # default to age of parent if split not found
                nd.age = nd.parent_node.age
            ## force parent nodes to be at least as old as their oldest child
            if collapse_negative_edges:
                for child in nd.child_nodes():
                    if child.age > nd.age:
                        nd.age = child.age
        if set_edge_lengths:
            tree.set_edge_lengths_from_node_ages(allow_negative_edges=allow_negative_edges)
        return tree

    def summarize_edge_lengths_on_tree(self,
            tree,
            split_distribution,
            summarization_func=None):
        """
        Sets the lengths of edges on `tree` (a `Tree` object) to the mean
        lengths of the corresponding edges on the input trees (in
        `split_distribution`, a `SplitDistribution` object) being
        summarized.
        `summarization_func` should take an iterable of floats, and return a float. If `None`, it
        defaults to calculating the mean (`lambda x: float(sum(x))/len(x)`).
        """
        if summarization_func is None:
            summarization_func = lambda x: float(sum(x))/len(x)
        if not hasattr(tree, "split_edges"):
            tree.update_splits()
        for split, edge in tree.split_edges.items():
            if (split in split_distribution.split_edge_lengths
                    and split_distribution.split_edge_lengths[split]):
                lengths = split_distribution.split_edge_lengths[split]
                edge.length = summarization_func(lengths)
            else:
                edge.length = 0.0
        return tree

        ## here we add the support values and/or edge lengths for the terminal taxa ##
        for node in leaves:
            if not is_rooted:
                split = con_tree.split_edges.normalize_key(node.edge.split_bitmask)
            else:
                split = node.edge.split_bitmask
            self.map_split_support_to_node(node, 1.0)
            if include_edge_lengths:
                elen = split_distribution.split_edge_lengths.get(split, [0.0])
                if len(elen) > 0:
                    mean, var = mean_and_sample_variance(elen)
                    node.edge.length = mean
                    if include_edge_length_var:
                        node.edge.length_var = var
                else:
                    node.edge.length = None
                    if include_edge_length_var:
                        node.edge.length_var = None
        #if include_edge_lengths:
            #self.map_edge_lengths_to_tree(tree=con_tree,
            #        split_distribution=split_distribution,
            #        summarization_func=summarization_func,
            #        include_edge_length_var=False)
        return con_tree

    def count_splits_on_trees(self, tree_iterator, split_distribution=None, trees_splits_encoded=False):
        """
        Given a list of trees file, a SplitsDistribution object (a new one, or,
        if passed as an argument) is returned collating the split data in the files.
        """
        if split_distribution is None:
            split_distribution = treesplit.SplitDistribution()
        taxon_set = split_distribution.taxon_set
        for tree_idx, tree in enumerate(tree_iterator):
            self.total_trees_counted += 1
            if taxon_set is None:
                assert(split_distribution.taxon_set is None)
                split_distribution.taxon_set = tree.taxon_set
                taxon_set = tree.taxon_set
            else:
                assert(taxon_set is tree.taxon_set)
            if not trees_splits_encoded:
                treesplit.encode_splits(tree)
            split_distribution.count_splits_on_tree(tree)
        return split_distribution

## TreeSummarizer
##############################################################################

##############################################################################
## TreeCounter

class TopologyCounter(object):
    """
    Tracks frequency of occurrences of topologies.
    """

    def split_sets_topology_hash(tree):
        """
        Set of all splits on tree: default topology hash.
        """
        return frozenset(tree.split_edges.keys())
    split_sets_topology_hash = staticmethod(split_sets_topology_hash)

    def normalized_newick_topology_hash(tree):
        """
        Returns simple and normalized newick string corresponding to topology.
        """
        tree.ladderize()
        newick_str = tree.as_string(
                'newick',
                edge_lengths=False,
                write_rooting=True,
                intenral_labels=False,
                annotations_as_comments=False,
                annotations_as_nhx=False,
                node_comments=False)
        newick_str.replace("\n", "")
        return newick_str
    normalized_newick_topology_hash = staticmethod(normalized_newick_topology_hash)

    def __init__(self,
            topology_hash_func=None,
            tree_store_func=None):
        """
        `topology_hash_func`, if not None, should be function that takes a DendroPy tree
        object as an argument and returns an object to be used to represent the topology
        in the count-tracking dictonary. If `None` (default), then the set of all
        splits on the tree is used.
        `tree_store_func` should be a function that takes a DendroPy tree
        object as an argument and returns an object to be used to represent the
        tree in the results. If `None` (default), the DendroPy tree object will be
        stored directly.
        The defaults for `topology_hash_func` and `tree_store_func` should be fine in most cases.
        However, when, for example, memory is a concern, but passing in
        normalized NEWICK string composition functions to both of these, references to trees can be
        garbage-collected after being processed, as well as avoiding the need to store
        large collections of splits.
        """
        self.topology_hash_map = {}
        self.total_trees_counted = 0
        if topology_hash_func is None:
            self.topology_hash_func = TopologyCounter.split_sets_topology_hash
        else:
            self.topology_hash_func = topology_hash_func
        if tree_store_func is None:
            self.tree_store_func = lambda x: x
        else:
            self.tree_store_func = tree_store_func

    def update_topology_hash_map(self,
            src_map):
        """
        Imports data from another counter.
        """
        trees_counted = 0
        for topology_hash in src_map:
            if topology_hash not in self.topology_hash_map:
                self.topology_hash_map[topology_hash] = list(src_map[topology_hash])
            else:
                self.topology_hash_map[topology_hash][0] = self.topology_hash_map[topology_hash][0] + src_map[topology_hash][0]
            self.total_trees_counted += src_map[topology_hash][0]

    def count(self,
            tree,
            tree_splits_encoded=False):
        """
        Logs/registers a tree.
        """
        if not tree_splits_encoded:
            treesplit.encode_splits(tree)
        topology = self.topology_hash_func(tree)
        if topology not in self.topology_hash_map:
            self.topology_hash_map[topology] = [1, self.tree_store_func(tree)]
        else:
            self.topology_hash_map[topology][0] = self.topology_hash_map[topology][0] + 1
        self.total_trees_counted += 1

    def calc_freqs(self):
        """
        returns an ordered dictionary (ordereddict) of topologies mapped to a
        tuple, (raw numbers of occurrences, proportion of total trees counted)
        in (descending) order of the proportion of occurrence.
        """
        freqs = OrderedDict()
        topology_counts = self.topology_hash_map.values()
        topology_counts.sort(reverse=True)
        for count, tree in topology_counts:
            freq = float(count) / self.total_trees_counted
            freqs[tree] = (count, freq)
        return freqs

## TreeCounter
##############################################################################

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
        if (p.edge.split_bitmask & common_mask) != l:
            break
        curr_n = p
        p = curr_n.parent_node

    without_lowest = common_mask^l

    taxa_mask = t.seed_node.edge.split_bitmask
    if (curr_n.edge.split_bitmask & common_mask) == l:
        # we did not make it to the root.  Make curr_n, the first_child of the root
        t.to_outgroup_position(curr_n, splits=True, delete_deg_two=True)
        avoid = curr_n
        nd_source = iter(t.seed_node.child_nodes())

        try:
            while True:
                curr_n = nd_source.next()
                if curr_n is not avoid:
                    cm = (curr_n.edge.split_bitmask & without_lowest)
                    if cm:
                        if cm == without_lowest:
                            r = t.seed_node
                            assert curr_n.parent_node is r
                            t.reroot_at(curr_n, splits=True, delete_deg_two=True)
                            t.to_outgroup_position(r, splits=True, delete_deg_two=True)
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
            cm = c.edge.split_bitmask
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
        t.to_outgroup_position(children[0], splits=True, delete_deg_two=True)
        t.to_outgroup_position(p, splits=True, delete_deg_two=True)
    else:
        # if the root first relevant, node then we just make the path leading
        #   to the lowest index node the first child of the root
        t.to_outgroup_position(lowest_on_path_to_l, splits=True, delete_deg_two=True)

def _collapse_paths_not_found(f, s, other_dict=None):
    to_del = []
    for masked_split, path in f.iteritems():
        if masked_split not in s:
            for edge in path:
                if other_dict:
                    del other_dict[edge.split_bitmask]
                collapse_edge(edge)
            to_del.append(masked_split)
    for k in to_del:
        del f[k]


def add_to_scm(to_modify, to_consume, rooted=False, gordons_supertree=False):
    """Adds the tree `to_consume` to the tree `to_modify` in a strict consensus
    merge operation.  Both trees must have had encode_splits called on them."""
    assert(to_modify.taxon_set is to_consume.taxon_set)
    taxon_set = to_consume.taxon_set
    if rooted:
        raise NotImplementedError("rooted form of add_to_scm not implemented")
    to_mod_root = to_modify.seed_node
    to_mod_split = to_mod_root.edge.split_bitmask

    to_consume_root = to_consume.seed_node
    to_consume_split = to_consume_root.edge.split_bitmask

    leaf_intersection = to_mod_split & to_consume_split
    if _IS_DEBUG_LOGGING:
        _LOG.debug("add_to_scm:\n  %s\n  + %s\n%s" % (str(to_modify), str(to_consume), format_split(leaf_intersection, taxon_set=taxon_set)))

    n_common_leaves = count_bits(leaf_intersection)
    if n_common_leaves < 2:
        _LOG.error('trees must have at least 2 common leaves')
        raise ValueError('trees must have at least 2 common leaves')
    if n_common_leaves == 2:
        # SCM with 2 leaves in common results in a polytomy
        collapse_clade(to_mod_root)
        collapse_clade(to_consume_root)
        leaves_to_steal = [c for c in to_consume_root.child_nodes() if not (leaf_intersection & c.edge.split_bitmask)]
        for leaf in leaves_to_steal:
            to_mod_root.add_child(leaf)
            to_mod_root.edge.split_bitmask |= leaf.edge.split_bitmask
        to_modify.split_edges = {to_mod_root.edge.split_bitmask : to_mod_root.edge}
        for child in to_mod_root.child_nodes():
            to_modify.split_edges[child.edge.split_bitmask] = child.edge
        return

    # at least 3 leaves in common
    tmse = to_modify.split_edges

    to_mod_relevant_splits = {}
    to_consume_relevant_splits = {}
    if not rooted:
        if _IS_DEBUG_LOGGING:
            to_modify.debug_check_tree(splits=True, logger_obj=_LOG)
            to_consume.debug_check_tree(splits=True, logger_obj=_LOG)

        reroot_on_lowest_common_index_path(to_modify, leaf_intersection)
        reroot_on_lowest_common_index_path(to_consume, leaf_intersection)

        if _IS_DEBUG_LOGGING:
            to_modify.debug_check_tree(splits=True, logger_obj=_LOG)
            to_consume.debug_check_tree(splits=True, logger_obj=_LOG)

        to_mod_root = to_modify.seed_node
        assert(to_mod_root.edge.split_bitmask == to_mod_split)
        to_consume_root = to_consume.seed_node
        assert(to_consume_root.edge.split_bitmask == to_consume_split)

    for s, e in tmse.iteritems():
        s = e.split_bitmask
        masked = s & leaf_intersection
        if masked and masked != leaf_intersection:
            e_list = to_mod_relevant_splits.setdefault(masked, [])
            e_list.append((s, e))

    for s, e in to_consume.split_edges.iteritems():
        s = e.split_bitmask
        masked = s & leaf_intersection
        if masked and masked != leaf_intersection:
            e_list = to_consume_relevant_splits.setdefault(masked, [])
            e_list.append((s, e))

    # Because each of these paths radiates away from the root (none of the paths
    #   cross the root), the split_bitmasks for deeper edges will be supersets
    #   of the split_bitmasks for shallower nodes.  Thus if we reverse sort we
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
    if _IS_DEBUG_LOGGING:
        to_modify.debug_check_tree(splits=True, logger_obj=_LOG)
        to_consume.debug_check_tree(splits=True, logger_obj=_LOG)


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
    to_steal = [i for i in to_consume_root.child_nodes() if (i.edge.split_bitmask & leaf_intersection) == 0]
    for child in to_steal:
        to_mod_root.add_child(child)
        to_mod_root.edge.split_bitmask |= child.edge.split_bitmask

    for masked_split, to_consume_path in to_consume_relevant_splits.iteritems():
        to_mod_path = to_mod_relevant_splits.get(masked_split)
        if _IS_DEBUG_LOGGING and to_mod_path is None: #to_mod_path is None:
            _LOG.debug("%s = mask" % format_split(leaf_intersection, taxon_set=taxon_set))
            _LOG.debug("%s = masked" % format_split(masked_split, taxon_set=taxon_set))
            _LOG.debug("%s = raw" % format_split(to_consume_path[-1].split_bitmask, taxon_set=taxon_set))
            for k, v in to_mod_relevant_splits.iteritems():
                _LOG.debug("%s in to_mod_relevant_splits" % format_split(k, taxon_set=taxon_set))

        assert to_mod_path is not None
        to_mod_head = to_mod_path[-1].head_node
        to_mod_head_edge = to_mod_head.edge
        to_consume_head = to_consume_path[-1].head_node
        for child in to_consume_head.child_nodes():
            if (child.edge.split_bitmask & leaf_intersection) == 0:
                # child is the root of a subtree that has no children in the leaf_intersection
                to_mod_head.add_child(child)
                to_mod_head_edge.split_bitmask |= child.edge.split_bitmask
        if len(to_consume_path) > 1:
            if len(to_mod_path) > 1:
                # collision
                if gordons_supertree:
                    for edge in to_mod_path[2:]:
                        p = edge.tail_node
                        c = edge.head_node
                        sibs = p.child_nodes()
                        for sib in sibs:
                            _LOG.debug("sib is %s" % (sib.compose_newick()))
                            if sib is not c:
                                if not sib.is_leaf():
                                    collapse_clade(sib)
                                    collapse_edge(sib.edge)
                        collapse_edge(p.edge)
                    mid_node = to_mod_path[0].head_node
                    for edge in to_consume_path[1:]:
                        p = edge.tail_node
                        avoid = edge.head_node
                        for child in p.child_nodes():
                            _LOG.debug("child is %s" % (child.compose_newick()))
                            if child is not avoid:
                                mid_node.add_child(child)
                                collapse_clade(child)
                                if not child.is_leaf():
                                    collapse_edge(child.edge)
                                mid_node.edge.split_bitmask |= child.edge.split_bitmask
                else:
                    for edge in to_mod_path[1:-1]:
                        collapse_edge(edge)
                    mid_node = to_mod_path[0].head_node
                    for edge in to_consume_path[1:]:
                        p = edge.tail_node
                        avoid = edge.head_node
                        for child in p.child_nodes():
                            if child is not avoid:
                                mid_node.add_child(child)
                                mid_node.edge.split_bitmask |= child.edge.split_bitmask
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
    encode_splits(to_modify)

def strict_consensus_merge(tree_list, rooted=False, gordons_supertree=False):
    tree_list = [copy.copy(i) for i in tree_list]
    return inplace_strict_consensus_merge(tree_list, rooted=rooted, gordons_supertree=gordons_supertree)

def inplace_strict_consensus_merge(trees_to_merge, rooted=False, gordons_supertree=False):
    """Returns a tree that is the strict consensus merger of the input trees.
    """
    tree_list = list(trees_to_merge)
    del trees_to_merge[1:]
    nTrees = len(tree_list)
    _LOG.debug('%d Trees to merge:\n%s\n' % (nTrees, '\n'.join([str(i) for i in tree_list])))
    if nTrees < 2:
        return tree_list[0]
    tree_iter = iter(tree_list)
    to_modify = tree_iter.next()

    if rooted:
        raise NotImplementedError("Rooted SCM is not implemented")
    else:
        to_modify.deroot()
    encode_splits(to_modify)
    if _IS_DEBUG_LOGGING:
        assert to_modify._debug_tree_is_valid(splits=False)
    for to_consume in tree_iter:
        if not rooted:
            to_consume.deroot()
        encode_splits(to_consume)
        if _IS_DEBUG_LOGGING:
            assert to_consume._debug_tree_is_valid(splits=True)
        add_to_scm(to_modify, to_consume, rooted, gordons_supertree=gordons_supertree)
        if _IS_DEBUG_LOGGING:
            assert to_modify._debug_tree_is_valid(splits=False)

    return to_modify
