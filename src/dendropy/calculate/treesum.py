#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
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

import math
import collections
import dendropy
from dendropy.datamodel import taxonmodel
from dendropy.calculate.statistics import mean_and_sample_variance

##############################################################################
## TreeSummarizer

class TreeSummarizer(object):
    "Summarizes a distribution of trees."

    def __init__(self, **kwargs):
        """
        __init__ kwargs:

            - ``support_as_labels`` (boolean)
            - ``support_as_edge_lengths`` (boolean)
            - ``support_as_percentages`` (boolean)
            - ``support_label_decimals`` (integer)
        """
        self.support_as_labels = kwargs.get("support_as_labels", True)
        self.support_as_edge_lengths = kwargs.get("support_as_edge_lengths", False)
        self.support_as_percentages = kwargs.get("support_as_percentages", False)
        self.add_node_metadata = kwargs.get("add_node_metadata", True)
        self.default_support_label_decimals = 4
        self.support_label_decimals = kwargs.get("support_label_decimals", self.default_support_label_decimals)
        self.weighted_splits = False

    def tree_from_splits(self,
            split_distribution,
            min_freq=0.5,
            rooted=None,
            include_edge_lengths=True):
        """Returns a consensus tree from splits in ``split_distribution``.

        If include_edge_length_var is True, then the sample variance of the
            edge length will also be calculated and will be stored as
            a length_var attribute.
        """
        taxon_namespace = split_distribution.taxon_namespace
        taxa_mask = taxon_namespace.all_taxa_bitmask()
        if self.weighted_splits:
            split_freqs = split_distribution.weighted_split_frequencies
        else:
            split_freqs = split_distribution.split_frequencies
        if rooted is None:
            if split_distribution.is_all_counted_trees_rooted():
                rooted = True
            elif split_distribution.is_all_counted_trees_strictly_unrooted:
                rooted = False
        #include_edge_lengths = self.support_as_labels and include_edge_lengths
        if self.support_as_edge_lengths and include_edge_lengths:
            raise Exception("Cannot map support as edge lengths if edge lengths are to be set on consensus tree")
        to_try_to_add = []
        _almost_one = lambda x: abs(x - 1.0) <= 0.0000001
        for s in split_freqs:
            freq = split_freqs[s]
            if (min_freq is None) or (freq >= min_freq) or (_almost_one(min_freq) and _almost_one(freq)):
                to_try_to_add.append((freq, s))
        to_try_to_add.sort(reverse=True)
        splits_for_tree = [i[1] for i in to_try_to_add]
        con_tree = dendropy.Tree.from_split_bitmasks(
                split_bitmasks=splits_for_tree,
                taxon_namespace=taxon_namespace,
                is_rooted=rooted)
        con_tree.encode_bipartitions()

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
            split = node.edge.bipartition.split_bitmask
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
        if self.support_as_edge_lengths:
            node.edge.length = support_value
        if self.add_node_metadata and attr_name:
            setattr(node, attr_name, support_value)
            node.annotations.drop(name=attr_name)
            node.annotations.add_bound_attribute(attr_name,
                    real_value_format_specifier=".{}f".format(self.support_label_decimals))
        return node

    def map_split_support_to_tree(self, tree, split_distribution):
        "Maps splits support to the given tree."
        if self.weighted_splits:
            split_freqs = split_distribution.weighted_split_frequencies
        else:
            split_freqs = split_distribution.split_frequencies
        tree.reindex_taxa(taxon_namespace=split_distribution.taxon_namespace)
        assert tree.taxon_namespace is split_distribution.taxon_namespace
        tree.encode_bipartitions()
        for edge in tree.postorder_edge_iter():
            split = edge.bipartition.split_bitmask
            if split in split_freqs:
                split_support = split_freqs[split]
            else:
                split_support = 0.0
            self.map_split_support_to_node(edge.head_node, split_support)
        return tree

    def annotate_nodes_and_edges(self,
            tree,
            split_distribution,
            is_bipartitions_updated=False,):
        """
        Summarizes edge length and age information in ``split_distribution`` for
        each node on target tree ``tree``.
        This will result in each node in ``tree`` being decorated with the following attributes:
            ``age_mean``,
            ``age_median``,
            ``age_sd``,
            ``age_hpd95``,
            ``age_range``,
        And each edge in ``tree`` being decorated with the following attributes:
            ``length_mean``,
            ``length_median``,
            ``length_sd``,
            ``length_hpd95``,
            ``length_range``,
        These attributes will be added to the annotations dictionary to be persisted.
        """
        assert tree.taxon_namespace is split_distribution.taxon_namespace
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
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
                        # clear annotations, which might be associated with either nodes
                        # or edges due to the way NEXUS/NEWICK node comments are parsed
                        nd.annotations.drop(name=attr_name)
                        edge.annotations.drop(name=attr_name)
                        summary_target.annotations.add_bound_attribute(attr_name)
                else:
                    for field in fields:
                        attr_name = summary_name + "_" + field
                        setattr(summary_target, attr_name, None)

    def summarize_node_ages_on_tree(self,
            tree,
            split_distribution,
            set_edge_lengths=True,
            collapse_negative_edges=False,
            allow_negative_edges=False,
            summarization_fn=None,
            is_bipartitions_updated=False):
        """
        Sets the ``age`` attribute of nodes on ``tree`` (a |Tree| object) to the
        result of ``summarization_fn`` applied to the vector of ages of the
        same node on the input trees (in ``split_distribution``, a
        `SplitDistribution` object) being summarized.
        ``summarization_fn`` should take an iterable of floats, and return a float. If |None|, it
        defaults to calculating the mean (``lambda x: float(sum(x))/len(x)``).
        If ``set_edge_lengths`` is |True|, then edge lengths will be set to so that the actual node ages
        correspond to the ``age`` attribute value.
        If ``collapse_negative_edges`` is True, then edge lengths with negative values will be set to 0.
        If ``allow_negative_edges`` is True, then no error will be raised if edges have negative lengths.
        """
        if summarization_fn is None:
            summarization_fn = lambda x: float(sum(x))/len(x)
        if is_bipartitions_updated:
            tree.encode_splits()
        #'height',
        #'height_median',
        #'height_95hpd',
        #'height_range',
        #'length',
        #'length_median',
        #'length_95hpd',
        #'length_range',
        for edge in tree.preorder_edge_iter():
            split = edge.bipartition.split_bitmask
            nd = edge.head_node
            if split in split_distribution.split_node_ages:
                ages = split_distribution.split_node_ages[split]
                nd.age = summarization_fn(ages)
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
            summarization_fn=None,
            is_bipartitions_updated=False):
        """
        Sets the lengths of edges on ``tree`` (a |Tree| object) to the mean
        lengths of the corresponding edges on the input trees (in
        ``split_distribution``, a `SplitDistribution` object) being
        summarized.
        ``summarization_fn`` should take an iterable of floats, and return a float. If |None|, it
        defaults to calculating the mean (``lambda x: float(sum(x))/len(x)``).
        """
        if summarization_fn is None:
            summarization_fn = lambda x: float(sum(x))/len(x)
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
        for edge in tree.postorder_edge_iter():
            split = edge.bipartition.split_bitmask
            if (split in split_distribution.split_edge_lengths
                    and split_distribution.split_edge_lengths[split]):
                lengths = split_distribution.split_edge_lengths[split]
                edge.length = summarization_fn(lengths)
            elif (split in split_distribution.split_edge_lengths
                    and not split_distribution.split_edge_lengths[split]):
                # no input trees had any edge lengths for this split
                edge.length = None
            else:
                # split on target tree that was not found in any of the input
                # trees
                edge.length = 0.0
        return tree

    def count_splits_on_trees(self, tree_iterator, split_distribution=None, is_bipartitions_updated=False):
        """
        Given a list of trees file, a SplitsDistribution object (a new one, or,
        if passed as an argument) is returned collating the split data in the files.
        """
        if split_distribution is None:
            split_distribution = dendropy.SplitDistribution()
        taxon_namespace = split_distribution.taxon_namespace
        for tree_idx, tree in enumerate(tree_iterator):
            if taxon_namespace is None:
                assert(split_distribution.taxon_namespace is None)
                split_distribution.taxon_namespace = tree.taxon_namespace
                taxon_namespace = tree.taxon_namespace
            else:
                assert(taxon_namespace is tree.taxon_namespace)
            split_distribution.count_splits_on_tree(tree,
                    is_bipartitions_updated=is_bipartitions_updated)
        return split_distribution

    def consensus_tree(self, trees, min_freq=0.5, is_bipartitions_updated=False):
        """
        Returns a consensus tree of all trees in ``trees``, with minumum frequency
        of split to be added to the consensus tree given by ``min_freq``.
        """
        taxon_namespace = trees[0].taxon_namespace
        split_distribution = dendropy.SplitDistribution(taxon_namespace=taxon_namespace)
        self.count_splits_on_trees(trees,
                split_distribution=split_distribution,
                is_bipartitions_updated=is_bipartitions_updated)
        tree = self.tree_from_splits(split_distribution, min_freq=min_freq)
        return tree

##############################################################################
## Convenience Wrappers

def consensus_tree(trees, min_freq=0.5, is_bipartitions_updated=False, **kwargs):
    """
    Returns a consensus tree of all trees in ``trees``, with minumum frequency
    of split to be added to the consensus tree given by ``min_freq``.
    """
    tsum = TreeSummarizer(**kwargs)
    return tsum.consensus_tree(trees,
            min_freq=min_freq,
            is_bipartitions_updated=is_bipartitions_updated)

##############################################################################
## TreeCounter

class TopologyCounter(object):
    """
    Tracks frequency of occurrences of topologies.
    """

    @staticmethod
    def hash_topology(tree):
        """
        Set of all splits on tree: default topology hash.
        """
        return frozenset(tree.bipartition_encoding)

    def __init__(self):
        self.topology_hash_map = {}
        self.total_trees_counted = 0

    def update_topology_hash_map(self,
            src_map):
        """
        Imports data from another counter.
        """
        trees_counted = 0
        for topology_hash in src_map:
            if topology_hash not in self.topology_hash_map:
                self.topology_hash_map[topology_hash] = src_map[topology_hash]
            else:
                self.topology_hash_map[topology_hash] = self.topology_hash_map[topology_hash] + src_map[topology_hash]
            self.total_trees_counted += src_map[topology_hash]

    def count(self,
            tree,
            is_bipartitions_updated=False):
        """
        Logs/registers a tree.
        """
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
        topology = self.hash_topology(tree)
        if topology not in self.topology_hash_map:
            self.topology_hash_map[topology] = 1
        else:
            self.topology_hash_map[topology] = self.topology_hash_map[topology] + 1
        self.total_trees_counted += 1

    def calc_hash_freqs(self):
        """
        Returns an ordered dictionary (collections.OrderedDict) of topology hashes mapped
        to a tuple, (raw numbers of occurrences, proportion of total trees
        counted) in (descending) order of the proportion of occurrence.
        """
        t_freqs = collections.OrderedDict()
        count_topology = [(v, k) for k, v in self.topology_hash_map.items()]
        count_topology.sort(reverse=True)
        for count, topology_hash in count_topology:
            freq = float(count) / self.total_trees_counted
            t_freqs[topology_hash] = (count, freq)
        return t_freqs

    def calc_tree_freqs(self, taxon_namespace, is_rooted=False):
        """
        Returns an ordered dictionary (collections.OrderedDict) of DendroPy trees mapped
        to a tuple, (raw numbers of occurrences, proportion of total trees
        counted) in (descending) order of the proportion of occurrence.
        """
        hash_freqs = self.calc_hash_freqs()
        tree_freqs = collections.OrderedDict()
        for topology_hash, (count, freq) in hash_freqs.items():
            tree = dendropy.Tree.from_bipartition_encoding(
                bipartition_encoding=topology_hash,
                taxon_namespace=taxon_namespace,
                is_rooted=is_rooted)
            tree_freqs[tree] = (count, freq)
        return tree_freqs

## TreeCounter
##############################################################################

