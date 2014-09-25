#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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

import math
import dendropy
import collections
from dendropy.calculate import treesplit
from dendropy.mathlib.statistics import mean_and_sample_variance

##############################################################################
## TreeSummarizer

class TreeSummarizer(object):
    "Summarizes a distribution of trees."

    def __init__(self, **kwargs):
        """
        __init__ kwargs:

            - `support_as_labels` (boolean)
            - `support_as_edge_lengths` (boolean)
            - `support_as_percentages` (boolean)
            - `support_label_decimals` (integer)
        """
        self.support_as_labels = kwargs.get("support_as_labels", True)
        self.support_as_edge_lengths = kwargs.get("support_as_edge_lengths", False)
        self.support_as_percentages = kwargs.get("support_as_percentages", False)
        self.add_node_metadata = kwargs.get("add_node_metadata", True)
        self.default_support_label_decimals = 4
        self.support_label_decimals = kwargs.get("support_label_decimals", self.default_support_label_decimals)
        self.total_trees_counted = 0
        self.weighted_splits = False

    def tree_from_splits(self,
            split_distribution,
            min_freq=0.5,
            rooted=None,
            include_edge_lengths=True):
        """Returns a consensus tree from splits in `split_distribution`.

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
        con_tree = treesplit.tree_from_splits(splits=splits_for_tree,
                taxon_namespace=taxon_namespace,
                is_rooted=rooted)
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
        treesplit.encode_splits(tree)
        for split in tree.split_edge_map:
            if split in split_freqs:
                split_support = split_freqs[split]
            else:
                split_support = 0.0
            self.map_split_support_to_node(tree.split_edge_map[split].head_node, split_support)
        return tree

    def annotate_nodes_and_edges(self,
            tree,
            split_distribution,
            recalculate_splits=False):
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
        assert tree.taxon_namespace is split_distribution.taxon_namespace
        if recalculate_splits or tree.split_edge_map is None:
            tree.encode_splits()
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
            summarization_func=None,
            recalculate_splits=False):
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
        if recalculate_splits or tree.split_edge_map is None:
            tree.encode_splits()
        #'height',
        #'height_median',
        #'height_95hpd',
        #'height_range',
        #'length',
        #'length_median',
        #'length_95hpd',
        #'length_range',
        #for split, edge in tree.split_edge_map.items():
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
            summarization_func=None,
            recalculate_splits=False):
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
        if recalculate_splits or tree.split_edge_map is None:
            tree.encode_splits()
        for split, edge in tree.split_edge_map.items():
            if (split in split_distribution.split_edge_lengths
                    and split_distribution.split_edge_lengths[split]):
                lengths = split_distribution.split_edge_lengths[split]
                #if len(lengths) != split_distribution.total_trees_counted:
                #    # not all input trees had edge lengths (at least, for this split)
                #    pass
                edge.length = summarization_func(lengths)
            elif (split in split_distribution.split_edge_lengths
                    and not split_distribution.split_edge_lengths[split]):
                # no input trees had any edge lengths for this split
                edge.length = None
            else:
                # split on target tree that was not found in any of the input
                # trees
                edge.length = 0.0
        return tree

        ## here we add the support values and/or edge lengths for the terminal taxa ##
        for node in leaves:
            if not is_rooted:
                split = con_tree.split_edge_map.normalize_key(node.edge.split_bitmask)
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
        taxon_namespace = split_distribution.taxon_namespace
        for tree_idx, tree in enumerate(tree_iterator):
            self.total_trees_counted += 1
            if taxon_namespace is None:
                assert(split_distribution.taxon_namespace is None)
                split_distribution.taxon_namespace = tree.taxon_namespace
                taxon_namespace = tree.taxon_namespace
            else:
                assert(taxon_namespace is tree.taxon_namespace)
            if not trees_splits_encoded:
                treesplit.encode_splits(tree)
            split_distribution.count_splits_on_tree(tree)
        return split_distribution

    def consensus_tree(self, trees, min_freq=0.5, trees_splits_encoded=False):
        """
        Returns a consensus tree of all trees in `trees`, with minumum frequency
        of split to be added to the consensus tree given by `min_freq`.
        """
        taxon_namespace = trees[0].taxon_namespace
        split_distribution = treesplit.SplitDistribution(taxon_namespace=taxon_namespace)
        self.count_splits_on_trees(trees,
                split_distribution=split_distribution,
                trees_splits_encoded=trees_splits_encoded)
        tree = self.tree_from_splits(split_distribution, min_freq=min_freq)
        return tree

    def calculate_tree_clade_credibilities(self,
            trees,
            split_distribution=None,
            burnin_offset=0,
            log_product_of_clade_posteriors_attr_name="log_product_of_clade_posteriors",
            sum_of_clade_posteriors_attr_name="sum_of_clade_posteriors",
            trees_splits_encoded=False):
        """
        Calculates the "clade credibility scores" of trees from a set of trees,
        storing the log product of clade posteriors in a new attribute of the tree,
        `log_product_of_clade_posteriors`, and the sum of clade posteriors in
        another attribute of the tree, `sum_of_clade_posteriors`. Returns a
        tuple, with the first element the tree with the highest product of
        clade posteriors and the second element the tree with highest sum of
        clade posteriors.

        From `Wikipedia <http://en.wikipedia.org/wiki/Maximum_clade_credibility_tree>`_ :

            A maximum clade credibility tree is a tree that summarises the
            results of a Bayesian phylogenetic inference. Whereas a
            majority-rule tree combines the most common clades, potentially
            resulting in a tree that was not sampled during the analysis, the
            maximum-credibility method evaluates each of the sampled posterior
            trees. Each clade within the tree is given a score based on the
            number of times that it appears in other sampled posterior trees,
            and these scores are added to give a total score for the tree. The
            tree with the highest score represents the maximum clade
            credibility tree.

            Since each clade's score is akin to a probability, it may be more
            appropriate to multiply, rather than add, each clade's score
            (expressed as a probability, or the fraction of posterior trees
            that contain the clade) to generate a total score for the
            tree. This would generate a Maximum credibility tree. However, both
            methods are used in various contexts.

        Parameters
        ----------
        trees : a :class:`TreeList` or iterable of trees
            The list of trees from which to compute the clade posterior
            probabilities.
        split_distribution : :class:`SplitDistribution`
            If split frequencies have already been computed, the
            :class:`SplitDistribution` used to manage the counts can be passed
            here. If `None`, and new one will be created and the splits will be counted.
        burnin_offset : int
            Number of trees to skip for burinin (only applies if a
            :class:`SplitDistribution` instance was not passed and the
            splits have to be counted.
        log_product_of_clade_posteriors_attr_name : str
            The attribute on each tree with which to store the log of the
            product of the clade posteriors.
        sum_of_clade_posteriors_attr_name : str
            The attribute on each tree with which to store the sum of the clade
            posteriors.
        tree_splits_encoded : bool
            If `True`, then the splits are assumed to have already been encoded
            and will not be updated on three trees.

        Returns
        -------
        mct_tree : :class:`Tree`
            The :class:`Tree` with the highest product of clade posterior probabilities.
        mcct_tree : :class:`Tree`
            The :class:`Tree` with the highest sum of clade posterior probabilities.

        Examples
        --------

            trees = dendropy.TreeList.get_from_path(
                    "data/mcmc1.rooted.nex",
                    "nexus")
            tsum = treesum.TreeSummarizer()
            mcct, mct = tsum.calculate_tree_clade_credibilities(trees)

            result_trees = (mcct, mct)
            tree_descs = ("Maximum Credibility Tree", "Maximum Clade Credibility Tree")

            for tree, tree_desc in zip(result_trees, tree_descs):
                print("{:>30}: {} '{}': {} {}".format(
                    tree_desc,
                    trees.index(mct)+1,
                    mct.label,
                    mct.log_product_of_clade_posteriors,
                    mct.sum_of_clade_posteriors))

        """
        if burnin_offset is not None:
            trees = trees[burnin_offset:]
        taxon_namespace = trees[0].taxon_namespace
        if split_distribution is None:
            split_distribution = treesplit.SplitDistribution(taxon_namespace=taxon_namespace)
            self.count_splits_on_trees(trees,
                    split_distribution=split_distribution,
                    trees_splits_encoded=trees_splits_encoded)
            tree_splits_encoded = True
        max_product_of_clade_posteriors_tree = None
        max_product_of_clade_posteriors = None
        max_sum_of_clade_posteriors_tree = None
        max_sum_of_clade_posteriors = None
        for tree in trees:
            log_product_of_clade_posteriors = 0
            sum_of_clade_posteriors = 0
            if not tree_splits_encoded:
                tree.encode_splits()
            for split in tree.split_edge_map:
                posterior = split_distribution[split]
                if posterior:
                    log_product_of_clade_posteriors += math.log(posterior)
                    sum_of_clade_posteriors += posterior
            setattr(tree, log_product_of_clade_posteriors_attr_name, log_product_of_clade_posteriors)
            setattr(tree, sum_of_clade_posteriors_attr_name, sum_of_clade_posteriors)
            if max_product_of_clade_posteriors is None or log_product_of_clade_posteriors > max_product_of_clade_posteriors:
                max_product_of_clade_posteriors = log_product_of_clade_posteriors
                max_product_of_clade_posteriors_tree = tree
            if max_sum_of_clade_posteriors is None or log_product_of_clade_posteriors > max_sum_of_clade_posteriors:
                max_sum_of_clade_posteriors = log_product_of_clade_posteriors
                max_sum_of_clade_posteriors_tree = tree
        return max_product_of_clade_posteriors_tree, max_sum_of_clade_posteriors_tree

##############################################################################
## Convenience Function to Get Consensus Tree

def consensus_tree(trees, min_freq=0.5, trees_splits_encoded=False, **kwargs):
    """
    Returns a consensus tree of all trees in `trees`, with minumum frequency
    of split to be added to the consensus tree given by `min_freq`.
    """
    tsum = TreeSummarizer(**kwargs)
    return tsum.consensus_tree(trees,
            min_freq=min_freq,
            trees_splits_encoded=trees_splits_encoded)

##############################################################################
## TreeCounter

class TopologyCounter(object):
    """
    Tracks frequency of occurrences of topologies.
    """

    def hash_topology(tree):
        """
        Set of all splits on tree: default topology hash.
        """
        return frozenset(tree.split_edge_map.keys())
    hash_topology = staticmethod(hash_topology)

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
            tree_splits_encoded=False):
        """
        Logs/registers a tree.
        """
        if not tree_splits_encoded:
            treesplit.encode_splits(tree)
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
            tree = treesplit.tree_from_splits(splits=topology_hash,
                taxon_namespace=taxon_namespace,
                is_rooted=is_rooted)
            tree_freqs[tree] = (count, freq)
        return tree_freqs

## TreeCounter
##############################################################################

