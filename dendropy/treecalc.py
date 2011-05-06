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
Tree metrics/statistics calculations.
"""
from itertools import izip
from math import sqrt
from dendropy import treesplit
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

class PatristicDistanceMatrix(object):
    """
    Calculates and maintains patristic distance information of taxa on a tree.
    `max_dist_taxa` and `max_dist_nodes` will return a tuple of taxon objects
    and corresponding nodes, respectively, that span the greatest path distance
    in the tree. The mid-point between the two is *guaranteed* to be on the
    closer to the first item of each pair.
    """

    def __init__(self, tree=None):
        self.tree = None
        self.taxon_set = None
        self._pat_dists = {}
        self.max_dist = None
        self.max_dist_taxa = None
        self.max_dist_nodes = None
        self._mrca = {}
        if tree is not None:
            self.calc(tree)

    def __call__(self, taxon1, taxon2):
        """
        Returns patristic distance between two taxon objects.
        """
        if taxon1 is taxon2:
            return 0.0
        try:
            return self._pat_dists[taxon1][taxon2]
        except KeyError:
            return self._pat_dists[taxon2][taxon1]

    def mrca(self, taxon1, taxon2):
        """
        Returns MRCA of two taxon objects.
        """
        if taxon1 is taxon2:
            return taxon1
        try:
            return self._mrca[taxon1][taxon2]
        except KeyError:
            return self._mrca[taxon2][taxon1]

    def calc(self, tree=None, create_midpoints=None):
        """
        Calculates the distances.
        """
        if tree is not None:
            self.tree = tree
        assert tree is not None
        if not hasattr(self.tree, "split_edges"):
            treesplit.encode_splits(self.tree)
        self.taxon_set = tree.taxon_set
        self._pat_dists = {}
        for i1, t1 in enumerate(self.taxon_set):
            self._pat_dists[t1] = {}
            self._mrca[t1] = {}
            self.max_dist = None
            self.max_dist_taxa = None
            self.max_dist_nodes = None

        for node in tree.postorder_node_iter():
            children = node.child_nodes()
            if len(children) == 0:
                node.desc_paths = {node : 0}
            else:
                node.desc_paths = {}
                for cidx1, c1 in enumerate(children):
                    for desc1, desc1_plen in c1.desc_paths.items():
                        node.desc_paths[desc1] = desc1_plen + c1.edge.length
                        for c2 in children[cidx1+1:]:
                            for desc2, desc2_plen in c2.desc_paths.items():
                                pat_dist = node.desc_paths[desc1] + desc2_plen + c2.edge.length
                                self._pat_dists[desc1.taxon][desc2.taxon] = pat_dist
                                self._mrca[desc1.taxon][desc2.taxon] = c1.parent_node
                                if pat_dist > self.max_dist:
                                    self.max_dist = pat_dist
                                    midpoint = float(pat_dist) / 2
                                    if midpoint - node.desc_paths[desc1] <= 0:
                                        self.max_dist_nodes = (desc1, desc2)
                                        self.max_dist_taxa = (desc1.taxon, desc2.taxon)
                                    else:
                                        self.max_dist_nodes = (desc2, desc1)
                                        self.max_dist_taxa = (desc2.taxon, desc1.taxon)
                    del(c1.desc_paths)

    def distances(self):
        """
        Returns list of patristic distances.
        """
        dists = []
        for dt in self._pat_dists.values():
            for d in dt.values():
                dists.append(d)
        return dists

    def sum_of_distances(self):
        """
        Returns sum of patristic distances on tree.
        """
        return sum(self.distances())

def patristic_distance(tree, taxon1, taxon2):
    """
    Given a tree with splits encoded, and two taxa on that tree, returns the
    patristic distance between the two. Much more inefficient than constructing
    a PatristicDistanceMatrix object.
    """
    if not hasattr(tree, "split_edges"):
        treesplit.encode_splits(tree)
    mrca = tree.mrca(taxa=[taxon1, taxon2])
    dist = 0
    n = tree.find_node(lambda x: x.taxon == taxon1)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    n = tree.find_node(lambda x: x.taxon == taxon2)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    return dist

def robinson_foulds_calc(length_diffs):
    """
    Given, `length_diffs`, a list of pairs of corresponding (length/weight) values
    of edges from two trees, this returns the Robinson-Foulds distance (sum of
    absolute differences) between the two trees.
    """
    return sum([abs(i[0] - i[1]) for i in length_diffs])

def brlen_scores_calc(length_diffs):
    """
    Given, `length_diffs`, a list of pairs of corresponding (length/weight) values
    of edges from two trees, this returns the branch length score (sum of
    squared differences) between the two trees. This is equivalent to the squared
    Euclidean distance between the two trees.
    """
    d = [pow(i[0] - i[1], 2) for i in length_diffs]
    return sum(d)

def brlen_dists_calc(length_diffs):
    """
    Given, `length_diffs`, a list of pairs of corresponding (length/weight) values
    of edges from two trees, this returns the branch length distance (square root of the
    sum of squared differences) between the two trees. This is equivalent to the Euclidean
    branch length distance between the two trees.
    """
    return sqrt(brlen_scores_calc(length_diffs))

def get_length_diffs(tree1,
        tree2,
        edge_length_attr="length",
        value_type=float,
        split_length_diff_map=False):
    """
    Returns a list of tuples, with the first element of each tuple representing
    the length of the branch subtending a particular split on ``tree1``, and
    the second element the length of the same branch on ``tree2``. If a
    particular split is found on one tree but not in the other, a value of zero
    is used for the missing split.
    """
    length_diffs = []
    split_length_diffs = {}
    if tree1.taxon_set is not tree2.taxon_set:
        raise TypeError("Trees have different TaxonSet objects: %s vs. %s" \
                % (hex(id(tree1.taxon_set)), hex(id(tree2.taxon_set))))
    if not hasattr(tree1, "split_edges"):
        treesplit.encode_splits(tree1)
    if not hasattr(tree2, "split_edges"):
        treesplit.encode_splits(tree2)
    split_edges2_copy = dict(tree2.split_edges) # O(n*(2*bind + dict_item_cost))
    split_edges1_ref = tree1.split_edges
    for split, edge in split_edges1_ref.iteritems(): # O n : 2*bind
        elen1 = getattr(edge, edge_length_attr) # attr + bind
        if elen1 is None:
            elen1 = 0 # worst-case: bind
        value1 = value_type(elen1) #  ctor + bind
        try:
            e2 = split_edges2_copy.pop(split) # attr + dict_lookup + bind
            elen2 = getattr(e2, edge_length_attr) # attr + bind
            if elen2 is None:
                # allow root edge to have split with no value: raise error if not root edge
                if e2.tail_node is None:
                    elen2 = 0.0
                else:
                    raise ValueError("Edge length attribute is 'None': Tree: %s ('%s'), Split: %s" % (tree2.oid, tree2.label, tree2.taxon_set.split_as_newick_string(split)))
        except KeyError: # excep
            elen2 = 0.0
        value2 = value_type(elen2) #  ctor + bind # best case
        length_diffs.append((value1,value2)) # ctor + listappend
        split_length_diffs[split] = length_diffs[-1]

    for split, edge in split_edges2_copy.iteritems(): # best-case not executed, worst case O(n) : 2*bind
        elen2 = getattr(edge, edge_length_attr) # attr +  bind
        if elen2 is None:
            elen2 = 0
        value2 = value_type(elen2) #  ctor + bind
        e1 = split_edges1_ref.get(split) # attr + dict_lookup + bind
        if e1 is None:
            elen1 = 0.0
        else:
            elen1 = getattr(e1, edge_length_attr) # attr  + bind
            if elen1 is None:
                # allow root edge to have split with no value: raise error if not root edge
                if e1.tail_node is None:
                    elen1 = 0.0
                else:
                    raise ValueError("Edge length attribute is 'None': Tree: %s ('%s'), Split: %s" % (tree1.oid, tree1.label, split))
                #elen1 = 0
        value1 = value_type(elen1)
        length_diffs.append((value1,value2)) # ctor + listappend
        split_length_diffs[split] = length_diffs[-1]
    # the numbers below do not reflect additions to the code to protect against
    #   edges with length None
    # loops
    #  best-case:
    #   O(n * (dict_lookup + 3*attr + 3*ctor + 7*bind + listappend))
    #  worst-case:
    #     separated: O(n * (2*dict_lookup + 4*attr + 3*ctor + 8*bind + listappend + excep) + n*(2*dict_lookup + 4*attr + 3*ctor + 8*bind + listappend))
    #   or:
    #     O(2n*(2*dict_lookup + 4*attr + 3*ctor + 8*bind + listappend + 0.5*excep))

    # total
    #  best-case:
    #       O(n * (dict_lookup + 3*attr + 3*ctor + 8*bind + listappend + dict_item_cost))
    #  worst-case:
    #     O(2n*(2*dict_lookup + 4*attr + 3*ctor + 9*bind + listappend + 0.5*(dict_item_cost + excep))
    if split_length_diff_map:
        return length_diffs, split_length_diffs
    else:
        return length_diffs

def splits_distance(tree1,
                    tree2,
                    dist_func=robinson_foulds_calc,
                    edge_length_attr="length",
                    value_type=float):
    """
    Returns distance between two trees, each represented by a dictionary of
    splits (as split_mask strings) to edges, using `dist_func` to calculate the
    distance based on `edge_length_attr` of the edges. `dist_func` is a function
    that takes a list of pairs of values, where the values correspond to the edge
    lengths of a given split on tree1 and tree2 respectively.
    """
    length_diffs = get_length_diffs(tree1, tree2, edge_length_attr=edge_length_attr, value_type=value_type)
    return dist_func(length_diffs)

def robinson_foulds_distance(tree1, tree2, edge_length_attr="length"):
    """
    Returns Robinson-Foulds distance between two trees based on `edge_length_attr`.
    Trees need to have been decorated with the `encode_splits` method of the splits
    module.
    """
    return splits_distance(tree1,
                           tree2,
                           dist_func=robinson_foulds_calc,
                           edge_length_attr=edge_length_attr,
                           value_type=float)

def euclidean_distance(tree1, tree2, edge_length_attr="length", value_type=float):
    """
    Returns Euclidean distance (a.k.a. Felsenstein's 2004 `branch length distance`)
    between two trees based on `edge_length_attr`.
    Trees need to have been decorated with the `encode_splits` method of the splits
    module.
    """
    return splits_distance(tree1,
                           tree2,
                           dist_func=brlen_dists_calc,
                           edge_length_attr=edge_length_attr,
                           value_type=value_type)

def symmetric_difference(tree1, tree2):
    "Returns the number of splits that are present in only 1 of the 2 trees."
    t = false_positives_and_negatives(tree1, tree2)
    return t[0] + t[1]

def false_positives_and_negatives(reference_tree, test_tree):
    """
    False pos = splits in test_tree NOT in reference_tree
    False neg = splits in reference_tree NOT in test_tree
    """
    sym_diff = 0
    false_positives = 0
    false_negatives = 0

    if reference_tree.taxon_set is not test_tree.taxon_set:
        raise TypeError("Trees have different TaxonSet objects: %s vs. %s" \
                % (hex(id(reference_tree.taxon_set)), hex(id(test_tree.taxon_set))))
    if not hasattr(reference_tree, "split_edges"):
        treesplit.encode_splits(reference_tree)
    if not hasattr(test_tree, "split_edges"):
        treesplit.encode_splits(test_tree)
    for split in reference_tree.split_edges:
        if split in test_tree.split_edges:
            pass
        else:
            false_negatives = false_negatives + 1
            sym_diff = sym_diff + 1

    for split in test_tree.split_edges:
        if split in reference_tree.split_edges:
            pass
        else:
            false_positives = false_positives + 1
            sym_diff = sym_diff + 1

    return false_positives, false_negatives


def fitch_down_pass(postorder_node_list, attr_name="state_sets", weight_list=None, taxa_to_state_set_map=None):
    """
    Reads `attr_name` attribute of leaves as an iterable of state sets, and
    sets that attribute for internal nodes using the "preliminary phase" of
    Fitch's (1971) unordered parsimony algorithm.
    Returns the parsimony score.

        `weight_list`
            is an optional vector of weights for each pattern.
        `taxa_to_state_set_map`
            if a child node does not have an attribute with name
        `attr_name`
            then the nodes.taxon will be used as a key in taxa_to_state_set_map
            to find the state set. This allows for the scoring of
            previously undecorated trees.

    Currently this requires a bifurcating tree (even at the root).
    """
    score = 0
    for nd in postorder_node_list:
        c = nd.child_nodes()
        if not c:
            try:
                ss = getattr(nd, attr_name)
            except AttributeError:
                if not taxa_to_state_set_map:
                    raise
                ss = taxa_to_state_set_map[nd.taxon]
                setattr(nd, attr_name, ss)
            continue
        left_c, right_c = c[:2]
        remaining = c[2:]
        try:
            left_ssl = getattr(left_c, attr_name)
        except AttributeError:
            if not taxa_to_state_set_map:
                raise
            left_ssl = taxa_to_state_set_map[left_c.taxon]

        while True:
            try:
                right_ssl = getattr(right_c, attr_name)
            except AttributeError:
                if not taxa_to_state_set_map:
                    raise
                right_ssl = taxa_to_state_set_map[right_c.taxon]
            result = []
            for n, ssp in enumerate(izip(left_ssl, right_ssl)):
                left_ss, right_ss = ssp
                inter = set.intersection(left_ss, right_ss)
                if inter:
                    result.append(inter)
                else:
                    if weight_list is None:
                        wt = 1
                    else:
                        wt = weight_list[n]
                    score += wt
                    result.append(set.union(left_ss, right_ss))
                #_LOG.debug("left = %s, right = %s, nd= %s" %
                #                (str(left_ss), str(right_ss), str(result)))
            if remaining:
                right_c = remaining.pop(0)
                left_ssl = result
            else:
                break

        setattr(nd, attr_name, result)
    return score

def fitch_up_pass(preorder_node_list, attr_name="state_sets", taxa_to_state_set_map=None):
    """
    Reads `attr_name` attribute of nodes as an iterable of state sets, and
    sets that attribute for internal nodes using the "final phase" of Fitch's
    (1971) unordered parsimony algorithm.

        `taxa_to_state_set_map`
            if a child node does not have an attribute with name
        `attr_name`
            then the nodes.taxon will be used as a key in taxa_to_state_set_map
            to find the state set. This allows for the scoring of
            previously undecorated trees.

    Currently this requires a bifurcating tree (even at the root).
    """
    for nd in preorder_node_list:
        c = nd.child_nodes()
        p = nd.parent_node
        if (not c) or (not p):
            continue
        assert(len(c) == 2)
        left_c, right_c = c
        try:
            left_ssl = getattr(left_c, attr_name)
        except AttributeError:
            if not taxa_to_state_set_map:
                raise
            left_ssl = taxa_to_state_set_map[left_c.taxon]
        try:
            right_ssl = getattr(right_c, attr_name)
        except AttributeError:
            if not taxa_to_state_set_map:
                raise
            right_ssl = taxa_to_state_set_map[right_c.taxon]
        par_ssl = getattr(p, attr_name)
        curr_ssl = getattr(nd, attr_name)
        result = []
        for n, ssp in enumerate(izip(par_ssl, curr_ssl, left_ssl, right_ssl)):
            par_ss, curr_ss, left_ss, right_ss = ssp

            down_parup_inter = set.intersection(par_ss, curr_ss)
            if down_parup_inter == par_ss:
                final_ss = down_parup_inter
            else:
                rl_inter = set.intersection(left_ss, right_ss)
                if not rl_inter:
                    final_ss = set.union(par_ss, curr_ss)
                else:
                    in_par_and_left = set.intersection(par_ss, left_ss)
                    in_par_and_right = set.intersection(par_ss, right_ss)
                    final_ss = set.union(in_par_and_left, in_par_and_right, curr_ss)
            #_LOG.debug("downpass = %s, par = %s, left = %s, right = %s, final_ss= %s" %
            #                    (str(curr_ss), str(par_ss), str(left_ss), str(right_ss), str(final_ss)))
            result.append(final_ss)
        setattr(nd, attr_name, result)
