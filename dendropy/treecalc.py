#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Tree metrics/statistics calculations.
"""

from math import sqrt
from dendropy import treesplit

class PatristicDistanceMatrix(object):
    """
    Calculates and maintains patristic distance information of taxa on a tree.
    """

    def __init__(self, tree=None):
        self.tree = None
        self.taxon_set = None
        self._pat_dists = {}
        if tree is not None:
            self.calc(tree)

    def __call__(self, taxon1, taxon2):
        """
        Returns patristic distance between two taxon objects.
        """
        try:
            return self._pat_dists[taxon1][taxon2]
        except KeyError, e:
            return self._pat_dists[taxon2][taxon1]

    def calc(self, tree):
        """
        Calculates the distances.
        """
        self.tree = tree
        if not hasattr(self.tree, "split_edges"):
            treesplit.encode_splits(self.tree)
        self.taxon_set = tree.taxon_set
        self._pat_dists = {}
        for i1, t1 in enumerate(self.taxon_set):
            self._pat_dists[t1] = {}

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
                                self._pat_dists[desc1.taxon][desc2.taxon] = \
                                    node.desc_paths[desc1] + desc2_plen + c2.edge.length
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

def find_mrca(tree, **kwargs):
    """Returns the MRCA of all taxa that:

        - are specified by the split bitmask given by the keyword argument
         `split_bitmask`
        - are in the list of Taxon objects given by the keyword argument
          'taxa'
        - have the labels specified by the list of strings given by the
          keyword argument 'taxon_labels'

    If `taxa` or `taxon_labels` are used, then a TaxonSet object must be
    given by the keyword `taxon_set`.
    Returns None if no appropriate node is found.
    Assumes that edges on tree have been decorated with treesplit.
    """
    return find_mrca_from_node(tree.seed_node, taxon_set=tree.taxon_set, **kwargs)

def find_mrca_from_node(start_node, **kwargs):
    """Returns the shallowest node in the tree (the node furthest from
    `start_node` in a direction away from the root) that has
    all of the taxa that:

        - are specified by the split bitmask given by the keyword argument
         `split_bitmask`
        - are in the list of Taxon objects given by the keyword argument
          'taxa'
        - have the labels specified by the list of strings given by the
          keyword argument 'taxon_labels'

    If `taxa` or `taxon_labels` are used, then a TaxonSet object must be
    given by the keyword `taxon_set`.
    Returns None if no appropriate node is found.
    Assumes that edges on tree have been decorated with treesplit.
    It is possible that split is not compatible with the subtree that is
        returned! (compatibility tests are not fully performed).
    This function is used to find the "insertion point" for a new split via a
        root to tip search.
    """

    split_bitmask = None

    if "split_bitmask" in kwargs:
        split_bitmask = kwargs["split_bitmask"]
    else:
        if "taxon_set" not in kwargs:
            raise TypeError("Must specify TaxonSet using keyword 'taxon_set' if" \
                + " specifying 'taxa' or 'taxon_labels' arguments")
        taxon_set = kwargs["taxon_set"]
        taxa = kwargs.get("taxa", None)
        if taxa is None:
            if "taxon_labels" in kwargs:
                taxa = taxon_set.get_taxa(labels=kwargs["taxon_labels"])
                if len(taxa) != len(kwargs["taxon_labels"]):
                    raise KeyError("Not all labels matched to taxa")
            else:
                raise TypeError("Must specify one of: 'split_bitmask', 'taxa' or 'taxon_labels'")
        if taxa is None:
            raise ValueError("No taxa matching criteria found")
        split_bitmask = taxon_set.get_taxa_bitmask(taxa=taxa)

    if split_bitmask is None or split_bitmask == 0:
        raise ValueError("Null split bitmask (0)")

    if not hasattr(start_node.edge, "split_bitmask"):
        raise ValueError("Splits not encoded")

    if (start_node.edge.split_bitmask & split_bitmask) != split_bitmask:
        return None

    curr_node = start_node
    last_match = start_node
    nd_source = iter(start_node.child_nodes())
    try:
        while True:
            cm = curr_node.edge.split_bitmask
            cms = (cm & split_bitmask)
            if cms:
                # for at least one taxon cm has 1 and split has 1
                if cms == split_bitmask:
                    # curr_node has all of the 1's that split has
                    if cm == split_bitmask:
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
        #   decorated with split_bitmask attributes, but there may be some hacky
        #   context in which we want to allow the function to be called with
        #   leaves that have not been encoded with split_bitmasks.
        return last_match

def patristic_distance(tree, taxon1, taxon2):
    """
    Given a tree with splits encoded, and two taxa on that tree, returns the
    patristic distance between the two.
    """
    if not hasattr(tree, "split_edges"):
        treesplit.encode_splits(tree)
    mrca = find_mrca_from_node(tree.seed_node, taxa=[taxon1, taxon2], taxon_set=tree.taxon_set)
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
    length_diffs = []
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
        except KeyError: # excep
            elen2 = 0
        if elen2 is None:
            elen2 = 0 # worst-case: bind
        value2 = value_type(elen2) #  ctor + bind # best case
        length_diffs.append((value1,value2)) # ctor + listappend

    for split, edge in split_edges2_copy.iteritems(): # best-case not executed, worst case O(n) : 2*bind
        elen2 = getattr(edge, edge_length_attr) # attr +  bind
        if elen2 is None:
            elen2 = 0
        value2 = value_type(elen2) #  ctor + bind
        e1 = split_edges1_ref.get(split) # attr + dict_lookup + bind
        if e1 is None:
            elen1 = 0
        else:
            elen1 = getattr(e1, edge_length_attr) # attr  + bind
        if elen1 is None:
            elen1 = 0
        value1 = value_type(elen1)
        length_diffs.append((value1,value2)) # ctor + listappend
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



