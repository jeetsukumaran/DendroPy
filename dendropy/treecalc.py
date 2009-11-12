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

def get_mrca(start_node, split, taxa_mask):
    """Returns the shallowest node in the tree (the node furthest from
    `start_node`) that has all of the taxa that are specified in `split` or
    None if no appropriate node is found.

    Assumes that edges on tree have been decorated with splitmask.

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

def patristic_distance(tree, taxon1, taxon2):
    """
    Given a tree with splits encoded, and two taxa on that tree, returns the
    patristic distance between the two.
    """
    split = tree.taxon_set.taxon_bitmask(taxon1) | tree.taxon_set.taxon_bitmask(taxon2)
    mrca = get_mrca(tree.seed_node, split, tree.taxon_set.all_taxa_bitmask())
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

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).

    If the lengths of different paths to the node differ by more than `prec`,
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.

    As a side effect a `depth` attribute is added to the nodes of the tree.

    Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
    molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
    (267). 2267-2272
    """
    # the equation is given by:
    #   T = \sum_{j=2}^n (jg_j)
    #   C = T \sqrt{\frac{1}{12(n-2)}}
    #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
    # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
    #   times between consecutive (in time, not along a branch) speciation times.
    node = None
    speciation_depths = []
    n = 0
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        n_ch = len(ch)
        if n_ch == 0:
            node.depth = 0.0
            n += 1
        elif n_ch > 2:
            raise ValueError("Polytomy encountered")
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    raise ValueError("Tree is not ultrametric")
            if n_ch == 2:
                speciation_depths.append(node.depth)
    if node is None:
        raise ValueError("Empty tree encountered")
    speciation_depths.sort(reverse=True)
    g = []
    older = speciation_depths[0]
    for age in speciation_depths[1:]:
        g.append(older - age)
        older = age
    g.append(older)
    if not g:
        raise ValueError("No internal nodes found (other than the root)")
    assert(len(g) == (n - 1))
    T = 0.0
    accum = 0.0
    for i in xrange(2, n):
        list_index = i - 2
        T += i * float(g[list_index])
        accum += T
    list_index = n - 2
    T += (n) * g[list_index]
    nmt = n - 2.0
    numerator = accum/nmt - T/2.0
    C = T*pow(1/(12*nmt), 0.5)
    return numerator/C

