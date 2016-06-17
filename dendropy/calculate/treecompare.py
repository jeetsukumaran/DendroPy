#! /usr/bin/env python

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
Statistics, metrics, measurements, and values calculated *between* *two* trees.
"""

import math
import collections
import itertools
import dendropy
from dendropy.utility import error

###############################################################################
## Public Functions

def symmetric_difference(tree1, tree2, is_bipartitions_updated=False):
    """
    Returns *unweighted* Robinson-Foulds distance between two trees.

    Trees need to share the same |TaxonNamespace| reference. The
    bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
    ``is_bipartitions_updated`` argument must be |False| to force recalculation
    of bipartitions.

    Parameters
    ----------
    tree1 : |Tree| object
        The first tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree2`` and must have
        bipartitions encoded.
    tree2 : |Tree| object
        The second tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree1`` and must have
        bipartitions encoded.
    is_bipartitions_updated : bool
        If |False|, then the bipartitions on *both* trees will be updated
        before comparison. If |True| then the bipartitions will only be
        calculated for a |Tree| object if they have not been calculated
        before, either explicitly or implicitly.

    Returns
    -------
    d : int
        The symmetric difference (a.k.a. the unweighted Robinson-Foulds
        distance) between ``tree1`` and ``tree2``.

    Examples
    --------

    ::

        import dendropy
        from dendropy.calculate import treecompare
        tns = dendropy.TaxonNamespace()
        tree1 = tree.get_from_path(
                "t1.nex",
                "nexus",
                taxon_namespace=tns)
        tree2 = tree.get_from_path(
                "t2.nex",
                "nexus",
                taxon_namespace=tns)
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
        print(treecompare.symmetric_difference(tree1, tree2))

    """
    t = false_positives_and_negatives(
            tree1,
            tree2,
            is_bipartitions_updated=is_bipartitions_updated)
    return t[0] + t[1]

def unweighted_robinson_foulds_distance(tree1, tree2, is_bipartitions_updated=False):
    """
    Alias for ``symmetric_difference()``.
    """
    return symmetric_difference(tree1, tree2, is_bipartitions_updated)

def weighted_robinson_foulds_distance(
        tree1,
        tree2,
        edge_weight_attr="length",
        is_bipartitions_updated=False):
    """
    Returns *weighted* Robinson-Foulds distance between two trees based on
    ``edge_weight_attr``.

    Trees need to share the same |TaxonNamespace| reference. The
    bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
    ``is_bipartitions_updated`` argument must be |False| to force recalculation of
    bipartitions.

    Parameters
    ----------
    tree1 : |Tree| object
        The first tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree2`` and must have
        bipartitions encoded.
    tree2 : |Tree| object
        The second tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree1`` and must have
        bipartitions encoded.
    edge_weight_attr : string
        Name of attribute on edges of trees to be used as the weight.
    is_bipartitions_updated : bool
        If |True|, then the bipartitions on *both* trees will be updated before
        comparison. If |False| (default) then the bipartitions will only be
        calculated for a |Tree| object if they have not been calculated
        before, either explicitly or implicitly.

    Returns
    -------
    d : float
        The edge-weighted Robinson-Foulds distance between ``tree1`` and ``tree2``.

    Examples
    --------

    ::

        import dendropy
        from dendropy.calculate import treecompare
        tns = dendropy.TaxonNamespace()
        tree1 = tree.get_from_path(
                "t1.nex",
                "nexus",
                taxon_namespace=tns)
        tree2 = tree.get_from_path(
                "t2.nex",
                "nexus",
                taxon_namespace=tns)
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
        print(treecompare.weighted_robinson_foulds_distance(tree1, tree2))

    """
    df = lambda length_diffs: sum([abs(i[0] - i[1]) for i in length_diffs])
    return _bipartition_difference(tree1,
                           tree2,
                           dist_fn=df,
                           edge_weight_attr=edge_weight_attr,
                           value_type=float,
                           is_bipartitions_updated=is_bipartitions_updated)

def false_positives_and_negatives(
        reference_tree,
        comparison_tree,
        is_bipartitions_updated=False):
    """
    Counts and returns number of false positive bipar (bipartitions found in
    ``comparison_tree`` but not in ``reference_tree``) and false negative
    bipartitions (bipartitions found in ``reference_tree`` but not in
    ``comparison_tree``).

    Trees need to share the same |TaxonNamespace| reference. The
    bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
    ``is_bipartitions_updated`` argument must be |False| to force recalculation of
    bipartitions.

    Parameters
    ----------
    reference_tree : |Tree| object
        The first tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree2`` and must have
        bipartitions encoded.
    comparison_tree : |Tree| object
        The second tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree1`` and must have
        bipartitions encoded.
    is_bipartitions_updated : bool
        If |True|, then the bipartitions on *both* trees will be updated
        before comparison. If |False| (default) then the bipartitions
        will only be calculated for a |Tree| object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    t : tuple(int)
        A pair of integers, with first integer being the number of false
        positives and the second being the number of false negatives.

    Examples
    --------

    ::

        import dendropy
        from dendropy.calculate import treecompare
        tns = dendropy.TaxonNamespace()
        tree1 = tree.get_from_path(
                "t1.nex",
                "nexus",
                taxon_namespace=tns)
        tree2 = tree.get_from_path(
                "t2.nex",
                "nexus",
                taxon_namespace=tns)
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
        print(treecompare.false_positives_and_negatives(tree1, tree2))

    """
    if reference_tree.taxon_namespace is not comparison_tree.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(reference_tree, comparison_tree)
    if not is_bipartitions_updated:
        reference_tree.encode_bipartitions()
        comparison_tree.encode_bipartitions()
    else:
        if reference_tree.bipartition_encoding is None:
            reference_tree.encode_bipartitions()
        if comparison_tree.bipartition_encoding is None:
            comparison_tree.encode_bipartitions()
    ref_bipartitions = set(reference_tree.bipartition_encoding)
    comparison_bipartitions = set(comparison_tree.bipartition_encoding)
    false_positives = comparison_bipartitions.difference(ref_bipartitions)
    false_negatives = ref_bipartitions.difference(comparison_bipartitions)
    return len(false_positives), len(false_negatives)

def euclidean_distance(
        tree1,
        tree2,
        edge_weight_attr="length",
        value_type=float,
        is_bipartitions_updated=False):
    """
    Returns the Euclidean distance (a.k.a. Felsenstein's 2004 "branch length
    distance") between two trees based on ``edge_weight_attr``.

    Trees need to share the same |TaxonNamespace| reference. The
    bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
    ``is_bipartitions_updated`` argument must be |False| to force recalculation of
    bipartitions.

    Parameters
    ----------
    tree1 : |Tree| object
        The first tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree2`` and must have
        bipartitions encoded.
    tree2 : |Tree| object
        The second tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree1`` and must have
        bipartitions encoded.
    edge_weight_attr : string
        Name of attribute on edges of trees to be used as the weight.
    is_bipartitions_updated : bool
        If |True|, then the bipartitions on *both* trees will be updated
        before comparison. If |False| (default) then the bipartitions
        will only be calculated for a |Tree| object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    d : int
        The Euclidean distance between ``tree1`` and ``tree2``.

    Examples
    --------

    ::

        import dendropy
        from dendropy.calculate import treecompare
        tns = dendropy.TaxonNamespace()
        tree1 = tree.get_from_path(
                "t1.nex",
                "nexus",
                taxon_namespace=tns)
        tree2 = tree.get_from_path(
                "t2.nex",
                "nexus",
                taxon_namespace=tns)
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
        print(treecompare.euclidean_distance(tree1, tree2))

    """
    df = lambda length_diffs: math.sqrt(sum([pow(i[0] - i[1], 2) for i in length_diffs]))
    return _bipartition_difference(tree1,
                           tree2,
                           dist_fn=df,
                           edge_weight_attr=edge_weight_attr,
                           value_type=value_type,
                           is_bipartitions_updated=is_bipartitions_updated)

def find_missing_bipartitions(reference_tree, comparison_tree, is_bipartitions_updated=False):
    """
    Returns a list of bipartitions that are in ``reference_tree``, but
    not in ``comparison_tree``.

    Trees need to share the same |TaxonNamespace| reference. The
    bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
    ``is_bipartitions_updated`` argument must be |False| to force recalculation of
    bipartitions.

    Parameters
    ----------
    reference_tree : |Tree| object
        The first tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree2`` and must have
        bipartitions encoded.
    comparison_tree : |Tree| object
        The second tree of the two trees being compared. This must share the
        same |TaxonNamespace| reference as ``tree1`` and must have
        bipartitions encoded.
    is_bipartitions_updated : bool
        If |True|, then the bipartitions on *both* trees will be updated
        before comparison. If |False| (default) then the bipartitions
        will only be calculated for a |Tree| object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    s : list[|Bipartition|]
        A list of bipartitions that are in the first tree but not in the second.

    """
    missing = []
    if reference_tree.taxon_namespace is not comparison_tree.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(reference_tree, comparison_tree)
    if not is_bipartitions_updated:
        reference_tree.encode_bipartitions()
        comparison_tree.encode_bipartitions()
    else:
        if reference_tree.bipartition_encoding is None:
            reference_tree.encode_bipartitions()
        if comparison_tree.bipartition_encoding is None:
            comparison_tree.encode_bipartitions()
    for bipartition in reference_tree.bipartition_encoding:
        if bipartition in comparison_tree.bipartition_encoding:
            pass
        else:
            missing.append(bipartition)
    return missing

##############################################################################
### TreeshapeKernel

class TreeShapeKernel(object):

    _TreeShapeKernelNodeCache = collections.namedtuple("_TreeShapeKernelNodeCache",
            ["production", "index", "edge_lengths", "sum_of_square_edge_lengths"])

    def __init__(self, **kwargs):
        """
        Calculator for tree shape kernel tricking.

        References
        ----------

        [1] Poon, A. F., Pond, S. L. K., Bennett, P., Richman, D. D., Brown, A.
        J. L., & Frost, S. D. (2007). Adaptation to human populations is
        revealed by within-host polymorphisms in HIV-1 and hepatitis C virus.
        PLoS Pathog, 3(3), e45.

        [2] Poon, A. F., Walker, L. W., Murray, H., McCloskey, R. M., Harrigan,
        P. R., & Liang, R. H. (2013). Mapping the shapes of phylogenetic trees
        from human and zoonotic RNA viruses. PLoS one, 8(11), e78122.

        [3] Poon, A. F., Walker, L. W., Murray, H., McCloskey, R. M., Harrigan,
        P. R., & Liang, R. H. (2013). Mapping the shapes of phylogenetic trees
        from human and zoonotic RNA viruses. PLoS one, 8(11), e78122.

        [4] Poon, A. F. (2015). Phylodynamic inference with kernel ABC and its
        application to HIV epidemiology. Molecular biology and evolution,
        msv123.

        """

        # kernel function
        # sigma=1,
        # gauss_factor=1,
        # decay_factor=0.1,
        self.sigma = kwargs.pop("sigma", 1)
        self.gauss_factor = kwargs.pop("gauss_factor", 1)
        self.decay_factor = kwargs.pop("decay_factor", 0.1)

        # cache management
        self._tree_cache = {}

    def remove_from_cache(self, tree):
        del self._tree_cache[tree]

    def update_cache(self, tree):
        """
        Pre-computes values needed for the kernel trick with this tree and
        caches them.
        """
        current_tree_cache = {}
        for nd in tree.leaf_node_iter():
            current_tree_cache[nd] = TreeShapeKernel._TreeShapeKernelNodeCache(
                    production=0,
                    index=0,
                    edge_lengths=None,
                    sum_of_square_edge_lengths=0)
            nd.production = 0
        for nd_idx, nd in enumerate(tree.postorder_internal_node_iter()):
            nterms = 0
            edge_lengths = []
            for ch in nd.child_node_iter():
                if current_tree_cache[ch].production == 0:
                    nterms += 1
                edge_lengths.append(ch.edge.length)
            production = nterms + 1
            index = nd_idx
            sum_of_square_edge_lengths = sum([elen**2 for elen in edge_lengths])
            current_tree_cache[nd] = TreeShapeKernel._TreeShapeKernelNodeCache(
                    production=production,
                    index=index,
                    edge_lengths=edge_lengths,
                    sum_of_square_edge_lengths=sum_of_square_edge_lengths)
        self._tree_cache[tree] = current_tree_cache
        return current_tree_cache

    def __call__(self,
            tree1,
            tree2,
            is_tree1_cache_updated=False,
            is_tree2_cache_updated=False,
            ):
        """
        Recursive function for computing tree convolution
        kernel.

        Parameters
        ----------
        tree1 : |Tree| instance
            First tree to be compared. If it has already been seen by self, its
            values will have been cached. If the tree has changed since it has been
            seen by self, it will need to be recached, either explicitly before
            the calculation by calling 'update_cache' or by specifying
            'is_tree1_cache_updated=False'
        tree2 : |Tree| instance
            Second tree to be compared. If it has already been seen by self, its
            values will have been cached. If the tree has changed since it has been
            seen by self, it will need to be recached, either explicitly before
            the calculation by calling 'update_cache' or by specifying
            'is_tree2_cache_updated=True'
        is_tree1_cache_updated : bool
            If ``tree1`` has not been seen before, then this is ignored as the
            cache will be updated regardless. If ``tree1`` has been seen, then
            the cached values representing it will be used unless
            ``is_tree1_cache_updated`` is ``False``.
        is_tree2_cache_updated : bool
            If ``tree1`` has not been seen before, then this is ignored as the
            cache will be updated regardless. If ``tree1`` has been seen, then
            the cached values representing it will be used unless
            ``is_tree1_cache_updated`` is ``False``.

        Acknowledgements
        ----------------

        Based in part on:

            KAMPHIR
            https://github.com/ArtPoon/kamphir.git
            KAMPHIR is written and maintained by: Art F.Y. Poon.
            With major contributions from: Rosemary McCloskey

        Copyright (c) 2015, Art Poon. All rights reserved.
        See https://github.com/ArtPoon/kamphir/blob/master/LICENSE.md for
        more license information.

        Original work adapted from Moschitti (2006) Making tree kernels
        practical for natural language learning. Proceedings of the
        11th Conference of the European Chapter of the Association
        for Computational Linguistics.
        """
        internal_nodes2 = list(tree2.postorder_internal_node_iter())
        k = 0
        if not is_tree1_cache_updated:
            tree1_cache = self.update_cache(tree1)
        else:
            try:
                tree1_cache = self._tree_cache[tree1]
            except KeyError:
                tree1_cache = self.update_cache(tree1)
        if not is_tree2_cache_updated:
            tree2_cache = self.update_cache(tree2)
        else:
            try:
                tree2_cache = self._tree_cache[tree2]
            except KeyError:
                tree2_cache = self.update_cache(tree2)
        dp_matrix = {}
        for ni, tree1_node in enumerate(tree1.postorder_internal_node_iter()):
            tree1_cache_node = tree1_cache[tree1_node]
            for tree2_node in internal_nodes2:
                tree2_cache_node = tree2_cache[tree2_node]
                if tree1_cache_node.production != tree2_cache_node.production:
                    continue
                res = self.decay_factor * math.exp( -1. / self.gauss_factor
                    * (tree1_cache_node.sum_of_square_edge_lengths + tree2_cache_node.sum_of_square_edge_lengths - 2*sum([(tree1_cache_node.edge_lengths[i]*tree2_cache_node.edge_lengths[i]) for i in range(len(tree1_cache_node.edge_lengths))])))
                ## TODO:
                ##  - (check and) handles cases where unequal number of children
                ##  - how to handle rotation mismatch problems? or do we assume
                ##    trees have equal rotations
                # for node_idx in range(2):
                #     c1 = tree1_node.clades[node_idx]
                #     c2 = tree2_node.clades[node_idx]
                for c1, c2 in zip(tree1_node.child_node_iter(), tree2_node.child_node_iter()):
                    if tree1_cache[c1].production != tree2_cache[c2].production:
                        continue
                    if tree1_cache[c1].production == 0:
                        # branches are terminal
                        res *= self.sigma + self.decay_factor
                    else:
                        try:
                            res *= self.sigma + dp_matrix[(tree1_cache[c1].index, tree2_cache[c2].index)]
                        except KeyError:
                            res *= self.sigma
                dp_matrix[(tree1_cache[tree1_node].index, tree2_cache[tree2_node].index)] = res
                k += res
        return k

##############################################################################
### AssemblageInducedTree

class AssemblageInducedTreeManager(object):

    def __init__(self, *args, **kwargs):
        self.is_exchangeable_assemblage_classifications = kwargs.pop("is_exchangeable_assemblage_classifications", True)
        self._num_assemblage_classifications = kwargs.pop("num_assemblages", None)
        self.induced_tree_factory = kwargs.pop("induced_tree_factory", None)
        self.induced_tree_node_factory = kwargs.pop("induced_tree_node_factory", None)
        self.skip_null_assemblages = kwargs.pop("skip_null_assemblages", False)
        self._tree_assemblage_induced_trees_map = {}

    def remove_from_cache(self, tree):
        del self._tree_assemblage_induced_trees_map[tree]

    def generate_induced_trees(self, tree, assemblage_leaf_sets):
        if assemblage_leaf_sets is None:
            if self._num_assemblage_classifications is not None:
                raise ValueError("Expecting {} assemblage leaf set classifications, but none provided".format(self._num_assemblage_classifications))
        else:
            if self._num_assemblage_classifications is None:
                self._num_assemblage_classifications = len(assemblage_leaf_sets)
            elif len(assemblage_leaf_sets) != self._num_assemblage_classifications:
                    raise ValueError("Expecting {} assemblage leaf set classifications, but only {} specified".format(
                        self._num_assemblage_classifications,
                        len(assemblage_leaf_sets)))
        induced_trees = []
        for idx, assemblage_leaf_set in enumerate(assemblage_leaf_sets):
            if len(assemblage_leaf_set) == 0:
                if self.skip_null_assemblages:
                    continue
                raise error.NullLeafSetException()
            node_filter_fn = lambda nd: nd in assemblage_leaf_set
            induced_tree = tree.extract_tree(
                               node_filter_fn=node_filter_fn,
                               is_apply_filter_to_leaf_nodes=True,
                               is_apply_filter_to_internal_nodes=False,
                               tree_factory=self.induced_tree_factory,
                               node_factory=self.induced_tree_node_factory)
            induced_trees.append(induced_tree)
        self._tree_assemblage_induced_trees_map[tree] = induced_trees
        return induced_trees

##############################################################################
### AssemblageInducedTreeShapeKernel

class AssemblageInducedTreeShapeKernel(TreeShapeKernel, AssemblageInducedTreeManager):

    @staticmethod
    def _euclidean_distance(v1, v2, is_weight_values_by_comparison_size=True):
        v1_size = len(v1)
        v2_size = len(v2)
        v1_idx = 0
        v2_idx = 0
        if v1_size > v2_size:
            v1_idx = v1_size - v2_size
            weight = float(v2_size)
        elif v2_size > v1_size:
            v2_idx = v2_size - v1_size
            weight = float(v1_size)
        else:
            weight = float(v1_size)
        if not is_weight_values_by_comparison_size:
            weight = 1.0
        ss = 0.0
        while v1_idx < v1_size and v2_idx < v2_size:
            ss += pow(v1[v1_idx]/weight - v2[v2_idx]/weight, 2)
            v1_idx += 1
            v2_idx += 1
        return math.sqrt(ss)

    def __init__(self, *args, **kwargs):
        self.exchangeable_assemblage_comparison_strategy = kwargs.pop("exchangeable_assemblage_comparison_strategy", "joint minimum")
        TreeShapeKernel.__init__(self, *args, **kwargs)
        AssemblageInducedTreeManager.__init__(self, *args, **kwargs)

    def remove_from_cache(self, tree):
        for induced_tree in self._tree_assemblage_induced_trees_map[tree]:
            TreeShapeKernel.remove_from_cache(self, induced_tree)
        TreeShapeKernel.remove_from_cache(self, tree)
        AssemblageInducedTreeManager.remove_from_cache(self, tree)

    def update_assemblage_induced_tree_cache(self,
            tree,
            assemblage_leaf_sets):
        self.update_cache(tree=tree)
        induced_trees = self.generate_induced_trees(tree=tree,
                assemblage_leaf_sets=assemblage_leaf_sets)
        for induced_tree in induced_trees:
            self.update_cache(tree=induced_tree)

    def __call__(self,
            tree1,
            tree2,
            tree1_assemblage_leaf_sets,
            tree2_assemblage_leaf_sets,
            is_tree1_cache_updated=False,
            is_tree2_cache_updated=False,
            ):
        main_trees_score = TreeShapeKernel.__call__(self,
                tree1=tree1,
                tree2=tree2,
                is_tree1_cache_updated=is_tree1_cache_updated,
                is_tree2_cache_updated=is_tree2_cache_updated,
                )
        if not is_tree1_cache_updated or tree1 not in self._tree_assemblage_induced_trees_map:
            if tree1_assemblage_leaf_sets is None:
                raise ValueError("Uncached tree requires specification of 'tree1_assemblage_leaf_sets'")
            self.update_assemblage_induced_tree_cache(
                    tree=tree1,
                    assemblage_leaf_sets=tree1_assemblage_leaf_sets)
        if not is_tree2_cache_updated or tree2 not in self._tree_assemblage_induced_trees_map:
            if tree2_assemblage_leaf_sets is None:
                raise ValueError("Uncached tree requires specification of 'tree2_assemblage_leaf_sets'")
            self.update_assemblage_induced_tree_cache(
                    tree=tree2,
                    assemblage_leaf_sets=tree2_assemblage_leaf_sets)
        ## ++ main tree score
        score_table = collections.OrderedDict()
        score_table["primary.tree.kernel.trick.distance"] = main_trees_score
        induced_trees1 = self._tree_assemblage_induced_trees_map[tree1]
        induced_trees2 = self._tree_assemblage_induced_trees_map[tree2]
        # assert len(induced_trees1) == len(induced_trees2) == self._num_assemblage_classifications
        if not self.is_exchangeable_assemblage_classifications:
            if len(induced_trees1) != len(induced_trees2):
                raise TypeError("Different numbers of induced trees not supported for non-exchangeable classifications: {} vs. {}".format(len(induced_trees1), len(induced_trees2)))
            for idx, (induced_tree1, induced_tree2) in enumerate(zip(induced_trees1, induced_trees2)):
                s = TreeShapeKernel.__call__(self,
                                tree1=induced_tree1,
                                tree2=induced_tree2,
                                is_tree1_cache_updated=True,
                                is_tree2_cache_updated=True,
                                )
                ## ++ raw scores direct comparisons of each of the induced trees
                score_table["induced.tree.{}.kernel.trick.distance".format(idx+1)] = s
        else:
            if self.exchangeable_assemblage_comparison_strategy == "joint minimum":
                # if lengths are different, we want to fix the smaller set
                if len(induced_trees1) > len(induced_trees2):
                    induced_trees2, induced_trees1 = induced_trees1, induced_trees2
                comparison_vector = [0.0] * len(induced_trees1)
                current_minimum_distance = None
                current_joint_minimum_vector = None
                for induced_trees_permutation in itertools.permutations(induced_trees2, len(induced_trees1)):
                    distances = []
                    for t2, t1 in zip(induced_trees_permutation, induced_trees1):
                        distances.append(TreeShapeKernel.__call__(self,
                                tree1=t1,
                                tree2=t2,
                                is_tree1_cache_updated=True,
                                is_tree2_cache_updated=True,))
                    euclidean_distance = self._euclidean_distance(distances, comparison_vector)
                    if current_minimum_distance is None or euclidean_distance < current_minimum_distance:
                        current_minimum_distance = euclidean_distance
                        current_joint_minimum_vector = distances
                for didx, d in enumerate(distances):
                    score_table["induced.tree.{}.kernel.trick.distance".format(didx+1)] = d
                for didx in range(didx+1, self._num_assemblage_classifications):
                    score_table["induced.tree.{}.kernel.trick.distance".format(didx+1)] = "NA"
            else:
                raise NotImplementedError()
        return score_table

###############################################################################
## Legacy

def robinson_foulds_distance(tree1, tree2, edge_weight_attr="length"):
    """
    DEPRECATED: Use :func:``symmetric_difference`` for the common
    unweighged Robinson-Fould's distance metric (i.e., the symmetric difference between two trees)
    :func:``weighted_robinson_foulds_distance`` or for the RF distance as defined by Felsenstein, 2004.
    """
    return weighted_robinson_foulds_distance(tree1, tree2, edge_weight_attr)

def mason_gamer_kellogg_score(tree1, tree2, is_bipartitions_updated=False):
    """
    Mason-Gamer and Kellogg. Testing for phylogenetic conflict among molecular
    data sets in the tribe Triticeae (Gramineae). Systematic Biology (1996)
    vol. 45 (4) pp. 524
    """
    if tree1.taxon_namespace is not tree2.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(tree1, tree2)
    if not is_bipartitions_updated:
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
    else:
        if tree1.bipartition_encoding is None:
            tree1.encode_bipartitions()
        if tree2.bipartition_encoding is None:
            tree2.encode_bipartitions()
    se1 = tree1.bipartition_encoding
    se2 = tree2.bipartition_encoding
    bipartitions = sorted(list(set(se1.keys() + se2.keys())))

###############################################################################
## Supporting

def _get_length_diffs(
        tree1,
        tree2,
        edge_weight_attr="length",
        value_type=float,
        is_bipartitions_updated=False,
        bipartition_length_diff_map=False):
    """
    Returns a list of tuples, with the first element of each tuple representing
    the length of the branch subtending a particular bipartition on ``tree1``, and
    the second element the length of the same branch on ``tree2``. If a
    particular bipartition is found on one tree but not in the other, a value of zero
    is used for the missing bipartition.
    """
    length_diffs = []
    bipartition_length_diffs = {}
    if tree1.taxon_namespace is not tree2.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(tree1, tree2)
    if not is_bipartitions_updated:
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
    else:
        if tree1.bipartition_encoding is None:
            tree1.encode_bipartitions()
        if tree2.bipartition_encoding is None:
            tree2.encode_bipartitions()

    tree1_bipartition_edge_map = dict(tree2.bipartition_edge_map) # O(n*(2*bind + dict_item_cost))
    tree2_bipartition_edge_map = tree1.bipartition_edge_map
    for bipartition in tree2_bipartition_edge_map: # O n : 2*bind
        edge = tree2_bipartition_edge_map[bipartition]
        elen1 = getattr(edge, edge_weight_attr) # attr + bind
        if elen1 is None:
            elen1 = 0 # worst-case: bind
        value1 = value_type(elen1) #  ctor + bind
        try:
            e2 = tree1_bipartition_edge_map.pop(bipartition) # attr + dict_lookup + bind
            elen2 = getattr(e2, edge_weight_attr) # attr + bind
            if elen2 is None:
                # allow root edge to have bipartition with no value: raise error if not root edge
                if e2.tail_node is None:
                    elen2 = 0.0
                else:
                    raise ValueError("Edge length attribute is 'None': Tree: %s ('%s'), Split: %s" % (id(tree2), tree2.label, bipartition.leafset_as_newick_string(tree2.taxon_namespace)))
        except KeyError: # excep
            elen2 = 0.0
        value2 = value_type(elen2) #  ctor + bind # best case
        # if abs(value2-value1) > 1e-5:
        #     print("{}: {}, {}".format(bipartition.leafset_as_newick_string(tree1.taxon_namespace), value2, value1))
        length_diffs.append((value1,value2)) # ctor + listappend
        bipartition_length_diffs[bipartition] = length_diffs[-1]

    for bipartition in tree1_bipartition_edge_map: # best-case not executed, worst case O(n) : 2*bind
        edge = tree1_bipartition_edge_map[bipartition]
        elen2 = getattr(edge, edge_weight_attr) # attr +  bind
        if elen2 is None:
            elen2 = 0
        value2 = value_type(elen2) #  ctor + bind
        e1 = tree2_bipartition_edge_map.get(bipartition) # attr + dict_lookup + bind
        if e1 is None:
            elen1 = 0.0
        else:
            elen1 = getattr(e1, edge_weight_attr) # attr  + bind
            if elen1 is None:
                # allow root edge to have bipartition with no value: raise error if not root edge
                if e1.tail_node is None:
                    elen1 = 0.0
                else:
                    raise ValueError("Edge length attribute is 'None': Tree: %s ('%s'), Split: %s" % (id(tree1), tree1.label, bipartition))
                #elen1 = 0
        value1 = value_type(elen1)
        length_diffs.append((value1,value2)) # ctor + listappend
        bipartition_length_diffs[bipartition] = length_diffs[-1]

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
    if bipartition_length_diff_map:
        return length_diffs, bipartition_length_diffs
    else:
        return length_diffs

def _bipartition_difference(
        tree1,
        tree2,
        dist_fn,
        edge_weight_attr="length",
        value_type=float,
        is_bipartitions_updated=False):
    """
    Returns distance between two trees, each represented by a dictionary of
    bipartitions (as bipartition_mask strings) to edges, using ``dist_fn`` to calculate the
    distance based on ``edge_weight_attr`` of the edges. ``dist_fn`` is a function
    that takes a list of pairs of values, where the values correspond to the edge
    lengths of a given bipartition on tree1 and tree2 respectively.
    """
    length_diffs = _get_length_diffs(
            tree1,
            tree2,
            edge_weight_attr=edge_weight_attr,
            value_type=value_type,
            is_bipartitions_updated=is_bipartitions_updated)
    return dist_fn(length_diffs)

