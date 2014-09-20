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
Statistics, metrics, measurements, and values calculated *between* *two* trees.
"""

from math import sqrt
from dendropy.calculate import treesplit
from dendropy.utility import error

###############################################################################
## Public Functions

def symmetric_difference(tree1, tree2, recalculate_splits=False):
    """
    Returns *unweighted* Robinson-Foulds distance between two trees.

    Trees need to share the same :class:`TaxonNamespace` reference. The splits
    hash bitmasks of the trees must be correct for the current tree structures
    (by calling :meth:`Tree.encode_splits()` method) or the
    `recalculate_splits` argument must be `True` to force recalculation of
    splits.

    Parameters
    ----------
    tree1 : :class:`dendropy.datamodel.Tree` object
        The first tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `tree2` and must have split
        bitmasks encoded.
    tree2 : :class:`dendropy.datamodel.Tree` object
        The second tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `tree1` and must have split
        bitmasks encoded.
    recalculate_splits : bool
        If `True`, then the split hash bitmasks on *both* trees will be updated
        before comparison. If `False` (default) then the split hash bitmasks
        will only be calculate for a :class:`Tree` object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    d : int
        The symmetric difference (a.k.a. the unweighted Robinson-Foulds
        distance) between `tree1` and `tree2`.

    Examples
    --------

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
        tree1.encode_splits()
        tree2.encode_splits()
        print(treecompare.symmetric_difference(tree1, tree2))

    """
    t = false_positives_and_negatives(
            tree1,
            tree2,
            recalculate_splits=recalculate_splits)
    return t[0] + t[1]

def unweighted_robinson_foulds_distance(tree1, tree2, recalculate_splits=False):
    """
    Alias for :func:`symmetric_difference()`.
    """
    return symmetric_difference(tree1, tree2, recalculate_splits)

def weighted_robinson_foulds_distance(
        tree1,
        tree2,
        edge_weight_attr="length",
        recalculate_splits=False):
    """
    Returns *weighted* Robinson-Foulds distance between two trees based on
    `edge_weight_attr`.

    Trees need to share the same :class:`TaxonNamespace` reference. The splits
    hash bitmasks of the trees must be correct for the current tree structures
    (by calling :meth:`Tree.encode_splits()` method) or the
    `recalculate_splits` argument must be `True` to force recalculation of
    splits.

    Parameters
    ----------
    tree1 : :class:`dendropy.datamodel.Tree` object
        The first tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `tree2` and must have split
        bitmasks encoded.
    tree2 : :class:`dendropy.datamodel.Tree` object
        The second tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `tree1` and must have split
        bitmasks encoded.
    edge_weight_attr : string
        Name of attribute on edges of trees to be used as the weight.
    recalculate_splits : bool
        If `True`, then the split hash bitmasks on *both* trees will be updated
        before comparison. If `False` (default) then the split hash bitmasks
        will only be calculate for a :class:`Tree` object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    d : float
        The edge-weighted Robinson-Foulds distance between `tree1` and `tree2`.

    Examples
    --------

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
        tree1.encode_splits()
        tree2.encode_splits()
        print(treecompare.weighted_robinson_foulds_distance(tree1, tree2))

    """
    df = lambda length_diffs: sum([abs(i[0] - i[1]) for i in length_diffs])
    return _splits_distance(tree1,
                           tree2,
                           dist_func=df,
                           edge_weight_attr=edge_weight_attr,
                           value_type=float,
                           recalculate_splits=recalculate_splits)

def false_positives_and_negatives(reference_tree, comparison_tree, recalculate_splits=False):
    """
    Counts and returns number of false positive splits (splits found in
    `comparison_tree` but not in `reference_tree`) and false negative splits (splits
    found in `reference_tree` but not in `comparison_tree`).

    Trees need to share the same :class:`TaxonNamespace` reference. The splits
    hash bitmasks of the trees must be correct for the current tree structures
    (by calling :meth:`Tree.encode_splits()` method) or the
    `recalculate_splits` argument must be `True` to force recalculation of
    splits.

    Parameters
    ----------
    reference_tree : :class:`dendropy.datamodel.Tree` object
        The first tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `comparison_tree` and must have split
        bitmasks encoded.
    comparison_tree : :class:`dendropy.datamodel.Tree` object
        The second tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `reference_tree` and must have split
        bitmasks encoded.
    recalculate_splits : bool
        If `True`, then the split hash bitmasks on *both* trees will be updated
        before comparison. If `False` (default) then the split hash bitmasks
        will only be calculate for a :class:`Tree` object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    t : tuple(int)
        A pair of integers, with first integer being the number of false
        positives and the second being the number of false negatives.

    Examples
    --------

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
        tree1.encode_splits()
        tree2.encode_splits()
        print(treecompare.false_positives_and_negatives(tree1, tree2))

    """
    sym_diff = 0
    false_positives = 0
    false_negatives = 0
    if reference_tree.taxon_namespace is not comparison_tree.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(reference_tree, comparison_tree)
    if recalculate_splits:
        treesplit.encode_splits(reference_tree)
        treesplit.encode_splits(comparison_tree)
    else:
        if reference_tree.split_edges is None:
            reference_tree.encode_splits()
        if comparison_tree.split_edges is None:
            comparison_tree.encode_splits()
    for split in reference_tree.split_edges:
        if split in comparison_tree.split_edges:
            pass
        else:
            false_negatives = false_negatives + 1
            sym_diff = sym_diff + 1

    for split in comparison_tree.split_edges:
        if split in reference_tree.split_edges:
            pass
        else:
            false_positives = false_positives + 1
            sym_diff = sym_diff + 1

    return false_positives, false_negatives

def euclidean_distance(
        tree1,
        tree2,
        edge_weight_attr="length",
        value_type=float,
        recalculate_splits=False):
    """
    Returns the Euclidean distance (a.k.a. Felsenstein's 2004 "branch length
    distance") between two trees based on `edge_weight_attr`.

    Trees need to share the same :class:`TaxonNamespace` reference. The splits
    hash bitmasks of the trees must be correct for the current tree structures
    (by calling :meth:`Tree.encode_splits()` method) or the
    `recalculate_splits` argument must be `True` to force recalculation of
    splits.

    Parameters
    ----------
    tree1 : :class:`dendropy.datamodel.Tree` object
        The first tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `tree2` and must have split
        bitmasks encoded.
    tree2 : :class:`dendropy.datamodel.Tree` object
        The second tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `tree1` and must have split
        bitmasks encoded.
    edge_weight_attr : string
        Name of attribute on edges of trees to be used as the weight.
    recalculate_splits : bool
        If `True`, then the split hash bitmasks on *both* trees will be updated
        before comparison. If `False` (default) then the split hash bitmasks
        will only be calculate for a :class:`Tree` object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    d : int
        The Euclidean distance between `tree1` and `tree2`.

    Examples
    --------

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
        tree1.encode_splits()
        tree2.encode_splits()
        print(treecompare.euclidean_distance(tree1, tree2))

    """
    df = lambda length_diffs: sqrt(sum([pow(i[0] - i[1], 2) for i in length_diffs]))
    return _splits_distance(tree1,
                           tree2,
                           dist_func=df,
                           edge_weight_attr=edge_weight_attr,
                           value_type=value_type,
                           recalculate_splits=recalculate_splits)

def find_missing_splits(reference_tree, comparison_tree, recalculate_splits=False):
    """
    Returns a list of splits that are in `reference_tree`, but
    not in `comparison_tree`.

    Trees need to share the same :class:`TaxonNamespace` reference. The splits
    hash bitmasks of the trees must be correct for the current tree structures
    (by calling :meth:`Tree.encode_splits()` method) or the
    `recalculate_splits` argument must be `True` to force recalculation of
    splits.

    Parameters
    ----------
    reference_tree : :class:`dendropy.datamodel.Tree` object
        The first tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `comparison_tree` and must have split
        bitmasks encoded.
    comparison_tree : :class:`dendropy.datamodel.Tree` object
        The second tree of the two trees being compared. This must share the
        same :class:`TaxonNamespace` reference as `reference_tree` and must have split
        bitmasks encoded.
    recalculate_splits : bool
        If `True`, then the split hash bitmasks on *both* trees will be updated
        before comparison. If `False` (default) then the split hash bitmasks
        will only be calculate for a :class:`Tree` object if they have not been
        calculated before, either explicitly or implicitly.

    Returns
    -------
    s : list[splits]
        A list of splits that are in the first tree but not in the second.

    """
    missing = []
    if reference_tree.taxon_namespace is not comparison_tree.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(reference_tree, comparison_tree)
    if recalculate_splits:
        treesplit.encode_splits(reference_tree)
        treesplit.encode_splits(comparison_tree)
    else:
        if reference_tree.split_edges is None:
            reference_tree.encode_splits()
        if comparison_tree.split_edges is None:
            comparison_tree.encode_splits()
    for split in reference_tree.split_edges:
        if split in comparison_tree.split_edges:
            pass
        else:
            missing.append(split)
    return missing

###############################################################################
## Legacy

def robinson_foulds_distance(tree1, tree2, edge_weight_attr="length"):
    """
    DEPRECATED: Use :func:`symmetric_difference` for the common
    unweighged Robinson-Fould's distance metric (i.e., the symmetric difference between two trees)
    :func:`weighted_robinson_foulds_distance` or for the RF distance as defined by Felsenstein, 2004.
    """
    return weighted_robinson_foulds_distance(tree1, tree2, edge_weight_attr)

def mason_gamer_kellogg_score(tree1, tree2, recalculate_splits=False):
    """
    Mason-Gamer and Kellogg. Testing for phylogenetic conflict among molecular
    data sets in the tribe Triticeae (Gramineae). Systematic Biology (1996)
    vol. 45 (4) pp. 524
    """
    if tree1.taxon_namespace is not tree2.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(tree1, tree2)
    if recalculate_splits:
        tree1.encode_splits()
        tree2.encode_splits()
    else:
        if tree1.split_edges is None:
            tree1.encode_splits()
        if tree2.split_edges is None:
            tree2.encode_splits()
    se1 = tree1.split_edges
    se2 = tree2.split_edges
    splits = sorted(list(set(se1.keys() + se2.keys())))

###############################################################################
## Supporting

def _get_length_diffs(
        tree1,
        tree2,
        edge_weight_attr="length",
        value_type=float,
        recalculate_splits=False,
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
    if tree1.taxon_namespace is not tree2.taxon_namespace:
        raise error.TaxonNamespaceIdentityError(tree1, tree2)
    if recalculate_splits:
        treesplit.encode_splits(tree1)
        treesplit.encode_splits(tree2)
    else:
        if tree1.split_edges is None:
            tree1.encode_splits()
        if tree2.split_edges is None:
            tree2.encode_splits()
    split_edges2_copy = dict(tree2.split_edges) # O(n*(2*bind + dict_item_cost))
    split_edges1_ref = tree1.split_edges
    for split in split_edges1_ref: # O n : 2*bind
        edge = split_edges1_ref[split]
        elen1 = getattr(edge, edge_weight_attr) # attr + bind
        if elen1 is None:
            elen1 = 0 # worst-case: bind
        value1 = value_type(elen1) #  ctor + bind
        try:
            e2 = split_edges2_copy.pop(split) # attr + dict_lookup + bind
            elen2 = getattr(e2, edge_weight_attr) # attr + bind
            if elen2 is None:
                # allow root edge to have split with no value: raise error if not root edge
                if e2.tail_node is None:
                    elen2 = 0.0
                else:
                    raise ValueError("Edge length attribute is 'None': Tree: %s ('%s'), Split: %s" % (tree2.oid, tree2.label, tree2.taxon_namespace.split_as_newick_string(split)))
        except KeyError: # excep
            elen2 = 0.0
        value2 = value_type(elen2) #  ctor + bind # best case
        length_diffs.append((value1,value2)) # ctor + listappend
        split_length_diffs[split] = length_diffs[-1]

    for split in split_edges2_copy: # best-case not executed, worst case O(n) : 2*bind
        edge = split_edges2_copy[split]
        elen2 = getattr(edge, edge_weight_attr) # attr +  bind
        if elen2 is None:
            elen2 = 0
        value2 = value_type(elen2) #  ctor + bind
        e1 = split_edges1_ref.get(split) # attr + dict_lookup + bind
        if e1 is None:
            elen1 = 0.0
        else:
            elen1 = getattr(e1, edge_weight_attr) # attr  + bind
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

def _splits_distance(
        tree1,
        tree2,
        dist_func,
        edge_weight_attr="length",
        value_type=float,
        recalculate_splits=False):
    """
    Returns distance between two trees, each represented by a dictionary of
    splits (as split_mask strings) to edges, using `dist_func` to calculate the
    distance based on `edge_weight_attr` of the edges. `dist_func` is a function
    that takes a list of pairs of values, where the values correspond to the edge
    lengths of a given split on tree1 and tree2 respectively.
    """
    length_diffs = _get_length_diffs(
            tree1,
            tree2,
            edge_weight_attr=edge_weight_attr,
            value_type=value_type,
            recalculate_splits=recalculate_splits)
    return dist_func(length_diffs)

