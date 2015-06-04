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
from dendropy.utility import error

###############################################################################
## Public Functions

def symmetric_difference(tree1, tree2, is_bipartitions_updated=False):
    """
    Returns *unweighted* Robinson-Foulds distance between two trees.

    Trees need to share the same |TaxonNamespace| reference. The
    bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
    ``is_bipartitions_updated`` argument must be `False` to force recalculation
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
        If `False`, then the bipartitions on *both* trees will be updated
        before comparison. If `True` then the bipartitions will only be
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
    Alias for `symmetric_difference()`.
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
    ``is_bipartitions_updated`` argument must be `False` to force recalculation of
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
        If `True`, then the bipartitions on *both* trees will be updated before
        comparison. If `False` (default) then the bipartitions will only be
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
    ``is_bipartitions_updated`` argument must be `False` to force recalculation of
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
        If `True`, then the bipartitions on *both* trees will be updated
        before comparison. If `False` (default) then the bipartitions
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
    false_positives = ref_bipartitions.difference(comparison_bipartitions)
    false_negatives = comparison_bipartitions.difference(ref_bipartitions)
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
    ``is_bipartitions_updated`` argument must be `False` to force recalculation of
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
        If `True`, then the bipartitions on *both* trees will be updated
        before comparison. If `False` (default) then the bipartitions
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
    df = lambda length_diffs: sqrt(sum([pow(i[0] - i[1], 2) for i in length_diffs]))
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
    ``is_bipartitions_updated`` argument must be `False` to force recalculation of
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
        If `True`, then the bipartitions on *both* trees will be updated
        before comparison. If `False` (default) then the bipartitions
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
        comparision_tree.encode_bipartitions()
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

