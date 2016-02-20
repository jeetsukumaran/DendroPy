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
Models, modeling and model-fitting of parsimony.
"""

from functools import reduce
import operator
import dendropy
from dendropy.utility.error import TaxonNamespaceIdentityError

class _NodeStateSetMap(dict):
    def __init__(self, taxon_state_sets_map=None):
        self.taxon_state_sets_map = taxon_state_sets_map
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            v = self.taxon_state_sets_map[key.taxon]
            self[key] = v
            return v

def _store_sets_as_attr(n, state_sets_attr_name, v):
    setattr(n, state_sets_attr_name, v)

def _retrieve_state_sets_from_attr(n, state_sets_attr_name, taxon_state_sets_map):
    try:
        return getattr(n, state_sets_attr_name)
    except AttributeError:
        v = taxon_state_sets_map[n.taxon]
        setattr(n, state_sets_attr_name, v)
        return v

def fitch_down_pass(
        postorder_nodes,
        state_sets_attr_name="state_sets",
        taxon_state_sets_map=None,
        weights=None,
        score_by_character_list=None,
        ):
    """
    Returns the parsimony score given a list of nodes in postorder and
    associated states, using Fitch's (1971) unordered parsimony algorithm.

    Parameters
    ----------
    postorder_nodes : iterable of/over |Node| objects
        An iterable of |Node| objects in in order of post-order
        traversal of the tree.
    state_sets_attr_name : str
        Name of attribute on |Node| objects in which state set lists
        will stored/accessed. If |None|, then state sets will not be stored on
        the tree.
    taxon_state_sets_map : dict[taxon] = state sets
        A dictionary that takes a taxon object as a key and returns a state set
        list as a value. This will be used to populate the state set of a node
        that has not yet had its state sets scored and recorded (typically,
        leaves of a tree that has not yet been processed).
    weights : iterable
        A list of weights for each pattern.
    score_by_character_list : None or list
        If not |None|, should be a reference to a list object.
        This list will be populated by the scores on a character-by-character
        basis.

    Returns
    -------
    s : int
        Parismony score of tree.

    Notes
    -----
    Currently this requires a bifurcating tree (even at the root).

    Examples
    --------

    Assume that we have a tree, ``tree``, and an associated data set, ``data``::

        import dendropy
        from dendropy.model.parsimony import fitch_down_pass

        taxa = dendropy.TaxonNamespace()
        data = dendropy.StandardCharacterMatrix.get_from_path(
                "apternodus.chars.nexus",
                "nexus",
                taxon_namespace=taxa)
        tree = dendropy.Tree.get_from_path(
                "apternodus.tre",
                "nexus",
                taxon_namespace=taxa)
        taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=True)

    The following will return the parsimony score of the ``tree`` with
    respect to the data in ``data``::

        score = fitch_down_pass(
                nodes=tree.postorder_node_iter(),
                taxon_state_sets_map=taxon_set_map)
        print(score)

    In the above, every |Node| object of ``tree`` will have an attribute
    added, "state_sets", that stores the list of state sets from the analysis::

        for nd in tree:
            print(nd.state_sets)

    If you want to store the list of state sets in a different attribute, e.g.,
    "analysis1_states"::

        score = fitch_down_pass(
                nodes=tree.postorder_node_iter(),
                state_sets_attr_name="analysis1_states",
                taxon_state_sets_map=taxon_set_map)
        print(score)
        for nd in tree:
            print(nd.analysis1_states)

    Or not to store these at all::

        score = fitch_down_pass(
                nodes=tree.postorder_node_iter(),
                state_sets_attr_name=None,
                taxon_state_sets_map=taxon_set_map)
        print(score)

    Scoring custom data can be done by something like the following::

        taxa = dendropy.TaxonNamespace()
        taxon_state_sets_map = {}
        t1 = taxa.require_taxon("A")
        t2 = taxa.require_taxon("B")
        t3 = taxa.require_taxon("C")
        t4 = taxa.require_taxon("D")
        t5 = taxa.require_taxon("E")
        taxon_state_sets_map[t1] = [ set([0,1]),  set([0,1]),  set([0]),     set([0]) ]
        taxon_state_sets_map[t2] = [ set([1]),    set([1]),    set([1]),     set([0]) ]
        taxon_state_sets_map[t3] = [ set([0]),    set([1]),    set([1]),     set([0]) ]
        taxon_state_sets_map[t4] = [ set([0]),    set([1]),    set([0,1]),   set([1]) ]
        taxon_state_sets_map[t5] = [ set([1]),    set([0]),    set([1]),     set([1]) ]
        tree = dendropy.Tree.get_from_string(
                "(A,(B,(C,(D,E))));", "newick",
                taxon_namespace=taxa)
        score = fitch_down_pass(tree.postorder_node_iter(),
                taxon_state_sets_map=taxon_state_sets_map)
        print(score)

    """
    if score_by_character_list is not None:
        assert len(score_by_character_list) == 0
        for idx in range(len(list(taxon_state_sets_map.values())[0])): # this is unacceptable!
            score_by_character_list.append(0)
    score = 0
    if state_sets_attr_name is None:
        node_state_set_map = _NodeStateSetMap(taxon_state_sets_map)
        get_node_state_sets = lambda node : node_state_set_map[node]
        set_node_state_sets = lambda node, v : node_state_set_map.__setitem__(node, v)
    else:
        get_node_state_sets = lambda node : _retrieve_state_sets_from_attr(node, state_sets_attr_name, taxon_state_sets_map)
        set_node_state_sets = lambda node, v : _store_sets_as_attr(node, state_sets_attr_name, v)
    for nd in postorder_nodes:
        c = nd.child_nodes()
        if not c:
            ss = get_node_state_sets(nd)
            continue
        left_c, right_c = c[:2]
        remaining = c[2:]
        left_ssl = get_node_state_sets(left_c)
        while True:
            right_ssl = get_node_state_sets(right_c)
            result = []
            for n, ssp in enumerate(zip(left_ssl, right_ssl)):
                left_ss, right_ss = ssp
                inter = left_ss.intersection(right_ss)
                if inter:
                    result.append(inter)
                else:
                    if weights is None:
                        wt = 1
                    else:
                        wt = weights[n]
                    score += wt
                    result.append(left_ss.union(left_ss, right_ss))
                    if score_by_character_list is not None:
                        # try:
                        #     score_by_character_list[n] += wt
                        # except IndexError:
                        #     score_by_character_list.append(wt)
                        score_by_character_list[n] += wt
            if remaining:
                right_c = remaining.pop(0)
                left_ssl = result
            else:
                break
        # setattr(nd, state_sets_attr_name, result)
        set_node_state_sets(nd, result)
    return score

def fitch_up_pass(
        preorder_node_list,
        state_sets_attr_name="state_sets",
        taxon_state_sets_map=None):
    """
    Finalizes the state set lists associated with each node using the "final
    phase" of Fitch's (1971) unordered parsimony algorithm.

    Parameters
    ----------
    postorder_nodes : iterable of/over |Node| objects
        An iterable of |Node| objects in in order of post-order
        traversal of the tree.
    state_sets_attr_name : str
        Name of attribute on |Node| objects in which state set lists
        will stored/accessed. If |None|, then state sets will not be stored on
        the tree.
    taxon_state_sets_map : dict[taxon] = state sets
        A dictionary that takes a taxon object as a key and returns a state set
        list as a value. This will be used to populate the state set of a node
        that has not yet had its state sets scored and recorded (typically,
        leaves of a tree that has not yet been processed).

    Notes
    -----
    Currently this requires a bifurcating tree (even at the root).

    Examples
    --------

    ::

        taxa = dendropy.TaxonNamespace()
        data = dendropy.StandardCharacterMatrix.get_from_path(
                "apternodus.chars.nexus",
                "nexus",
                taxon_namespace=taxa)
        tree = dendropy.Tree.get_from_path(
                "apternodus.tre",
                "nexus",
                taxon_namespace=taxa)
        taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=True)
        score = fitch_down_pass(tree.postorder_node_iter(),
                taxon_state_sets_map=taxon_state_sets_map)
        print(score)
        fitch_up_pass(tree.preorder_node_iter())
        for nd in tree:
            print(nd.state_sets)

    """
    node_state_sets_map = {}
    for nd in preorder_node_list:
        c = nd.child_nodes()
        p = nd.parent_node
        if (not c) or (not p):
            continue
        assert(len(c) == 2)
        left_c, right_c = c
        try:
            left_ssl = getattr(left_c, state_sets_attr_name)
        except AttributeError:
            if not taxon_state_sets_map:
                raise
            left_ssl = taxon_state_sets_map[left_c.taxon]
        try:
            right_ssl = getattr(right_c, state_sets_attr_name)
        except AttributeError:
            if not taxon_state_sets_map:
                raise
            right_ssl = taxon_state_sets_map[right_c.taxon]
        par_ssl = getattr(p, state_sets_attr_name)
        curr_ssl = getattr(nd, state_sets_attr_name)
        result = []
        for n, ssp in enumerate(zip(par_ssl, curr_ssl, left_ssl, right_ssl)):
            par_ss, curr_ss, left_ss, right_ss = ssp

            down_parup_inter = par_ss.intersection(curr_ss)
            if down_parup_inter == par_ss:
                final_ss = down_parup_inter
            else:
                rl_inter = left_ss.intersection(right_ss)
                if not rl_inter:
                    final_ss = par_ss.union(curr_ss)
                else:
                    in_par_and_left = par_ss.intersection(left_ss)
                    in_par_and_right = par_ss.intersection(right_ss)
                    final_ss = in_par_and_left.union(in_par_and_right, curr_ss)
            #_LOG.debug("downpass = %s, par = %s, left = %s, right = %s, final_ss= %s" %
            #                    (str(curr_ss), str(par_ss), str(left_ss), str(right_ss), str(final_ss)))
            result.append(final_ss)
        setattr(nd, state_sets_attr_name, result)


def parsimony_score(
        tree,
        chars,
        gaps_as_missing=True,
        weights=None,
        score_by_character_list=None,
        ):
    """
    Calculates the score of a tree, ``tree``, given some character data,
    ``chars``, under the parsimony model using the Fitch algorithm.

    Parameters
    ----------
    tree : a |Tree| instance
        A |Tree| to be scored. Must reference the same |TaxonNamespace| as
        ``chars``.
    chars : a |CharacterMatrix| instance
        A |CharacterMatrix|-derived object with data to be scored. Must have
        the same |TaxonNamespace| as ``tree``.
    gap_as_missing : bool
        If |True| [default], then gaps will be treated as missing data.
        If |False|, then gaps will be treated as a new/additional state.
    weights : iterable
        A list of weights for each pattern/column in the matrix.
    score_by_character_list : None or list
        If not |None|, should be a reference to a list object.
        This list will be populated by the scores on a character-by-character
        basis.

    Returns
    -------
    pscore : int
        The parsimony score of the tree given the data.

    Examples
    --------

    ::

        import dendropy
        from dendropy.calculate import treescore

        # establish common taxon namespace
        taxon_namespace = dendropy.TaxonNamespace()

        # Read data; if data is, e.g., "standard", use StandardCharacterMatrix.
        # If unsure of data type, can do:
        #       dataset = dendropy.DataSet.get(
        #               path="path/to/file.nex",
        #               schema="nexus",
        #               taxon_namespace=tns,)
        #       chars = dataset.char_matrices[0]
        chars = dendropy.DnaCharacterMatrix.get(
                path="pythonidae.chars.nexus",
                schema="nexus",
                taxon_namespace=taxon_namespace)
        tree = dendropy.Tree.get(
                path="pythonidae.mle.newick",
                schema="newick",
                taxon_namespace=taxon_namespace)

        # We store the site-specific scores here
        # This is optional; if we do not want to
        # use the per-site scores, just pass in |None|
        # for the ``score_by_character_list`` argument
        # or do not specify this argument at all.
        score_by_character_list = []

        score = treescore.parsimony_score(
                tree,
                chars,
                gaps_as_missing=False,
                score_by_character_list=score_by_character_list)

        # Print the results: the score
        print("Score: {}".format(score))

        # Print the results: the per-site scores
        for idx, x in enumerate(score_by_character_list):
            print("{}: {}".format(idx+1, x))

    Notes
    -----

    If the same data is going to be used to score multiple trees or multiple times,
    it is probably better to generate the 'taxon_state_sets_map' once and call
    "fitch_down_pass" directly yourself, as this function generates a new map
    each time.

    """
    if tree.taxon_namespace is not chars.taxon_namespace:
        raise TaxonNamespaceIdentityError(tree, data)
    taxon_state_sets_map = chars.taxon_state_sets_map(gaps_as_missing=gaps_as_missing)
    nodes = tree.postorder_node_iter()
    pscore = fitch_down_pass(nodes,
            taxon_state_sets_map=taxon_state_sets_map,
            weights=weights,
            score_by_character_list=score_by_character_list)
    return pscore

