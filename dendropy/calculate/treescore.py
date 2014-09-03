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
Statistics, metrics, measurements, and values calculated on a tree with
reference to external data of some kind under various criteria.
"""

from functools import reduce
import operator
import dendropy

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
            for n, ssp in enumerate(zip(left_ssl, right_ssl)):
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
        for n, ssp in enumerate(zip(par_ssl, curr_ssl, left_ssl, right_ssl)):
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

