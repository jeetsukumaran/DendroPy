#! /usr/bin/env python

##############################################################################
##  DendroPy Phylocoalescenttic Computing Library.
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
##     for phylocoalescenttic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

import math
import dendropy
from dendropy.model import reconcile
from dendropy.model import coalescent
from dendropy.utility import constants
# from dendropy.calculate import combinatorics

class StructuredCoalescent(object):

    """
    Provides methods to work with the "Multispecies Coalescent", a.k.a. the
    "Truncated" or "Censored" Coalescent. That is, the coalescent process
    conditioned on a structuring process such as population subdivision or
    speciation. This is typically represented as gene or coalescent trees
    embedded within a population or species tree.
    """

    def __init__(self,
            structure_tree,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION):
        self._structure_tree = None
        self.ultrametricity_precision = ultrametricity_precision
        if structure_tree is not None:
            self._set_structure_tree(structure_tree=structure_tree)

    def _get_structure_tree(self):
        return self._structure_tree
    def _set_structure_tree(self, structure_tree):
        self._structure_tree = structure_tree
        if self._structure_tree.seed_node.age is None:
            self._structure_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        # for edge in self._structure_tree.postorder_edge_iter():
        #     edge.head_coalescent_edges = {}
        #     edge.tail_coalescent_edges = {}
    structure_tree = property(_get_structure_tree, _set_structure_tree)

    def score_coalescent_tree(self,
            coalescent_tree,
            coalescent_to_structure_node_map_fn=None,
            default_haploid_pop_size=1.0,
            ):
        """
        Returns the log-probability of a coalescent (or gene) tree conditioned
        on the structure (species or population) tree.

        Parameters
        ----------
        coalescent_tree : |Tree|
            The tree instance to be scored.
        coalescent_to_structure_node_map_fn : function object
            A function that takes a (leaf) node on the coalescent tree as an
            argument and returns the corresponding node on the structure tree
            as a value.
        default_haploid_pop_size : numeric
            Population size for each edge of the coalescent tree.

        Returns
        -------
        p : numeric
            Log probability of ``coalescent_tree`` given structuring imposed by
            ``self._structure_tree``.
        """
        edge_head_coalescent_edges, edge_tail_coalescent_edges = self._fit_coalescent_tree(coalescent_tree=coalescent_tree,
                coalescent_to_structure_node_map_fn=coalescent_to_structure_node_map_fn)
        logP = 0.0

        for structure_tree_edge in edge_head_coalescent_edges:
            coalescing_edges = edge_head_coalescent_edges[structure_tree_edge] - edge_tail_coalescent_edges[structure_tree_edge]
            k = len(edge_head_coalescent_edges[structure_tree_edge])
            t0 = structure_tree_edge.head_node.age
            t1 = 0.0
            # probability of coalescences within this edge
            for ce in coalescing_edges:
                if k == 1:
                    # t1 = structure_tree_edge.tail_node.age
                    break
                if ce.head_node is coalescent_tree.seed_node:
                    ### Whut??
                    break
                t1 = ce.tail_node.age
                wt = t1 - t0
                k2N = (float(k * (k-1)) / 2) / default_haploid_pop_size
                # k2N = float(combinatorics.choose(k, 2)) / default_haploid_pop_size
                logP = logP + math.log(k2N) - (k2N * wt)
                k -= 1
                t0 = t1
            # probability of non-coalescences within this edge
            remaining_lineages = len(edge_tail_coalescent_edges[structure_tree_edge])
            if remaining_lineages > 1 and structure_tree_edge.tail_node is not None:
                remaining_time = structure_tree_edge.tail_node.age - t1
                logP += -1 * remaining_lineages / default_haploid_pop_size * remaining_time
        return logP

    def _fit_coalescent_tree(self,
            coalescent_tree,
            coalescent_to_structure_node_map_fn):
        """
        Map edges of coalescent tree into structure tree (i.e., self).
        """
        # if self.fit_structure_edge_lengths:
        #     self.fit_edge_lengths(self.coalescent_trees)
        if coalescent_tree.seed_node.age is None:
            coalescent_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        coalescent_leaves = coalescent_tree.leaf_nodes()
        structure_to_coalescent = {}
        for coalescent_nd in coalescent_leaves:
            structure_leaf = coalescent_to_structure_node_map_fn(coalescent_nd)
            x = structure_to_coalescent.setdefault(structure_leaf, set())
            x.add(coalescent_nd.edge)
        edge_head_coalescent_edges = {}
        edge_tail_coalescent_edges = {}
        for structure_edge in self._structure_tree.postorder_edge_iter():
            if structure_edge.is_terminal():
                edge_head_coalescent_edges[structure_edge] = structure_to_coalescent[structure_edge.head_node]
            else:
                edge_head_coalescent_edges[structure_edge] = set()
                for nd in structure_edge.head_node.child_nodes():
                    # edge_head_coalescent_edges[structure_edge].update(nd.edge.tail_coalescent_edges[coalescent_tree])
                    edge_head_coalescent_edges[structure_edge].update(edge_tail_coalescent_edges[nd.edge])

            if structure_edge.tail_node is None:
                if structure_edge.length is not None:
                    target_age =  structure_edge.head_node.age + structure_edge.length
                else:
                    # assume all coalesce?
                    # structure_edge.tail_coalescent_edges[coalescent_tree] = set([coalescent_tree.seed_node.edge])
                    edge_tail_coalescent_edges[structure_edge] = set([coalescent_tree.seed_node.edge])
                    continue
            else:
                target_age = structure_edge.tail_node.age

            # structure_edge.tail_coalescent_edges[coalescent_tree] = set()
            edge_tail_coalescent_edges[structure_edge] = set()
            for coalescent_edge in edge_head_coalescent_edges[structure_edge]:
                if coalescent_edge.tail_node is not None:
                    remaining = target_age - coalescent_edge.tail_node.age
                elif coalescent_edge.length is not None:
                    remaining = target_age - (coalescent_edge.head_node.age + coalescent_edge.length)
                else:
                    continue
                while remaining > 0:
                    if coalescent_edge.tail_node is not None:
                        coalescent_edge = coalescent_edge.tail_node.edge
                    else:
                        if coalescent_edge.length is not None and (remaining - coalescent_edge.length) <= 0:
                            coalescent_edge = None
                            remaining = 0
                            break
                        else:
                            remaining = 0
                            break
                    if coalescent_edge and remaining > 0:
                        try:
                            remaining -= coalescent_edge.length
                        except TypeError:
                            pass
                if coalescent_edge is not None:
                    # structure_edge.tail_coalescent_edges[coalescent_tree].add(coalescent_edge)
                    edge_tail_coalescent_edges[structure_edge].add(coalescent_edge)
        return edge_head_coalescent_edges, edge_tail_coalescent_edges


