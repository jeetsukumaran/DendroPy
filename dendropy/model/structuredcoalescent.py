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
            coalescent_to_structure_map_fn,
            default_haploid_pop_size=1.0,
            is_coalescent_to_structure_map_by_node=False,
            ):
        """
        Returns the log-probability of a coalescent (or gene) tree conditioned
        on the structure (species or population) tree.

        Parameters
        ----------
        coalescent_tree : |Tree|
            The tree instance to be scored.
        coalescent_to_structure_map_fn : function object
            A function that takes either a |Taxon| instance (if
            ``is_coalescent_to_structure_map_by_node`` is False) or
            |Node| instance (if ``is_coalescent_to_structure_map_by_node`` is
            True) representing a lineage on the coalescent or gene tree
            (specified by ``coalescent_tree``), and returns the |Taxon|
            instance (if ``is_coalescent_to_structure_map_by_node`` is False)
            or |Node| instance (if ``is_coalescent_to_structure_map_by_node``
            is True) corresponding to the species or population on the
            structure tree with which it is associated.
        default_haploid_pop_size : numeric
            Population size for each edge of the coalescent tree.
        is_coalescent_to_structure_map_by_node : str
            Specifies the expected type of argument and return value of the
            mapping function, ``coalescent_to_structure_map_fn``. By default
            this is |False|, and the mapping function is thus expected to take
            a |Taxon| instance representing a lineage on the coalescent or gene
            tree and return a |Taxon| instance representing a species or
            population on the species tree with which that lineage is
            associated. This is more efficient if you have many
            moderately-sized gene trees that share the same taxa: you only need
            to construct and provide a single mapping. On the other hand, if
            you are dealing with HUGE trees, it might be optimum to skip
            processing the taxon namespace, i.e. deserializing the tip labels
            into rich |Taxon| objects, and just deal with them as plain string
            labels. In this case, you would want to map the coalescent trees to
            the structure trees by nodes based on the labels, and will species
            ``True`` for the ``is_coalescent_to_structure_map_by_node``
            argument. Then the mapping function
            ``coalescent_to_structure_map_fn`` will be expected to take a
            |Node| instance representing a lineage on the coalescent or gene
            tree and return a |Node| instance representing a species or
            population on the species tree with which that lineage is
            associated.

        Returns
        -------
        p : numeric
            Log probability of ``coalescent_tree`` given structuring imposed by
            ``self._structure_tree``.
        """
        edge_head_coalescent_edges, edge_tail_coalescent_edges, edge_coalescent_nodes = self._fit_coalescent_tree(
                coalescent_tree=coalescent_tree,
                coalescent_to_structure_map_fn=coalescent_to_structure_map_fn,
                is_coalescent_to_structure_map_by_node=is_coalescent_to_structure_map_by_node)
        logP = 0.0

        # for nd in self._structure_tree.postorder_node_iter():
        #     if not nd._child_nodes:
        #         nd._haploid_pop_size = default_haploid_pop_size
        #     else:
        #         nd._haploid_pop_size = sum(ch._haploid_pop_size for ch in nd._child_nodes)

        # for nd in self._structure_tree.preorder_node_iter():
        #     if nd is self._structure_tree.seed_node:
        #         nd._haploid_pop_size = default_haploid_pop_size
        #     else:
        #         nd._haploid_pop_size = nd.parent_node._haploid_pop_size / 2.0

        for structure_tree_edge in self._structure_tree.postorder_edge_iter():
        # for structure_tree_edge in edge_coalescent_nodes:
            haploid_pop_size = default_haploid_pop_size
            # haploid_pop_size = structure_tree_edge.head_node._haploid_pop_size
            k = len(edge_head_coalescent_edges[structure_tree_edge])
            t0 = structure_tree_edge.head_node.age
            t1 = structure_tree_edge.head_node.age
            oldest_coalescent_event_age = t1
            # probability of coalescences within this edge
            coalescing_nodes = sorted(edge_coalescent_nodes[structure_tree_edge], key=lambda nd: nd.age)
            print("\n{}".format(edge_tail_coalescent_edges[structure_tree_edge]))
            print("\n{} => {}: {}".format(
                len(edge_head_coalescent_edges[structure_tree_edge]),
                len(edge_tail_coalescent_edges[structure_tree_edge]),
                [ce.age for ce in coalescing_nodes],
                ))
            for cnd in coalescing_nodes:
                print(k)
                if k == 1:
                    # t1 = structure_tree_edge.tail_node.age
                    break
                t1 = cnd.age
                wt = t1 - t0
                k2N = (float(k * (k-1)) / 2) / haploid_pop_size
                # k2N = float(combinatorics.choose(k, 2)) / default_haploid_pop_size
                logP = logP + math.log(k2N) - (k2N * wt)
                k -= 1
                t0 = t1
                if t1 > oldest_coalescent_event_age:
                    oldest_coalescent_event_age = t1
            # assert k == len(edge_tail_coalescent_edges[structure_tree_edge]), "{} vs {}: {}".format(k, len(edge_tail_coalescent_edges[structure_tree_edge]), edge_tail_coalescent_edges[structure_tree_edge])
            # probability of non-coalescences within this edge
            remaining_lineages = len(edge_tail_coalescent_edges[structure_tree_edge])
            if remaining_lineages > 1 and structure_tree_edge.tail_node is not None:
                remaining_time = structure_tree_edge.tail_node.age - oldest_coalescent_event_age
                # logP += -1 * remaining_lineages / default_haploid_pop_size * remaining_time
                logP += -1 * remaining_lineages / haploid_pop_size * remaining_time
        return logP

    def _fit_coalescent_tree(self,
            coalescent_tree,
            coalescent_to_structure_map_fn,
            is_coalescent_to_structure_map_by_node=False):
        """
        Map edges of coalescent tree into structure tree (i.e., self).
        """
        # if self.fit_structure_edge_lengths:
        #     self.fit_edge_lengths(self.coalescent_trees)
        if coalescent_tree.seed_node.age is None:
            coalescent_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        coalescent_leaves = coalescent_tree.leaf_nodes()
        structure_to_coalescent = {}
        if is_coalescent_to_structure_map_by_node:
            for coalescent_nd in coalescent_leaves:
                structure_leaf = coalescent_to_structure_map_fn(coalescent_nd)
                x = structure_to_coalescent.setdefault(structure_leaf, set())
                x.add(coalescent_nd.edge)
        else:
            for coalescent_nd in coalescent_leaves:
                structure_taxon = coalescent_to_structure_map_fn(coalescent_nd.taxon)
                x = structure_to_coalescent.setdefault(structure_taxon, set())
                x.add(coalescent_nd.edge)
        edge_head_coalescent_edges = {}
        edge_tail_coalescent_edges = {}
        edge_coalescent_nodes = {}
        for structure_edge in self._structure_tree.postorder_edge_iter():
            if structure_edge.is_terminal():
                if is_coalescent_to_structure_map_by_node:
                    edge_head_coalescent_edges[structure_edge] = structure_to_coalescent[structure_edge.head_node]
                else:
                    edge_head_coalescent_edges[structure_edge] = structure_to_coalescent[structure_edge.head_node.taxon]
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
                    # edge_tail_coalescent_edges[structure_edge] = set([coalescent_tree.seed_node.edge])
                    edge_tail_coalescent_edges[structure_edge] = set([])
                    # print(">>>>>>>>>>> {}".format(len(edge_head_coalescent_edges[structure_edge])))
                    # print(">>>>>>>>>>> {}".format(len([e.tail_node for e in edge_head_coalescent_edges[structure_edge]])))
                    # print(">>>>>>>>>>> {}".format(["{}<={}".format(e.head_node.label, e.tail_node.label) for e in edge_head_coalescent_edges[structure_edge]]))
                    # print(">>>>>>>>>>> {}".format(len(set(e.tail_node for e in edge_head_coalescent_edges[structure_edge]))))
                    edge_coalescent_nodes[structure_edge] = set()
                    for ex in edge_head_coalescent_edges[structure_edge]:
                        edge_coalescent_nodes[structure_edge].add(ex.tail_node)
                        if ex.tail_node is not None:
                            parent_node = ex.tail_node.parent_node
                            while parent_node is not None and parent_node not in edge_coalescent_nodes[structure_edge]:
                                edge_coalescent_nodes[structure_edge].add(parent_node)
                                parent_node = parent_node.parent_node
                    # print("==== {}".format([nd.label for nd in edge_coalescent_nodes[structure_edge]]))
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
                edge_coalescent_nodes[structure_edge] = set([e.tail_node for e in (edge_head_coalescent_edges[structure_edge] - edge_tail_coalescent_edges[structure_edge])])
        return edge_head_coalescent_edges, edge_tail_coalescent_edges, edge_coalescent_nodes


