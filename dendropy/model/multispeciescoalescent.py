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
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

from __future__ import division
import math
import dendropy
import heapq
from dendropy.model import reconcile
from dendropy.model import coalescent
from dendropy.utility import constants
from dendropy.utility import error
# from dendropy.calculate import combinatorics

class MultispeciesCoalescent(object):

    """
    Provides methods to work with the "Multispecies Coalescent", a.k.a. the
    "Truncated" or "Censored" Coalescent. That is, the coalescent process
    conditioned on a structuring process such as population subdivision or
    speciation. This is typically represented as gene or coalescent trees
    embedded within a population or species tree.
    """

    def __init__(self,
            species_tree,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION):
        self._species_tree = None
        self.ultrametricity_precision = ultrametricity_precision
        self.is_enforce_structure_integrity = True
        if species_tree is not None:
            self._set_species_tree(species_tree=species_tree)

    def _get_species_tree(self):
        return self._species_tree
    def _set_species_tree(self, species_tree):
        self._species_tree = species_tree
        self._species_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
    species_tree = property(_get_species_tree, _set_species_tree)

    def score_coalescent_tree(self,
            coalescent_tree,
            coalescent_species_lineage_map_fn,
            population_theta_fn=None,
            is_coalescent_species_lineage_map_by_node=False,
            ):
        """
        Returns the log-probability of a coalescent (or gene) tree conditioned
        on the structure (species or population) tree.

        Parameters
        ----------
        coalescent_tree : |Tree|
            The tree instance to be scored.
        coalescent_species_lineage_map_fn : function object
            A function that takes either a |Taxon| instance (if
            ``is_coalescent_species_lineage_map_by_node`` is False) or
            |Node| instance (if ``is_coalescent_species_lineage_map_by_node`` is
            True) representing a lineage on the coalescent or gene tree
            (specified by ``coalescent_tree``), and returns the |Taxon|
            instance (if ``is_coalescent_species_lineage_map_by_node`` is False)
            or |Node| instance (if ``is_coalescent_species_lineage_map_by_node``
            is True) corresponding to the species or population on the
            species tree with which it is associated.
        population_theta_fn : function object
            Function that takes an edge on the species structure tree as an argument
            and returns the population parameter (theta) for that population or
            species. If not specified, all edges are assumed to have a theta
            value of 1.0.
        is_coalescent_species_lineage_map_by_node : bool
            Specifies the expected type of argument and return value of the
            mapping function, ``coalescent_species_lineage_map_fn``. By default
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
            the species tree by nodes based on the labels, and will species
            ``True`` for the ``is_coalescent_species_lineage_map_by_node``
            argument. Then the mapping function
            ``coalescent_species_lineage_map_fn`` will be expected to take a
            |Node| instance representing a lineage on the coalescent or gene
            tree and return a |Node| instance representing a species or
            population on the species tree with which that lineage is
            associated.

        Returns
        -------
        p : numeric
            Log probability of ``coalescent_tree`` given structuring imposed by
            ``self._species_tree``.
        """
        edge_head_coalescent_edges, edge_tail_coalescent_edges, edge_coalescent_nodes = self._fit_coalescent_tree(
                coalescent_tree=coalescent_tree,
                coalescent_species_lineage_map_fn=coalescent_species_lineage_map_fn,
                is_coalescent_species_lineage_map_by_node=is_coalescent_species_lineage_map_by_node)
        if population_theta_fn is None:
            population_theta_fn = lambda e: 1.0
        logP = 0.0

        # def _debug_log(x):
        #     print(x)

        for species_tree_edge in self._species_tree.postorder_edge_iter():
        # for species_tree_edge in edge_coalescent_nodes:
            theta = population_theta_fn(species_tree_edge)
            coalescing_nodes = sorted(edge_coalescent_nodes[species_tree_edge], key=lambda nd: nd.age if nd else float("inf"))

            # if species_tree_edge.tail_node is None:
            #     _debug_log("-- Structure Edge {}: from {} to infinity (root)".format(self._compose_edge_desc(species_tree_edge), species_tree_edge.head_node.age, ))
            # else:
            #     _debug_log("-- Structure Edge {}: from {} to {}".format(self._compose_edge_desc(species_tree_edge), species_tree_edge.head_node.age, species_tree_edge.tail_node.age))
            # _debug_log("          In: {} edges: {}".format(len(edge_head_coalescent_edges[species_tree_edge]), [self._compose_edge_desc(e) for e in edge_head_coalescent_edges[species_tree_edge]]))
            # _debug_log("         Out: {} edges: {}".format(len(edge_tail_coalescent_edges[species_tree_edge]), [self._compose_edge_desc(e) for e in edge_tail_coalescent_edges[species_tree_edge]]))
            # _debug_log("Coalescences: {} nodes: {}".format(len(coalescing_nodes),
            #     ["{} ^ {} (= {} at {})".format(self._compose_edge_desc(c._child_nodes[0].edge), self._compose_edge_desc(c._child_nodes[1].edge), self._compose_edge_desc(c.edge), c.age) for c in coalescing_nodes if c is not None]))

            j = len(edge_head_coalescent_edges[species_tree_edge])
            t0 = species_tree_edge.head_node.age
            t1 = species_tree_edge.head_node.age
            oldest_coalescent_event_age = None
            subP = 0.0
            for cnd in coalescing_nodes:
                if j == 1:
                    break
                t1 = cnd.age
                wt = t1 - t0
                # q = math.log( (j*(j-1.0))  /theta) + (-j * (j-1.0) * theta * wt)
                q = math.log(2.0/theta) + (-j * (j-1) * theta * wt)
                subP += q
                j -= 1
                t0 = t1
                if oldest_coalescent_event_age is None or t1 > oldest_coalescent_event_age:
                    oldest_coalescent_event_age = t1

            remaining_lineages = j
            if remaining_lineages > 1:
                if oldest_coalescent_event_age is None:
                    remaining_time = species_tree_edge.tail_node.age - species_tree_edge.head_node.age
                else:
                    remaining_time = species_tree_edge.tail_node.age - oldest_coalescent_event_age
                q = -1 * (remaining_lineages*(remaining_lineages-1))/theta * remaining_time
                subP += q
            logP += subP
        return logP

    def _fit_coalescent_tree(self,
            coalescent_tree,
            coalescent_species_lineage_map_fn,
            is_coalescent_species_lineage_map_by_node=False):
        """
        Map edges of coalescent tree into species tree (i.e., self).
        """
        # if self.fit_species_edge_lengths:
        #     self.fit_edge_lengths(self.coalescent_trees)
        coalescent_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        coalescent_leaves = coalescent_tree.leaf_nodes()
        structure_to_coalescent = {}
        if is_coalescent_species_lineage_map_by_node:
            for coalescent_nd in coalescent_leaves:
                structure_leaf = coalescent_species_lineage_map_fn(coalescent_nd)
                x = structure_to_coalescent.setdefault(structure_leaf, set())
                x.add(coalescent_nd.edge)
        else:
            for coalescent_nd in coalescent_leaves:
                structure_taxon = coalescent_species_lineage_map_fn(coalescent_nd.taxon)
                x = structure_to_coalescent.setdefault(structure_taxon, set())
                x.add(coalescent_nd.edge)
        edge_head_coalescent_edges = {}
        edge_tail_coalescent_edges = {}
        edge_coalescent_nodes = {}

        # def _debug_log(x):
        #     print(x)
        # def _dump_edge_set(ee):
        #     return ", ".join(self._compose_edge_desc(e) for e in ee)

        for species_edge in self._species_tree.postorder_edge_iter():

            # if species_edge.tail_node is None:
            #     _debug_log("\nStructure Edge {}: from {} to infinity (root)".format(self._compose_edge_desc(species_edge), species_edge.head_node.age, ))
            # else:
            #     _debug_log("\nStructure Edge {}: from {} to {}".format(self._compose_edge_desc(species_edge), species_edge.head_node.age, species_edge.tail_node.age))

            ## add initial/inherited coalescent edges to structure edges
            if species_edge.is_terminal():
                if is_coalescent_species_lineage_map_by_node:
                    edge_head_coalescent_edges[species_edge] = structure_to_coalescent[species_edge.head_node]
                else:
                    edge_head_coalescent_edges[species_edge] = structure_to_coalescent[species_edge.head_node.taxon]
                # _debug_log("    Initializing terminal with edges: {}".format(_dump_edge_set(edge_head_coalescent_edges[species_edge])))
            else:
                edge_head_coalescent_edges[species_edge] = set()
                # _debug_log("    Initializing internal")
                for nd in species_edge.head_node.child_node_iter():
                    # edge_head_coalescent_edges[species_edge].update(nd.edge.tail_coalescent_edges[coalescent_tree])
                    # _debug_log("        Adding from {}: {}".format(self._compose_edge_desc(nd.edge), _dump_edge_set(edge_tail_coalescent_edges[nd.edge])))
                    edge_head_coalescent_edges[species_edge].update(edge_tail_coalescent_edges[nd.edge])
                # _debug_log("        Internal now has edges: {}".format(_dump_edge_set(edge_head_coalescent_edges[species_edge])))

            ## initialize data containers
            edge_coalescent_nodes[species_edge] = set()
            if len(edge_head_coalescent_edges[species_edge]) == 1:
                edge_tail_coalescent_edges[species_edge] = set(edge_head_coalescent_edges[species_edge])
                continue
            edge_tail_coalescent_edges[species_edge] = set([])

            if species_edge.tail_node is None:
                ## root edge
                # _debug_log("    Structure root edge: coalesce all")
                for ex in edge_head_coalescent_edges[species_edge]:
                    edge_coalescent_nodes[species_edge].add(ex.tail_node)
                    if ex.tail_node is not None:
                        parent_node = ex.tail_node.parent_node
                        while parent_node is not None and parent_node not in edge_coalescent_nodes[species_edge]:
                            edge_coalescent_nodes[species_edge].add(parent_node)
                            parent_node = parent_node.parent_node
                # _debug_log("        Structure edge now has edges: {}".format(_dump_edge_set(edge_head_coalescent_edges[species_edge])))
                continue

            structure_end_time = species_edge.tail_node.age
            # _debug_log("    Structure end time: {}".format(structure_end_time))
            current_lineages = list( (e.tail_node.age, e) for e in edge_head_coalescent_edges[species_edge] )
            heapq.heapify(current_lineages)
            if self.is_enforce_structure_integrity:
                valid_coalescing_lineages = set(edge_head_coalescent_edges[species_edge])
            while len(current_lineages) > 1:
                coalescent_age, coalescent_edge = current_lineages[0]
                if coalescent_age > structure_end_time:
                    # _debug_log("    Time exceeded with {} lineages remaining".format(len(current_lineages)))
                    break
                # _debug_log("    {} lineages remaining".format(len(current_lineages)))
                heapq.heappop(current_lineages)
                if coalescent_edge.tail_node is not None:
                    if coalescent_edge.tail_node in edge_coalescent_nodes[species_edge]:
                        continue
                    # _debug_log("        Adding node {} (age={}) to set of coalescing nodes".format(self._compose_edge_desc(coalescent_edge.tail_node.edge), coalescent_edge.tail_node.age))
                    if self.is_enforce_structure_integrity:
                        for chnd in coalescent_edge.tail_node._child_nodes:
                            if chnd.edge is not coalescent_edge and chnd.edge not in valid_coalescing_lineages:
                                msg = "Invalid coalescence within structure tree edge {}: coalescent tree lineage {} cannot coalesce with lineage {} because the latter is not in the same population at this time".format(
                                        self._compose_edge_desc(species_edge), self._compose_edge_desc(coalescent_edge), self._compose_edge_desc(chnd.edge), )
                                raise error.InvalidMultispeciesCoalescentStructureError(msg)
                        valid_coalescing_lineages.add(coalescent_edge.tail_node.edge)
                    edge_coalescent_nodes[species_edge].add(coalescent_edge.tail_node)
                    new_edge = coalescent_edge.tail_node.edge
                    heapq.heappush(current_lineages, (new_edge.tail_node.age, new_edge))
                else:
                    assert False
            if current_lineages:
                edge_tail_coalescent_edges[species_edge] = set(x[1] for x in current_lineages)
        return edge_head_coalescent_edges, edge_tail_coalescent_edges, edge_coalescent_nodes

    def _compose_edge_desc(self, e):
        return "+".join(x.taxon.label for x in e.head_node.leaf_iter())
