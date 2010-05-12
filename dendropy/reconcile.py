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
Classes and Methods for working with tree reconciliation, fitting, embedding,
contained/containing etc.
"""

from dendropy import dataobject

def fit_contained_tree(contained_tree, containing_tree, contained_taxon_to_containing_taxon_map, optimize_containing_edge_lengths=True):
    """
    Fits a contained (gene) tree into a containing (population or species) tree.
    Adjusts the edge lengths / node ages of population/species tree,
    ``containing_tree``, to best explain the contained ``contained_tree``.

        ``contained_tree``
            A DendroPy ``Tree`` object representing a genealogy, where the
            ``taxon`` attribute of each leaf node is a sampled sequence.

        ``containing_tree``
            A Dendropy ``Tree`` or ``ContainingTree`` object representing a
            containing tree (e.g., species or population tree).

        ``contained_taxon_to_containing_taxon_map``
            A dictionary with Taxon objects from the gene tree as keys and the
            corresponding containing tree Taxon object (i.e., the population or
            species to which the gene belongs) as values.
    """

    # Map nodes of gene tree to containing tree, based on ``contained_taxon_to_containing_taxon_map``
    containing_taxon_to_node_map = {}
    for nd in containing_tree.leaf_iter():
        containing_taxon_to_node_map[nd.taxon] = nd
    contained_taxon_to_node_map = {}
    for nd in contained_tree.leaf_iter():
        contained_taxon_to_node_map[nd.taxon] = nd
    for gt, st in contained_taxon_to_containing_taxon_map.items():
        containing_nd = containing_taxon_to_node_map[contained_taxon_to_containing_taxon_map[gt]]
        contained_nd = contained_taxon_to_node_map[gt]
        try:
            containing_nd.contained_tree_nodes[contained_tree].append(contained_nd)
        except AttributeError:
            containing_nd.contained_tree_nodes = {contained_tree: [contained_nd]}
        except KeyError:
            containing_nd.contained_tree_nodes[contained_tree] = [contained_nd]
    if optimize_containing_edge_lengths:
        # Pre-calculate MRCA's.
        contained_taxa_mrca = {}
        contained_tree_leaf_nodes = contained_tree.leaf_nodes()
        for gi1, gnd1 in enumerate(contained_tree_leaf_nodes[:-1]):
            for gi2, gnd2 in enumerate(contained_tree_leaf_nodes[gi1:]):
                contained_taxa_mrca[(gnd1.taxon, gnd2.taxon)] = contained_tree.ancestor(gnd1, gnd2)
                contained_taxa_mrca[(gnd2.taxon, gnd1.taxon)] = contained_taxa_mrca[(gnd1.taxon, gnd2.taxon)]

        # For each split in the containing tree, find the youngest coalescent age
        # of genes between each of the daughter species/populations, and set the
        # containing node to that age.  e.g., if containing node 'A' has daughters
        # 'a1' and 'a2', find the youngest coalescent age between any gene in 'a1'
        # and 'a2', and set that as the age of A.
        if not hasattr(contained_tree.seed_node, 'age'):
            contained_tree.add_ages_to_nodes(check_prec=False)
        for containing_node in containing_tree.preorder_node_iter():
            containing_node_children = containing_node.child_nodes()
            if not containing_node_children:
                containing_node.age = 0
                continue
            containing_node_subtree_groups = [ cnd.leaf_nodes() for cnd in containing_node_children ]
            youngest_coalescence_node = None
            for xi, x in enumerate(containing_node_subtree_groups[:-1]):               # for each group of leaf nodes
                for yi, y in enumerate(containing_node_subtree_groups[xi+1:]):         # for each other group of leaf nodes
                    contained_leaves1 = sum((i.contained_tree_nodes[contained_tree] for i in x), [])      # collect gene leaves in group 1
                    contained_leaves2 = sum((i.contained_tree_nodes[contained_tree] for i in y), [])      # collect gene leaves in group 2
                    for g1 in contained_leaves1:                                     # for each leaf in group 1
                        for g2 in contained_leaves2:                                 # for each leaf in group 2
                            mrca_node = contained_taxa_mrca[(g1.taxon, g2.taxon)]
                            if youngest_coalescence_node is None or mrca_node.age < youngest_coalescence_node.age:
                                youngest_coalescence_node = mrca_node
            containing_node.age = youngest_coalescence_node.age
            #if not hasattr(containing_node, "contained_tree_nodes"):
            #    containing_node.contained_tree_nodes = []
            #containing_node.contained_tree_nodes.append(youngest_coalescence_node)

            # Make sure that the ancestors of this containing node are not younger
            # than it.
            parent_node = containing_node.parent_node
            while parent_node is not None:
                if parent_node.age < containing_node.age:
                    parent_node.age = containing_node.age
                parent_node = parent_node.parent_node

        # Set the edge lengths.
        for nd in containing_tree.preorder_node_iter():
            if nd.parent_node is not None:
                nd.edge.length = nd.parent_node.age - nd.age
                if nd.edge.length < 0:
                    nd.edge.length = 0
    else:
        if not hasattr(contained_tree.seed_node, 'age'):
            contained_tree.add_ages_to_nodes(check_prec=False)
        if not hasattr(containing_tree.seed_node, 'age'):
            containing_tree.add_ages_to_nodes(check_prec=False)

    # Map the edges
    num_deep_coalescences = 0
    for containing_edge in containing_tree.postorder_edge_iter():
        if not hasattr(containing_edge, 'tail_contained_edges'):
            containing_edge.tail_contained_edges = {}
        if not hasattr(containing_edge, 'head_contained_edges'):
            containing_edge.head_contained_edges = {}
        if not hasattr(containing_edge, 'uncoalesced_edges'):
            containing_edge.num_uncoalesced_edges = {}
        child_nodes = containing_edge.head_node.child_nodes()
        if not child_nodes:
            containing_edge.head_contained_edges[contained_tree] = set([n.edge for n in containing_edge.head_node.contained_tree_nodes[contained_tree]])
        else:
            containing_edge.head_contained_edges[contained_tree] = set()
            for n in child_nodes:
                containing_edge.head_contained_edges[contained_tree].update(n.edge.tail_contained_edges[contained_tree])
            containing_edge.num_uncoalesced_edges[contained_tree] = len(containing_edge.head_contained_edges[contained_tree]) - len(child_nodes)
            num_deep_coalescences += containing_edge.num_uncoalesced_edges[contained_tree]
        if containing_edge.tail_node is None:
            containing_edge.tail_contained_edges[contained_tree] = containing_edge.head_contained_edges[contained_tree]
            continue
        containing_edge.tail_contained_edges[contained_tree] = set()
        for contained_edge in containing_edge.head_contained_edges[contained_tree]:
            #if contained_edge.tail_node is None:
            #    break
            remaining = containing_edge.tail_node.age - contained_edge.tail_node.age
            while remaining > 0:
                try:
                    contained_edge = contained_edge.tail_node.edge
                except AttributeError:
                    #contained_edge = None
                    break
                remaining -= contained_edge.length
            if contained_edge is not None:
                containing_edge.tail_contained_edges[contained_tree].add(contained_edge)

    if not hasattr(containing_tree, 'num_deep_coalescences'):
        containing_tree.num_deep_coalescences = {}
    containing_tree.num_deep_coalescences[contained_tree] = num_deep_coalescences

    return containing_tree

def reconciliation_discordance(gene_tree, species_tree):
    """
    Given two trees (with splits encoded), this returns the number of gene
    duplications implied by the gene tree reconciled on the species tree, based
    on the algorithm described here:

        Goodman, M. J. Czelnusiniak, G. W. Moore, A. E. Romero-Herrera, and
        G. Matsuda. 1979. Fitting the gene lineage into its species lineage,
        a parsimony strategy illustrated bu cladograms constructed from globin
        sequences. Syst. Zool. 19: 99-113.

        Maddison, W. P. 1997. Gene trees in species dataobject. Syst. Biol. 46:
        523-536.

    This function requires that the gene tree and species tree *have the same
    leaf set*. Note that for correct results,

        (a) trees must be rooted (i.e., is_rooted = True)
        (b) split masks must have been added as rooted (i.e., when
            encode_splits was called, is_rooted must have been set to True)

    """
    taxa_mask = species_tree.taxon_set.all_taxa_bitmask()

    species_node_gene_nodes = {}
    gene_node_species_nodes = {}

    for gnd in gene_tree.postorder_node_iter():
        gn_children = gnd.child_nodes()
        if len(gn_children) > 0:
            ssplit = 0
            for gn_child in gn_children:
                ssplit = ssplit | gene_node_species_nodes[gn_child].edge.split_bitmask
            sanc = species_tree.mrca(start_node=species_tree.seed_node, split_bitmask=ssplit)
            gene_node_species_nodes[gnd] = sanc
            if sanc not in species_node_gene_nodes:
                species_node_gene_nodes[sanc] = []
            species_node_gene_nodes[sanc].append(gnd)
        else:
            gene_node_species_nodes[gnd] = species_tree.find_node(lambda x : x.taxon == gnd.taxon)

    contained_gene_lineages = {}
    for snd in species_tree.postorder_node_iter():
        if snd in species_node_gene_nodes:
            for gnd in species_node_gene_nodes[snd]:
                for gnd_child in gnd.child_nodes():
                    sanc = gene_node_species_nodes[gnd_child]
                    p = sanc
                    while p is not None and p != snd:
                        if p.edge not in contained_gene_lineages:
                            contained_gene_lineages[p.edge] = 0
                        contained_gene_lineages[p.edge] += 1
                        p = p.parent_node

    dc = 0
    for v in contained_gene_lineages.values():
        dc += v - 1

    return dc

def monophyletic_partition_discordance(tree, taxon_set_partition):
    """
    Returns the number of deep coalescences on tree `tree` that would result
    if the taxa in `tax_sets` formed K mutually-exclusive monophyletic groups,
    where K = len(tax_sets)
    `taxon_set_partition` == TaxonSetPartition
    """

    tax_sets = taxon_set_partition.subsets()
    dc_tree = dataobject.Tree()
    dc_tree.taxon_set = dataobject.TaxonSet()

    for t in range(len(tax_sets)):
        dc_tree.taxon_set.append(dataobject.Taxon(label=str(t)))

    def _get_dc_taxon(nd):
        for idx, tax_set in enumerate(tax_sets):
            if nd.taxon in tax_set:
                return dc_tree.taxon_set[idx]
        assert "taxon not found in partition: '%s'" % nd.taxon.label

    src_dc_map = {}
    for snd in tree.postorder_node_iter():
        nnd = dataobject.Node()
        src_dc_map[snd] = nnd
        children = snd.child_nodes()
        if len(children) == 0:
            nnd.taxon = _get_dc_taxon(snd)
        else:
            taxa_set = []
            for cnd in children:
                dc_node = src_dc_map[cnd]
                if len(dc_node.child_nodes()) > 1:
                    nnd.add_child(dc_node)
                else:
                    ctax = dc_node.taxon
                    if ctax is not None and ctax not in taxa_set:
                        taxa_set.append(ctax)
                    del src_dc_map[cnd]
            if len(taxa_set) > 1:
                for t in taxa_set:
                    cnd = dataobject.Node()
                    cnd.taxon = t
                    nnd.add_child(cnd)
            else:
                if len(nnd.child_nodes()) == 0:
                    nnd.taxon = taxa_set[0]
                elif len(taxa_set) == 1:
                    cnd = dataobject.Node()
                    cnd.taxon = taxa_set[0]
                    nnd.add_child(cnd)
    dc_tree.seed_node = nnd
    return len(dc_tree.leaf_nodes()) - len(tax_sets)

