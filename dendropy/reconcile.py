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

def fit_gene_tree(gene_tree, containing_tree, gene_taxon_to_containing_taxon_map):
    """
    Fits a gene tree into a containing (population or species) tree.
    Adjusts the node ages of population/species tree, ``containing_tree``, to
    best explain the contained ``gene_tree``.

        ``gene_tree``
            A DendroPy ``Tree`` object representing a genealogy, where the
            ``taxon`` attribute of each leaf node is a sampled sequence.

        ``containing_tree``
            A DendroPy ``Tree`` object representing a containing tree (e.g., a
            species or population tree).

        ``gene_taxon_to_containing_taxon_map``
            A dictionary with Taxon objects from the gene tree as keys and the
            corresponding containing tree Taxon object (i.e., the population or
            species to which the gene belongs) as values.

    In addition to edge lengths being set, the nodes of ``containing_tree`` will have
    two new or reassigned attributes:

        ``age``
            Age of the node, in terms of time units back from present.

        ``gene_tree_nodes``
            List of nodes of the gene tree that coalesce at this node.

    """

    # Map nodes of gene tree to containing tree, based on ``gene_taxon_to_containing_taxon_map``
    containing_taxon_to_node_map = {}
    for nd in containing_tree.leaf_iter():
        containing_taxon_to_node_map[nd.taxon] = nd
    gene_taxon_to_node_map = {}
    for nd in gene_tree.leaf_iter():
        gene_taxon_to_node_map[nd.taxon] = nd
    for gt, st in gene_taxon_to_containing_taxon_map.items():
        containing_nd = containing_taxon_to_node_map[gene_taxon_to_containing_taxon_map[gt]]
        gene_nd = gene_taxon_to_node_map[gt]
        try:
            containing_nd.gene_tree_nodes.append(gene_nd)
        except AttributeError:
            containing_nd.gene_tree_nodes = [gene_nd]

    # Pre-calculate MRCA's.
    gene_taxa_mrca = {}
    gene_tree_leaf_nodes = gene_tree.leaf_nodes()
    for gi1, gnd1 in enumerate(gene_tree_leaf_nodes[:-1]):
        for gi2, gnd2 in enumerate(gene_tree_leaf_nodes[gi1:]):
            gene_taxa_mrca[(gnd1.taxon, gnd2.taxon)] = gene_tree.ancestor(gnd1, gnd2)
            gene_taxa_mrca[(gnd2.taxon, gnd1.taxon)] = gene_taxa_mrca[(gnd1.taxon, gnd2.taxon)]

    # For each split in the containing tree, find the youngest coalescent age
    # of genes between each of the daughter species/populations, and set the
    # containing node to that age.  e.g., if containing node 'A' has daughters
    # 'a1' and 'a2', find the youngest coalescent age between any gene in 'a1'
    # and 'a2', and set that as the age of A.
    gene_tree.add_ages_to_nodes(check_prec=False)
    for containing_node in containing_tree.preorder_node_iter():
        containing_node_children = containing_node.child_nodes()
        if not containing_node_children:
            containing_node.age = 0
            continue
        containing_node_subtree_groups = [ cnd.leaf_nodes() for cnd in containing_node_children ]
        youngest_coalescence_node = None
        for xi, x in enumerate(containing_node_subtree_groups[:-1]):               # for each group of leaf nodes
            for yi, y in enumerate(containing_node_subtree_groups[xi+1:]):         # for each other group of leaf nodes
                gene_leaves1 = sum((i.gene_tree_nodes for i in x), [])      # collect gene leaves in group 1
                gene_leaves2 = sum((i.gene_tree_nodes for i in y), [])      # collect gene leaves in group 2
                for g1 in gene_leaves1:                                     # for each leaf in group 1
                    for g2 in gene_leaves2:                                 # for each leaf in group 2
                        mrca_node = gene_taxa_mrca[(g1.taxon, g2.taxon)]
                        if youngest_coalescence_node is None or mrca_node.age < youngest_coalescence_node.age:
                            youngest_coalescence_node = mrca_node
        containing_node.age = youngest_coalescence_node.age
        if not hasattr(containing_node, "gene_tree_nodes"):
            containing_node.gene_tree_nodes = []
        containing_node.gene_tree_nodes.append(youngest_coalescence_node)

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

