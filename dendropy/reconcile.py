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
Classes and Methods for working with tree reconciliation, fitting, embedding,
contained/containing etc.
"""

from dendropy import dataobject
from dendropy import coalescent
from dendropy import treemanip

class ContainingTree(dataobject.Tree):
    """
    A "containing tree" is a (usually rooted) tree data structure within which
    other trees are "embedded". For example, species trees and their embedded
    gene trees; host trees and their embedded parasite trees; biogeographical
    "area" trees and their embedded species or taxon trees.
    """

    def __init__(self,
            containing_tree,
            embedded_taxon_set,
            embedded_to_containing_taxon_map,
            embedded_trees=None,
            fit_containing_edge_lengths=True,
            collapse_empty_edges=True,
            ultrametricity_check_prec=False,
            ignore_root_deep_coalescences=True,
            **kwargs):
        """
        __init__ converts `self` to ContainingTree class, embedding the trees
        given in the list, ``embedded_trees.``


        Mandatory Arguments:

            ``containing_tree``
                A ``Tree`` or ``Tree``-like object that describes the topological
                constraints or conditions of the containing tree (e.g., species,
                host, or biogeographical area trees).

            ``embedded_taxon_set``
                A ``TaxonSet`` object that will be used to manage the taxa of
                the embedded trees.

            ``embedded_to_containing_taxon_map``
                A ``TaxonSetMapping`` object mapping ``Taxon`` objects in the
                embedded ``TaxonSet`` to corresponding ``Taxon`` objects in the
                containing tree.

        Optional Arguments:

            ``embedded_trees``
                An iterable container of ``Tree`` or ``Tree``-like objects that
                will be embedded into ``containing_tree``; e.g. gene or
                parasite trees.

            ``fit_containing_edge_lengths``
                If ``True`` [default], then the branch lengths of
                ``containing_tree`` will be adjusted to fit the embedded tree
                as they are added. Otherwise, the containing tree edge lengths
                will not be changed.

            ``collapse_empty_edges``
                If ``True`` [default], after edge lengths are adjusted,
                zero-length branches will be collapsed.

            ``ultrametricity_check_prec``
                If ``False`` [default], then trees will not be checked for
                ultrametricity. Otherwise this is the threshold within which
                all node to tip distances for sister nodes must be equal.

            ``ignore_root_deep_coalescences``
                If ``True`` [default], then deep coalescences in the root will
                not be counted.

        Other Keyword Arguments: Will be passed to Tree().

    """
        if "taxon_set" not in kwargs:
            kwargs["taxon_set"] = containing_tree.taxon_set
        dataobject.Tree.__init__(self, containing_tree, **kwargs)
        self.original_tree = containing_tree
        for edge in self.postorder_edge_iter():
            edge.head_embedded_edges = {}
            edge.tail_embedded_edges = {}
            edge.containing_taxa = set()
            edge.embedded_taxa = set()
        self._embedded_taxon_set = embedded_taxon_set
        self._embedded_to_containing_taxon_map = None
        self._embedded_trees = None
        self._set_embedded_to_containing_taxon_map(embedded_to_containing_taxon_map)
        self.fit_containing_edge_lengths = fit_containing_edge_lengths
        self.collapse_empty_edges = collapse_empty_edges
        self.ultrametricity_check_prec = ultrametricity_check_prec
        self.ignore_root_deep_coalescences = ignore_root_deep_coalescences
        if embedded_trees:
            self._set_embedded_trees(embedded_trees)
        if self.embedded_trees:
            self.rebuild(rebuild_taxa=False)

    def _set_embedded_taxon_set(self, taxon_set):
        self._embedded_taxon_set = taxon_set

    def _get_embedded_taxon_set(self):
        if self._embedded_taxon_set is None:
            self._embedded_taxon_set = dataobject.TaxonSet()
        return self._embedded_taxon_set

    embedded_taxon_set = property(_get_embedded_taxon_set)

    def _set_embedded_to_containing_taxon_map(self, embedded_to_containing_taxon_map):
        """
        Sets mapping of ``Taxon`` objects of the genes/parasite/etc. to that of
        the population/species/host/etc.
        Creates mapping (e.g., species to genes) and decorates edges of self
        with sets of both containing ``Taxon`` objects and the embedded
        ``Taxon`` objects that map to them.
        """
        if isinstance(embedded_to_containing_taxon_map, dataobject.TaxonSetMapping):
            if self._embedded_taxon_set is not embedded_to_containing_taxon_map.domain_taxon_set:
                raise ValueError("Domain TaxonSet of TaxonSetMapping ('domain_taxon_set') not the same as 'embedded_taxon_set' TaxonSet")
            self._embedded_to_containing_taxon_map = embedded_to_containing_taxon_map
        else:
            self._embedded_to_containing_taxon_map = dataobject.TaxonSetMapping(
                    mapping_dict=embedded_to_containing_taxon_map,
                    domain_taxon_set=self.embedded_taxon_set,
                    range_taxon_set=self.taxon_set)
        self.build_edge_taxa_sets()

    def _get_embedded_to_containing_taxon_map(self):
        return self._embedded_to_containing_taxon_map

    embedded_to_containing_taxon_map = property(_get_embedded_to_containing_taxon_map)

    def _set_embedded_trees(self, trees):
        if hasattr(trees, 'taxon_set'):
            if self._embedded_taxon_set is None:
                self._embedded_taxon_set = trees.taxon_set
            elif self._embedded_taxon_set is not trees.taxon_set:
                raise ValueError("'embedded_taxon_set' of ContainingTree is not the same TaxonSet object of 'embedded_trees'")
        self._embedded_trees = dataobject.TreeList(trees, taxon_set=self._embedded_taxon_set)
        if self._embedded_taxon_set is None:
            self._embedded_taxon_set = self._embedded_trees.taxon_set

    def _get_embedded_trees(self):
        if self._embedded_trees is None:
            self._embedded_trees = dataobject.TreeList(taxon_set=self._embedded_taxon_set)
        return self._embedded_trees

    embedded_trees = property(_get_embedded_trees)

    def _get_containing_to_embedded_taxa_map(self):
        return self._embedded_to_containing_taxon_map.reverse

    containing_to_embedded_taxa_map = property(_get_containing_to_embedded_taxa_map)

    def clear(self):
        """
        Clears all embedded trees and mapped edges.
        """
        self.embedded_trees = dataobject.TreeList(taxon_set=self._embedded_to_containing_taxon_map.domain_taxa)
        self.clear_embedded_edges()

    def clear_embedded_edges(self):
        """
        Clears all embedded mapped edges.
        """
        for edge in self.postorder_edge_iter():
            edge.head_embedded_edges = {}
            edge.tail_embedded_edges = {}

    def fit_edge_lengths(self, embedded_trees):
        """
        Recalculate node ages / edge lengths of containing tree to accomodate
        embedded trees.
        """

        # set the ages
        for node in self.postorder_node_iter():
            if node.is_internal():
                disjunct_leaf_set_list_split_bitmasks = []
                for i in node.child_nodes():
                    disjunct_leaf_set_list_split_bitmasks.append(self.taxon_set.get_taxa_bitmask(taxa=i.edge.containing_taxa))
                min_age = float('inf')
                for et in embedded_trees:
                    min_age = self._find_youngest_intergroup_age(et, disjunct_leaf_set_list_split_bitmasks, min_age)
                node.age = max( [min_age] + [cn.age for cn in node.child_nodes()] )
            else:
                node.age = 0

        # set the corresponding edge lengths
        self.set_edge_lengths_from_node_ages()

        # collapse 0-length branches
        if self.collapse_empty_edges:
           self.collapse_unweighted_edges()

    def rebuild(self, rebuild_taxa=True):
        """
        Recalculate edge taxa sets, node ages / edge lengths of containing
        tree, and embed edges of embedded trees.
        """
        if rebuild_taxa:
            self.build_edge_taxa_sets()
        if self.fit_containing_edge_lengths:
            self.fit_edge_lengths(self.embedded_trees)
        self.clear_embedded_edges()
        for et in self.embedded_trees:
            self.embed_tree(et)

    def embed_tree(self, embedded_tree):
        """
        Map edges of embedded tree into containing tree (i.e., self).
        """
        if self.seed_node.age is None:
            self.calc_node_ages(check_prec=self.ultrametricity_check_prec)
        if embedded_tree not in self.embedded_trees:
            self.embedded_trees.append(embedded_tree)
        if embedded_tree.seed_node.age is None:
            embedded_tree.calc_node_ages(check_prec=self.ultrametricity_check_prec)
        contained_leaves = embedded_tree.leaf_nodes()
        taxon_to_contained = {}
        for nd in contained_leaves:
            containing_taxon = self.embedded_to_containing_taxon_map[nd.taxon]
            x = taxon_to_contained.setdefault(containing_taxon, set())
            x.add(nd.edge)
        for containing_edge in self.postorder_edge_iter():
            if containing_edge.is_terminal():
                containing_edge.head_embedded_edges[embedded_tree] = taxon_to_contained[containing_edge.head_node.taxon]
            else:
                containing_edge.head_embedded_edges[embedded_tree] = set()
                for nd in containing_edge.head_node.child_nodes():
                    containing_edge.head_embedded_edges[embedded_tree].update(nd.edge.tail_embedded_edges[embedded_tree])

            if containing_edge.tail_node is None:
                if containing_edge.length is not None:
                    target_age =  containing_edge.head_node.age + containing_edge.length
                else:
                    # assume all coalesce?
                    containing_edge.tail_embedded_edges[embedded_tree] = set(embedded_tree.seed_node.edge)
                    continue
            else:
                target_age = containing_edge.tail_node.age

            containing_edge.tail_embedded_edges[embedded_tree] = set()
            for embedded_edge in containing_edge.head_embedded_edges[embedded_tree]:
                if embedded_edge.tail_node is not None:
                    remaining = target_age - embedded_edge.tail_node.age
                elif embedded_edge.length is not None:
                    remaining = target_age - (embedded_edge.head_node.age + embedded_age.length)
                else:
                    continue
                while remaining > 0:
                    if embedded_edge.tail_node is not None:
                        embedded_edge = embedded_edge.tail_node.edge
                    else:
                        if embedded_edge.length is not None and (remaining - embedded_edge.length) <= 0:
                            embedded_edge = None
                            remaining = 0
                            break
                        else:
                            remaining = 0
                            break
                    if embedded_edge and remaining > 0:
                        remaining -= embedded_edge.length
                if embedded_edge is not None:
                    containing_edge.tail_embedded_edges[embedded_tree].add(embedded_edge)

    def build_edge_taxa_sets(self):
        """
        Rebuilds sets of containing and corresponding embedded taxa at each
        edge.
        """
        for edge in self.postorder_edge_iter():
            if edge.is_terminal():
                edge.containing_taxa = set([edge.head_node.taxon])
            else:
                edge.containing_taxa = set()
                for i in edge.head_node.child_nodes():
                    edge.containing_taxa.update(i.edge.containing_taxa)
            edge.embedded_taxa = set()
            for t in edge.containing_taxa:
                edge.embedded_taxa.update(self.containing_to_embedded_taxa_map[t])

    def num_deep_coalescences(self):
        """
        Returns total number of deep coalescences of the embedded trees.
        """
        return sum(self.deep_coalescences().values())

    def deep_coalescences(self):
        """
        Returns dictionary where the embedded trees are keys, and the number of
        deep coalescences corresponding to the tree are values.
        """
        dc = {}
        for tree in self.embedded_trees:
            for edge in self.postorder_edge_iter():
                if edge.tail_node is None and self.ignore_root_deep_coalescences:
                    continue
                try:
                    dc[tree] += len(edge.tail_embedded_edges[tree]) - 1
                except KeyError:
                    dc[tree] = len(edge.tail_embedded_edges[tree]) - 1
        return dc

    def embed_contained_kingman(self, rng=None, pop_size_attr='pop_size', label=None):
        """
        Simulates, *embeds*, and returns a "censored" (Kingman) neutral coalescence tree
        conditional on self.

            ``rng``
                Random number generator to use. If ``None``, the default will
                be used.

            ``pop_size_attr``
                Name of attribute of self's edges that specify the population
                size. If this attribute does not exist, then the population
                size is taken to be 1.

        Note that all edge-associated taxon sets must be up-to-date (otherwise,
        ``build_edge_taxa_sets()`` should be called).
        """
        et = self.simulate_contained_kingman(rng=rng,
                pop_size_attr=pop_size_attr,
                label=label)
        self.embed_tree(et)
        return et

    def simulate_contained_kingman(self, rng=None, pop_size_attr='pop_size', label=None):
        """
        Simulates and returns a "censored" (Kingman) neutral coalescence tree
        conditional on self.

            ``rng``
                Random number generator to use. If ``None``, the default will
                be used.

            ``pop_size_attr``
                Name of attribute of self's edges that specify the population
                size. If this attribute does not exist, then the population
                size is taken to be 1.

        Note that all edge-associated taxon sets must be up-to-date (otherwise,
        ``build_edge_taxa_sets()`` should be called), and that the tree
        is *not* added to the set of embedded trees. For the latter, call
        ``embed_contained_kingman``.
        """

        # Dictionary that maps nodes of containing tree to list of
        # corresponding nodes on gene tree, initially populated with leaf
        # nodes.
        embedded_nodes = {}
        for nd in self.leaf_iter():
            embedded_nodes[nd] = []
            for gt in nd.edge.embedded_taxa:
                gn = dataobject.Node(taxon=gt)
                embedded_nodes[nd].append(gn)

        # Generate the tree structure
        for edge in self.postorder_edge_iter():
            if edge.head_node.parent_node is None:
                # root: run unconstrained coalescence until just one gene node
                # remaining
                if len(embedded_nodes[edge.head_node]) > 1:
                    final = coalescent.coalesce(nodes=embedded_nodes[edge.head_node],
                            pop_size=pop_size,
                            period=None,
                            rng=rng)
                else:
                    final = embedded_nodes[edge.head_node]
            else:
                # run until next coalescence event, as determined by this edge
                # size.
                if hasattr(edge, pop_size_attr):
                    pop_size = getattr(edge, pop_size_attr)
                else:
                    pop_size = 1
                remaining = coalescent.coalesce(nodes=embedded_nodes[edge.head_node],
                        pop_size=pop_size,
                        period=edge.length,
                        rng=rng)
                try:
                    embedded_nodes[edge.tail_node].extend(remaining)
                except KeyError:
                    embedded_nodes[edge.tail_node] = remaining

        # Create and return the full tree
        embedded_tree = dataobject.Tree(taxon_set=self.embedded_taxon_set, label=label)
        embedded_tree.seed_node = final[0]
        embedded_tree.is_rooted = True
        return embedded_tree

    def _find_youngest_intergroup_age(self, embedded_tree, disjunct_leaf_set_list_split_bitmasks, starting_min_age=None):
        """
        Find the age of the youngest MRCA of disjunct leaf sets.
        """
        if starting_min_age is None:
            starting_min_age = float('inf')
        if embedded_tree.seed_node.age is None:
            embedded_tree.calc_node_ages(check_prec=self.ultrametricity_check_prec)
        if not hasattr(embedded_tree, 'split_edges'):
            embedded_tree.update_splits()
        for nd in embedded_tree.age_order_node_iter(include_leaves=False):
            if nd.age > starting_min_age:
                break
            prev_intersections = False
            for bm in disjunct_leaf_set_list_split_bitmasks:
                if bm & nd.edge.split_bitmask:
                    if prev_intersections:
                        return nd.age
                    prev_intersections = True
        return starting_min_age

    def write_as_mesquite(self, out, **kwargs):
        """
        For debugging purposes, write out a Mesquite-format file.
        """
        from dendropy.dataio import nexuswriter
        nw = nexuswriter.NexusWriter(**kwargs)
        nw.is_write_block_titles = True
        out.write("#NEXUS\n\n")
        nw.write_taxa_block(self.taxon_set, out)
        out.write('\n')
        nw.write_taxa_block(self.embedded_trees.taxon_set, out)
        if self.embedded_trees.taxon_set.label:
            domain_title = self.embedded_trees.taxon_set.label
        else:
            domain_title = self.embedded_trees.taxon_set.oid
        embedded_taxon_set = self.embedded_trees.taxon_set
        embedded_label = self.embedded_trees.label
        out.write('\n')
        self._embedded_to_containing_taxon_map.write_mesquite_association_block(out)
        out.write('\n')
        nw.write_trees_block(dataobject.TreeList([self], taxon_set=self.taxon_set), out)
        out.write('\n')
        nw.write_trees_block(dataobject.TreeList(self.embedded_trees, taxon_set=embedded_taxon_set, label=embedded_label), out)
        out.write('\n')

def __RETIRED__fit_contained_tree(contained_tree, containing_tree, contained_taxon_to_containing_taxon_map, optimize_containing_edge_lengths=True):
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

    # Optimize containing branch lengths to minimize numbers of deep coalescences.
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
        if contained_tree.seed_node.age is None:
            contained_tree.calc_node_ages(check_prec=False)
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
        if contained_tree.seed_node.age is None:
            contained_tree.calc_node_ages(check_prec=False)
        if containing_tree.seed_node.age is None:
            containing_tree.calc_node_ages(check_prec=False)

    # Map the edges
    num_deep_coalescences = 0
    for containing_edge in containing_tree.postorder_edge_iter():
        if not hasattr(containing_edge, 'tail_embedded_edges'):
            containing_edge.tail_embedded_edges = {}
        if not hasattr(containing_edge, 'head_embedded_edges'):
            containing_edge.head_embedded_edges = {}
        if not hasattr(containing_edge, 'uncoalesced_edges'):
            containing_edge.num_uncoalesced_edges = {}
        child_nodes = containing_edge.head_node.child_nodes()
        if not child_nodes:
            containing_edge.head_embedded_edges[contained_tree] = set([n.edge for n in containing_edge.head_node.contained_tree_nodes[contained_tree]])
        else:
            containing_edge.head_embedded_edges[contained_tree] = set()
            for n in child_nodes:
                containing_edge.head_embedded_edges[contained_tree].update(n.edge.tail_embedded_edges[contained_tree])
            containing_edge.num_uncoalesced_edges[contained_tree] = len(containing_edge.head_embedded_edges[contained_tree]) - len(child_nodes)
            num_deep_coalescences += containing_edge.num_uncoalesced_edges[contained_tree]
        if containing_edge.tail_node is None:
            containing_edge.tail_embedded_edges[contained_tree] = containing_edge.head_embedded_edges[contained_tree]
            continue
        containing_edge.tail_embedded_edges[contained_tree] = set()
        for embedded_edge in containing_edge.head_embedded_edges[contained_tree]:
            #if embedded_edge.tail_node is None:
            #    break
            remaining = containing_edge.tail_node.age - embedded_edge.tail_node.age
            while remaining > 0:
                try:
                    embedded_edge = embedded_edge.tail_node.edge
                except AttributeError:
                    #embedded_edge = None
                    break
                remaining -= embedded_edge.length
            if embedded_edge is not None:
                containing_edge.tail_embedded_edges[contained_tree].add(embedded_edge)

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
        a parsimony strategy illustrated by cladograms constructed from globin
        sequences. Syst. Zool. 19: 99-113.

        Maddison, W. P. 1997. Gene trees in species trees. Syst. Biol. 46:
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

