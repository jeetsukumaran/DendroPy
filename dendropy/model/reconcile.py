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
Classes and Methods for working with tree reconciliation, fitting, embedding,
contained/containing etc.
"""

import dendropy
from dendropy.model import coalescent

class ContainingTree(dendropy.Tree):
    """
    A "containing tree" is a (usually rooted) tree data structure within which
    other trees are "contained". For example, species trees and their contained
    gene trees; host trees and their contained parasite trees; biogeographical
    "area" trees and their contained species or taxon trees.
    """

    def __init__(self,
            containing_tree,
            contained_taxon_namespace,
            contained_to_containing_taxon_map,
            contained_trees=None,
            fit_containing_edge_lengths=True,
            collapse_empty_edges=True,
            ultrametricity_precision=False,
            ignore_root_deep_coalescences=True,
            **kwargs):
        """
        __init__ converts ``self`` to ContainingTree class, embedding the trees
        given in the list, ``contained_trees.``


        Mandatory Arguments:

            ``containing_tree``
                A |Tree| or |Tree|-like object that describes the topological
                constraints or conditions of the containing tree (e.g., species,
                host, or biogeographical area trees).

            ``contained_taxon_namespace``
                A |TaxonNamespace| object that will be used to manage the taxa of
                the contained trees.

            ``contained_to_containing_taxon_map``
                A |TaxonNamespaceMapping| object mapping |Taxon| objects in the
                contained |TaxonNamespace| to corresponding |Taxon| objects in the
                containing tree.

        Optional Arguments:

            ``contained_trees``
                An iterable container of |Tree| or |Tree|-like objects that
                will be contained into ``containing_tree``; e.g. gene or
                parasite trees.

            ``fit_containing_edge_lengths``
                If |True| [default], then the branch lengths of
                ``containing_tree`` will be adjusted to fit the contained tree
                as they are added. Otherwise, the containing tree edge lengths
                will not be changed.

            ``collapse_empty_edges``
                If |True| [default], after edge lengths are adjusted,
                zero-length branches will be collapsed.

            ``ultrametricity_precision``
                If |False| [default], then trees will not be checked for
                ultrametricity. Otherwise this is the threshold within which
                all node to tip distances for sister nodes must be equal.

            ``ignore_root_deep_coalescences``
                If |True| [default], then deep coalescences in the root will
                not be counted.

        Other Keyword Arguments: Will be passed to Tree().

    """
        if "taxon_namespace" not in kwargs:
            kwargs["taxon_namespace"] = containing_tree.taxon_namespace
        dendropy.Tree.__init__(self,
                containing_tree,
                taxon_namespace=containing_tree.taxon_namespace)
        self.original_tree = containing_tree
        for edge in self.postorder_edge_iter():
            edge.head_contained_edges = {}
            edge.tail_contained_edges = {}
            edge.containing_taxa = set()
            edge.contained_taxa = set()
        self._contained_taxon_namespace = contained_taxon_namespace
        self._contained_to_containing_taxon_map = None
        self._contained_trees = None
        self._set_contained_to_containing_taxon_map(contained_to_containing_taxon_map)
        self.fit_containing_edge_lengths = fit_containing_edge_lengths
        self.collapse_empty_edges = collapse_empty_edges
        self.ultrametricity_precision = ultrametricity_precision
        self.ignore_root_deep_coalescences = ignore_root_deep_coalescences
        if contained_trees:
            self._set_contained_trees(contained_trees)
        if self.contained_trees:
            self.rebuild(rebuild_taxa=False)

    def _set_contained_taxon_namespace(self, taxon_namespace):
        self._contained_taxon_namespace = taxon_namespace

    def _get_contained_taxon_namespace(self):
        if self._contained_taxon_namespace is None:
            self._contained_taxon_namespace = dendropy.TaxonNamespace()
        return self._contained_taxon_namespace

    contained_taxon_namespace = property(_get_contained_taxon_namespace)

    def _set_contained_to_containing_taxon_map(self, contained_to_containing_taxon_map):
        """
        Sets mapping of |Taxon| objects of the genes/parasite/etc. to that of
        the population/species/host/etc.
        Creates mapping (e.g., species to genes) and decorates edges of self
        with sets of both containing |Taxon| objects and the contained
        |Taxon| objects that map to them.
        """
        if isinstance(contained_to_containing_taxon_map, dendropy.TaxonNamespaceMapping):
            if self._contained_taxon_namespace is not contained_to_containing_taxon_map.domain_taxon_namespace:
                raise ValueError("Domain TaxonNamespace of TaxonNamespaceMapping ('domain_taxon_namespace') not the same as 'contained_taxon_namespace' TaxonNamespace")
            self._contained_to_containing_taxon_map = contained_to_containing_taxon_map
        else:
            self._contained_to_containing_taxon_map = dendropy.TaxonNamespaceMapping(
                    mapping_dict=contained_to_containing_taxon_map,
                    domain_taxon_namespace=self.contained_taxon_namespace,
                    range_taxon_namespace=self.taxon_namespace)
        self.build_edge_taxa_sets()

    def _get_contained_to_containing_taxon_map(self):
        return self._contained_to_containing_taxon_map

    contained_to_containing_taxon_map = property(_get_contained_to_containing_taxon_map)

    def _set_contained_trees(self, trees):
        if hasattr(trees, 'taxon_namespace'):
            if self._contained_taxon_namespace is None:
                self._contained_taxon_namespace = trees.taxon_namespace
            elif self._contained_taxon_namespace is not trees.taxon_namespace:
                raise ValueError("'contained_taxon_namespace' of ContainingTree is not the same TaxonNamespace object of 'contained_trees'")
        self._contained_trees = dendropy.TreeList(trees, taxon_namespace=self._contained_taxon_namespace)
        if self._contained_taxon_namespace is None:
            self._contained_taxon_namespace = self._contained_trees.taxon_namespace

    def _get_contained_trees(self):
        if self._contained_trees is None:
            self._contained_trees = dendropy.TreeList(taxon_namespace=self._contained_taxon_namespace)
        return self._contained_trees

    contained_trees = property(_get_contained_trees)

    def _get_containing_to_contained_taxa_map(self):
        return self._contained_to_containing_taxon_map.reverse

    containing_to_contained_taxa_map = property(_get_containing_to_contained_taxa_map)

    def clear(self):
        """
        Clears all contained trees and mapped edges.
        """
        self.contained_trees = dendropy.TreeList(taxon_namespace=self._contained_to_containing_taxon_map.domain_taxa)
        self.clear_contained_edges()

    def clear_contained_edges(self):
        """
        Clears all contained mapped edges.
        """
        for edge in self.postorder_edge_iter():
            edge.head_contained_edges = {}
            edge.tail_contained_edges = {}

    def fit_edge_lengths(self, contained_trees):
        """
        Recalculate node ages / edge lengths of containing tree to accomodate
        contained trees.
        """

        # set the ages
        for node in self.postorder_node_iter():
            if node.is_internal():
                disjunct_leaf_set_list_split_bitmasks = []
                for i in node.child_nodes():
                    disjunct_leaf_set_list_split_bitmasks.append(self.taxon_namespace.taxa_bitmask(taxa=i.edge.containing_taxa))
                min_age = float('inf')
                for et in contained_trees:
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
        tree, and embed edges of contained trees.
        """
        if rebuild_taxa:
            self.build_edge_taxa_sets()
        if self.fit_containing_edge_lengths:
            self.fit_edge_lengths(self.contained_trees)
        self.clear_contained_edges()
        for et in self.contained_trees:
            self.embed_tree(et)

    def embed_tree(self, contained_tree):
        """
        Map edges of contained tree into containing tree (i.e., self).
        """
        if self.seed_node.age is None:
            self.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        if contained_tree not in self.contained_trees:
            self.contained_trees.append(contained_tree)
        if self.fit_containing_edge_lengths:
            self.fit_edge_lengths(self.contained_trees)
        if contained_tree.seed_node.age is None:
            contained_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        contained_leaves = contained_tree.leaf_nodes()
        taxon_to_contained = {}
        for nd in contained_leaves:
            containing_taxon = self.contained_to_containing_taxon_map[nd.taxon]
            x = taxon_to_contained.setdefault(containing_taxon, set())
            x.add(nd.edge)
        for containing_edge in self.postorder_edge_iter():
            if containing_edge.is_terminal():
                containing_edge.head_contained_edges[contained_tree] = taxon_to_contained[containing_edge.head_node.taxon]
            else:
                containing_edge.head_contained_edges[contained_tree] = set()
                for nd in containing_edge.head_node.child_nodes():
                    containing_edge.head_contained_edges[contained_tree].update(nd.edge.tail_contained_edges[contained_tree])

            if containing_edge.tail_node is None:
                if containing_edge.length is not None:
                    target_age =  containing_edge.head_node.age + containing_edge.length
                else:
                    # assume all coalesce?
                    containing_edge.tail_contained_edges[contained_tree] = set([contained_tree.seed_node.edge])
                    continue
            else:
                target_age = containing_edge.tail_node.age

            containing_edge.tail_contained_edges[contained_tree] = set()
            for contained_edge in containing_edge.head_contained_edges[contained_tree]:
                if contained_edge.tail_node is not None:
                    remaining = target_age - contained_edge.tail_node.age
                elif contained_edge.length is not None:
                    remaining = target_age - (contained_edge.head_node.age + contained_edge.length)
                else:
                    continue
                while remaining > 0:
                    if contained_edge.tail_node is not None:
                        contained_edge = contained_edge.tail_node.edge
                    else:
                        if contained_edge.length is not None and (remaining - contained_edge.length) <= 0:
                            contained_edge = None
                            remaining = 0
                            break
                        else:
                            remaining = 0
                            break
                    if contained_edge and remaining > 0:
                        remaining -= contained_edge.length
                if contained_edge is not None:
                    containing_edge.tail_contained_edges[contained_tree].add(contained_edge)

    def build_edge_taxa_sets(self):
        """
        Rebuilds sets of containing and corresponding contained taxa at each
        edge.
        """
        for edge in self.postorder_edge_iter():
            if edge.is_terminal():
                edge.containing_taxa = set([edge.head_node.taxon])
            else:
                edge.containing_taxa = set()
                for i in edge.head_node.child_nodes():
                    edge.containing_taxa.update(i.edge.containing_taxa)
            edge.contained_taxa = set()
            for t in edge.containing_taxa:
                edge.contained_taxa.update(self.containing_to_contained_taxa_map[t])

    def num_deep_coalescences(self):
        """
        Returns total number of deep coalescences of the contained trees.
        """
        return sum(self.deep_coalescences().values())

    def deep_coalescences(self):
        """
        Returns dictionary where the contained trees are keys, and the number of
        deep coalescences corresponding to the tree are values.
        """
        dc = {}
        for tree in self.contained_trees:
            for edge in self.postorder_edge_iter():
                if edge.tail_node is None and self.ignore_root_deep_coalescences:
                    continue
                try:
                    dc[tree] += len(edge.tail_contained_edges[tree]) - 1
                except KeyError:
                    dc[tree] = len(edge.tail_contained_edges[tree]) - 1
        return dc

    def embed_contained_kingman(self,
            edge_pop_size_attr='pop_size',
            default_pop_size=1,
            label=None,
            rng=None,
            use_expected_tmrca=False):
        """
        Simulates, *embeds*, and returns a "censored" (Kingman) neutral coalescence tree
        conditional on self.

            ``rng``
                Random number generator to use. If |None|, the default will
                be used.

            ``edge_pop_size_attr``
                Name of attribute of self's edges that specify the population
                size. If this attribute does not exist, then the population
                size is taken to be 1.

        Note that all edge-associated taxon sets must be up-to-date (otherwise,
        ``build_edge_taxa_sets()`` should be called).
        """
        et = self.simulate_contained_kingman(
                edge_pop_size_attr=edge_pop_size_attr,
                default_pop_size=default_pop_size,
                label=label,
                rng=rng,
                use_expected_tmrca=use_expected_tmrca)
        self.embed_tree(et)
        return et

    def simulate_contained_kingman(self,
            edge_pop_size_attr='pop_size',
            default_pop_size=1,
            label=None,
            rng=None,
            use_expected_tmrca=False):
        """
        Simulates and returns a "censored" (Kingman) neutral coalescence tree
        conditional on self.

            ``rng``
                Random number generator to use. If |None|, the default will
                be used.

            ``edge_pop_size_attr``
                Name of attribute of self's edges that specify the population
                size. If this attribute does not exist, then the population
                size is taken to be 1.

        Note that all edge-associated taxon sets must be up-to-date (otherwise,
        ``build_edge_taxa_sets()`` should be called), and that the tree
        is *not* added to the set of contained trees. For the latter, call
        ``embed_contained_kingman``.
        """

        # Dictionary that maps nodes of containing tree to list of
        # corresponding nodes on gene tree, initially populated with leaf
        # nodes.
        contained_nodes = {}
        for nd in self.leaf_node_iter():
            contained_nodes[nd] = []
            for gt in nd.edge.contained_taxa:
                gn = dendropy.Node(taxon=gt)
                contained_nodes[nd].append(gn)

        # Generate the tree structure
        for edge in self.postorder_edge_iter():
            if edge.head_node.parent_node is None:
                # root: run unconstrained coalescence until just one gene node
                # remaining
                if hasattr(edge, edge_pop_size_attr):
                    pop_size = getattr(edge, edge_pop_size_attr)
                else:
                    pop_size = default_pop_size
                if len(contained_nodes[edge.head_node]) > 1:
                    final = coalescent.coalesce_nodes(nodes=contained_nodes[edge.head_node],
                            pop_size=pop_size,
                            period=None,
                            rng=rng,
                            use_expected_tmrca=use_expected_tmrca)
                else:
                    final = contained_nodes[edge.head_node]
            else:
                # run until next coalescence event, as determined by this edge
                # size.
                if hasattr(edge, edge_pop_size_attr):
                    pop_size = getattr(edge, edge_pop_size_attr)
                else:
                    pop_size = default_pop_size
                remaining = coalescent.coalesce_nodes(nodes=contained_nodes[edge.head_node],
                        pop_size=pop_size,
                        period=edge.length,
                        rng=rng,
                        use_expected_tmrca=use_expected_tmrca)
                try:
                    contained_nodes[edge.tail_node].extend(remaining)
                except KeyError:
                    contained_nodes[edge.tail_node] = remaining

        # Create and return the full tree
        contained_tree = dendropy.Tree(taxon_namespace=self.contained_taxon_namespace, label=label)
        contained_tree.seed_node = final[0]
        contained_tree.is_rooted = True
        return contained_tree

    def _find_youngest_intergroup_age(self, contained_tree, disjunct_leaf_set_list_split_bitmasks, starting_min_age=None):
        """
        Find the age of the youngest MRCA of disjunct leaf sets.
        """
        if starting_min_age is None:
            starting_min_age = float('inf')
        if contained_tree.seed_node.age is None:
            contained_tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        for nd in contained_tree.ageorder_node_iter(include_leaves=False):
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
        nw._write_taxa_block(out, self.taxon_namespace)
        out.write('\n')
        nw._write_taxa_block(out, self.contained_trees.taxon_namespace)
        if self.contained_trees.taxon_namespace.label:
            domain_title = self.contained_trees.taxon_namespace.label
        else:
            domain_title = self.contained_trees.taxon_namespace.oid
        contained_taxon_namespace = self.contained_trees.taxon_namespace
        contained_label = self.contained_trees.label
        out.write('\n')
        self._contained_to_containing_taxon_map.write_mesquite_association_block(out)
        out.write('\n')
        nw._write_trees_block(out, dendropy.TreeList([self], taxon_namespace=self.taxon_namespace))
        out.write('\n')
        nw._write_trees_block(out, dendropy.TreeList(self.contained_trees, taxon_namespace=contained_taxon_namespace, label=contained_label))
        out.write('\n')

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
    taxa_mask = species_tree.taxon_namespace.all_taxa_bitmask()
    species_node_gene_nodes = {}
    gene_node_species_nodes = {}
    for gnd in gene_tree.postorder_node_iter():
        gn_children = gnd.child_nodes()
        if len(gn_children) > 0:
            ssplit = 0
            for gn_child in gn_children:
                ssplit = ssplit | gene_node_species_nodes[gn_child].edge.leafset_bitmask
            sanc = species_tree.mrca(start_node=species_tree.seed_node, leafset_bitmask=ssplit)
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

def monophyletic_partition_discordance(tree, taxon_namespace_partition):
    """
    Returns the number of deep coalescences on tree ``tree`` that would result
    if the taxa in ``tax_sets`` formed K mutually-exclusive monophyletic groups,
    where K = len(tax_sets)
    ``taxon_namespace_partition`` == TaxonNamespacePartition
    """

    tax_sets = taxon_namespace_partition.subsets()

    # from dendropy.model import parsimony
    # taxon_state_sets_map = {}
    # assert tree.taxon_namespace is taxon_namespace_partition.taxon_namespace
    # for taxon in tree.taxon_namespace:
    #     taxon_state_sets_map[taxon] = [0 for i in range(len(tax_sets))]
    # for idx, ts in enumerate(tax_sets):
    #     for taxon in ts:
    #         taxon_state_sets_map[taxon][idx] = 1
    # for taxon in tree.taxon_namespace:
    #     taxon_state_sets_map[taxon] = [set([i]) for i in taxon_state_sets_map[taxon]]
    # return parsimony.fitch_down_pass(
    #         postorder_nodes=tree.postorder_node_iter(),
    #         taxon_state_sets_map=taxon_state_sets_map
    #         )

    dc_tree = dendropy.Tree()
    dc_tree.taxon_namespace = dendropy.TaxonNamespace()
    for t in range(len(tax_sets)):
        dc_tree.taxon_namespace.add_taxon(dendropy.Taxon(label=str(t)))
    def _get_dc_taxon(nd):
        for idx, tax_set in enumerate(tax_sets):
            if nd.taxon in tax_set:
                return dc_tree.taxon_namespace[idx]
        assert "taxon not found in partition: '%s'" % nd.taxon.label
    src_dc_map = {}
    for snd in tree.postorder_node_iter():
        nnd = dendropy.Node()
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
                    cnd = dendropy.Node()
                    cnd.taxon = t
                    nnd.add_child(cnd)
            else:
                if len(nnd.child_nodes()) == 0:
                    nnd.taxon = taxa_set[0]
                elif len(taxa_set) == 1:
                    cnd = dendropy.Node()
                    cnd.taxon = taxa_set[0]
                    nnd.add_child(cnd)
    dc_tree.seed_node = nnd
    return len(dc_tree.leaf_nodes()) - len(tax_sets)

