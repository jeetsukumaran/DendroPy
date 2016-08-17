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
Taxon-to-taxon phylogenetic distances.
"""

import math
import collections
import csv
from dendropy.calculate import statistics
from dendropy.utility import GLOBAL_RNG
from dendropy.utility import container
from dendropy.utility import error
import dendropy

class PhylogeneticDistanceMatrix(object):
    """
    Calculates and maintains patristic distance information of taxa on a tree.
    """

    @classmethod
    def from_tree(cls, tree):
        """
        Creates and returns a |PhylogeneticDistanceMatrix| based
        on the given tree.

        Note that this creates a "snapshot" of the current state of the tree.
        Subsequent changes to the tree will not be reflected in
        |PhylogeneticDistanceMatrix| instances previously created.

        Also note that syntactically you may prefer to use::

            pdm = tree.phylogenetic_distance_matrix()

        instead of::

            pdm = PhylogeneticDistanceMatrix.from_tree(tree)

        Parameters
        ----------
        tree : a |Tree| instance
            The |Tree| from which to get the phylogenetic distances.

        Returns
        -------
        pdm : A |PhylogeneticDistanceMatrix| instance

        Examples
        --------

        ::

            import dendropy
            tree = dendropy.Tree.get(path="tree.nex",
                    schema="nexus")
            pdm1 = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)

            # following is equivalent to above and probably preferred:
            pdm2 = tree.phylogenetic_distance_matrix()

        """
        pdm = cls()
        pdm.compile_from_tree(tree=tree)
        return pdm

    @classmethod
    def from_csv(cls,
            src,
            taxon_namespace=None,
            is_allow_new_taxa=None,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            default_data_type=float,
            label_transform_fn=None,
            **csv_reader_kwargs
            ):
        """
        Instantiates a new PhylogeneticDistanceMatrix instance with data
        from an external source.

        Parameters
        ----------
        src : file or file-like
            Source of data. This is a token delimited-file (e.g., a
            comma-delimited or tab-delimited file) providing a table which
            lists taxon labels in both rows and columns. The cells of the table
            are numeric (typically real) values that indicate the distance
            between the taxa of the current row and column. Note that *only*
            the upper right section of the table is considered. The diagonals
            values are typically zeroes and, in either case, ignored along with
            the lower diagonal. Despite being ignored by the
            PhylogeneticDistanceMatrix object, the values are parsed by the
            underlying reader and thus have to be valid numerical values.
        taxon_namespace : |TaxonNamespace| instance
            The taxon namespace with which to manage taxa. If this has
            not already been pre-populated with the taxon names, then
            ``is_allow_new_taxa`` should be set to |True|.
        is_allow_new_taxa : bool
            If |False|: we do *not* expect to encounter any new taxa in the
            data file, and it is an error if we do. If |True|: we do expect to
            encounter new taxa in the data file. The default value of this
            depends on the value passed to ``taxon_namespace``. If
            ``taxon_namespace`` is ``None`` or an empty |TaxonNamespace|
            instance, then unless explicitly set to |False|,
            ``is_allow_new_taxa`` will default to |True|: allowing of creation
            of new taxa corresponding to labels found in the data source. On
            the other hand, if ``taxon_namespace`` is not None and its value is
            a |TaxonNamespace| instance with at least one taxon, unless
            explicitly set to |True|, ``is_allow_new_taxa`` will default to
            |False|, and it will be an error if taxon labels are found in the
            data source that do not correspond (exactly) to |Taxon| objects
            defined in the taxon namespace.  This is to err on the side of
            caution, to avoid (or rather, highlight) problems due to incorrect
            or mismatching labels between the data source and the current taxon
            namespace.
        is_first_row_column_names : bool
            By default |True|: assumes that first row lists the taxon names.
            Set to |False| if there is no header row.
        is_first_column_row_names : bool
            By default |True|: assumes that first column lists the taxon names.
            Set to |False| if there is now row name column.
        label_transform_fn : function object
            If not None, this should be a function object that takes a string
            as an argument and returns another string. This function will be
            applied to row and column labels before they are matched to taxon
            labels in the |TaxonNamespace| instance given by
            ``taxon_namespace``.
        \*\*csv_reader_kwargs : keyword arguments
            This arguments will be passed to the underlying CSV reader.
            The most important one is probably 'delimiter'.

        Returns
        -------
        pdm : A |PhylogeneticDistanceMatrix| instance

        Examples
        --------

        ::

            import dendropy
            pdm1 = dendropy.PhylogeneticDistanceMatrix.from_csv(
                    src=open("data.csv"),
                    delimiter=",")
            pdm2 = dendropy.PhylogeneticDistanceMatrix.from_csv(
                    src=open("data.tsv"),
                    delimiter="\t")

        """
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        if len(taxon_namespace) == 0 and is_allow_new_taxa is None:
            is_allow_new_taxa = True
        old_taxon_namespace_mutability = taxon_namespace.is_mutable
        taxon_namespace.is_mutable = is_allow_new_taxa
        data_table = container.DataTable.from_csv(
                src,
                is_first_row_column_names=is_first_row_column_names,
                is_first_column_row_names=is_first_column_row_names,
                default_data_type=default_data_type,
                label_transform_fn=label_transform_fn,
                **csv_reader_kwargs
                )
        if not is_first_row_column_names and not is_first_column_row_names:
            distances = {}
            seen_row_labels = set()
            seen_row_taxa_labels = set()
            for i1, t1_label in enumerate(data_table.row_name_iter()):
                if len(taxon_namespace) <= i1:
                    t1 = taxon_namespace.require_taxon(label=t1_label)
                else:
                    t1 = taxon_namespace[i1]
                seen_row_labels.add(t1_label)
                assert t1.label not in seen_row_taxa_labels
                seen_row_taxa_labels.add(t1.label)
                distances[t1] = {}
                for i2, t2_label in enumerate(data_table.column_name_iter()):
                    if t2_label in seen_row_labels:
                        continue
                    if len(taxon_namespace) <= i2:
                        t2 = taxon_namespace.require_taxon(label=t2_label)
                    else:
                        t2 = taxon_namespace[i2]
                    distances[t1][t2] = data_table[t1_label, t2_label]
        else:
            distances = {}
            seen_taxa = set()
            taxa = []
            if is_first_column_row_names:
                name_iter = data_table.row_name_iter()
            else:
                name_iter = data_table.column_name_iter()
            for label in name_iter:
                t1 = taxon_namespace.require_taxon(label=label)
                assert t1 not in seen_taxa
                seen_taxa.add(t1)
                taxa.append(t1)
            seen_row_taxa = set()
            for i1, t1 in enumerate(taxa):
                assert t1 in seen_taxa
                seen_row_taxa.add(t1)
                distances[t1] = {}
                for i2, t2 in enumerate(taxa):
                    if t2 in seen_row_taxa:
                        continue
                    distances[t1][t2] = data_table[i1, i2]
        # else:
            # raise NotImplementedError()

        taxon_namespace.is_mutable = old_taxon_namespace_mutability
        pdm = cls()
        pdm.compile_from_dict(
                distances=distances,
                taxon_namespace=taxon_namespace)
        return pdm

    def __init__(self):
        self.clear()

    def clear(self):
        self.taxon_namespace = None
        self._mapped_taxa = set()
        self._all_distinct_mapped_taxa_pairs = set()
        self._tree_length = None
        self._num_edges = None
        self._taxon_phylogenetic_distances = {}
        self._taxon_phylogenetic_path_steps = {}
        self._mrca = {}

    def compile_from_tree(self, tree):
        """
        Calculates the distances. Note that the path length (in number of
        steps) between taxa that span the root will be off by one if
        the tree is unrooted.
        """
        self.clear()
        self.taxon_namespace = tree.taxon_namespace
        # for i1, t1 in enumerate(self.taxon_namespace):
        #     self._taxon_phylogenetic_distances[t1] = {}
        #     self._taxon_phylogenetic_path_steps[t1] = {}
        #     self._mrca[t1] = {}
        self._tree_length = 0.0
        self._num_edges = 0
        for node in tree.postorder_node_iter():
            try:
                self._tree_length += node.edge.length
            except TypeError: # None for edge length
                pass
            self._num_edges += 1
            children = node.child_nodes()
            if len(children) == 0:
                node.desc_paths = {node : (0,0)}
            else:
                node.desc_paths = {}
                for cidx1, c1 in enumerate(children):
                    for desc1, (desc1_plen, desc1_psteps) in c1.desc_paths.items():
                        if c1.edge_length is None:
                            c1_edge_length = 0.0
                        else:
                            c1_edge_length = c1.edge.length
                        node.desc_paths[desc1] = (desc1_plen + c1_edge_length, desc1_psteps + 1)
                        assert desc1.taxon is not None
                        if desc1.taxon not in self._taxon_phylogenetic_distances:
                            self._mapped_taxa.add(desc1.taxon)
                            self._taxon_phylogenetic_distances[desc1.taxon] = {}
                            self._taxon_phylogenetic_distances[desc1.taxon][desc1.taxon] = 0.0
                            self._taxon_phylogenetic_path_steps[desc1.taxon] = {}
                            self._taxon_phylogenetic_path_steps[desc1.taxon][desc1.taxon] = 0
                            self._mrca[desc1.taxon] = {desc1.taxon: desc1}
                        for c2 in children[cidx1+1:]:
                            for desc2, (desc2_plen, desc2_psteps) in c2.desc_paths.items():
                                self._mapped_taxa.add(desc2.taxon)
                                self._mrca[desc1.taxon][desc2.taxon] = c1.parent_node
                                # self._all_distinct_mapped_taxa_pairs.add( tuple([desc1.taxon, desc2.taxon]) )
                                self._all_distinct_mapped_taxa_pairs.add( frozenset([desc1.taxon, desc2.taxon]) )
                                if c2.edge_length is None:
                                    c2_edge_length = 0.0
                                else:
                                    c2_edge_length = c2.edge.length
                                pat_dist = node.desc_paths[desc1][0] + desc2_plen + c2_edge_length
                                self._taxon_phylogenetic_distances[desc1.taxon][desc2.taxon] = pat_dist
                                path_steps = node.desc_paths[desc1][1] + desc2_psteps + 1
                                self._taxon_phylogenetic_path_steps[desc1.taxon][desc2.taxon] = path_steps
                    del(c1.desc_paths)
        self._mirror_lookups()
        # assert self._tree_length == tree.length()

    def compile_from_dict(self, distances, taxon_namespace):
        self.clear()
        self.taxon_namespace = taxon_namespace
        for t1 in distances:
            self._mapped_taxa.add(t1)
            self._taxon_phylogenetic_distances[t1] = {}
            for t2 in distances[t1]:
                self._taxon_phylogenetic_distances[t1][t2] = distances[t1][t2]
        self._mirror_lookups()

    def _mirror_lookups(self):
        for ddata in (
                self._taxon_phylogenetic_distances,
                self._taxon_phylogenetic_path_steps,
                self._mrca,
                ):
            for taxon1 in ddata:
                for taxon2 in ddata[taxon1]:
                    # assert taxon1 is not taxon2
                    if taxon2 not in ddata:
                        ddata[taxon2] = {}
                    ddata[taxon2][taxon1] = ddata[taxon1][taxon2]

    def __eq__(self, o):
        if self.taxon_namespace is not o.taxon_namespace:
            return False
        return (True
                and (self._mapped_taxa == o._mapped_taxa)
                and (self._all_distinct_mapped_taxa_pairs == o._all_distinct_mapped_taxa_pairs)
                and (self._taxon_phylogenetic_distances == o._taxon_phylogenetic_distances)
                and (self._taxon_phylogenetic_path_steps == o._taxon_phylogenetic_path_steps)
                and (self._mrca == o._mrca)
                and (self._tree_length == o._tree_length)
                and (self._num_edges == o._num_edges)
                )

    def __hash__(self):
        return id(self)

    def __call__(self, taxon1, taxon2):
        return self.patristic_distance(taxon1, taxon2)

    def __copy__(self):
        return self.clone()

    def __iter__(self):
        for taxon in self._taxon_phylogenetic_distances:
            yield taxon

    def clone(self):
        o = self.__class__()
        o.taxon_namespace = self.taxon_namespace
        o._mapped_taxa = set(self._mapped_taxa)
        o._all_distinct_mapped_taxa_pairs = set(self._all_distinct_mapped_taxa_pairs)
        o._tree_length = self._tree_length
        o._num_edges = self._num_edges
        for src, dest in (
                (self._taxon_phylogenetic_distances, o._taxon_phylogenetic_distances,),
                (self._taxon_phylogenetic_path_steps, o._taxon_phylogenetic_path_steps,),
                (self._mrca, o._mrca,),
                ):
            for t1 in src:
                dest[t1] = {}
                for t2 in src[t1]:
                    dest[t1][t2] = src[t1][t2]
        return o

    def mrca(self, taxon1, taxon2):
        """
        Returns MRCA of two taxon objects.
        """
        return self._mrca[taxon1][taxon2]

    def distance(self,
            taxon1,
            taxon2,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Returns distance between taxon1 and taxon2.
        """
        if is_weighted_edge_distances:
            return self.patristic_distance(taxon1, taxon2, is_normalize_by_tree_size=is_normalize_by_tree_size)
        else:
            return self.path_edge_count(taxon1, taxon2, is_normalize_by_tree_size=is_normalize_by_tree_size)

    def patristic_distance(self, taxon1, taxon2, is_normalize_by_tree_size=False):
        """
        Returns patristic distance between two taxon objects.
        """
        if taxon1 is taxon2:
            return 0.0
        d = self._taxon_phylogenetic_distances[taxon1][taxon2]
        if is_normalize_by_tree_size:
            return d / self._tree_length
        else:
            return d

    def path_edge_count(self, taxon1, taxon2, is_normalize_by_tree_size=False):
        """
        Returns the number of edges between two taxon objects.
        """
        if taxon1 is taxon2:
            return 0
        d = self._taxon_phylogenetic_path_steps[taxon1][taxon2]
        if is_normalize_by_tree_size:
            return float(d) / self._num_edges
        else:
            return d

    def distances(self,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Returns list of patristic distances.
        """
        dmatrix, normalization_factor = self._get_distance_matrix_and_normalization_factor(
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,
                )
        results = []
        for t1, t2 in self._all_distinct_mapped_taxa_pairs:
            results.append(dmatrix[t1][t2]/normalization_factor)
        return results

    def max_pairwise_distance_taxa(self,
            is_weighted_edge_distances=True):
        if is_weighted_edge_distances:
            dists = self._taxon_phylogenetic_distances
        else:
            dists = self._taxon_phylogenetic_path_steps
        max_dist = None
        max_dist_taxa = None
        for t1, t2 in self._all_distinct_mapped_taxa_pairs:
            pat_dist = dists[t1][t2]
            if max_dist is None or pat_dist > max_dist:
                max_dist = pat_dist
                max_dist_taxa = (t1, t2)
        return max_dist_taxa

    def sum_of_distances(self,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Returns sum of patristic distances on tree.
        """
        return sum(self.distances(is_weighted_edge_distances=is_weighted_edge_distances,is_normalize_by_tree_size=is_normalize_by_tree_size))

    def taxon_iter(self, filter_fn=None):
        """
        Iterates over taxa in matrix. Note that this could be a subset of the taxa in
        the associated taxon namespace.
        """
        for t1 in self._mapped_taxa:
            if not filter_fn or filter_fn(t1):
                yield t1

    def distinct_taxon_pair_iter(self, filter_fn=None):
        """
        Iterates over all distinct pairs of taxa in matrix.
        """
        for t1, t2 in self._all_distinct_mapped_taxa_pairs:
            if not filter_fn or (filter_fn(t1) and filter_fn(t2)):
                yield t1, t2

    def mean_pairwise_distance(self,
            filter_fn=None,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Calculates the phylogenetic ecology statistic "MPD"[1,2] for the tree
        (only considering taxa for which ``filter_fn`` returns True when
        applied if ``filter_fn`` is specified).

        The mean pairwise distance (mpd) is given by:

            .. math::
                mpd = \\frac{ \\sum_{i}^{n} \\sum_{j}^{n} \\delta_{i,j} }{n \\choose 2},

        where :math:`i \\neq j`, :math:`\\delta_{i,j}` is the phylogenetic
        distance between species :math:`i` and :math:`j`, and :math:`n` is the number
        of species in the sample.

        Parameters
        ----------
        filter_fn : function object or None
            If |None|, then all leaves will be considered. Otherwise should
            be a function object that takes a Taxon instance as an argument and
            returns |True| if it is to be included in the calculation or
            |False| otherwise.
            In trees sampled from multiple communites, ``filter_fn`` can be
            used to restrict the calculation to only one community based on
            some criteria.
        is_weighted_edge_distances : bool
            If |True| then the edge-weighted distance, i.e., considering edge
            lengths, is returned. Otherwise the the path steps or the number of
            edges rather then the sum of is_weighted_edge_distances edges, connecting two
            taxa is considered.
        is_normalize_by_tree_size : bool
            If |True| then the results are normalized by the total tree length
            or steps/edges (depending on whether edge-weighted or unweighted
            distances are used, respectively). Otherwise, raw distances are
            used.

        Returns
        -------
        mpd : float
            The Mean Pairwise Distance (MPD) statistic for the daata.

        Examples
        --------

        ::

            import dendropy
            tree = dendropy.Tree.get(path="data.nex",
                    schema="nexus")
            pdm = dendropy.PhylogeneticDistanceMatrix(tree)

            # consider all tips
            mpd1 = pdm.mean_pairwise_distance()

            # only tips within the same community, based on the node annotation
            # "community"
            mpds_by_community = {}
            for community_label in ("1", "2", "3",):
                filter_fn = lambda x: x.annotations["community"] == community_label
                mpd = pdm.mean_pairwise_distance(filter_fn=filter_fn)
                mpds_by_community[community_label] = mpd

        References
        ----------

        [1] Webb, C.O. 2000. Exploring the phylogenetic structure of
        ecological communities: An example for rainforest trees. The
        American Naturalist 156: 145-155.

        [2] Swenson, N.G. Functional and Phylogenetic Ecology in R.


        """
        comparison_regime = self.distinct_taxon_pair_iter(filter_fn=filter_fn)
        return self._calculate_mean_pairwise_distance(
                comparison_regime=comparison_regime,
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,
                )

    def mean_nearest_taxon_distance(self,
            filter_fn=None,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Calculates the phylogenetic ecology statistic "MNTD"[1,2] for the tree
        (only considering taxa for which ``filter_fn`` returns True when
        applied if ``filter_fn`` is specified).

        The mean nearest taxon distance (mntd) is given by:

            .. math::
                mntd = \\frac{ \\sum_{i}^{n} min(\\delta_{i,j}) }{n},

        where :math:`i \\neq j`, :math:`\\delta_{i,j}` is the phylogenetic
        distance between species :math:`i` and :math:`j`, and :math:`n` is the number
        of species in the sample.

        Parameters
        ----------
        filter_fn : function object or None
            If |None|, then all leaves will be considered. Otherwise should
            be a function object that takes a Taxon instance as an argument and
            returns |True| if it is to be included in the calculation or
            |False| otherwise.
            In trees sampled from multiple communites, ``filter_fn`` can be
            used to restrict the calculation to only one community based on
            some criteria.
        is_weighted_edge_distances : bool
            If |True| then the edge-weighted distance, i.e., considering edge
            lengths, is returned. Otherwise the the path steps or the number of
            edges rather then the sum of is_weighted_edge_distances edges, connecting two
            taxa is considered.
        is_normalize_by_tree_size : bool
            If |True| then the results are normalized by the total tree length
            or steps/edges (depending on whether edge-weighted or unweighted
            distances are used, respectively). Otherwise, raw distances are
            used.

        Returns
        -------
        mntd : float
            The Mean Nearest Taxon Distance (MNTD) statistic for the daata.

        Examples
        --------

        ::

            import dendropy
            tree = dendropy.Tree.get(path="data.nex",
                    schema="nexus")
            pdm = dendropy.PhylogeneticDistanceMatrix(tree)

            # consider all tips
            mntd = pdm.mean_nearest_taxon_distance()

            # only tips within the same community, based on the node annotation
            # "community"
            mntds_by_community = {}
            for community_label in ("1", "2", "3",):
                filter_fn = lambda x: x.annotations["community"] == community_label
                mntd = pdm.mean_pairwise_distance(filter_fn=filter_fn)
                mntds_by_community[community_label] = mntd

        References
        ----------

        [1] Webb, C.O. 2000. Exploring the phylogenetic structure of
        ecological communities: An example for rainforest trees. The
        American Naturalist 156: 145-155.

        [2] Swenson, N.G. Functional and Phylogenetic Ecology in R.

        """
        comparison_regime = self._get_taxon_to_all_other_taxa_comparisons(filter_fn=filter_fn)
        return self._calculate_mean_nearest_taxon_distance(
                comparison_regime=comparison_regime,
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,)

    def standardized_effect_size_mean_pairwise_distance(self,
            assemblage_memberships,
            num_randomization_replicates=1000,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False,
            is_skip_single_taxon_assemblages=False,
            null_model_type="taxa.label",
            rng=None):
        """
        Returns the standardized effect size value for the MPD statistic under
        a null model under various community compositions.

        The S.E.S. is given by:

            .. math::
                SES(statistic) = \\frac{observed - mean(model_{null})}{sd(model_{null})}

        This removes any bias associated with the decrease in variance in the
        MPD statistic value as species richness increases to the point where
        communities become saturated. Equivalent to -1 times the Nearest
        Relative Index (NRI) when using phylogenetic distances.

        In contrast to the function calculating the non-standardized effect
        size version of this statistic, which uses filter function to specify
        the subset of taxa to be considerd, here a collection of (multiple)
        sets of taxa constituting a community is specified. This to allow
        calculation of the null model statistic across all community sets for
        each randomization replicate.

        Parameters
        ----------
        assemblage_memberships : iterable of iterable of |Taxon| objects
            A collection of collections, e.g. a list of sets, with the elements
            of each set being |Taxon| instances. Each set specifies the
            composition of a community. The standardized effect size of this
            statistic will be calculated for each community as specified by a
            set of |Taxon| instances.
        num_randomization_replicates : int
            Number of randomization replicates.
        is_weighted_edge_distances: bool
            If ``True`` then edge lengths will be considered for distances.
            Otherwise, just the number of edges.

        Returns
        -------
        r : list of results
            A list of results, with each result corresponding to a community
            set given in ``assemblage_memberships``. Each result consists of a named
            tuple with the following elements:

                -   obs       : the observed value of the statistic
                -   null_model_mean : the mean value of the statistic under the null
                                model
                -   null_model_sd   : the standard deviation of the statistic under
                                the null model
                -   z         : the standardized effect value of the statistic
                                (= SES as defined in [1] above)
                -   p         : the p-value of the observed value of the

        Examples
        --------

        ::

            import dendropy
            tree = dendropy.Tree.get_from_path(
                    src="data/community.tree.newick",
                    schema="newick",
                    rooting="force-rooted")
            pdm = tree.phylogenetic_distance_matrix()
            assemblage_membership_definitions = pdm.assemblage_membership_definitions_from_csv("data/comm1.csv")
            results = pdm.standardized_effect_size_mean_pairwise_distance(assemblage_memberships=assemblage_membership_definitions.values())
            print(results)

        """
        if assemblage_memberships is None:
            assemblage_memberships = [ set(self._mapped_taxa) ]
        comparison_regimes = []
        for idx, assemblage_membership in enumerate(assemblage_memberships):
            if len(assemblage_membership) == 1:
                if is_skip_single_taxon_assemblages:
                    continue
                else:
                    raise error.SingleTaxonAssemblageException("{}: {}".format(idx, assemblage_membership))
            filter_fn = lambda taxon: taxon in assemblage_membership
            comparison_regime = list(self.distinct_taxon_pair_iter(filter_fn=filter_fn))
            comparison_regimes.append(comparison_regime)
        results = self._calculate_standardized_effect_size(
                statisticf_name="_calculate_mean_pairwise_distance",
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,
                comparison_regimes=comparison_regimes,
                null_model_type=null_model_type,
                num_randomization_replicates=num_randomization_replicates,
                rng=rng)
        return results

    def standardized_effect_size_mean_nearest_taxon_distance(self,
            assemblage_memberships,
            num_randomization_replicates=1000,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False,
            is_skip_single_taxon_assemblages=False,
            null_model_type="taxa.label",
            rng=None):
        """
        Returns the standardized effect size value for the MNTD statistic under
        a null model under various community compositions.

        The S.E.S. is given by:

            .. math::
                SES(statistic) = \\frac{observed - mean(model_{null})}{sd(model_{null})}

        This removes any bias associated with the decrease in variance in the
        MPD statistic value as species richness increases to the point where
        communities become saturated. Equivalent to -1 times the Nearest Taxon
        Index when using phylogenetic distances.

        In contrast to the function calculating the non-standardized effect
        size version of this statistic, which uses filter function to specify
        the subset of taxa to be considerd, here a collection of (multiple)
        sets of taxa constituting a community is specified. This to allow
        calculation of the null model statistic across all community sets for
        each randomization replicate.

        Parameters
        ----------
        assemblage_memberships : iterable of iterable of |Taxon| objects
            A collection of collections, e.g. a list of sets, with the elements
            of each set being |Taxon| instances. Each set specifies the
            composition of a community. The standardized effect size of this
            statistic will be calculated for each community as specified by a
            set of |Taxon| instances.
        num_randomization_replicates : int
            Number of randomization replicates.
        is_weighted_edge_distances: bool
            If ``True`` then edge lengths will be considered for distances.
            Otherwise, just the number of edges.

        Returns
        -------
        r : list of results
            A list of results, with each result corresponding to a community
            set given in ``assemblage_memberships``. Each result consists of a named
            tuple with the following elements:

                -   obs       : the observed value of the statistic
                -   null_model_mean : the mean value of the statistic under the null
                                model
                -   null_model_sd   : the standard deviation of the statistic under
                                the null model
                -   z         : the standardized effect value of the statistic
                                (= SES as defined in [1] above)
                -   p         : the p-value of the observed value of the
                                statistic under the null model.
        Examples
        --------

        ::

            import dendropy
            tree = dendropy.Tree.get_from_path(
                    src="data/community.tree.newick",
                    schema="newick",
                    rooting="force-rooted")
            pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)
            assemblage_memberships = pdm.assemblage_membership_definitions_from_csv("data/comm1.csv")
            results = pdm.standardized_effect_size_mean_nearest_taxon_distance(assemblage_memberships=assemblage_memberships)
            print(results)

        """
        if assemblage_memberships is None:
            assemblage_memberships = [ set(self._mapped_taxa) ]
        comparison_regimes = []
        for idx, assemblage_membership in enumerate(assemblage_memberships):
            if len(assemblage_membership) == 1:
                if is_skip_single_taxon_assemblages:
                    continue
                else:
                    raise error.SingleTaxonAssemblageException("{}: {}".format(idx, assemblage_membership))
            filter_fn = lambda taxon: taxon in assemblage_membership
            comparison_regime = self._get_taxon_to_all_other_taxa_comparisons(filter_fn=filter_fn)
            comparison_regimes.append(comparison_regime)
        results = self._calculate_standardized_effect_size(
                statisticf_name="_calculate_mean_nearest_taxon_distance",
                comparison_regimes=comparison_regimes,
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,
                null_model_type=null_model_type,
                num_randomization_replicates=num_randomization_replicates,
                rng=rng)
        return results

    def shuffle_taxa(self,
            is_shuffle_phylogenetic_distances=True,
            is_shuffle_phylogenetic_path_steps=True,
            is_shuffle_mrca=True,
            rng=None):
        """
        Randomly shuffles taxa in-situ.
        """
        if rng is None:
            rng = GLOBAL_RNG
        reordered_taxa = list(self._mapped_taxa)
        rng.shuffle(reordered_taxa)
        current_to_shuffled_taxon_map = dict(zip(self._mapped_taxa, reordered_taxa))
        to_shuffle = []
        if is_shuffle_phylogenetic_distances:
            to_shuffle.append("_taxon_phylogenetic_distances")
        if is_shuffle_phylogenetic_path_steps:
            to_shuffle.append("_taxon_phylogenetic_path_steps")
        if is_shuffle_mrca:
            to_shuffle.append("_mrca")
        for attr_name in to_shuffle:
            src = getattr(self, attr_name)
            dest = {}

            ## 5m8.076s
            # for t1, t2 in self._all_distinct_mapped_taxa_pairs:
            #     x1 = current_to_shuffled_taxon_map[t1]
            #     x2 = current_to_shuffled_taxon_map[t2]
            #     d = src[t1][t2]
            #     try:
            #         dest[x1][x2] = d
            #     except KeyError:
            #         dest[x1] = {x2: d}
            #         if t1 in src[t1]:
            #             dest[x1][x1] = src[t1][t1]
            #     try:
            #         dest[x2][x1] = d
            #     except KeyError:
            #         dest[x2] = {x1: d}
            #         if t2 in src[t2]:
            #             dest[x2][x2] = src[t2][t2]
            # setattr(self, attr_name, dest)

            # 4m48.025s
            for t1 in src:
                x1 = current_to_shuffled_taxon_map[t1]
                dest[x1] = {}
                for t2 in src[t1]:
                    x2 = current_to_shuffled_taxon_map[t2]
                    dest[x1][x2] = src[t1][t2]
            setattr(self, attr_name, dest)

        return current_to_shuffled_taxon_map

    def nj_tree(self,
            is_weighted_edge_distances=True,
            tree_factory=None,
            ):
        """
        Returns an Neighbor-Joining (NJ) tree based on the distances in the matrix.

        Calculates and returns a tree under the Neighbor-Joining algorithm of
        Saitou and Nei (1987) for the data in the matrix.

        Parameters
        ----------
        is_weighted_edge_distances: bool
            If ``True`` then edge lengths will be considered for distances.
            Otherwise, just the number of edges.

        Returns
        -------
        t : |Tree|
            A |Tree| instance corresponding to the Neighbor-Joining (NJ) tree
            for this data.

        Examples
        --------

        ::

            import dendropy

            # Read data from a CSV file into a PhylogeneticDistanceMatrix
            # object
            with open("distance_matrix.csv") as src:
                pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                        src,
                        is_first_row_column_names=True,
                        is_first_column_row_names=True,
                        is_allow_new_taxa=True,
                        delimiter=",",
                        )

            # Calculate the tree
            nj_tree = pdm.nj_tree()

            # Print it
            print(nj_tree.as_string("nexus"))


        References
        ----------
        Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method
        for reconstructing phylogenetic trees. Molecular Biology and Evolution,
        4: 406-425.

        """

        if is_weighted_edge_distances:
            original_dmatrix = self._taxon_phylogenetic_distances
        else:
            original_dmatrix = self._taxon_phylogenetic_path_steps
        if tree_factory is None:
            tree_factory = dendropy.Tree
        tree = tree_factory(taxon_namespace=self.taxon_namespace)
        tree.is_rooted = False

        # initialize node pool
        node_pool = []
        for t1 in self._mapped_taxa:
            nd = tree.node_factory()
            nd.taxon = t1
            nd._nj_distances = {}
            node_pool.append(nd)

        # initialize factor
        n = len(self._mapped_taxa)

        # cache calculations
        for nd1 in node_pool:
            nd1._nj_xsub = 0.0
            for nd2 in node_pool:
                if nd1 is nd2:
                    continue
                d = original_dmatrix[nd1.taxon][nd2.taxon]
                nd1._nj_distances[nd2] = d
                nd1._nj_xsub += d

        while n > 1:

            # calculate the Q-matrix
            min_q = None
            nodes_to_join = None
            for idx1, nd1 in enumerate(node_pool[:-1]):
                for idx2, nd2 in enumerate(node_pool[idx1+1:]):
                    v1 = (n - 2) * nd1._nj_distances[nd2]
                    qvalue = v1 - nd1._nj_xsub - nd2._nj_xsub
                    if min_q is None or qvalue < min_q:
                        min_q = qvalue
                        nodes_to_join = (nd1, nd2)

            # create the new node
            new_node = tree.node_factory()

            # attach it to the tree
            for node_to_join in nodes_to_join:
                new_node.add_child(node_to_join)
                node_pool.remove(node_to_join)

            # calculate the distances for the new node
            new_node._nj_distances = {}
            new_node._nj_xsub = 0.0
            for node in node_pool:
                # actual node-to-node distances
                v1 = 0.0
                for node_to_join in nodes_to_join:
                    v1 += node._nj_distances[node_to_join]
                v3 = nodes_to_join[0]._nj_distances[nodes_to_join[1]]
                dist = 0.5 * (v1 - v3)
                new_node._nj_distances[node] = dist
                node._nj_distances[new_node] = dist

                # Adjust/recalculate the values needed for the Q-matrix
                # calculations
                new_node._nj_xsub += dist
                node._nj_xsub += dist
                for node_to_join in nodes_to_join:
                    node._nj_xsub -= node_to_join._nj_distances[node]

            # calculate the branch lengths
            if n > 2:
                v1 = 0.5 * nodes_to_join[0]._nj_distances[nodes_to_join[1]]
                v4  = 1.0/(2*(n-2)) * (nodes_to_join[0]._nj_xsub - nodes_to_join[1]._nj_xsub)
                delta_f = v1 + v4
                delta_g = nodes_to_join[0]._nj_distances[nodes_to_join[1]] - delta_f
                nodes_to_join[0].edge.length = delta_f
                nodes_to_join[1].edge.length = delta_g
            else:
                d = nodes_to_join[0]._nj_distances[nodes_to_join[1]]
                nodes_to_join[0].edge.length = d / 2
                nodes_to_join[1].edge.length = d / 2

            # clean up
            for node_to_join in nodes_to_join:
                del node_to_join._nj_distances
                del node_to_join._nj_xsub

            # add the new node to the pool of nodes
            node_pool.append(new_node)

            # adjust count
            n -= 1

        tree.seed_node = node_pool[0]
        del tree.seed_node._nj_distances
        del tree.seed_node._nj_xsub
        return tree

    def upgma_tree(self,
            is_weighted_edge_distances=True,
            tree_factory=None,
            ):
        """
        Returns an Unweighted Pair Group Method with Arithmetic Mean (UPGMA) tree
        based on the distances in the matrix.

        Parameters
        ----------
        is_weighted_edge_distances: bool
            If ``True`` then edge lengths will be considered for distances.
            Otherwise, just the number of edges.

        Returns
        -------
        t : |Tree|
            A |Tree| instance corresponding to the UPGMA tree
            for this data.

        Examples
        --------

        ::

            import dendropy

            # Read data from a CSV file into a PhylogeneticDistanceMatrix
            # object
            with open("distance_matrix.csv") as src:
                pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                        src,
                        is_first_row_column_names=True,
                        is_first_column_row_names=True,
                        is_allow_new_taxa=True,
                        delimiter=",",
                        )

            # Calculate the tree
            upgma_tree = pdm.upgma_tree()

            # Print it
            print(upgma_tree.as_string("nexus"))

        """

        if is_weighted_edge_distances:
            original_dmatrix = self._taxon_phylogenetic_distances
        else:
            original_dmatrix = self._taxon_phylogenetic_path_steps
        if tree_factory is None:
            tree_factory = dendropy.Tree
        tree = tree_factory(taxon_namespace=self.taxon_namespace)
        tree.is_rooted = True
        node_pool = []
        for t1 in self._mapped_taxa:
            nd = tree.node_factory()
            nd.taxon = t1
            nd._upgma_cluster = set([nd])
            nd._upgma_distance_from_tip = 0.0
            nd._upgma_distances = {}
            node_pool.append(nd)
        for idx1, nd1 in enumerate(node_pool[:-1]):
            for idx2, nd2 in enumerate(node_pool[idx1+1:]):
                d = original_dmatrix[nd1.taxon][nd2.taxon]
                nd1._upgma_distances[nd2] = d
                nd2._upgma_distances[nd1] = d
        while len(node_pool) > 1:
            min_distance = None
            nodes_to_join = None
            for idx1, nd1 in enumerate(node_pool[:-1]):
                for idx2, nd2 in enumerate(node_pool[idx1+1:]):
                    d = nd1._upgma_distances[nd2]
                    if min_distance is None or d < min_distance:
                        nodes_to_join = (nd1, nd2)
                        min_distance = d
            new_node = tree.node_factory()
            new_node._upgma_cluster = set()
            new_node._upgma_distances = {}
            elen = min_distance / 2.0
            for node_to_join in nodes_to_join:
                new_node.add_child(node_to_join)
                new_node._upgma_cluster.update(node_to_join._upgma_cluster)
                node_to_join.edge.length = elen - node_to_join._upgma_distance_from_tip
                node_pool.remove(node_to_join)
            new_node._upgma_distance_from_tip = nodes_to_join[0].edge.length + nodes_to_join[0]._upgma_distance_from_tip
            for idx1, nd1 in enumerate(node_pool):
                d1 = 0.0
                count = 0.0
                for node_to_join in nodes_to_join:
                    d2 = node_to_join._upgma_distances[nd1]
                    xc = len(node_to_join._upgma_cluster)
                    d1 += (d2 * xc)
                    count += xc
                d = d1 / count
                nd1._upgma_distances[new_node] = d
                new_node._upgma_distances[nd1] = d
            for node_to_join in nodes_to_join:
                del node_to_join._upgma_cluster
                del node_to_join._upgma_distance_from_tip
                del node_to_join._upgma_distances
            node_pool.append(new_node)
        tree.seed_node = node_pool[0]
        del tree.seed_node._upgma_cluster
        del tree.seed_node._upgma_distance_from_tip
        del tree.seed_node._upgma_distances
        return tree

    def as_data_table(self, is_weighted_edge_distances=True):
        """
        Returns this as a table.
        """
        if is_weighted_edge_distances:
            df = self.patristic_distance
        else:
            df = self.path_edge_count
        dt = container.DataTable()
        for t1 in self._mapped_taxa:
            dt.add_row(row_name=t1.label)
            dt.add_column(column_name=t1.label)
        for t1 in self._mapped_taxa:
            for t2 in self._mapped_taxa:
                dt[t1.label, t2.label] = df(t1, t2)
        return dt

    def write_csv(self,
            out,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=True,
            label_transform_fn=None,
            **csv_writer_kwargs
            ):
        if isinstance(out, str):
            dest = open(out, "w")
        else:
            dest = out
        if label_transform_fn is None:
            label_transform_fn = lambda x: x
        dmatrix, normalization_factor = self._get_distance_matrix_and_normalization_factor(
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,
                )
        if "delimiter" not in csv_writer_kwargs:
            csv_writer_kwargs["delimiter"] = ","
        writer = csv.writer(dest, csv_writer_kwargs)
        if is_first_row_column_names:
            row = []
            if is_first_column_row_names:
                row.append("")
            for taxon in self._mapped_taxa:
                row.append(label_transform_fn(taxon.label))
            writer.writerow(row)
            # dest.write(delimiter.join(row))
            # dest.write("\n")
        for taxon1 in self._mapped_taxa:
            row = []
            if is_first_column_row_names:
                row.append(label_transform_fn(taxon1.label))
            for taxon2 in self._mapped_taxa:
                d = dmatrix[taxon1][taxon2] / normalization_factor
                row.append("{}".format(d))
            writer.writerow(row)
            # dest.write(delimiter.join(row))
            # dest.write("\n")

    def assemblage_membership_definitions_from_csv(
            self,
            src,
            default_data_type=float,
            **csv_reader_kwargs):
        """
        Convenience method to return list of community sets from a delimited
        file that lists taxon (labels) in columns and community
        presence/absences or abundances in rows.
        """
        if isinstance(src, str):
            with open(src) as srcf:
                data_table = container.DataTable.from_csv(
                        src,
                        default_data_type=default_data_type,
                        **csv_reader_kwargs
                        )
        else:
            data_table = container.DataTable.from_csv(
                    src,
                    default_data_type=default_data_type,
                    **csv_reader_kwargs
                    )
        mapped_taxon_labels = set([taxon.label for taxon in self.taxon_iter()])
        for column_name in data_table.column_name_iter():
            assert column_name in mapped_taxon_labels
        assemblage_memberships = collections.OrderedDict()
        for row_name in data_table.row_name_iter():
            assemblage_membership = set()
            for taxon in self.taxon_iter():
                if data_table[row_name, taxon.label] > 0:
                    assemblage_membership.add(taxon)
            assemblage_memberships[row_name] = assemblage_membership
        return assemblage_memberships

    def _get_taxon_to_all_other_taxa_comparisons(self, filter_fn=None):
        permutations = collections.defaultdict(list)
        for taxon1 in self._mapped_taxa:
            # permutations[taxon1] = []
            if filter_fn and not filter_fn(taxon1):
                continue
            for taxon2 in self._mapped_taxa:
                if taxon1 is taxon2:
                    continue
                if filter_fn and not filter_fn(taxon2):
                    continue
                permutations[taxon1].append(taxon2)
        return permutations

    def _get_distance_matrix_and_normalization_factor(self,
            is_weighted_edge_distances,
            is_normalize_by_tree_size):
        if is_weighted_edge_distances:
            dmatrix = self._taxon_phylogenetic_distances
            if is_normalize_by_tree_size:
                normalization_factor = self._tree_length
            else:
                normalization_factor = 1.0
        else:
            dmatrix = self._taxon_phylogenetic_path_steps
            if is_normalize_by_tree_size:
                normalization_factor = float(self._num_edges)
            else:
                normalization_factor = 1.0
        return dmatrix, normalization_factor

    def _calculate_mean_pairwise_distance(self,
            comparison_regime,
            is_weighted_edge_distances,
            is_normalize_by_tree_size):
        dmatrix, normalization_factor = self._get_distance_matrix_and_normalization_factor(
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,)
        distances = []
        for taxon1, taxon2 in comparison_regime:
            distances.append(dmatrix[taxon1][taxon2])
        if distances:
            return (sum(distances) / normalization_factor) / (len(distances) * 1.0)
        else:
            raise error.NullAssemblageException("No taxa in assemblage")

    def _calculate_mean_nearest_taxon_distance(self,
            comparison_regime,
            is_weighted_edge_distances,
            is_normalize_by_tree_size):
        dmatrix, normalization_factor = self._get_distance_matrix_and_normalization_factor(
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,)
        distances = []
        for taxon1 in comparison_regime:
            # subdistances = [dmatrix[taxon1][taxon2] for taxon2 in comparison_regime[taxon1]]
            # distances.append(min(subdistances))
            min_distance = dmatrix[taxon1][comparison_regime[taxon1][0]]
            for taxon2 in comparison_regime[taxon1][1:]:
                d = dmatrix[taxon1][taxon2]
                if d < min_distance:
                    min_distance = d
            distances.append(min_distance)
        if distances:
            return (sum(distances) / normalization_factor) / (len(distances) * 1.0)
        else:
            raise error.NullAssemblageException("No taxa in assemblage")

    def _calculate_standardized_effect_size(self,
            statisticf_name,
            comparison_regimes,
            is_weighted_edge_distances,
            is_normalize_by_tree_size,
            null_model_type="taxa.label",
            num_randomization_replicates=1000,
            rng=None):
        result_type = collections.namedtuple("PhylogeneticCommunityStandardizedEffectSizeStatisticCalculationResult",
                ["obs", "null_model_mean", "null_model_sd", "z", "rank", "p",])
        statisticf_kwargs={
            "is_weighted_edge_distances": is_weighted_edge_distances,
            "is_normalize_by_tree_size": is_normalize_by_tree_size
        }
        observed_stat_values = {}
        null_model_stat_values = {}
        null_model_matrix = self.clone()
        assert null_model_matrix == self
        if is_weighted_edge_distances:
            is_shuffle_phylogenetic_distances = True
            is_shuffle_phylogenetic_path_steps = False
        else:
            is_shuffle_phylogenetic_distances = False
            is_shuffle_phylogenetic_path_steps = True
        for rep_idx in range(num_randomization_replicates):
            null_model_matrix.shuffle_taxa(
                    is_shuffle_phylogenetic_distances=is_shuffle_phylogenetic_distances,
                    is_shuffle_phylogenetic_path_steps=is_shuffle_phylogenetic_distances,
                    is_shuffle_mrca=False,
                    rng=rng)
            for comparison_regime_idx, comparison_regime in enumerate(comparison_regimes):
                statisticf_kwargs["comparison_regime"] = comparison_regime
                if rep_idx == 0:
                    observed_stat_values[comparison_regime_idx] = getattr(self, statisticf_name)(**statisticf_kwargs)
                    null_model_stat_values[comparison_regime_idx] = []
                stat_value = getattr(null_model_matrix, statisticf_name)(**statisticf_kwargs)
                null_model_stat_values[comparison_regime_idx].append(stat_value)
        results = []
        for comparison_regime_idx, comparison_regime in enumerate(comparison_regimes):
            obs_value = observed_stat_values[comparison_regime_idx]
            stat_values = null_model_stat_values[comparison_regime_idx]
            null_model_mean, null_model_var = statistics.mean_and_sample_variance(stat_values)
            rank = statistics.rank(
                    value_to_be_ranked=obs_value,
                    value_providing_rank=stat_values)
            if null_model_var > 0:
                null_model_sd = math.sqrt(null_model_var)
                z = (obs_value - null_model_mean) / null_model_sd
            else:
                null_model_sd = 0.0
                z = None
            p = float(rank) / len(stat_values)
            result = result_type(
                    obs=obs_value,
                    null_model_mean=null_model_mean,
                    null_model_sd=null_model_sd,
                    z=z,
                    rank=rank,
                    p=p)
            results.append(result)
        return results

class NodeDistanceMatrix(object):

    @classmethod
    def from_tree(cls, tree):
        ndm = cls()
        ndm.compile_from_tree(tree=tree)
        return ndm

    def __init__(self):
        self.clear()

    def clear(self):
        self._tree_length = None
        self._num_edges = None
        self._node_phylogenetic_distances = {}
        self._node_phylogenetic_path_steps = {}
        self._mrca = {}

    def compile_from_tree(self, tree):
        self.clear()
        self._tree_length = 0.0
        self._num_edges = 0
        for node1 in tree.postorder_node_iter():
            try:
                self._tree_length += node1.edge.length
            except TypeError: # None for edge length
                pass
            self._num_edges += 1
            if node1 not in self._node_phylogenetic_distances:
                self._node_phylogenetic_distances[node1] = {node1: 0.0}
                self._node_phylogenetic_path_steps[node1] = {node1: 0}
                self._mrca[node1] = {node1: node1}
            children = node1.child_nodes()
            for ch_idx, ch1 in enumerate(children):
                ch1_elen = ch1.edge.length if ch1.edge.length is not None else 0.0
                for ch1_subtree_node in list(self._node_phylogenetic_distances[ch1].keys()):
                    if ch1_subtree_node not in self._node_phylogenetic_distances[node1]:
                        d = self._node_phylogenetic_distances[ch1][ch1_subtree_node] + ch1_elen
                        d2 = self._node_phylogenetic_path_steps[ch1][ch1_subtree_node] + 1
                        self._node_phylogenetic_distances[node1][ch1_subtree_node] = d
                        self._node_phylogenetic_distances[ch1_subtree_node][node1] = d
                        self._node_phylogenetic_path_steps[node1][ch1_subtree_node] = d2
                        self._node_phylogenetic_path_steps[ch1_subtree_node][node1] = d2
                self._node_phylogenetic_distances[node1][ch1] = ch1_elen
                self._node_phylogenetic_distances[ch1][node1] = ch1_elen
                self._node_phylogenetic_path_steps[node1][ch1] = 1
                self._node_phylogenetic_path_steps[ch1][node1] = 1
                for ch2 in children[ch_idx+1:]:
                    self._mrca[ch1][ch2] = node1
                    self._mrca[ch2][ch1] = node1
                    ch2_elen = ch2.edge.length if ch2.edge.length is not None else 0.0
                    d = ch1_elen + ch2_elen
                    self._node_phylogenetic_distances[ch1][ch2] = d
                    self._node_phylogenetic_distances[ch2][ch1] = d
                    self._node_phylogenetic_path_steps[ch1][ch2] = 2
                    self._node_phylogenetic_path_steps[ch2][ch1] = 2
            # Below is ugly, ugly, ugly. Basic idea is to link nodes of each
            # the subtrees of each of the child nodes of node1. Assumes
            # that any pairwise comparison of nodes descending from node1
            # (as given by nodes in a pairwise comparison with node1) not
            # already made have their MRCA at node1.
            for snd1 in self._node_phylogenetic_distances[node1]:
                for snd2 in self._node_phylogenetic_distances[node1]:
                    if snd1 is snd2:
                        continue
                    if snd1 not in self._node_phylogenetic_distances:
                        self._node_phylogenetic_distances[snd1] = {}
                        self._node_phylogenetic_path_steps[snd1] = {}
                    if snd2 not in self._node_phylogenetic_distances:
                        self._node_phylogenetic_distances[snd2] = {}
                        self._node_phylogenetic_path_steps[snd2] = {}
                    if snd2 not in self._node_phylogenetic_distances[snd1]:
                        self._node_phylogenetic_distances[snd1][snd2] = self._node_phylogenetic_distances[node1][snd1] + self._node_phylogenetic_distances[node1][snd2]
                        self._node_phylogenetic_path_steps[snd1][snd2] = self._node_phylogenetic_path_steps[node1][snd1] + self._node_phylogenetic_path_steps[node1][snd2]
                    if snd1 not in self._node_phylogenetic_distances[snd2]:
                        self._node_phylogenetic_distances[snd2][snd1] = self._node_phylogenetic_distances[node1][snd1] + self._node_phylogenetic_distances[node1][snd2]
                        self._node_phylogenetic_path_steps[snd2][snd1] = self._node_phylogenetic_path_steps[node1][snd1] + self._node_phylogenetic_path_steps[node1][snd2]
                    if snd1 not in self._mrca:
                        self._mrca[snd1] = {}
                    if snd2 not in self._mrca:
                        self._mrca[snd2] = {}
                    if snd2 not in self._mrca[snd1]:
                        self._mrca[snd1][snd2] = node1
                        self._mrca[snd2][snd1] = node1

    def __eq__(self, o):
        if self.node_namespace is not o.node_namespace:
            return False
        return (True
                and (self._node_phylogenetic_distances == o._node_phylogenetic_distances)
                and (self._node_phylogenetic_path_steps == o._node_phylogenetic_path_steps)
                and (self._mrca == o._mrca)
                and (self._tree_length == o._tree_length)
                and (self._num_edges == o._num_edges)
                )

    def __iter__(self):
        for node in self._node_phylogenetic_distances:
            yield node

    def __hash__(self):
        return id(self)

    def __call__(self, node1, node2):
        return self.patristic_distance(node1, node2)

    def __copy__(self):
        return self.clone()

    def clone(self):
        o = self.__class__()
        o._tree_length = self._tree_length
        o._num_edges = self._num_edges
        for src, dest in (
                (self._node_phylogenetic_distances, o._node_phylogenetic_distances,),
                (self._node_phylogenetic_path_steps, o._node_phylogenetic_path_steps,),
                (self._mrca, o._mrca,),
                ):
            for t1 in src:
                dest[t1] = {}
                for t2 in src[t1]:
                    dest[t1][t2] = src[t1][t2]
        return o

    def mrca(self, node1, node2):
        """
        Returns MRCA of two node objects.
        """
        return self._mrca[node1][node2]

    def distance(self,
            node1,
            node2,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Returns distance between node1 and node2.
        """
        if is_weighted_edge_distances:
            return self.patristic_distance(node1, node2, is_normalize_by_tree_size=is_normalize_by_tree_size)
        else:
            return self.path_edge_count(node1, node2, is_normalize_by_tree_size=is_normalize_by_tree_size)


    def patristic_distance(self, node1, node2, is_normalize_by_tree_size=False):
        """
        Returns patristic distance between two node objects.
        """
        if node1 is node2:
            return 0.0
        d = self._node_phylogenetic_distances[node1][node2]
        if is_normalize_by_tree_size:
            return d / self._tree_length
        else:
            return d

    def path_edge_count(self, node1, node2, is_normalize_by_tree_size=False):
        """
        Returns the number of edges between two node objects.
        """
        if node1 is node2:
            return 0
        d = self._node_phylogenetic_path_steps[node1][node2]
        if is_normalize_by_tree_size:
            return float(d) / self._num_edges
        else:
            return d

    def distances(self,
            is_weighted_edge_distances=True,
            is_normalize_by_tree_size=False):
        """
        Returns list of patristic distances.
        """
        dmatrix, normalization_factor = self._get_distance_matrix_and_normalization_factor(
                is_weighted_edge_distances=is_weighted_edge_distances,
                is_normalize_by_tree_size=is_normalize_by_tree_size,
                )
        results = []
        nodes = list(dmatrix.keys())
        for node_idx1, node1 in enumerate(nodes[:-1]):
            for node_idx2, node2 in enumerate(nodes[node_idx1+1:]):
                results.append(dmatrix[node1][node2]/normalization_factor)
        return results

    def _get_distance_matrix_and_normalization_factor(self,
            is_weighted_edge_distances,
            is_normalize_by_tree_size):
        if is_weighted_edge_distances:
            dmatrix = self._node_phylogenetic_distances
            if is_normalize_by_tree_size:
                normalization_factor = self._tree_length
            else:
                normalization_factor = 1.0
        else:
            dmatrix = self._node_phylogenetic_path_steps
            if is_normalize_by_tree_size:
                normalization_factor = float(self._num_edges)
            else:
                normalization_factor = 1.0
        return dmatrix, normalization_factor

