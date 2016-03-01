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
Statistics, metrics, measurements, and values calculated on (single) trees.
"""

import math
import collections
import csv
import heapq
from dendropy.calculate import statistics
from dendropy.utility import GLOBAL_RNG
from dendropy.utility import container
from dendropy.utility import error
import dendropy

EULERS_CONSTANT = 0.5772156649015328606065120900824024310421

class PhylogeneticDistanceMatrix(object):
    """
    Calculates and maintains patristic distance information of taxa on a tree.
    ``max_dist_taxa`` and ``max_dist_nodes`` will return a tuple of taxon objects
    and corresponding nodes, respectively, that span the greatest path distance
    in the tree. The mid-point between the two is *guaranteed* to be on the
    closer to the first item of each pair.
    """

    @classmethod
    def from_tree(cls, tree, is_bipartitions_updated=False):
        pdm = cls()
        pdm.compile_from_tree(tree=tree,
                is_bipartitions_updated=is_bipartitions_updated)
        return pdm

    @classmethod
    def from_csv(cls,
            src,
            taxon_namespace=None,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            default_data_type=float,
            is_allow_new_taxa=None,
            label_transform_fn=None,
            **csv_reader_kwargs
            ):
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
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

    def compile_from_tree(self, tree, is_bipartitions_updated=False):
        """
        Calculates the distances. Note that the path length (in number of
        steps) between taxa that span the root will be off by one if
        the tree is unrooted.
        """
        self.clear()
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
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
                        node.desc_paths[desc1] = (desc1_plen + c1.edge.length, desc1_psteps + 1)
                        assert desc1.taxon is not None
                        if desc1.taxon not in self._taxon_phylogenetic_distances:
                            self._mapped_taxa.add(desc1.taxon)
                            self._taxon_phylogenetic_distances[desc1.taxon] = {}
                            self._taxon_phylogenetic_distances[desc1.taxon][desc1.taxon] = 0.0
                            self._taxon_phylogenetic_path_steps[desc1.taxon] = {}
                            self._taxon_phylogenetic_path_steps[desc1.taxon][desc1.taxon] = 0
                            self._mrca[desc1.taxon] = {}
                        for c2 in children[cidx1+1:]:
                            for desc2, (desc2_plen, desc2_psteps) in c2.desc_paths.items():
                                self._mapped_taxa.add(desc2.taxon)
                                self._mrca[desc1.taxon][desc2.taxon] = c1.parent_node
                                # self._all_distinct_mapped_taxa_pairs.add( tuple([desc1.taxon, desc2.taxon]) )
                                self._all_distinct_mapped_taxa_pairs.add( frozenset([desc1.taxon, desc2.taxon]) )
                                pat_dist = node.desc_paths[desc1][0] + desc2_plen + c2.edge.length
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
        if taxon1 is taxon2:
            return taxon1
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
                mpd = \\frac{ \\sum_{i}^{n} \\sum_{j}^{n} \\delta_{i,j} }{\\choose{n,2}},

        where :math:`i \\neq j`, :math:`\\delta_{i,j}` is the phylogenetic
        distance between species :math:`i` and :math:`j`, and $n$ is the number
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
            from dendropy.calculate import treemeasure
            tree = dendropy.Tree.get(path="data.nex",
                    schema="nexus")
            pdm = treemeasure.PhylogeneticDistanceMatrix(tree)

            # consider all tip
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
        distance between species :math:`i` and :math:`j`, and $n$ is the number
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
            from dendropy.calculate import treemeasure
            tree = dendropy.Tree.get(path="data.nex",
                    schema="nexus")
            pdm = treemeasure.PhylogeneticDistanceMatrix(tree)

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
            null_model_type="taxa.label",
            rng=None):
        """
        Returns the standardized effect size value for the MPD statistic under
        a null model under various community compositions.

        The S.E.S. is given by:

            .. math::
                SES(statistic) = \\frac{observed - mean(null)}{standard_deviation(null)} [1]

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

            tree = dendropy.Tree.get_from_path(
                    src="data/community.tree.newick",
                    schema="newick",
                    rooting="force-rooted")
            pdm = treemeasure.PhylogeneticDistanceMatrix.from_tree(tree)
            assemblage_memberships = pdm.read_assemblage_memberships_from_delimited_source("data/comm1.csv")
            results = pdm.standardized_effect_size_mean_pairwise_distance(assemblage_memberships=assemblage_memberships)
            print(results)
                                statistic under the null model.
        """
        if assemblage_memberships is None:
            assemblage_memberships = [ set(self._mapped_taxa) ]
        comparison_regimes = []
        for assemblage_membership in assemblage_memberships:
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
            null_model_type="taxa.label",
            rng=None):
        """
        Returns the standardized effect size value for the MNTD statistic under
        a null model under various community compositions.

        The S.E.S. is given by:

            .. math::
                SES(statistic) = \\frac{observed - mean(null)}{standard_deviation(null)} [1]

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

            tree = dendropy.Tree.get_from_path(
                    src="data/community.tree.newick",
                    schema="newick",
                    rooting="force-rooted")
            pdm = treemeasure.PhylogeneticDistanceMatrix.from_tree(tree)
            assemblage_memberships = pdm.read_assemblage_memberships_from_delimited_source("data/comm1.csv")
            results = pdm.standardized_effect_size_mean_nearest_taxon_distance(assemblage_memberships=assemblage_memberships)
            print(results)

        """
        if assemblage_memberships is None:
            assemblage_memberships = [ set(self._mapped_taxa) ]
        comparison_regimes = []
        for assemblage_membership in assemblage_memberships:
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

            # Read data from a CSV file into a PhylogeneticDistanceMatrix
            # object
            with open("distance_matrix.csv") as src:
                pdm = treemeasure.PhylogeneticDistanceMatrix.from_csv(
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
        node_pool = []
        for t1 in self._mapped_taxa:
            nd = tree.node_factory()
            nd.taxon = t1
            node_pool.append(nd)
        n = len(self._mapped_taxa)

        working_dmatrix = {}
        for nd1 in node_pool:
            working_dmatrix[nd1] = {}
            for nd2 in node_pool:
                if nd1 is nd2:
                    continue
                working_dmatrix[nd1][nd2] = original_dmatrix[nd1.taxon][nd2.taxon]
        # for idx1, nd1 in enumerate(node_pool[:-1]):
        #     working_dmatrix[nd1] = {}
        #     for idx2, nd2 in enumerate(node_pool[idx1+1:]):
        #         working_dmatrix[nd1][nd2] = original_dmatrix[nd1.taxon][nd2.taxon]
        while len(node_pool) > 1:
            # _dump_d(working_dmatrix, node_pool)
            qmatrix = []
            xsub_values = {}
            for idx1, nd1 in enumerate(node_pool[:-1]):
                for idx2, nd2 in enumerate(node_pool[idx1+1:]):
                    v1 = (n - 2) * working_dmatrix[nd1][nd2]
                    xsub1 = []
                    xsub2 = []
                    for ndx in node_pool:
                        if ndx is not nd1:
                            try:
                                xsub1.append( working_dmatrix[nd1][ndx] )
                            except KeyError:
                                xsub1.append( working_dmatrix[ndx][nd1] )
                                pass
                        if ndx is not nd2:
                            try:
                                xsub2.append( working_dmatrix[nd2][ndx] )
                            except KeyError:
                                xsub2.append( working_dmatrix[ndx][nd2] )
                                pass
                    xsubv_nd1 = sum(xsub1)
                    xsubv_nd2 = sum(xsub2)
                    xsub_values[nd1] = xsubv_nd1
                    xsub_values[nd2] = xsubv_nd2
                    qvalue = v1 - xsubv_nd1 - xsubv_nd2
                    heapq.heappush(qmatrix, (qvalue, (nd1, nd2)))
            # _dump_q(qmatrix)
            to_join = heapq.heappop(qmatrix)
            node_to_join1 = to_join[1][0]
            node_to_join2 = to_join[1][1]
            node_pool.remove(node_to_join1)
            node_pool.remove(node_to_join2)

            new_node = tree.node_factory()
            new_node.add_child(node_to_join1)
            new_node.add_child(node_to_join2)

            working_dmatrix[new_node] = {}
            for node in node_pool:
                try:
                    v1 = working_dmatrix[node][node_to_join1]
                except KeyError:
                    v1 = working_dmatrix[node_to_join1][node]
                try:
                    v2 = working_dmatrix[node][node_to_join2]
                except KeyError:
                    v2 = working_dmatrix[node_to_join2][node]
                try:
                    v3 = working_dmatrix[node_to_join1][node_to_join2]
                except KeyError:
                    v3 = working_dmatrix[node_to_join2][node_to_join1]
                dist = 0.5 * (v1 + v2 - v3)
                # dist = 0.5 * (v1 - node_to_join1.edge.length) + 0.5 * (v2 - node_to_join2.edge.length)
                working_dmatrix[node][new_node] = dist
                working_dmatrix[new_node][node] = dist
            node_pool.append(new_node)

            if n > 2:
                v1 = 0.5 * working_dmatrix[node_to_join1][node_to_join2]
                v4  = 1.0/(2*(n-2)) * (xsub_values[node_to_join1] - xsub_values[node_to_join2])
                delta_f = v1 + v4
                delta_g = working_dmatrix[node_to_join1][node_to_join2] - delta_f
                node_to_join1.edge.length = delta_f
                node_to_join2.edge.length = delta_g
            else:
                d = working_dmatrix[node_to_join1][node_to_join2]
                node_to_join1.edge.length = d / 2
                node_to_join2.edge.length = d / 2

            n = len(node_pool)
                # for nd in node_pool:
                #     for ch_nd in nd.child_node_iter():
                #         if ch_nd.edge.length is None:
                #             try:
                #                 ch_nd.edge.length = working_dmatrix[nd][ch_nd]
                #             except KeyError:
                #                 ch_nd.edge.length = working_dmatrix[ch_nd][nd]

        tree.seed_node = node_pool[0]
        return tree

    def upgma_tree(self,
            is_weighted_edge_distances=True,
            tree_factory=None,
            ):
        if is_weighted_edge_distances:
            original_dmatrix = self._taxon_phylogenetic_distances
        else:
            original_dmatrix = self._taxon_phylogenetic_path_steps
        if tree_factory is None:
            tree_factory = dendropy.Tree
        tree = tree_factory(taxon_namespace=self.taxon_namespace)
        node_pool = []
        for t1 in self._mapped_taxa:
            nd = tree.node_factory()
            nd.taxon = t1
            nd._upgma_cluster = set([nd])
            nd._upgma_distance_from_tip = 0.0
            node_pool.append(nd)
        working_dmatrix = {}
        for idx1, nd1 in enumerate(node_pool[:-1]):
            working_dmatrix[nd1] = {}
            for idx2, nd2 in enumerate(node_pool[idx1+1:]):
                d = original_dmatrix[nd1.taxon][nd2.taxon]
                working_dmatrix[nd1][nd2] = d
        while len(node_pool) > 1:
            min_distance = None
            nodes_to_join = None
            for idx1, nd1 in enumerate(node_pool[:-1]):
                for idx2, nd2 in enumerate(node_pool[idx1+1:]):
                    d = working_dmatrix[nd1][nd2]
                    if min_distance is None or d < min_distance:
                        nodes_to_join = (nd1, nd2)
                        min_distance = d
            new_node = tree.node_factory()
            new_node._upgma_cluster = set()
            elen = min_distance / 2.0
            for node_to_join in nodes_to_join:
                new_node.add_child(node_to_join)
                new_node._upgma_cluster.update(node_to_join._upgma_cluster)
                node_to_join.edge.length = elen - node_to_join._upgma_distance_from_tip
                node_pool.remove(node_to_join)
            new_node._upgma_distance_from_tip = nodes_to_join[0].edge.length + nodes_to_join[0]._upgma_distance_from_tip

            # diff1 =  nodes_to_join[0]._upgma_distance_from_tip - nodes_to_join[1]._upgma_distance_from_tip
            # nodes_to_join[0].edge.length = min_distance/2.0 - diff1
            # nodes_to_join[1].edge.length = min_distance/2.0 + diff1
            # new_node._upgma_distance_from_tip = nodes_to_join[0].edge.length + nodes_to_join[0]._upgma_distance_from_tip
            # print("--")
            # for node_to_join in nodes_to_join:
            #     print("{}: {}".format(
            #         "+".join(x.taxon.label for x in node_to_join._upgma_cluster),
            #         node_to_join.edge.length))
            # print("--")

            working_dmatrix[new_node] = {}
            for idx1, nd1 in enumerate(node_pool):
                d1 = 0.0
                count = 0.0
                # assert nodes_to_join[0] is not nodes_to_join[1]
                for node_to_join in nodes_to_join:
                    try:
                        d2 = working_dmatrix[node_to_join][nd1]
                    except KeyError:
                        d2 = working_dmatrix[nd1][node_to_join]
                    xc = len(node_to_join._upgma_cluster)
                    d1 += (d2 * xc)
                    count += xc
                d = d1 / count
                try:
                    working_dmatrix[nd1][new_node] = d
                except KeyError:
                    working_dmatrix[nd1] = {new_node: d}
                working_dmatrix[new_node][nd1] = d
                # print("{} vs {} = {}".format(
                #     "+".join(x.taxon.label for x in new_node._upgma_cluster),
                #     "+".join(x.taxon.label for x in nd1._upgma_cluster),
                #     d))
            for node_to_join in nodes_to_join:
                del node_to_join._upgma_cluster
                del node_to_join._upgma_distance_from_tip
            node_pool.append(new_node)
        tree.seed_node = node_pool[0]
        del tree.seed_node._upgma_cluster
        del tree.seed_node._upgma_distance_from_tip
        # print(tree.as_string("newick"))
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

    def read_assemblage_memberships_from_delimited_source(
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
        assemblage_memberships = []
        for row_name in data_table.row_name_iter():
            assemblage_membership = set()
            for taxon in self.taxon_iter():
                if data_table[row_name, taxon.label] > 0:
                    assemblage_membership.add(taxon)
            assemblage_memberships.append(assemblage_membership)
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

## legacy: will soon be deprecated
class PatrisiticDistanceMatrix(PhylogeneticDistanceMatrix):

    def __init__(self, tree):
        PhylogeneticDistanceMatrix.__init__(self)
        self.compile_from_tree(tree=tree)

def patristic_distance(tree, taxon1, taxon2, is_bipartitions_updated=False):
    """
    Given a tree with bipartitions encoded, and two taxa on that tree, returns the
    patristic distance between the two. Much more inefficient than constructing
    a PhylogeneticDistanceMatrix object.
    """
    mrca = tree.mrca(taxa=[taxon1, taxon2], is_bipartitions_updated=is_bipartitions_updated)
    dist = 0
    n = tree.find_node(lambda x: x.taxon == taxon1)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    n = tree.find_node(lambda x: x.taxon == taxon2)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    return dist

###########################################################################
### Metrics -- Unary

def B1(tree):
    """
    Returns the B1 statistic: the reciprocal of the sum of the maximum
    number of nodes between each interior node and tip over all internal
    nodes excluding root.
    """
    b1 = 0.0
    nd_mi = {}
    for nd in tree.postorder_node_iter():
        if nd._parent_node is None:
            continue
        child_nodes = nd._child_nodes
        if len(child_nodes) == 0:
            nd_mi[nd] = 0.0
            continue
        mi = max(nd_mi[ch] for ch in child_nodes)
        mi += 1
        nd_mi[nd] = mi
        b1 += 1.0/mi
    return b1

def colless_tree_imbalance(tree, normalize="max"):
    """
    Returns Colless' tree imbalance or I statistic: the sum of differences
    of numbers of children in left and right subtrees over all internal
    nodes. ``normalize`` specifies the normalization:

        * "max" or True [DEFAULT]
            normalized to maximum value for tree of
            this size
        * "yule"
            normalized to the Yule model
        * "pda"
            normalized to the PDA (Proportional to Distinguishable
            Arrangements) model
        * None or False
            no normalization

    """
    colless = 0.0
    num_leaves = 0
    subtree_leaves = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            subtree_leaves[nd] = 1
            num_leaves += 1
        else:
            total_leaves = 0
            if len(nd._child_nodes) > 2:
                raise TypeError("Colless' tree imbalance statistic requires strictly bifurcating trees")
            left = subtree_leaves[nd._child_nodes[0]]
            right = subtree_leaves[nd._child_nodes[1]]
            colless += abs(right-left)
            subtree_leaves[nd] = right + left
    if normalize == "yule":
        colless = float(colless - (num_leaves * math.log(num_leaves)) - (num_leaves * (EULERS_CONSTANT - 1.0 - math.log(2))))/num_leaves
    elif normalize == "pda":
        colless = colless / pow(num_leaves, 3.0/2)
    elif normalize is True or normalize == "max":
        ## note that Mooers 1995 (Evolution 49(2):379-384)
        ## remarks that the correct normalization factor is
        ## 2/((num_leaves - 1) * (num_leaves -2))
        colless = colless * (2.0/(num_leaves * (num_leaves-3) + 2))
    elif normalize is not None and normalize is not False:
        raise TypeError("``normalization`` accepts only None, True, False, 'yule' or 'pda' as argument values")
    return colless

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).

    If the lengths of different paths to the node differ by more than ``prec``,
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.

    As a side effect a ``age`` attribute is added to the nodes of the tree.

    Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
    molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
    (267). 2267-2272
    """
    # the equation is given by:
    #   T = \sum_{j=2}^n (jg_j)
    #   C = T \sqrt{\frac{1}{12(n-2)}}
    #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
    # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
    #   times between consecutive (in time, not along a branch) speciation times.
    node = None
    speciation_ages = []
    n = 0
    if tree.seed_node.age is None:
        tree.calc_node_ages(ultrametricity_precision=prec)
    for node in tree.postorder_node_iter():
        if len(node.child_nodes()) == 2:
            speciation_ages.append(node.age)
        else:
            n += 1
    if node is None:
        raise ValueError("Empty tree encountered")
    speciation_ages.sort(reverse=True)
    g = []
    older = speciation_ages[0]
    for age in speciation_ages[1:]:
        g.append(older - age)
        older = age
    g.append(older)
    if not g:
        raise ValueError("No internal nodes found (other than the root)")
    assert(len(g) == (n - 1))
    T = 0.0
    accum = 0.0
    for i in range(2, n):
        list_index = i - 2
        T += i * float(g[list_index])
        accum += T
    list_index = n - 2
    T += (n) * g[list_index]
    nmt = n - 2.0
    numerator = accum/nmt - T/2.0
    C = T*pow(1/(12*nmt), 0.5)
    return numerator/C

def N_bar(tree):
    """
    Returns the $\bar{N}$ statistic: the average number of nodes above a
    terminal node.
    """
    leaf_count = 0
    nbar = 0
    for leaf_node in tree.leaf_node_iter():
        leaf_count += 1
        for parent in leaf_node.ancestor_iter(inclusive=False):
            nbar += 1
    return float(nbar) / leaf_count

def sackin_index(tree, normalize=True):
    """
    Returns the Sackin's index: the sum of the number of ancestors for each
    tip of the tree. The larger the Sackin's index, the less balanced the
    tree. ``normalize`` specifies the normalization:

        * True [DEFAULT]
            normalized to number of leaves; this results in a value
            equivalent to that given by Tree.N_bar()
        * "yule"
            normalized to the Yule model
        * "pda"
            normalized to the PDA (Proportional to Distinguishable
            Arrangements) model
        * None or False
            no normalization

    """
    leaf_count = 0
    num_anc = 0
    for leaf_node in tree.leaf_node_iter():
        leaf_count += 1
        for parent in leaf_node.ancestor_iter(inclusive=False):
            num_anc += 1
    if normalize == "yule":
        x = sum(1.0/j for j in range(2, leaf_count+1))
        s = float(num_anc - (2 * leaf_count * x))/leaf_count
    elif normalize == "pda":
        s = float(num_anc)/(pow(leaf_count, 3.0/2))
    elif normalize is True:
        s = float(num_anc)/leaf_count
    elif normalize is None or normalize is False:
        s = float(num_anc)
    elif normalize is not None and normalize is not False:
        raise TypeError("``normalization`` accepts only None, True, False, 'yule' or 'pda' as argument values")
    return s

def treeness(tree):
    """
    Returns the proportion of total tree length that is taken up by
    internal branches.
    """
    internal = 0.0
    external = 0.0
    for nd in tree.postorder_node_iter():
        if not nd._parent_node:
            continue
        if nd.is_leaf():
            external += nd.edge.length
        else:
            internal += nd.edge.length
    return internal/(external + internal)

