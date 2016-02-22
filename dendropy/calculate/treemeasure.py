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
from dendropy.calculate import statistics
from dendropy.utility import GLOBAL_RNG
from dendropy.utility import container

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
        pdc = cls()
        pdc.compile(tree=tree,
                is_bipartitions_updated=is_bipartitions_updated)
        return pdc

    def __init__(self):
        self.clear()

    def clear(self):
        self.taxon_namespace = None
        self._mapped_taxa = set()
        self._taxon_phylogenetic_distances = {}
        self._taxon_phylogenetic_path_steps = {}
        self._mrca = {}

    def compile(self, tree, is_bipartitions_updated=False):
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
        for node in tree.postorder_node_iter():
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
                            self._taxon_phylogenetic_path_steps[desc1.taxon] = {}
                            self._mrca[desc1.taxon] = {}
                        for c2 in children[cidx1+1:]:
                            for desc2, (desc2_plen, desc2_psteps) in c2.desc_paths.items():
                                self._mapped_taxa.add(desc2.taxon)
                                self._mrca[desc1.taxon][desc2.taxon] = c1.parent_node
                                pat_dist = node.desc_paths[desc1][0] + desc2_plen + c2.edge.length
                                self._taxon_phylogenetic_distances[desc1.taxon][desc2.taxon] = pat_dist
                                path_steps = node.desc_paths[desc1][1] + desc2_psteps + 1
                                self._taxon_phylogenetic_path_steps[desc1.taxon][desc2.taxon] = path_steps
                    del(c1.desc_paths)

    def __eq__(self, o):
        if self.taxon_namespace is not o.taxon_namespace:
            return False
        return (True
                and (self._mapped_taxa == o._mapped_taxa)
                and (self._taxon_phylogenetic_distances == o._taxon_phylogenetic_distances)
                and (self._taxon_phylogenetic_path_steps == o._taxon_phylogenetic_path_steps)
                and (self._mrca == o._mrca)
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

    def patristic_distance(self, taxon1, taxon2):
        """
        Returns patristic distance between two taxon objects.
        """
        if taxon1 is taxon2:
            return 0.0
        try:
            return self._taxon_phylogenetic_distances[taxon1][taxon2]
        except KeyError:
            return self._taxon_phylogenetic_distances[taxon2][taxon1]

    def mrca(self, taxon1, taxon2):
        """
        Returns MRCA of two taxon objects.
        """
        if taxon1 is taxon2:
            return taxon1
        try:
            return self._mrca[taxon1][taxon2]
        except KeyError:
            return self._mrca[taxon2][taxon1]

    def path_edge_count(self, taxon1, taxon2):
        """
        Returns the number of edges between two taxon objects.
        """
        if taxon1 is taxon2:
            return 0
        try:
            return self._taxon_phylogenetic_path_steps[taxon1][taxon2]
        except KeyError:
            return self._taxon_phylogenetic_path_steps[taxon2][taxon1]

    def max_pairwise_distance_taxa(self, is_weighted_edge_distances=True):
        if is_weighted_edge_distances:
            dists = self._taxon_phylogenetic_distances
        else:
            dists = self._taxon_phylogenetic_path_steps
        max_dist = None
        max_dist_taxa = None
        for t1 in dists:
            for t2 in dists[t1]:
                pat_dist = dists[t1][t2]
                if max_dist is None or pat_dist > max_dist:
                    max_dist = pat_dist
                    max_dist_taxa = (t1, t2)
        return max_dist_taxa

    def distances(self):
        """
        Returns list of patristic distances.
        """
        dists = []
        for dt in self._taxon_phylogenetic_distances.values():
            for d in dt.values():
                dists.append(d)
        return dists

    def sum_of_distances(self):
        """
        Returns sum of patristic distances on tree.
        """
        return sum(self.distances())

    def iter_taxa(self):
        """
        Iterates over taxa in matrix. Note that this could be a subset of the taxa in
        the associated taxon namespace.
        """
        for t in self._mapped_taxa:
            yield t

    def mean_pairwise_distance(self, filter_fn=None, is_weighted_edge_distances=True):
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
        if is_weighted_edge_distances:
            dmatrix = self._taxon_phylogenetic_distances
        else:
            dmatrix = self._taxon_phylogenetic_path_steps
        seen_comps = set()
        distances = []
        for taxon1 in dmatrix:
            if filter_fn and not filter_fn(taxon1):
                continue
            for taxon2 in dmatrix:
                if taxon1 is taxon2:
                    continue
                if filter_fn and not filter_fn(taxon2):
                    continue
                comp_hash = frozenset([taxon1, taxon2])
                if comp_hash in seen_comps:
                    continue
                distances.append(self.__call__(taxon1, taxon2))
                seen_comps.add( comp_hash )
        if distances:
            return sum(distances) / (len(distances) * 1.0)
        else:
            return 0

    def standardized_effect_size_mean_pairwise_distance(self,
            assemblage_memberships,
            num_randomization_replicates=1000,
            is_weighted_edge_distances=True,
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
        """
        results = self._calculate_standardized_effect_size(
                statisticf_name="mean_pairwise_distance",
                statisticf_kwargs={"is_weighted_edge_distances": is_weighted_edge_distances},
                assemblage_memberships=assemblage_memberships,
                null_model_type=null_model_type,
                num_randomization_replicates=num_randomization_replicates,
                rng=rng)
        return results

    def mean_nearest_taxon_distance(self, filter_fn=None, is_weighted_edge_distances=True):
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
        if is_weighted_edge_distances:
            dmatrix = self._taxon_phylogenetic_distances
        else:
            dmatrix = self._taxon_phylogenetic_path_steps
        distances = []
        for taxon1 in dmatrix:
            if filter_fn and not filter_fn(taxon1):
                continue
            subdistances = []
            for taxon2 in dmatrix:
                if taxon1 is taxon2:
                    continue
                if filter_fn and not filter_fn(taxon2):
                    continue
                subdistances.append(self.__call__(taxon1, taxon2))
            if subdistances:
                distances.append(min(subdistances))
        if distances:
            return sum(distances) / (len(distances) * 1.0)
        else:
            return 0

    def standardized_effect_size_mean_nearest_taxon_distance(self,
            assemblage_memberships,
            num_randomization_replicates=1000,
            is_weighted_edge_distances=True,
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
        """
        results = self._calculate_standardized_effect_size(
                statisticf_name="mean_nearest_taxon_distance",
                statisticf_kwargs={"is_weighted_edge_distances": is_weighted_edge_distances},
                assemblage_memberships=assemblage_memberships,
                null_model_type=null_model_type,
                num_randomization_replicates=num_randomization_replicates,
                rng=rng)
        return results

    def shuffle_taxa(self, rng=None):
        """
        Randomly shuffles taxa in-situ.
        """
        if rng is None:
            rng = GLOBAL_RNG
        # mapped_taxa = list(mapped_taxa)
        # rng.shuffle(mapped_taxa)
        reordered_taxa = list(self._mapped_taxa)
        rng.shuffle(reordered_taxa)
        current_to_shuffled_taxon_map = dict(zip(self._mapped_taxa, reordered_taxa))
        # assert len(current_to_shuffled_taxon_map) == len(self._mapped_taxa), "{} != {}".format(len(current_to_shuffled_taxon_map), len(self._mapped_taxa))
        for attr_name in (
                "_taxon_phylogenetic_distances",
                "_taxon_phylogenetic_path_steps",
                "_mrca",
                ):
            src = getattr(self, attr_name)
            dest = {}
            for t1 in src:
                x1 = current_to_shuffled_taxon_map[t1]
                dest[x1] = {}
                for t2 in src[t1]:
                    x2 = current_to_shuffled_taxon_map[t2]
                    dest[x1][x2] = src[t1][t2]
            setattr(self, attr_name, dest)
        return current_to_shuffled_taxon_map

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

    def read_assemblage_memberships_from_delimited_source(
            self,
            src,
            delimiter=",",
            default_data_type=int,):
        """
        Convenience method to return list of community sets from a delimited
        file that lists taxon (labels) in columns and community
        presence/absences or abundances in rows.
        """
        if isinstance(src, str):
            with open(src) as srcf:
                reader = csv.reader(srcf, delimiter=delimiter)
                data_table = container.DataTable.get_from_csv_reader(reader, default_data_type=default_data_type)
        else:
            reader = csv.reader(src, delimiter=delimiter)
            data_table = container.DataTable.get_from_csv_reader(reader, default_data_type=default_data_type)
        mapped_taxon_labels = set([taxon.label for taxon in self.iter_taxa()])
        for column_name in data_table.iter_column_names():
            assert column_name in mapped_taxon_labels
        assemblage_memberships = []
        for row_name in data_table.iter_row_names():
            assemblage_membership = set()
            for taxon in self.iter_taxa():
                if data_table[row_name, taxon.label] > 0:
                    assemblage_membership.add(taxon)
            assemblage_memberships.append(assemblage_membership)
        return assemblage_memberships

    def _calculate_standardized_effect_size(self,
            statisticf_name,
            statisticf_kwargs=None,
            assemblage_memberships=None,
            null_model_type="taxa.label",
            num_randomization_replicates=1000,
            rng=None):
        result_type = collections.namedtuple("PhylogeneticCommunityStandardizedEffectSizeStatisticCalculationResult",
                ["obs", "null_model_mean", "null_model_sd", "z", "rank", "p",])
        if assemblage_memberships is None:
            assemblage_memberships = [ set(self._mapped_taxa) ]
        if statisticf_kwargs is None:
            statisticf_kwargs = {}
        observed_stat_values = {}
        null_model_stat_values = {}
        null_model_matrix = self.clone()
        assert null_model_matrix == self
        for rep_idx in range(num_randomization_replicates):
            null_model_matrix.shuffle_taxa(rng=rng)
            for community_idx, assemblage_membership in enumerate(assemblage_memberships):
                filter_fn = lambda taxon: taxon in assemblage_membership
                statisticf_kwargs["filter_fn"] = filter_fn
                if rep_idx == 0:
                    observed_stat_values[community_idx] = getattr(self, statisticf_name)(**statisticf_kwargs)
                    null_model_stat_values[community_idx] = []
                stat_value = getattr(null_model_matrix, statisticf_name)(**statisticf_kwargs)
                null_model_stat_values[community_idx].append(stat_value)
        results = []
        for community_idx, assemblage_membership in enumerate(assemblage_memberships):
            obs_value = observed_stat_values[community_idx]
            stat_values = null_model_stat_values[community_idx]
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
        self.compile(tree=tree)

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

