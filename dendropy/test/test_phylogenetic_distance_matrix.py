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

import unittest
import dendropy
import csv
from dendropy.utility import container
from dendropy.utility.textprocessing import StringIO
from dendropy.test.support import pathmap
from dendropy.calculate import treemeasure
from dendropy.calculate import probability
from dendropy.calculate import combinatorics

class PhylogeneticDistanceMatrixCloneTest(unittest.TestCase):

    def setUp(self):
        self.tree = dendropy.Tree.get_from_string("(((a:1, b:1):1, c:2):1, (d:2, (e:1,f:1):1):1):0;", schema="newick")

    def test_clone(self):
        pdm0 = self.tree.phylogenetic_distance_matrix()
        pdm1 = pdm0.clone()
        self.assertIsNot(pdm0, pdm1)
        self.assertIs(pdm0.taxon_namespace, pdm1.taxon_namespace)
        self.assertEqual(len(pdm0.taxon_namespace), len(pdm0._mapped_taxa))
        self.assertEqual(len(pdm1.taxon_namespace), len(pdm1._mapped_taxa))
        for src, dest in (
                    (pdm0._taxon_phylogenetic_distances, pdm1._taxon_phylogenetic_distances,),
                    (pdm0._taxon_phylogenetic_path_steps, pdm1._taxon_phylogenetic_path_steps,),
                    (pdm0._mrca, pdm1._mrca,),
                ):
            self.assertIsNot(src, dest)
            for t1 in src:
                self.assertIn(t1, dest)
                self.assertIsNot(src[t1], dest[t1])
        for t1 in self.tree.taxon_namespace:
            for t2 in self.tree.taxon_namespace:
                self.assertEqual(pdm0.patristic_distance(t1, t2), pdm1.patristic_distance(t1, t2))
                self.assertEqual(pdm0.mrca(t1, t2), pdm1.mrca(t1, t2))
                self.assertEqual(pdm0.path_edge_count(t1, t2), pdm1.path_edge_count(t1, t2))
        self.assertEqual(set(pdm0.distances()), set(pdm1.distances()))
        self.assertEqual(pdm0.sum_of_distances(), pdm1.sum_of_distances())
        self.assertEqual(pdm0, pdm1)

class PhylogeneticDistanceMatrixCompileTest(unittest.TestCase):

        def setUp(self):
            # library(ape)
            # tree = read.nexus("data/pythonidae.mle.nex")
            # pdm = cophenetic.phylo(tree)
            # write.csv(format(pdm,digits=22), "pythonidae.mle.weighted.pdm.csv")
            with open(pathmap.other_source_path("pythonidae.mle.weighted.pdm.csv")) as src:
                self.reference_pdm_weighted_table = container.DataTable.from_csv(src, default_data_type=float, delimiter=",")
            with open(pathmap.other_source_path("pythonidae.mle.unweighted.pdm.csv")) as src:
                self.reference_pdm_unweighted_table = container.DataTable.from_csv(src, default_data_type=float, delimiter=",")
            self.tree = dendropy.Tree.get(path=pathmap.tree_source_path(
                "pythonidae.mle.nex"),
                schema="nexus",
                preserve_underscores=True)
            self.pdm = self.tree.phylogenetic_distance_matrix()

        def test_mapped_taxa(self):
            n1 = len(self.tree.taxon_namespace)
            self.assertEqual(self.pdm.taxon_namespace, self.tree.taxon_namespace)
            self.assertEqual(n1, len(self.tree.taxon_namespace))
            for taxon1 in self.tree.taxon_namespace:
                self.assertIn(taxon1, self.pdm._mapped_taxa)
            for taxon in self.pdm.taxon_iter():
                self.assertIn(taxon1, self.pdm._mapped_taxa)
                self.assertIn(taxon1, self.tree.taxon_namespace)

        def test_all_distinct_mapped_taxa_pairs(self):
            n1 = len(self.tree.taxon_namespace)
            taxon_pair_iter1 = iter(self.pdm._all_distinct_mapped_taxa_pairs)
            taxon_pair_iter2 = self.pdm.distinct_taxon_pair_iter()
            for tpi in (taxon_pair_iter1, taxon_pair_iter2):
                seen_pairs = set()
                visited_taxa = set()
                for taxon1, taxon2 in tpi:
                    s = frozenset([taxon1, taxon2])
                    self.assertIn(taxon1, self.pdm._mapped_taxa)
                    self.assertIn(taxon1, self.tree.taxon_namespace)
                    self.assertIn(taxon2, self.pdm._mapped_taxa)
                    self.assertIn(taxon2, self.tree.taxon_namespace)
                    self.assertNotIn(s, seen_pairs)
                    seen_pairs.add(s)
                    visited_taxa.add(taxon1)
                    visited_taxa.add(taxon2)
                self.assertEqual(len(visited_taxa), n1)
                self.assertEqual(len(seen_pairs), combinatorics.choose(n1, 2))

        def test_tree_length(self):
            self.assertEqual(self.pdm._tree_length, self.tree.length())

        def test_tree_num_edges(self):
            self.assertEqual(self.pdm._num_edges, combinatorics.num_edges_on_tree(
                num_leaves=len(self.tree.taxon_namespace), is_rooted=True))

        def test_tree_weighted_unnormalized_pairwise_distance_query(self):
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    exp = self.reference_pdm_weighted_table[taxon1.label, taxon2.label]
                    obs1 = self.pdm._taxon_phylogenetic_distances[taxon1][taxon2]
                    self.assertAlmostEqual(obs1, exp, 6)
                    obs2 = self.pdm.patristic_distance(taxon1, taxon2)
                    self.assertAlmostEqual(obs2, exp, 6)
                    obs3 = self.pdm.distance(
                            taxon1,
                            taxon2,
                            is_weighted_edge_distances=True,
                            is_normalize_by_tree_size=False)
                    self.assertAlmostEqual(obs3, exp, 6)

        def test_tree_weighted_normalized_pairwise_distance_query(self):
            tree_length = self.tree.length()
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    exp = self.reference_pdm_weighted_table[taxon1.label, taxon2.label] / tree_length
                    obs2 = self.pdm.patristic_distance(taxon1, taxon2, is_normalize_by_tree_size=True)
                    self.assertAlmostEqual(obs2, exp, 6)
                    obs3 = self.pdm.distance(
                            taxon1,
                            taxon2,
                            is_weighted_edge_distances=True,
                            is_normalize_by_tree_size=True)
                    self.assertAlmostEqual(obs3, exp, 6)

        def test_tree_unweighted_unnormalized_pairwise_distance_query(self):
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    exp = self.reference_pdm_unweighted_table[taxon1.label, taxon2.label]
                    obs1 = self.pdm._taxon_phylogenetic_path_steps[taxon1][taxon2]
                    self.assertAlmostEqual(obs1, exp, 6)
                    obs2 = self.pdm.path_edge_count(taxon1, taxon2)
                    self.assertAlmostEqual(obs2, exp, 6)
                    obs3 = self.pdm.distance(
                            taxon1,
                            taxon2,
                            is_weighted_edge_distances=False,
                            is_normalize_by_tree_size=False)
                    self.assertAlmostEqual(obs3, exp, 6)

        def test_tree_unweighted_normalized_pairwise_distance_query(self):
            num_edges_on_tree = combinatorics.num_edges_on_tree(
                    num_leaves=len(self.tree.taxon_namespace), is_rooted=True)
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    exp = self.reference_pdm_unweighted_table[taxon1.label, taxon2.label] / num_edges_on_tree
                    obs2 = self.pdm.path_edge_count(taxon1, taxon2, is_normalize_by_tree_size=True)
                    self.assertAlmostEqual(obs2, exp, 6)
                    obs3 = self.pdm.distance(
                            taxon1,
                            taxon2,
                            is_weighted_edge_distances=False,
                            is_normalize_by_tree_size=True)
                    self.assertAlmostEqual(obs3, exp, 6)

        def test_tree_weighted_unnormalized_pairwise_distance_collection(self):
            # this test relies on populating the expected result list
            # in the *same order* as the function it is testing
            exp_list = []
            for taxon1, taxon2 in self.pdm.distinct_taxon_pair_iter():
                exp = self.reference_pdm_weighted_table[taxon1.label, taxon2.label]
                exp_list.append(exp)
            obs = self.pdm.distances(
                    is_weighted_edge_distances=True,
                    is_normalize_by_tree_size=False,
                    )
            self.assertEqual(len(obs), len(exp_list))
            for obs, exp in zip(obs, exp_list):
                self.assertAlmostEqual(obs, exp, 6)
            self.assertAlmostEqual(self.pdm.sum_of_distances(
                    is_weighted_edge_distances=True,
                    is_normalize_by_tree_size=False,
                    ), sum(exp_list), 6)

        def test_tree_weighted_normalized_pairwise_distance_collection(self):
            # this test relies on populating the expected result list
            # in the *same order* as the function it is testing
            tree_length = self.tree.length()
            exp_list = []
            for taxon1, taxon2 in self.pdm.distinct_taxon_pair_iter():
                exp = self.reference_pdm_weighted_table[taxon1.label, taxon2.label] / tree_length
                exp_list.append(exp)
            obs = self.pdm.distances(
                    is_weighted_edge_distances=True,
                    is_normalize_by_tree_size=True,
                    )
            self.assertEqual(len(obs), len(exp_list))
            for obs, exp in zip(obs, exp_list):
                self.assertAlmostEqual(obs, exp, 6)
            self.assertAlmostEqual(self.pdm.sum_of_distances(
                    is_weighted_edge_distances=True,
                    is_normalize_by_tree_size=True,
                    ), sum(exp_list), 6)

        def test_tree_unweighted_unnormalized_pairwise_distance_collection(self):
            # this test relies on populating the expected result list
            # in the *same order* as the function it is testing
            exp_list = []
            for taxon1, taxon2 in self.pdm.distinct_taxon_pair_iter():
                exp = self.reference_pdm_unweighted_table[taxon1.label, taxon2.label]
                exp_list.append(exp)
            obs = self.pdm.distances(
                    is_weighted_edge_distances=False,
                    is_normalize_by_tree_size=False,
                    )
            self.assertEqual(len(obs), len(exp_list))
            for obs, exp in zip(obs, exp_list):
                self.assertAlmostEqual(obs, exp, 6)
            self.assertAlmostEqual(self.pdm.sum_of_distances(
                    is_weighted_edge_distances=False,
                    is_normalize_by_tree_size=False,
                    ), sum(exp_list), 6)

        def test_tree_unweighted_normalized_pairwise_distance_collection(self):
            # this test relies on populating the expected result list
            # in the *same order* as the function it is testing
            num_edges_on_tree = combinatorics.num_edges_on_tree(
                    num_leaves=len(self.tree.taxon_namespace), is_rooted=True)
            exp_list = []
            for taxon1, taxon2 in self.pdm.distinct_taxon_pair_iter():
                exp = self.reference_pdm_unweighted_table[taxon1.label, taxon2.label] / num_edges_on_tree
                exp_list.append(exp)
            obs = self.pdm.distances(
                    is_weighted_edge_distances=False,
                    is_normalize_by_tree_size=True,
                    )
            self.assertEqual(len(obs), len(exp_list))
            for obs, exp in zip(obs, exp_list):
                self.assertAlmostEqual(obs, exp, 6)
            self.assertAlmostEqual(self.pdm.sum_of_distances(
                    is_weighted_edge_distances=False,
                    is_normalize_by_tree_size=True,
                    ), sum(exp_list), 6)

        def test_max_pairwise_weighted_distance_taxa(self):
            maxd = None
            maxt = None
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    d = self.reference_pdm_weighted_table[taxon1.label, taxon2.label]
                    if maxd is None or d > maxd:
                        maxd = d
                        maxt = set([taxon1, taxon2])
            obs_maxt = self.pdm.max_pairwise_distance_taxa(is_weighted_edge_distances=True)
            self.assertEqual(set(obs_maxt), maxt)

        def test_max_pairwise_weighted_distance_taxa(self):
            maxd = None
            maxt = None
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    d = self.reference_pdm_weighted_table[taxon1.label, taxon2.label]
                    if maxd is None or d > maxd:
                        maxd = d
                        maxt = set([taxon1, taxon2])
            obs_maxt = self.pdm.max_pairwise_distance_taxa(is_weighted_edge_distances=True)
            self.assertEqual(set(obs_maxt), maxt)

        def test_max_pairwise_unweighted_distance_taxa(self):
            maxd = None
            maxts = []
            for taxon1 in self.tree.taxon_namespace:
                for taxon2 in self.tree.taxon_namespace:
                    d = self.reference_pdm_unweighted_table[taxon1.label, taxon2.label]
                    if maxd is None or d > maxd:
                        maxd = d
                        maxts.append(set([taxon1, taxon2]))
                    elif d == maxd:
                        maxts.append(set([taxon1, taxon2]))
            obs_maxt = self.pdm.max_pairwise_distance_taxa(is_weighted_edge_distances=False)
            found_match = False
            for maxt in maxts:
                if set(obs_maxt) == maxt:
                    found_match = True
                    break
            else:
                self.fail()

class PhylogeneticDistanceMatrixShuffleTest(unittest.TestCase):

    def test_shuffle(self):
        tree = dendropy.Tree.get_from_path(
                    src=pathmap.tree_source_path("community.tree.newick"),
                    schema="newick",
                    rooting="force-rooted")
        pdc0 = tree.phylogenetic_distance_matrix()
        pdc1 = tree.phylogenetic_distance_matrix()
        current_to_shuffled_taxon_map = pdc1.shuffle_taxa()
        keys = set(current_to_shuffled_taxon_map.keys())
        values = set(current_to_shuffled_taxon_map.values())
        self.assertEqual(len(keys), len(values), "\n\n({}): {}\n\n({}): {}".format(len(keys), keys, len(values), values))
        self.assertEqual(keys, values)
        for taxon in tree.taxon_namespace:
            self.assertIn(taxon, current_to_shuffled_taxon_map)
        for nd in tree.leaf_node_iter():
            self.assertIn(current_to_shuffled_taxon_map[nd.taxon], tree.taxon_namespace)
            nd.taxon = current_to_shuffled_taxon_map[nd.taxon]
        pdc2 = tree.phylogenetic_distance_matrix()
        same_as_before = []
        different = []
        for t1 in tree.taxon_namespace:
            for t2 in tree.taxon_namespace:
                d2 = pdc2.patristic_distance(t1, t2)
                d1 = pdc1.patristic_distance(t1, t2)
                self.assertEqual(d1, d2)
                if t1 is not t2:
                    d0 = pdc0.patristic_distance(t1, t2)
                    if d0 == d1:
                        same_as_before.append( (t1, t2) )
                    else:
                        different.append( (t1, t2) )
                else:
                    self.assertEqual(d1, 0)
        self.assertTrue(len(different) > 0)
        self.assertEqual(pdc1, pdc2)
        self.assertNotEqual(pdc0, pdc1)

class TreePatristicDistTest(unittest.TestCase):

    def setUp(self):
        self.tree = dendropy.Tree.get_from_string("(((a:1, b:1):1, c:2):1, (d:2, (e:1,f:1):1):1):0;", schema="newick")

    def testPatDistMatrix(self):
        pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(self.tree)
        def _chk_distance(pdm, t1, t2, exp_distance):
            tax1 = self.tree.taxon_namespace.require_taxon(label=t1)
            tax2 = self.tree.taxon_namespace.require_taxon(label=t2)
            pd = pdm(tax1, tax2)
            self.assertEqual(pd, exp_distance, "{}: {} <-> {}: {} instead of {}".format(self.tree, t1, t2, pd, exp_distance))
        _chk_distance(pdm, "a", "b", 2)
        _chk_distance(pdm, "a", "c", 4)
        _chk_distance(pdm, "b", "c", 4)
        _chk_distance(pdm, "a", "d", 6)
        _chk_distance(pdm, "f", "d", 4)
        _chk_distance(pdm, "c", "d", 6)

    def testPatDistFunc(self):
        self.tree.encode_bipartitions()
        def _chk_distance(t1, t2, exp_distance):
            tax1 = self.tree.taxon_namespace.get_taxon(label=t1)
            tax2 = self.tree.taxon_namespace.get_taxon(label=t2)
            pd = treemeasure.patristic_distance(self.tree, tax1, tax2)
            self.assertEqual(pd, exp_distance)
        _chk_distance("a", "b", 2)
        _chk_distance("a", "c", 4)
        _chk_distance("b", "c", 4)
        _chk_distance("a", "d", 6)
        _chk_distance("f", "d", 4)
        _chk_distance("c", "d", 6)

class PhylogeneticEcologyStatsTests(unittest.TestCase):

    def setUp(self):
        self.tree = dendropy.Tree.get_from_path(
                src=pathmap.tree_source_path("community.tree.newick"),
                schema="newick",
                rooting="force-rooted")
        self.pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(self.tree)
        assemblage_data_filepath = pathmap.other_source_path("community.data.tsv")
        with open(assemblage_data_filepath) as src:
            self.data_table = container.DataTable.from_csv(src, default_data_type=int, delimiter="\t")
        self.assemblage_membership_definitions = self.pdm.assemblage_membership_definitions_from_csv(
                assemblage_data_filepath,
                delimiter="\t")
        self.assemblage_memberships = self.assemblage_membership_definitions.values()

    def _low_precision_equality(self, v1, v2, error=1.0):
        return abs(v1-v2) <= error

    def test_nonabundance_edgeweighted_unnormalized_mpd(self):
        # my.sample = read.table("data/PD.example.sample.txt", sep = "\t", row.names = 1, header = T)
        # my.phylo = read.tree("data/PD.example.phylo.txt")
        # pd.matrix = cophenetic(my.phylo)
        # mpd(my.sample, cophenetic(my.phylo), abundance.weighted=F)
        #    [1] 3.2225706087019050372078 1.9156605943056665974922 1.9156605943290001548007
        #    [4] 1.9395923093204667786438 0.1934132401466666650869
        expected_results = {
                "C1": 3.222570608701905037208,
                "C2": 1.915660594305666597492,
                "C3": 1.915660594329000154801,
                "C4": 1.939592309320466778644,
                "C5": 0.1934132401466666650869,
        }
        for row_name in self.data_table.row_name_iter():
            filter_fn = lambda taxon: self.data_table[row_name, taxon.label] > 0
            d = self.pdm.mean_pairwise_distance(filter_fn=filter_fn)
            self.assertAlmostEqual(d, expected_results[row_name])

    def test_nonabundance_edgeweighted_unnormalized_mntd(self):
        # my.sample = read.table("data/PD.example.sample.txt", sep = "\t", row.names = 1, header = T)
        # my.phylo = read.tree("data/PD.example.phylo.txt")
        # pd.matrix = cophenetic(my.phylo)
        # mntd(my.sample, cophenetic(my.phylo), abundance.weighted=F)
        #   [1] 1.6347319809428570991372 1.0891173926393333815099 1.0891173926543333827510
        #   [4] 0.1180230301583333335502 0.1426761318733333339104
        expected_results = {
                "C1": 1.6347319809428570991372,
                "C2": 1.0891173926393333815099,
                "C3": 1.0891173926543333827510,
                "C4": 0.1180230301583333335502,
                "C5": 0.1426761318733333339104,
        }
        for row_name in self.data_table.row_name_iter():
            filter_fn = lambda taxon: self.data_table[row_name, taxon.label] > 0
            d = self.pdm.mean_nearest_taxon_distance(filter_fn=filter_fn)
            self.assertAlmostEqual(d, expected_results[row_name])

    def test_nonabundance_edgeweighted_unnormalized_ses_mpd(self):
        # suppressMessages(library(picante))
        # dists = as.matrix(read.csv("data/dist.csv",header=T,row.names=1))
        # comm = as.matrix(read.csv("data/community.data.tsv",sep="\t",header=T,row.names=1))
        # results.mpd = ses.mpd(comm, dists, null.model="taxa.labels",abundance.weighted=F,runs=100000)
        # write.csv(format(results.mpd, digits=22), quote=F, "community.data.weighted.unnormalized.ses.mpd.csv")
        # results.mntd = ses.mntd(comm, dists, null.model="taxa.labels",abundance.weighted=F,runs=100000)
        # write.csv(format(results.mntd,digits=22), quote=F, "community.data.weighted.unnormalized.ses.mntd.csv")
        with open(pathmap.other_source_path("community.data.weighted.unnormalized.ses.mpd.csv")) as src:
            expected_results_data_table = container.DataTable.from_csv(src, default_data_type=float)
        # for row_name in expected_results_data_table.row_name_iter():
        #     for column_name in expected_results_data_table.column_name_iter():
        #         v = expected_results_data_table[row_name, column_name]
        #         print("{}, {}: {} ({})".format(row_name, column_name, v, type(v)))
        obs_results = self.pdm.standardized_effect_size_mean_pairwise_distance(
                assemblage_memberships=self.assemblage_memberships,
                num_randomization_replicates=100,
                is_weighted_edge_distances=True,
                is_normalize_by_tree_size=False,
                )
        self.assertEqual(len(obs_results), expected_results_data_table.num_rows())
        for obs_result, expected_result_row_name in zip(obs_results, expected_results_data_table.row_name_iter()):
            self.assertTrue(self._low_precision_equality(
                    obs_result.obs,
                    expected_results_data_table[expected_result_row_name, "mpd.obs"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.null_model_mean,
                    expected_results_data_table[expected_result_row_name, "mpd.rand.mean"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.null_model_sd,
                    expected_results_data_table[expected_result_row_name, "mpd.rand.sd"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.z,
                    expected_results_data_table[expected_result_row_name, "mpd.obs.z"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.p,
                    expected_results_data_table[expected_result_row_name, "mpd.obs.p"],
                    ))

    def test_nonabundance_edgeweighted_unnormalized_ses_mntd(self):
        # suppressMessages(library(picante))
        # dists = as.matrix(read.csv("data/dist.csv",header=T,row.names=1))
        # comm = as.matrix(read.csv("data/community.data.tsv",sep="\t",header=T,row.names=1))
        # results.mntd = ses.mntd(comm, dists, null.model="taxa.labels",abundance.weighted=F,runs=100000)
        # write.csv(format(results.mntd, digits=22), quote=F, "community.data.weighted.unnormalized.ses.mntd.csv")
        # results.mntd = ses.mntd(comm, dists, null.model="taxa.labels",abundance.weighted=F,runs=100000)
        # write.csv(format(results.mntd,digits=22), quote=F, "community.data.weighted.unnormalized.ses.mntd.csv")
        with open(pathmap.other_source_path("community.data.weighted.unnormalized.ses.mntd.csv")) as src:
            expected_results_data_table = container.DataTable.from_csv(src, default_data_type=float, delimiter=",")
        # for row_name in expected_results_data_table.row_name_iter():
        #     for column_name in expected_results_data_table.column_name_iter():
        #         v = expected_results_data_table[row_name, column_name]
        #         print("{}, {}: {} ({})".format(row_name, column_name, v, type(v)))
        obs_results = self.pdm.standardized_effect_size_mean_nearest_taxon_distance(
                assemblage_memberships=self.assemblage_memberships,
                num_randomization_replicates=100,
                is_weighted_edge_distances=True,
                is_normalize_by_tree_size=False,
                )
        self.assertEqual(len(obs_results), expected_results_data_table.num_rows())
        for obs_result, expected_result_row_name in zip(obs_results, expected_results_data_table.row_name_iter()):
            self.assertTrue(self._low_precision_equality(
                    obs_result.obs,
                    expected_results_data_table[expected_result_row_name, "mntd.obs"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.null_model_mean,
                    expected_results_data_table[expected_result_row_name, "mntd.rand.mean"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.null_model_sd,
                    expected_results_data_table[expected_result_row_name, "mntd.rand.sd"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.z,
                    expected_results_data_table[expected_result_row_name, "mntd.obs.z"],
                    ))
            self.assertTrue(self._low_precision_equality(
                    obs_result.p,
                    expected_results_data_table[expected_result_row_name, "mntd.obs.p"],
                    ))

class PhylogeneticDistanceMatrixReader(unittest.TestCase):

    def setUp(self):
        self.s00 = """\
         0  , 1  , 2  , 3  , 4  , 5
         6  , 0  , 7  , 8  , 9  , 10
         11 , 12 , 0  , 13 , 14 , 15
         16 , 17 , 18 ,  0 , 19 , 20
         21 , 22 , 23 , 24 , 0  , 25
         26 , 27 , 28 , 29 , 30 , 0
        """
        self.s11 = """\
        .  , a  , b  , c  , d  , e  , f
        a  , 0  , 1  , 2  , 3  , 4  , 5
        b  , 6  , 0  , 7  , 8  , 9  , 10
        c  , 11 , 12 , 0  , 13 , 14 , 15
        d  , 16 , 17 , 18 ,  0 , 19 , 20
        e  , 21 , 22 , 23 , 24 , 0  , 25
        f  , 26 , 27 , 28 , 29 , 30 , 0
        """
        self.s10 = """\
        a  , 0  , 1  , 2  , 3  , 4  , 5
        b  , 6  , 0  , 7  , 8  , 9  , 10
        c  , 11 , 12 , 0  , 13 , 14 , 15
        d  , 16 , 17 , 18 ,  0 , 19 , 20
        e  , 21 , 22 , 23 , 24 , 0  , 25
        f  , 26 , 27 , 28 , 29 , 30 , 0
        """
        self.s01 = """\
        a  , b  , c  , d  , e  , f
        0  , 1  , 2  , 3  , 4  , 5
        6  , 0  , 7  , 8  , 9  , 10
        11 , 12 , 0  , 13 , 14 , 15
        16 , 17 , 18 ,  0 , 19 , 20
        21 , 22 , 23 , 24 , 0  , 25
        26 , 27 , 28 , 29 , 30 , 0
        """
        self.expected_labels = {
            0 : "a",
            1 : "b",
            2 : "c",
            3 : "d",
            4 : "e",
            5 : "f",
        }
        self.expected_distances = {
            frozenset([0,1]): 1,
            frozenset([0,2]): 2,
            frozenset([0,3]): 3,
            frozenset([0,4]): 4,
            frozenset([0,5]): 5,
            frozenset([1,2]): 7,
            frozenset([1,3]): 8,
            frozenset([1,4]): 9,
            frozenset([1,5]): 10,
            frozenset([2,2]): 0,
            frozenset([2,3]): 13,
            frozenset([2,4]): 14,
            frozenset([2,5]): 15,
            frozenset([3,4]): 19,
            frozenset([3,5]): 20,
            frozenset([4,5]): 25,
        }

    def check_pdm(self, pdm, is_check_labels=True):
        self.assertEqual(len(pdm.taxon_namespace), 6)
        for i1 in range(6):
            t1 = pdm.taxon_namespace[i1]
            if is_check_labels:
                self.assertEqual(t1.label, self.expected_labels[i1])
            for i2 in range(6):
                t2 = pdm.taxon_namespace[i2]
                if is_check_labels:
                    self.assertEqual(t2.label, self.expected_labels[i2])
                if i1 == i2:
                    self.assertIs(t1, t2)
                    self.assertEqual(pdm.patristic_distance(t1, t2), 0)
                else:
                    self.assertEqual(pdm.patristic_distance(t1, t2), self.expected_distances[frozenset([i1, i2])])

    def test_read_new_taxon_namespace_with_no_row_and_no_column_names(self):
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s00),
                is_first_row_column_names=False,
                is_first_column_row_names=False,
                is_allow_new_taxa=True)
        self.check_pdm(pdm, is_check_labels=False)

    def test_read_new_taxon_namespace_with_row_and_no_column_names(self):
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s10),
                is_first_row_column_names=False,
                is_first_column_row_names=True,
                is_allow_new_taxa=True)
        self.check_pdm(pdm, is_check_labels=True)

    def test_read_new_taxon_namespace_with_no_row_and_column_names(self):
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s01),
                is_first_row_column_names=True,
                is_first_column_row_names=False,
                is_allow_new_taxa=True)
        self.check_pdm(pdm, is_check_labels=True)

    def test_read_new_taxon_namespace_with_row_and_column_names(self):
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s11),
                is_first_row_column_names=True,
                is_first_column_row_names=True,
                is_allow_new_taxa=True)
        self.check_pdm(pdm, is_check_labels=False)

    def test_read_existing_taxon_namespace_with_no_row_and_no_column_names(self):
        taxon_namespace = dendropy.TaxonNamespace()
        for label in self.expected_labels.values():
            t1 = taxon_namespace.require_taxon(label=label)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s00),
                taxon_namespace=taxon_namespace,
                is_first_row_column_names=False,
                is_first_column_row_names=False,
                is_allow_new_taxa=True)
        self.assertIs(taxon_namespace, pdm.taxon_namespace)
        self.check_pdm(pdm, is_check_labels=False)

    def test_read_existing_taxon_namespace_with_no_row_and_column_names(self):
        taxon_namespace = dendropy.TaxonNamespace()
        for label in self.expected_labels.values():
            t1 = taxon_namespace.require_taxon(label=label)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s01),
                taxon_namespace=taxon_namespace,
                is_first_row_column_names=True,
                is_first_column_row_names=False,
                is_allow_new_taxa=True)
        self.assertIs(taxon_namespace, pdm.taxon_namespace)
        self.check_pdm(pdm, is_check_labels=True)

    def test_read_existing_taxon_namespace_with_row_and_no_column_names(self):
        taxon_namespace = dendropy.TaxonNamespace()
        for label in self.expected_labels.values():
            t1 = taxon_namespace.require_taxon(label=label)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s10),
                taxon_namespace=taxon_namespace,
                is_first_row_column_names=False,
                is_first_column_row_names=True,
                is_allow_new_taxa=True)
        self.assertIs(taxon_namespace, pdm.taxon_namespace)
        self.check_pdm(pdm, is_check_labels=True)

    def test_read_existing_taxon_namespace_with_row_and_column_names(self):
        taxon_namespace = dendropy.TaxonNamespace()
        for label in self.expected_labels.values():
            t1 = taxon_namespace.require_taxon(label=label)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                StringIO(self.s11),
                taxon_namespace=taxon_namespace,
                is_first_row_column_names=True,
                is_first_column_row_names=True,
                is_allow_new_taxa=True)
        self.assertIs(taxon_namespace, pdm.taxon_namespace)
        self.check_pdm(pdm, is_check_labels=True)

class PdmTreeChecker(object):

    def check_tree(self, obs_tree, expected_tree):
        bipartitions1 = obs_tree.encode_bipartitions()
        bipartitions2 = expected_tree.encode_bipartitions()
        self.assertEqual(len(bipartitions1), len(bipartitions2))
        b1 = [b.split_as_bitstring() for b in bipartitions1]
        b2 = [b.split_as_bitstring() for b in bipartitions2]
        self.assertEqual(set(b1), set(b2))
        self.assertEqual(set(bipartitions1), set(bipartitions2))
        # try:
        #     self.assertEqual(set(bipartitions1), set(bipartitions2))
        # except AssertionError as e:
        #     print(e)
        #     print(e.__dict__)
        #     raise
        for b1 in expected_tree.bipartition_edge_map:
            self.assertIn(b1, obs_tree.bipartition_edge_map)
            self.assertAlmostEqual(
                    expected_tree.bipartition_edge_map[b1].length,
                    obs_tree.bipartition_edge_map[b1].length,
                    7,
                    "{}: {} != {}".format(b1.leafset_as_newick_string(obs_tree.taxon_namespace), expected_tree.bipartition_edge_map[b1].length, obs_tree.bipartition_edge_map[b1].length,)
                    )

class PdmNeighborJoiningTree(PdmTreeChecker, unittest.TestCase):

    def test_njtree_from_distance_matrices(self):

        # library(ape)
        # ---
        # z = matrix( c(0,5,9,9,8, 5,0,10,10,9, 9,10,0,8,7, 9,10,8,0,3, 8,9,7,3,0), byrow=T, nrow=5)
        # rownames(z)  <- c("a", "b", "c", "d", "e")
        # colnames(z)  <- c("a", "b", "c", "d", "e")
        # t = nj(z)
        # write.tree(t)
        # ---
        # p1 = read.csv("pythonidae.mle.unweighted.pdm.csv", header=T, row.names=1)
        # m1 = as.matrix(p1)
        # nj(m1)
        # t = nj(m1)
        # write.tree(t)
        test_runs = [
                ("wpnjex.csv", "(e:1,d:2,((a:2,b:3):3,c:4):2);"),
                ("saitou_and_nei_1987_table1.csv", "(h:6,g:2,((((a:5,b:2):2,c:1):1,d:3):2,(e:1,f:4):2):1);"),
                ("pythonidae.mle.unweighted.pdm.csv", "(Morelia_spilota:1,Morelia_bredli:1,((((((Morelia_kinghorni:1,Morelia_nauta:1):1,Morelia_clastolepis:1):1,Morelia_amethistina:1):1,Morelia_tracyae:1):1,Morelia_oenpelliensis:1):1,(((((Liasis_albertisii:1,Bothrochilus_boa:1):1,((Antaresia_melanocephalus:1,Antaresia_ramsayi:1):1,((Liasis_fuscus:1,Liasis_mackloti:1):1,(Apodora_papuana:1,Liasis_olivaceus:1):1):1):1):1,Morelia_boeleni:1):1,((Python_timoriensis:1,Python_reticulatus:1):1,((((Python_sebae:1,Python_molurus:1):1,Python_curtus:1):1,Python_regius:1):1,((Xenopeltis_unicolor:1,Candoia_aspera:1):1,Loxocemus_bicolor:1):1):1):1):1,((((Antaresia_stimsoni:1,Antaresia_childreni:1):1,Antaresia_perthensis:1):1,Antaresia_maculosa:1):1,((Morelia_viridisN:1,Morelia_viridisS:1):1,Morelia_carinata:1):1):1):1):1);"),
                ("pythonidae.mle.weighted.pdm.csv", "((Liasis_albertisii:0.0542142498,Bothrochilus_boa:0.0638595214):0.038444,(((Apodora_papuana:0.0670782319,Liasis_olivaceus:0.0430801028):0.010168,(Liasis_fuscus:0.0194903208,Liasis_mackloti:0.0141916418):0.048505):0.013422,(Antaresia_melanocephalus:0.0380695554,Antaresia_ramsayi:0.0325474267):0.043626):0.007734,(((((((Antaresia_stimsoni:0.0152390165,Antaresia_childreni:0.023141749):0.032397,Antaresia_perthensis:0.0760812159):0.012848,Antaresia_maculosa:0.0679212061):0.011617,((Morelia_viridisN:0.0377499268,Morelia_viridisS:0.0473589755):0.027329,Morelia_carinata:0.0660356718):0.013482):0.015469,((((((Morelia_kinghorni:0.0075825724,Morelia_nauta:0.0086155842):0.004182,Morelia_clastolepis:0.0045446653):0.018597,Morelia_amethistina:0.0227641045):0.007181,Morelia_tracyae:0.0377936102):0.024796,Morelia_oenpelliensis:0.0579745143):0.004283,(Morelia_bredli:0.0274921037,Morelia_spilota:0.0241663426):0.026356):0.031732):0.006602,(((((Python_sebae:0.0629755585,Python_molurus:0.0335903967):0.02165,Python_curtus:0.1067094932):0.016163,Python_regius:0.1058922755):0.032743,((Xenopeltis_unicolor:0.1983677797,Candoia_aspera:0.4092923305):0.048508,Loxocemus_bicolor:0.2627888765):0.060789):0.030952,(Python_timoriensis:0.074479767,Python_reticulatus:0.0562613055):0.06004):0.027099):0.002859,Morelia_boeleni:0.0843874314):0.002713);"),
                ]
        for data_filename, expected_tree_str in test_runs:
            with open(pathmap.other_source_path(data_filename)) as src:
                pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                        src,
                        is_first_row_column_names=True,
                        is_first_column_row_names=True,
                        is_allow_new_taxa=True,
                        delimiter=",")
            obs_tree = pdm.nj_tree()
            # print(obs_tree.as_string("newick"))
            # print(obs_tree.as_ascii_plot(plot_metric="length"))
            expected_tree = dendropy.Tree.get(
                    data=expected_tree_str,
                    schema="newick",
                    rooting="force-unrooted",
                    taxon_namespace=pdm.taxon_namespace,
                    preserve_underscores=True)
            self.check_tree(obs_tree=obs_tree,
                    expected_tree=expected_tree)

    def test_njtree_from_weighted_and_unweighted_distances(self):

        tree = dendropy.Tree.get(path=pathmap.tree_source_path(
            "pythonidae.mle.nex"),
            schema="nexus",
            preserve_underscores=True)
        pdm = tree.phylogenetic_distance_matrix()
        test_runs = [
                (False,  "(Morelia_spilota:1,Morelia_bredli:1,((((((Morelia_kinghorni:1,Morelia_nauta:1):1,Morelia_clastolepis:1):1,Morelia_amethistina:1):1,Morelia_tracyae:1):1,Morelia_oenpelliensis:1):1,(((((Liasis_albertisii:1,Bothrochilus_boa:1):1,((Antaresia_melanocephalus:1,Antaresia_ramsayi:1):1,((Liasis_fuscus:1,Liasis_mackloti:1):1,(Apodora_papuana:1,Liasis_olivaceus:1):1):1):1):1,Morelia_boeleni:1):1,((Python_timoriensis:1,Python_reticulatus:1):1,((((Python_sebae:1,Python_molurus:1):1,Python_curtus:1):1,Python_regius:1):1,((Xenopeltis_unicolor:1,Candoia_aspera:1):1,Loxocemus_bicolor:1):1):1):1):1,((((Antaresia_stimsoni:1,Antaresia_childreni:1):1,Antaresia_perthensis:1):1,Antaresia_maculosa:1):1,((Morelia_viridisN:1,Morelia_viridisS:1):1,Morelia_carinata:1):1):1):1):1);"),
                (True,   "((Liasis_albertisii:0.0542142498,Bothrochilus_boa:0.0638595214):0.038444,(((Apodora_papuana:0.0670782319,Liasis_olivaceus:0.0430801028):0.010168,(Liasis_fuscus:0.0194903208,Liasis_mackloti:0.0141916418):0.048505):0.013422,(Antaresia_melanocephalus:0.0380695554,Antaresia_ramsayi:0.0325474267):0.043626):0.007734,(((((((Antaresia_stimsoni:0.0152390165,Antaresia_childreni:0.023141749):0.032397,Antaresia_perthensis:0.0760812159):0.012848,Antaresia_maculosa:0.0679212061):0.011617,((Morelia_viridisN:0.0377499268,Morelia_viridisS:0.0473589755):0.027329,Morelia_carinata:0.0660356718):0.013482):0.015469,((((((Morelia_kinghorni:0.0075825724,Morelia_nauta:0.0086155842):0.004182,Morelia_clastolepis:0.0045446653):0.018597,Morelia_amethistina:0.0227641045):0.007181,Morelia_tracyae:0.0377936102):0.024796,Morelia_oenpelliensis:0.0579745143):0.004283,(Morelia_bredli:0.0274921037,Morelia_spilota:0.0241663426):0.026356):0.031732):0.006602,(((((Python_sebae:0.0629755585,Python_molurus:0.0335903967):0.02165,Python_curtus:0.1067094932):0.016163,Python_regius:0.1058922755):0.032743,((Xenopeltis_unicolor:0.1983677797,Candoia_aspera:0.4092923305):0.048508,Loxocemus_bicolor:0.2627888765):0.060789):0.030952,(Python_timoriensis:0.074479767,Python_reticulatus:0.0562613055):0.06004):0.027099):0.002859,Morelia_boeleni:0.0843874314):0.002713);"),
                ]
        for is_weighted_edge_distances, expected_tree_str in test_runs:
            obs_tree = pdm.nj_tree(is_weighted_edge_distances=is_weighted_edge_distances)
            expected_tree = dendropy.Tree.get(
                    data=expected_tree_str,
                    schema="newick",
                    rooting="force-unrooted",
                    taxon_namespace=pdm.taxon_namespace,
                    preserve_underscores=True)
            self.check_tree(obs_tree=obs_tree,
                    expected_tree=expected_tree)

class PdmUpgmaTree(PdmTreeChecker, unittest.TestCase):

    def test_upgma_average_from_distance_matrices(self):
        # library(phangorn)
        # d = read.csv("wpupgmaex.csv", header=T, row.names=1)
        # upgma(d)
        test_runs = [
                ("wpupgmaex.csv", "((e:11,(a:8.5,b:8.5):2.5):5.5,(c:14,d:14):2.5);"),
                ("pythonidae.mle.weighted.pdm.csv", "(Candoia_aspera:0.3358679339,(Loxocemus_bicolor:0.239139024,(Xenopeltis_unicolor:0.2306593655,((Python_regius:0.1021235458,(Python_curtus:0.0883212354,(Python_sebae:0.0482829776,Python_molurus:0.0482829776):0.0400382578):0.01380231042):0.04278465647,((Python_timoriensis:0.06537053625,Python_reticulatus:0.06537053625):0.06029549731,(((Liasis_albertisii:0.0590368856,Bothrochilus_boa:0.0590368856):0.03300696511,(Morelia_boeleni:0.08681248898,((Antaresia_melanocephalus:0.03530849105,Antaresia_ramsayi:0.03530849105):0.04351804164,((Liasis_fuscus:0.0168409813,Liasis_mackloti:0.0168409813):0.04845559302,(Apodora_papuana:0.05507916735,Liasis_olivaceus:0.05507916735):0.01021740697):0.01352995836):0.007985956296):0.005231361731):0.005812651357,(((Morelia_tracyae:0.0359450459,(Morelia_amethistina:0.02553168923,(Morelia_clastolepis:0.0084128718,(Morelia_kinghorni:0.0080990783,Morelia_nauta:0.0080990783):0.0003137935):0.01711881743):0.01041335667):0.02235606786,(Morelia_oenpelliensis:0.05722136872,(Morelia_bredli:0.02582922315,Morelia_spilota:0.02582922315):0.03139214558):0.001079745035):0.03700402529,((Morelia_carinata:0.06795956147,(Morelia_viridisN:0.04255445115,Morelia_viridisS:0.04255445115):0.02540511032):0.01460551598,(Antaresia_maculosa:0.07026059995,(Antaresia_perthensis:0.06383429933,(Antaresia_stimsoni:0.01919038275,Antaresia_childreni:0.01919038275):0.04464391658):0.006426300625):0.0123044775):0.01274006159):0.002551363025):0.02780953149):0.01924216872):0.08575116319):0.008479658536):0.09672890989);"),
                ("laurasiatherian.distances.ml.csv", "(Platypus:0.115554082,((Opposum:0.04554785647,(Bandicoot:0.03760589713,(Wallaroo:0.02994721074,Possum:0.02994721074):0.00765868639):0.007941959338):0.05561264108,(Elephant:0.09767314505,(Tenrec:0.09197201988,(Hedghog:0.08893681608,((Cebus:0.07429923126,(Baboon:0.06398456889,Human:0.06398456889):0.01031466238):0.01347573751,((Mouse:0.04952979734,Vole:0.04952979734):0.03752885067,(Gymnure:0.08423117145,((GuineaPig:0.0740040685,CaneRat:0.0740040685):0.009033085883,(Armadillo:0.0815558576,((Squirrel:0.0719024199,Dormouse:0.0719024199):0.006011912967,(Loris:0.07535539083,((Rabbit:0.06082665549,Pika:0.06082665549):0.01312116941,(LongTBat:0.07253285264,(Aardvark:0.07207716472,(FruitBat:0.06375844723,((((FurSeal:0.03311262951,(HarbSeal:0.004907721762,GraySeal:0.004907721762):0.02820490774):0.01108478546,(Cat:0.04153849221,Dog:0.04153849221):0.002658922751):0.01112227127,((Mole:0.04702531576,Shrew:0.04702531576):0.006612877503,((Rbat:0.0483680062,(FlyingFox:0.01426606479,RyFlyFox:0.01426606479):0.03410194142):0.002903547772,(Pig:0.04899261297,((Horse:0.009073986173,Donkey:0.009073986173):0.02905609806,(WhiteRhino:0.02269245743,IndianRhin:0.02269245743):0.0154376268):0.01086252874):0.002278941003):0.002366639287):0.001681492971):0.004794138467,(Alpaca:0.05736084145,(Hippo:0.05532140845,((Cow:0.02807794212,Sheep:0.02807794212):0.02602799145,(SpermWhale:0.03397215079,(FinWhale:0.01328705847,BlueWhale:0.01328705847):0.02068509233):0.02013378278):0.001215474875):0.002039433007):0.002752983249):0.003644622527):0.008318717494):0.0004556879157):0.001414972262):0.001407565924):0.00255894204):0.003641524737):0.001481296783):0.001194017065):0.002827476566):0.0007163207591):0.001161847302):0.003035203803):0.005701125173):0.0034873525):0.01439358448);"),
                ## note:following fails, probably due to different arbitrary resolutions of equal distances
                # ("pythonidae.mle.unweighted.pdm.csv", "((((Morelia_carinata:1.5,(Morelia_viridisN:1,Morelia_viridisS:1):0.5):1.458333333,((Antaresia_stimsoni:1,Antaresia_childreni:1):0.75,(Antaresia_maculosa:1.5,Antaresia_perthensis:1.5):0.25):1.208333333):1.416666667,((Morelia_bredli:1,Morelia_spilota:1):2.166666667,((Morelia_clastolepis:1.5,(Morelia_kinghorni:1,Morelia_nauta:1):0.5):0.8333333333,(Morelia_oenpelliensis:1.75,(Morelia_tracyae:1.5,Morelia_amethistina:1.5):0.25):0.5833333333):0.8333333333):1.208333333):0.6861111111,(((Morelia_boeleni:2,(Liasis_albertisii:1,Bothrochilus_boa:1):1):0.8333333333,((Antaresia_melanocephalus:1,Antaresia_ramsayi:1):1.5,((Apodora_papuana:1,Liasis_olivaceus:1):1,(Liasis_fuscus:1,Liasis_mackloti:1):1):0.5):0.3333333333):1.888888889,(((Python_sebae:1,Python_molurus:1):0.75,(Python_regius:1.5,Python_curtus:1.5):0.25):1.275,((Python_timoriensis:1,Python_reticulatus:1):1.833333333,(Loxocemus_bicolor:1.5,(Xenopeltis_unicolor:1,Candoia_aspera:1):0.5):1.333333333):0.1916666667):1.697222222):0.3388888889);"),
                ]
        for data_filename, expected_tree_str in test_runs:
            with open(pathmap.other_source_path(data_filename)) as src:
                pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                        src,
                        is_first_row_column_names=True,
                        is_first_column_row_names=True,
                        is_allow_new_taxa=True,
                        delimiter=",")
            obs_tree = pdm.upgma_tree()
            expected_tree = dendropy.Tree.get(
                    data=expected_tree_str,
                    schema="newick",
                    rooting="force-rooted",
                    taxon_namespace=pdm.taxon_namespace,
                    preserve_underscores=True)
            self.check_tree(obs_tree=obs_tree,
                    expected_tree=expected_tree)

class NodeToNodeDistancesTest(unittest.TestCase):

    def test_distances(self):
        ## get distances from ape
        # library(ape)
        # tr = read.nexus("pythonidae.mle.nex")
        # tr$node.label <- (Ntip(tr)+1):(nrow(tr$edge)+1)
        # tr$tip.label <- (1:Ntip(tr))
        # write.tree(tr)
        # d = dist.nodes(tr)
        # write.csv(d, "file.csv")
        test_runs = [
                ("hiv1.newick", True, "hiv1.node-to-node-dists.csv"),
                ("pythonidae.mle.numbered-nodes.newick", True, "pythonidae.mle.node-to-node-dists.csv"),
                ("hiv1.newick", False, "hiv1.unweighted.node-to-node-dists.csv"),
                ("pythonidae.mle.numbered-nodes.newick", False, "pythonidae.mle.unweighted.node-to-node-dists.csv"),
                ]
        for tree_filename, is_weighted, distances_filename in test_runs:
            tree = dendropy.Tree.get_from_path(
                    src=pathmap.tree_source_path(tree_filename),
                    schema='newick',
                    suppress_leaf_node_taxa=True)
            ndm = tree.node_distance_matrix()
            reference_table = container.DataTable.from_csv(
                    src=open(pathmap.other_source_path(distances_filename)),
                    default_data_type=float,
                    delimiter=",")
            for nd1 in tree.postorder_node_iter():
                for nd2 in tree.postorder_node_iter():
                    d = ndm.distance(nd1, nd2, is_weighted_edge_distances=is_weighted)
                    e = reference_table[nd1.label, nd2.label]
                    self.assertAlmostEqual(d, e)

    def test_mrca(self):
        test_runs = [
                "hiv1.newick",
                "pythonidae.mle.numbered-nodes.newick",
                ]
        for tree_filename in test_runs:
            tree = dendropy.Tree.get_from_path(
                    src=pathmap.tree_source_path(tree_filename),
                    schema='newick',
                    rooting="force-rooted")
            tree.encode_bipartitions()
            ndm = tree.node_distance_matrix()
            for nd1 in tree.postorder_node_iter():
                for nd2 in tree.postorder_node_iter():
                    leafset_bitmask = nd1.leafset_bitmask | nd2.leafset_bitmask
                    exp_mrca = tree.mrca(leafset_bitmask=leafset_bitmask)
                    obs_mrca = ndm.mrca(nd1, nd2)
                    # print("{} | {} = {} ({})".format(
                    #     nd1.leafset_bitmask,
                    #     nd2.leafset_bitmask,
                    #     leafset_bitmask,
                    #     obs_mrca.edge.bipartition.leafset_bitmask))
                    self.assertIs(exp_mrca, obs_mrca)

if __name__ == "__main__":
    unittest.main()

