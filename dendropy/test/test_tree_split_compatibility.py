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
import json
from dendropy.test.support import pathmap
from dendropy.test.support import dendropytest

def generate_bipartitions(
        trees,
        generation_mode="from-trees",
        is_rooted=None):
    if is_rooted is None:
        is_rooted = trees[0].is_rooted
    if generation_mode == "all":
        bipartitions = []
        all_taxa_bitmask = trees.taxon_namespace.all_taxa_bitmask()
        for s in range(1, all_taxa_bitmask):
            b = dendropy.Bipartition(
                    leafset_bitmask=s,
                    tree_leafset_bitmask=all_taxa_bitmask,
                    is_rooted=is_rooted,
                    compile_bipartition=True)
            bipartitions.append(b)
    elif generation_mode == "from-trees":
        bipartitions = []
        for tree in trees:
            bipartitions.extend(tree.encode_bipartitions())
    else:
        raise ValueError(generation_mode)
    unique_bipartitions = {}
    for b in bipartitions:
        if is_rooted:
            unique_bipartitions[b.leafset_bitmask] = b
        else:
            unique_bipartitions[b.split_bitmask] = b
    unique_bipartitions_keys = sorted(unique_bipartitions.keys())
    bipartitions = []
    for key in unique_bipartitions_keys:
        unique_bipartitions[key].is_mutable = False
        bipartitions.append(unique_bipartitions[key])
    return bipartitions

class TestSplitCompatibility(dendropytest.ExtendedTestCase):

    def test_compatibility(self):
        regimes = (
            ("dendropy-test-trees-n12-x2.nexus", "all"),
            ("dendropy-test-trees-n33-unrooted-x100a.nexus", "from-trees"),
            ("dendropy-test-trees-n10-rooted-treeshapes.nexus", "all"),
        )
        for trees_filename_idx, (trees_filename, bipartition_generation_mode) in enumerate(regimes):
            trees_filepath = pathmap.tree_source_path(trees_filename)
            trees = dendropy.TreeList.get_from_path(
                    trees_filepath,
                    "nexus",)
            bipartitions = generate_bipartitions(trees, bipartition_generation_mode, is_rooted=trees[0].is_rooted)
            # for bipartition1_idx, bipartition1 in enumerate(bipartitions):
            for bipartition1_idx, bipartition1 in enumerate(bipartitions):
                for tree_idx, tree in enumerate(trees):
                    compatible_bipartitions = set()
                    incompatible_bipartitions = set()
                    bipartition_encoding = tree.encode_bipartitions()
                    for biparition2_idx, bipartition2 in enumerate(bipartition_encoding):
                        if bipartition2.is_compatible_with(bipartition1):
                            self.assertTrue(bipartition1.is_compatible_with(bipartition2))
                            compatible_bipartitions.add(bipartition2)
                        else:
                            self.assertFalse(bipartition1.is_compatible_with(bipartition2))
                            incompatible_bipartitions.add(bipartition2)
                    is_compatible = tree.is_compatible_with_bipartition(bipartition1)
                    self.assertEqual(len(compatible_bipartitions) + len(incompatible_bipartitions), len(bipartition_encoding))
                    if is_compatible:
                        self.assertEqual(len(incompatible_bipartitions), 0,
                                "Tree {} of '{}': bipartition {} (leafset = {}, index = {}) found compatible with tree, but is incompatible with following bipartitions on tree: {}".
                                format(
                                    tree_idx,
                                    trees_filename,
                                    bipartition1.split_as_bitstring(),
                                    bipartition1.leafset_as_bitstring(),
                                    bipartition1_idx,
                                    [b.split_as_bitstring() for b in incompatible_bipartitions],
                                    ))
                        self.assertEqual(len(compatible_bipartitions), len(bipartition_encoding))
                    else:
                        self.assertTrue(len(incompatible_bipartitions) > 0,
                                "Tree {} of '{}': bipartition {} (leafset = {}, index = {}) found incompatible with tree, but is compatible with all bipartitions on tree: {}".
                                format(
                                    tree_idx,
                                    trees_filename,
                                    bipartition1.split_as_bitstring(),
                                    bipartition1.leafset_as_bitstring(),
                                    bipartition1_idx,
                                    [b.split_as_bitstring() for b in compatible_bipartitions],
                                    ))



    # def test_compatibility(self):
    #     regime_source_names = (
    #             ("dendropy-test-trees-n10-rooted-treeshapes.nexus", "dendropy-test-trees-n10-rooted-treeshapes.split-compatibilities-all.json", ),
    #     )
    #     for tree_filename, ref_filename in regime_source_names:

    #         with open(pathmap.tree_source_path(ref_filename), "r") as src:
    #             regimed = json.load(src)

    #         bipartitions = {}
    #         for bidx in regimed["bipartitions"]:
    #             bdesc = regimed["bipartitions"][bidx]
    #             b = dendropy.Bipartition(
    #                     leafset_bitmask=bdesc["leafset_bitmask"],
    #                     tree_leafset_bitmask=bdesc["tree_leafset_bitmask"],
    #                     is_rooted=bdesc["is_rooted"])
    #             bipartitions[int(bidx)] = b
    #         bipartition_indexes = sorted(bipartitions.keys())

    #         tree_data = {}
    #         for tree_idx in regimed["trees"]:
    #             i = int(tree_idx)
    #             tree_data[i] = {}
    #             tree_data[i]["compatible_bipartitions"] = set([int(j) for j in regimed["trees"][tree_idx]["compatible_bipartitions"]])
    #             tree_data[i]["incompatible_bipartitions"] = set([int(j) for j in regimed["trees"][tree_idx]["incompatible_bipartitions"]])

    #         trees = dendropy.TreeList.get_from_path(
    #                 pathmap.tree_source_path(tree_filename),
    #                 "nexus")
    #         for tree_idx, tree in enumerate(trees):
    #             tree_bipartitions = tree.encode_bipartitions()
    #             for bipartition_idx in bipartition_indexes:
    #                 is_compatible = tree.is_compatible_with_bipartition(bipartition)

    #                 bipartition = bipartitions[bipartition_idx]
    #                 if bipartition_idx in tree_data[tree_idx]["compatible_bipartitions"]:
    #                     expected_is_compatible = True
    #                 elif bipartition_idx in tree_data[tree_idx]["incompatible_bipartitions"]:
    #                     expected_is_compatible = False
    #                 else:
    #                     raise Exception("Bipartition {} ('{}') not reported for tree {}".format(bipartition_idx, bipartition.split_as_bitstring(), tree_idx))
    #                 self.assertIs(is_compatible, expected_is_compatible,
    #                         "Bipartition {} ('{}') for tree {} in '{}': expected '{}' for compatibility but instead observed '{}'".format(
    #                             bipartition_idx,
    #                             bipartition.split_as_bitstring(),
    #                             tree_idx,
    #                             tree_filename,
    #                             expected_is_compatible,
    #                             is_compatible,
    #                             ))

if __name__ == "__main__":
    unittest.main()
