#! /usr/bin/env python

import unittest
import dendropy
from dendropy.test.support import pathmap

class LineagesThroughTimeTest(unittest.TestCase):

    def check_tree(self, tree, dist, num_lineages):
        working_tree = dendropy.Tree(tree)
        working_tree.truncate_from_root(dist)
        num_wt_leaves = len(working_tree.leaf_nodes())
        self.assertEqual(num_lineages, num_wt_leaves)

    def test_ultrametric(self):
        trees = dendropy.TreeList.get_from_stream(pathmap.tree_source_stream("pythonidae.reference-trees.nexus"), "nexus")
        for tree in trees:
            dists = tree.calc_node_root_distances()
            min_dist, max_dist = tree.minmax_leaf_distance_from_root()
            trunc_dists = [(max_dist * f) for f in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
            for td in trunc_dists:
                num_lineages = tree.num_lineages_at(td)
                self.check_tree(tree, td, num_lineages)

if __name__ == "__main__":
    unittest.main()



