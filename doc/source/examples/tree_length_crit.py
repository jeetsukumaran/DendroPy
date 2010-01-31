#! /usr/bin/env python

import dendropy

trees = dendropy.TreeList.get_from_path("pythonidae.random.bd0301.tre", "nexus")
tree_lengths = [tree.length() for tree in trees]
tree_lengths.sort()
crit_index_95 = int(0.95 * len(tree_lengths))
crit_index_99 = int(0.99 * len(tree_lengths))

print("95%% critical value: %s" % tree_lengths[crit_index_95])
print("99%% critical value: %s" % tree_lengths[crit_index_99])
