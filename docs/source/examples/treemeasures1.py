import collections
import dendropy
from dendropy.calculate import treemeasure
from dendropy.calculate import statistics

# Since we do not want to waste memory by keeping the actual trees around
# after we are done calculating the statistics, we use the tree yielder
# instead of:
#       dendropy.TreeList.get(
#           path="pythonidae.beast-mcmc.trees",
#           schema="nexus",
#           tree_offset=200)

tree_stats = collections.defaultdict(list)
for tree_idx, tree in enumerate(dendropy.Tree.yield_from_files(
            files=["pythonidae.beast-mcmc.trees"],
            schema="nexus")):
    if tree_idx < 200:
        continue # burnin
    tree_stats["B1"].append(treemeasure.B1(tree))
    tree_stats["colless"].append(treemeasure.colless_tree_imbalance(tree))
    tree_stats["PBH"].append(treemeasure.pybus_harvey_gamma(tree))
    tree_stats["sackin"].append(treemeasure.sackin_index(tree))
    tree_stats["treeness"].append(treemeasure.treeness(tree))

for key in tree_stats:
    values = tree_stats[key]
    mean, var = statistics.mean_and_sample_variance(values)
    hpd = statistics.empirical_hpd(values)
    print("{:15}: mean = {}, variance = {}, hpd = ({}, {})".format(key, mean, var, hpd[0], hpd[1]))
