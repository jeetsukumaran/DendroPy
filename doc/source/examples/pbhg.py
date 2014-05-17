#! /usr/bin/env python

import dendropy

trees = dendropy.TreeList.get_from_path("pythonidae.beast-mcmc.trees",
        "nexus",
        tree_offset=200)
pbhg = []
for idx, tree in enumerate(trees):
    pbhg.append(tree.pybus_harvey_gamma())
print("Mean Pybus-Harvey-Gamma: %s" \
    % (float(sum(pbhg))/len(pbhg)))


