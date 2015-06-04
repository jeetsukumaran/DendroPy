#! /usr/bin/env python

import dendropy
from dendropy.calculate import treemeasure

trees = dendropy.TreeList.get(
        path="pythonidae.beast-mcmc.trees",
        schema="nexus",
        tree_offset=200)
pbhg = []
for idx, tree in enumerate(trees):
    pbhg.append(treemeasure.pybus_harvey_gamma(tree))
print("Mean Pybus-Harvey-Gamma: %s" \
    % (float(sum(pbhg))/len(pbhg)))


