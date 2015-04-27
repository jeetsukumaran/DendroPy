#! /usr/bin/env python

import dendropy
from dendropy.calculate import treecompare

distances = []
taxa = dendropy.TaxonNamespace()
mle_tree = dendropy.Tree.get(
    path='pythonidae.mle.nex',
    schema='nexus',
    taxon_namespace=taxa)
burnin = 20
tree_yielder = dendropy.Tree.yield_from_files(
        files=[open('pythonidae.mcmc.nex', 'r')],
        schema='nexus',
        taxon_namespace=taxa,
        )
for tree_idx, mcmc_tree in enumerate(tree_yielder):
    if tree_idx < burnin:
        continue
    distances.append(treecompare.symmetric_difference(mle_tree, mcmc_tree))
print("Mean symmetric distance between MLE and MCMC trees: %d"
        % float(sum(distances)/len(distances)))
