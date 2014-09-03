#! /usr/bin/env python

import dendropy
from dendropy import tree_source_iter
from dendropy import treecompare

distances = []
taxa = dendropy.TaxonNamespace()
mle_tree = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus', taxon_namespace=taxa)
for mcmc_tree in tree_source_iter(
        stream=open('pythonidae.mcmc.nex', 'rU'),
        schema='nexus',
        taxon_namespace=taxa,
        tree_offset=200):
    distances.append(treecompare.symmetric_difference(mle_tree, mcmc_tree))
print("Mean symmetric distance between MLE and MCMC trees: %d"
        % float(sum(distances)/len(distances)))
