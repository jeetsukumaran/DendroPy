#! /usr/bin/env python

import dendropy
from dendropy.calculate import treecompare

distances = []
taxa = dendropy.TaxonNamespace()
mle_tree = dendropy.Tree.get(
    path='pythonidae.mle.nex',
    schema='nexus',
    taxon_namespace=taxa)
mcmc_tree_file_paths = ['pythonidae.mb.run1.t',
        'pythonidae.mb.run2.t',
        'pythonidae.mb.run3.t',
        'pythonidae.mb.run4.t']
for mcmc_tree in dendropy.Tree.yield_from_files(
        files=mcmc_tree_file_paths,
        schema='nexus',
        taxon_namespace=taxa):
    distances.append(treecompare.symmetric_difference(mle_tree, mcmc_tree))
print("Mean symmetric distance between MLE and MCMC trees: %d"
        % float(sum(distances)/len(distances)))
