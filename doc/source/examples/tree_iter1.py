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
source_files = [
        open("pythonidae.mcmc1.nex", "r"), # Note: for 'Tree.yield_from_files',
        open("pythonidae.mcmc2.nex", "r"), # sources can be specified as file
        "pythonidae.mcmc3.nex", "r",       # objects or strings, with strings
        "pythonidae.mcmc4.nex", "r",       # assumed to specify file paths
        ]
tree_yielder = dendropy.Tree.yield_from_files(
        files=source_files,
        schema='nexus',
        taxon_namespace=taxa,
        )
for tree_idx, mcmc_tree in enumerate(tree_yielder):
    if tree_idx < burnin:
        # skip burnin
        continue
    distances.append(treecompare.symmetric_difference(mle_tree, mcmc_tree))
print("Mean symmetric distance between MLE and MCMC trees: %d"
        % float(sum(distances)/len(distances)))
