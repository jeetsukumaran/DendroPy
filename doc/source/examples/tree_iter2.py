#! /usr/bin/env python

import dendropy
from dendropy import multi_tree_source_iter
from dendropy import treecalc

distances = []
taxa = dendropy.TaxonSet()
mle_tree = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus', taxon_set=taxa)
mcmc_tree_file_paths = ['pythonidae.mb.run1.t',
        'pythonidae.mb.run2.t',
        'pythonidae.mb.run3.t',
        'pythonidae.mb.run4.t']
for mcmc_tree in multi_tree_source_iter(
        mcmc_tree_file_paths,
        schema='nexus',
        taxon_set=taxa):
    distances.append(treecalc.symmetric_difference(mle_tree, mcmc_tree))
print("Mean symmetric distance between MLE and MCMC trees: %d"
        % float(sum(distances)/len(distances)))
