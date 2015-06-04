#! /usr/bin/env python

import dendropy
from dendropy.calculate import treemeasure

trees = dendropy.TreeList.get(
        path="pythonidae.beast-mcmc.trees",
        schema="nexus",
        tree_offset=200)
maculosa_childreni_ages = []
for idx, tree in enumerate(trees):
    tree.calc_node_ages()
    taxon_labels = ["Antaresia maculosa","Antaresia childreni"]
    mrca = tree.mrca(taxon_labels=taxon_labels)
    maculosa_childreni_ages.append(mrca.age)
print("Mean age of MRCA of 'Antaresia maculosa' and 'Antaresia childreni': %s" \
    % (float(sum(maculosa_childreni_ages))/len(maculosa_childreni_ages)))


