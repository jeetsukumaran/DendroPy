#! /usr/bin/env python

import dendropy

trees = dendropy.TreeList.get_from_path("pythonidae.beast-mcmc.trees",
        "nexus",
        tree_offset=200)
maculosa_childreni_ages = []
for idx, tree in enumerate(trees):
    tree.calc_node_ages()
    node1 = tree.find_node_with_taxon_label(label="Antaresia maculosa")
    node2 = tree.find_node_with_taxon_label(label="Antaresia childreni")
    mrca = dendropy.Tree.mrca(node1, node2)
    maculosa_childreni_ages.append(mrca.age)
print("Mean age of MRCA of 'Antaresia maculosa' and 'Antaresia childreni': %s" \
    % (float(sum(maculosa_childreni_ages))/len(maculosa_childreni_ages)))


