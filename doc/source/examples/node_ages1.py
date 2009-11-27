#! /usr/bin/env python

import dendropy

trees = dendropy.TreeList.get_from_path("pythonidae.beast-mcmc.trees",
        "nexus",
        from_index=200)
for idx, tree in enumerate(d.trees_blocks[0]):
    tree.add_ages_to_nodes()
    node1 = tree.find_node_with_taxon_label(label="Antaresia maculosa")
    node2 = tree.find_node_with_taxon_label(label="Antaresia childreni")
    mrca = dendropy.Tree.ancestor(node1, node2)
    age =
