#! /usr/bin/env python

import dendropy

tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus")
node1 = tree.find_node_with_taxon_label('Antaresia maculosa')
node2 = tree.find_node_with_taxon_label("Antaresia childreni")
mrca = tree.ancestor(node1, node2)
print(mrca.description())
