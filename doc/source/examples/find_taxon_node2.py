#! /usr/bin/env python

import dendropy

tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus")
node = tree.find_node_with_taxon_label('Antaresia maculosa')
print(node.description())
