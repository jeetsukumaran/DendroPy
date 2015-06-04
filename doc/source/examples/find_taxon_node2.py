#! /usr/bin/env python

import dendropy

tree = dendropy.Tree.get(path="pythonidae.mle.nex", schema="nexus")
node = tree.find_node_with_taxon_label('Antaresia maculosa')
print(node.description())
