#! /usr/bin/env python

import dendropy

tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus")
filter = lambda taxon: True if taxon.label=='Antaresia maculosa' else False
node = tree.find_node_with_taxon(filter)
print(node.description())
