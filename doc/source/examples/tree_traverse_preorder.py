#! /usr/bin/env python

import dendropy

mle = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
for node in mle.preorder_node_iter():
    print("%s: %d child node(s)" % (node.description(0), len(node.child_nodes())))
