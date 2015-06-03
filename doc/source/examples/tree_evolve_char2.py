#! /usr/bin/env python

import random
import dendropy

def evolve_char(tree, start=1.0):
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.value = 1.0
        else:
            node.value = random.gauss(node.parent_node.value, node.edge.length)
    return tree

mle = dendropy.Tree.get(
        path='pythonidae.mle.nex',
        schema='nexus')
evolve_char(mle)
for node in mle.leaf_iter():
    print("%s : %s" % (node.taxon, node.value))
