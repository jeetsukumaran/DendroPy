#! /usr/bin/env python

import random
import dendropy

def process_node(node, start=1.0):
    if node.parent_node is None:
        node.value = start
    else:
        node.value = random.gauss(node.parent_node.value, node.edge.length)
    for child in node.child_nodes():
        process_node(child)
    if node.taxon is not None:
        print("%s : %s" % (node.taxon, node.value))

mle = dendropy.Tree.get(
    path='pythonidae.mle.nex',
    schema='nexus')
process_node(mle.seed_node)

