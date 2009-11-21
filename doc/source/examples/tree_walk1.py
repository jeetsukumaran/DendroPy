#! /usr/bin/env python

import dendropy

def visit_node(node):
    children = node.child_nodes()
    print("%s: %d child node(s)" % (node.description(0), len(children)))
    for child in children:
        visit_node(child)

mle = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
visit_node(mle.seed_node)
