#! /usr/bin/env python

import dendropy

mle = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
short = lambda edge: True if edge.length < 0.01 else False
for edge in mle.postorder_edge_iter(short):
    print(edge.length)
