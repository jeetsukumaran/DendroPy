#! /usr/bin/env python

import dendropy

mle = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
multifurcating = lambda x: True if len(x.child_nodes()) > 2 else False
for nd in mle.postorder_node_iter(multifurcating):
    print(nd.description(0))
