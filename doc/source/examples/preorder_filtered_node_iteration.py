#! /usr/bin/env python
# -*- coding: utf-8 -*-

import dendropy

mle = dendropy.Tree.get(
    path='pythonidae.mle.nex',
    schema='nexus')
multifurcating = lambda x: True if len(x.child_nodes()) > 2 else False
for nd in mle.postorder_node_iter(multifurcating):
    print(nd.description(0))
