#! /usr/bin/env python

import dendropy

mle = dendropy.Tree.get(
    path='pythonidae.mle.nex',
    schema='nexus')
mle_len = mle.length()
for edge in mle.postorder_edge_iter():
    edge.length = None
print(mle.as_string(schema="newick"))
