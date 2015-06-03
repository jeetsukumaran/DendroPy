#! /usr/bin/env python

import dendropy

mle = dendropy.Tree.get(path='pythonidae.mle.nex', schema='nexus')
mle_len = mle.length()
for edge in mle.postorder_edge_iter():
    if edge.length is None:
        edge.length = 0
    else:
        edge.length = float(edge.length)/mle_len
print(mle.as_string(schema="newick"))
