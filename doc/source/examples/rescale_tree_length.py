#! /usr/bin/env python

import dendropy

mle = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
mle_len = mle.length()
for edge in mle.postorder_edge_iter():
    edge.length = float(edge.length)/mle_len
print(mle.as_string("newick"))
