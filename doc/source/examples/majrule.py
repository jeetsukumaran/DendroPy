#! /usr/bin/env python

import dendropy

trees = dendropy.TreeList()
for tree_file in ['pythonidae.mb.run1.t',
        'pythonidae.mb.run2.t',
        'pythonidae.mb.run3.t',
        'pythonidae.mb.run4.t']:
    trees.read(
            path=tree_file,
            schema='nexus',
            tree_offset=20)
con_tree = trees.consensus(min_freq=0.95)
print(con_tree.as_string(schema='newick'))
