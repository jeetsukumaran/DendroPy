#! /usr/bin/env python

import dendropy

trees = dendropy.TreeList()
for tree_file in ['pythonidae_cytb.mb.run1.t',
        'pythonidae_cytb.mb.run2.t',
        'pythonidae_cytb.mb.run3.t',
        'pythonidae_cytb.mb.run4.t']:
    trees.read_from_path(
            tree_file,
            'nexus',
            tree_offset=20)
con_tree = trees.consensus(min_freq=0.95)
print(con_tree.as_string('newick'))
