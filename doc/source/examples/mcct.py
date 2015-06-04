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

mcct = trees.maximum_product_of_split_support_tree()
print("\nTree {} maximizes the product of split support (log product = {}): {}".format(trees.index(mcct), mcct.log_product_of_split_support, mcct))

msst = trees.maximum_sum_of_split_support_tree()
print("\nTree {} maximizes the sum of split support (sum = {}): {}".format(trees.index(msst), msst.sum_of_split_support, msst))
