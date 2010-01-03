#! /usr/bin/env python

import dendropy
trees = dendropy.TreeList()
trees.read_from_path('pythonidae_cytb.mb.run1.t', 'nexus', tree_offset=20)
trees.read_from_path('pythonidae_cytb.mb.run2.t', 'nexus', tree_offset=20)
trees.read_from_path('pythonidae_cytb.mb.run3.t', 'nexus', tree_offset=20)
trees.read_from_path('pythonidae_cytb.mb.run4.t', 'nexus', tree_offset=20)
ds = dendropy.DataSet(trees)
ds.read_from_path('pythonidae_cytb.fasta', 'dnafasta', taxon_set=ds.taxon_sets[0])
ds.write_to_path('pythonidae_combined.nex', 'nexus')
