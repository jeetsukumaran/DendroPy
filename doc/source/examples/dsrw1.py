#! /usr/bin/env python

import dendropy
trees = dendropy.TreeList()
trees.read(path='pythonidae.mb.run1.t', schema='nexus', tree_offset=10)
trees.read(path='pythonidae.mb.run2.t', schema='nexus', tree_offset=10)
trees.read(path='pythonidae.mb.run3.t', schema='nexus', tree_offset=10)
trees.read(path='pythonidae.mb.run4.t', schema='nexus', tree_offset=10)
ds = dendropy.DataSet([trees])
ds.read(path='pythonidae_cytb.fasta',
        schema='fasta',
        data_type='dna',
        taxon_namespace=ds.taxon_namespaces[0])
ds.write(path='pythonidae_combined.nex', schema='nexus')
