import dendropy
taxa = dendropy.TaxonNamespace()
trees = dendropy.TreeList(taxon_namespace=taxa)
trees.read(path='pythonidae.mb.run1.t', schema='nexus', tree_offset=10)
trees.read(path='pythonidae.mb.run2.t', schema='nexus', tree_offset=10)
trees.read(path='pythonidae.mb.run3.t', schema='nexus', tree_offset=10)
trees.read(path='pythonidae.mb.run4.t', schema='nexus', tree_offset=10)
ds = dendropy.DataSet([trees])
ds.read(path='pythonidae_cytb.fasta',
        schema='fasta',
        data_type='dna',
        )
ds.write(path='pythonidae_combined.nex', schema='nexus')
