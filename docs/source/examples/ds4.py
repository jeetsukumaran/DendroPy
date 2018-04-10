import dendropy
treelist1 = dendropy.TreeList.get(
        path='pythonidae.mle.nex',
        schema='nexus')
cytb = dendropy.DnaCharacterMatrix.get(
    path='pythonidae_cytb.fasta',
    schema='fasta')
ds = dendropy.DataSet([cytb, treelist1])
ds.unify_taxon_namespaces()

