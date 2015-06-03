import dendropy
ds = dendropy.DataSet()
treelist1 = dendropy.TreeList.get(
        path='pythonidae.mle.nex',
        schema='nexus')
cytb = dendropy.DnaCharacterMatrix.get(
    path='pythonidae_cytb.fasta',
    schema='fasta')
ds.add(treelist1)
ds.add(cytb)
ds.unify_taxon_namespaces()
