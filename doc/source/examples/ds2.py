import dendropy

ds1 = dendropy.DataSet()
taxon_namespace = dendropy.TaxonNamespace()
ds1.attach_taxon_namespace(taxon_namespace)
ds1.read(
    path='pythonidae.mle.nex',
    schema='nexus',)
ds1.read(
    path='pythonidae.chars.nexus',
    schema='nexus',)

ds2 = dendropy.DataSet(ds1)

# for tx1, tx2 in zip(ds1, ds2):
#     print(tx1 is tx2)
